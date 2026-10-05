"""Asynchronous GUI adapter for the provider-neutral molecular assistant."""

from __future__ import annotations

from collections.abc import Callable
from typing import Any

from PyQt6.QtCore import QObject, QThread, pyqtSignal, pyqtSlot

from chemuson.molecular_assistant import (
    MolecularAssistant,
    MolecularAssistantRequest,
    MolecularAssistantResult,
    MolecularAssistantStatus,
    OpenAICompatibleConfig,
    OPENAI_COMPATIBLE_PROFILES,
    OpenAICompatibleProvider,
    get_openai_compatible_profile,
)


ResultGenerator = Callable[[str, OpenAICompatibleConfig], MolecularAssistantResult]


def _generate_structure(
    description: str,
    config: OpenAICompatibleConfig,
) -> MolecularAssistantResult:
    provider = OpenAICompatibleProvider(config)
    return MolecularAssistant(provider).generate(MolecularAssistantRequest(description))


def _provider_failure(config: OpenAICompatibleConfig) -> MolecularAssistantResult:
    return MolecularAssistantResult(
        status=MolecularAssistantStatus.PROVIDER_ERROR,
        provider_id=config.provider_id,
        model_id=config.model,
        reason_code="provider_error",
    )


class _MolecularAssistantWorker(QObject):
    """Execute provider I/O and isolated ChemIO validation on a QThread."""

    finished = pyqtSignal(int, object)

    def __init__(
        self,
        job_id: int,
        description: str,
        config: OpenAICompatibleConfig,
        generator: ResultGenerator,
    ) -> None:
        super().__init__()
        self._job_id = int(job_id)
        self._description = description
        self._config = config
        self._generator = generator

    @pyqtSlot()
    def run(self) -> None:
        try:
            result = self._generator(self._description, self._config)
        except Exception:
            result = _provider_failure(self._config)
        if not isinstance(result, MolecularAssistantResult):
            result = _provider_failure(self._config)
        self.finished.emit(self._job_id, result)


class MolecularAssistantController(QObject):
    """Own worker threads and relay only typed, stable M23 results to the GUI."""

    job_started = pyqtSignal(int)
    job_finished = pyqtSignal(int, object)

    def __init__(
        self,
        parent: QObject | None = None,
        *,
        generator: ResultGenerator | None = None,
    ) -> None:
        super().__init__(parent)
        self._generator = generator or _generate_structure
        self._next_job_id = 1
        self._jobs: dict[int, tuple[QThread, _MolecularAssistantWorker]] = {}
        self._pending_results: dict[int, MolecularAssistantResult] = {}
        self._abandoned_jobs: set[int] = set()
        self._shutting_down = False

    @property
    def provider_profiles(self):
        """Expose immutable profile metadata to the UI without a direct M23 import."""
        return OPENAI_COMPATIBLE_PROFILES

    def start_job(
        self,
        description: str,
        *,
        base_url: str,
        model: str,
        api_key: str | None = None,
        supports_json_output: bool = False,
        provider_id: str = "openai-compatible",
    ) -> int | None:
        """Validate explicit per-request configuration and start background work."""
        if self._shutting_down:
            return None
        if not isinstance(description, str) or not description.strip():
            return None
        profile = get_openai_compatible_profile(provider_id)
        if profile is None:
            return None
        if profile.api_key_required and (
            not isinstance(api_key, str) or not api_key.strip()
        ):
            return None
        try:
            config = OpenAICompatibleConfig(
                base_url=base_url.strip() if isinstance(base_url, str) else base_url,
                model=model.strip() if isinstance(model, str) else model,
                api_key=api_key if api_key else None,
                supports_json_output=supports_json_output,
                provider_id=profile.profile_id,
            )
        except (TypeError, ValueError):
            return None

        job_id = self._next_job_id
        self._next_job_id += 1
        thread = QThread(self)
        thread.setObjectName(f"MolecularAssistant-{job_id}")
        worker = _MolecularAssistantWorker(
            job_id,
            description,
            config,
            self._generator,
        )
        worker.moveToThread(thread)
        thread.started.connect(worker.run)
        worker.finished.connect(self._record_worker_finished)
        worker.finished.connect(thread.quit)
        worker.finished.connect(worker.deleteLater)
        thread.finished.connect(lambda job_id=job_id: self._on_thread_finished(job_id))
        thread.finished.connect(thread.deleteLater)
        self._jobs[job_id] = (thread, worker)
        self.job_started.emit(job_id)
        thread.start()
        return job_id

    def active_jobs(self) -> tuple[int, ...]:
        """Return job identifiers whose QThreads have not completed."""
        return tuple(sorted(self._jobs))

    def has_active_jobs(self) -> bool:
        """Whether any owned QThread is still running."""
        return any(thread.isRunning() for thread, _worker in self._jobs.values())

    @property
    def shutdown_complete(self) -> bool:
        """Whether shutdown began and every owned QThread has stopped."""
        return self._shutting_down and not self.has_active_jobs()

    def begin_shutdown(self) -> None:
        """Suppress late results and request cooperative worker interruption."""
        if self._shutting_down:
            return
        self._shutting_down = True
        for job_id, (thread, _worker) in self._jobs.items():
            self._abandoned_jobs.add(job_id)
            thread.requestInterruption()

    def abandon_job(self, job_id: int) -> None:
        """Ignore a late result without claiming to abort an in-flight HTTP call."""
        job_id = int(job_id)
        if job_id in self._jobs:
            self._abandoned_jobs.add(job_id)

    @pyqtSlot(int, object)
    def _record_worker_finished(self, job_id: int, result: Any) -> None:
        if self._shutting_down:
            return
        if isinstance(result, MolecularAssistantResult):
            self._pending_results[int(job_id)] = result

    def _on_thread_finished(self, job_id: int) -> None:
        job_id = int(job_id)
        self._jobs.pop(job_id, None)
        result = self._pending_results.pop(job_id, None)
        if self._shutting_down or job_id in self._abandoned_jobs:
            self._abandoned_jobs.discard(job_id)
            return
        if result is not None:
            self.job_finished.emit(job_id, result)
