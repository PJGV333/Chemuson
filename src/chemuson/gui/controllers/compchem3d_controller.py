from __future__ import annotations

"""Async controller for the 3D / CompChem dock."""

from dataclasses import dataclass
from typing import Any

from PyQt6.QtCore import QObject, QThread, pyqtSignal, pyqtSlot

from chemuson.core.model import MolGraph
from chemuson.geometry3d import (
    CoordinateSet3D,
    OptimizationFrame,
    OptimizationResult,
    OptimizationSettings,
    conformer_3d_for_graph,
)


@dataclass(frozen=True)
class CompChemJobSpec:
    """Description of a 3D/compchem worker job."""

    operation: str
    backend: str
    settings: OptimizationSettings
    force: bool = False


class CompChem3DWorker(QObject):
    """Runs 3D generation/optimization away from the UI thread."""

    frame_ready = pyqtSignal(int, object)
    finished = pyqtSignal(int, object)

    def __init__(
        self,
        job_id: int,
        graph: MolGraph,
        spec: CompChemJobSpec,
        coordset: CoordinateSet3D | None = None,
    ) -> None:
        super().__init__()
        self._job_id = int(job_id)
        self._graph = graph
        self._spec = spec
        self._coordset = coordset

    @pyqtSlot()
    def run(self) -> None:
        try:
            if self._spec.operation == "generate":
                result = self._generate()
            else:
                result = self._optimize()
        except Exception as exc:
            result = OptimizationResult(
                None,
                method=self._spec.backend,
                message=f"Fallo de backend: {exc}",
            )
        self.finished.emit(self._job_id, result)

    def _generate(self) -> OptimizationResult:
        if self._spec.backend != "rdkit":
            return OptimizationResult(
                None,
                method=self._spec.backend,
                message="La generación de conformero usa RDKit aislado; Open Babel se usa para optimizar si está disponible.",
            )
        conformer = conformer_3d_for_graph(
            self._graph,
            timeout_s=self._spec.settings.timeout_s,
            force=self._spec.force,
            settings=self._spec.settings,
        )
        coordset = conformer.coordinate_set
        if coordset is None:
            return OptimizationResult(
                None,
                method=conformer.method or "rdkit",
                message=conformer.message or "RDKit no devolvió conformero.",
                metadata={"cache_state": "miss"},
            )
        cache_state = "hit" if conformer.from_cache else "miss"
        return OptimizationResult(
            coordset,
            converged=True,
            energy=coordset.energy,
            method=coordset.method,
            message=conformer.message,
            metadata={"cache_state": cache_state, "cache_key": conformer.cache_key},
        )

    def _optimize(self) -> OptimizationResult:
        if self._spec.backend == "openbabel":
            from chemuson.geometry3d import obabel_backend

            result = obabel_backend.optimize(self._graph, self._coordset, self._spec.settings)
            metadata = dict(getattr(result, "metadata", {}) or {})
            metadata.setdefault("cache_state", "miss")
            return OptimizationResult(
                result.coordinates,
                converged=result.converged,
                energy=result.energy,
                method=result.method,
                message=result.message,
                frames=result.frames,
                metadata=metadata,
            )

        from chemuson.geometry3d import rdkit_backend

        frames: list[OptimizationFrame] = []
        for frame in rdkit_backend.optimize_iter(self._graph, self._coordset, self._spec.settings):
            frames.append(frame)
            self.frame_ready.emit(self._job_id, frame)
        if frames:
            last = frames[-1]
            return OptimizationResult(
                last.coordinates,
                converged=last.converged,
                energy=last.energy,
                method=last.coordinates.method,
                message=last.message,
                frames=tuple(frames),
                metadata={"cache_state": "miss"},
            )
        result = rdkit_backend.optimize(self._graph, self._coordset, self._spec.settings)
        metadata = dict(getattr(result, "metadata", {}) or {})
        metadata.setdefault("cache_state", "miss")
        return OptimizationResult(
            result.coordinates,
            converged=result.converged,
            energy=result.energy,
            method=result.method,
            message=result.message,
            frames=result.frames,
            metadata=metadata,
        )


class CompChem3DController(QObject):
    """Owns async compchem jobs and relays Qt signals."""

    job_started = pyqtSignal(int)
    frame_ready = pyqtSignal(int, object)
    job_finished = pyqtSignal(int, object)

    def __init__(self, parent: QObject | None = None) -> None:
        super().__init__(parent)
        self._next_job_id = 1
        self._jobs: dict[int, tuple[QThread, CompChem3DWorker]] = {}
        self._pending_results: dict[int, Any] = {}

    def start_job(
        self,
        graph: MolGraph,
        spec: CompChemJobSpec,
        coordset: CoordinateSet3D | None = None,
    ) -> int:
        job_id = self._next_job_id
        self._next_job_id += 1
        thread = QThread(self)
        worker = CompChem3DWorker(job_id, graph, spec, coordset)
        worker.moveToThread(thread)
        thread.started.connect(worker.run)
        worker.frame_ready.connect(self.frame_ready)
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
        return tuple(sorted(self._jobs))

    def _record_worker_finished(self, job_id: int, result: Any) -> None:
        self._pending_results[int(job_id)] = result

    def _on_thread_finished(self, job_id: int) -> None:
        job_id = int(job_id)
        self._jobs.pop(job_id, None)
        result = self._pending_results.pop(job_id, None)
        if result is not None:
            self.job_finished.emit(job_id, result)
