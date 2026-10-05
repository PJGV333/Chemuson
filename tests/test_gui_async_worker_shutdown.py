from __future__ import annotations

import threading
import time

from PyQt6.QtCore import QObject, QPointF, QThread, pyqtSignal, pyqtSlot
from PyQt6.QtTest import QSignalSpy
from PyQt6.QtWidgets import QApplication

from chemuson.core.model import MolGraph
from chemuson.gui import main_window as main_window_module
from chemuson.gui.canvas import canvas_structure
from chemuson.gui.controllers import compchem3d_controller, template_controller
from chemuson.gui.controllers.compchem3d_controller import CompChemJobSpec
from chemuson.gui.controllers.template_controller import TemplateControllerContext
from chemuson.gui.main_window import ChemusonWindow
from chemuson.geometry3d import OptimizationSettings
from chemuson.molecular_assistant import (
    MolecularAssistantResult,
    MolecularAssistantStatus,
)


def _wait_for(predicate, timeout_s: float = 8.0) -> bool:
    deadline = time.monotonic() + timeout_s
    while time.monotonic() < deadline:
        QApplication.processEvents()
        if predicate():
            return True
        time.sleep(0.002)
    QApplication.processEvents()
    return bool(predicate())


def _single_carbon_graph() -> MolGraph:
    graph = MolGraph()
    graph.add_atom("C", 0.0, 0.0)
    return graph


def _allow_close(window: ChemusonWindow) -> None:
    window._confirm_discard_changes = lambda _canvas: True
    window._apply_pending_portable_update_on_exit = lambda: True
    window._apply_pending_windows_update_on_exit = lambda: True


def test_close_cancel_does_not_begin_worker_shutdown(monkeypatch):
    entered = threading.Event()
    release = threading.Event()
    property_updates = []
    thread_finished = threading.Event()

    class HeldDescriptorWorker(QObject):
        finished = pyqtSignal(int, dict, str)

        def __init__(self, job_id, graph):
            super().__init__()
            self.job_id = job_id

        @pyqtSlot()
        def run(self):
            entered.set()
            release.wait(5.0)
            self.finished.emit(self.job_id, {"logp": 1.0}, "")

    monkeypatch.setattr(main_window_module, "_DescriptorWorker", HeldDescriptorWorker)
    window = ChemusonWindow()
    window._properties_update_timer.stop()
    _allow_close(window)
    window._confirm_discard_changes = lambda _canvas: False
    monkeypatch.setattr(
        window.chemical_properties_dock,
        "update_properties",
        lambda rows: property_updates.append(rows),
    )
    window.show()
    window._start_descriptor_job(window.canvas.model, [])
    thread = window._descriptor_jobs[1][0]
    thread.finished.connect(thread_finished.set)

    try:
        assert entered.wait(3.0)
        assert not window.close()
        assert not window._shutdown_started
        assert not thread.isInterruptionRequested()

        release.set()
        assert _wait_for(lambda: thread_finished.is_set() and bool(property_updates))
        assert property_updates
        assert window.isVisible()
        _allow_close(window)
        assert window.close()
    finally:
        release.set()
        if window._shutdown_started and not window._async_shutdown_complete:
            _wait_for(lambda: window._async_shutdown_complete)
        window.close()
        window.deleteLater()
        QApplication.processEvents()


def test_close_waits_for_all_descendant_worker_families_and_suppresses_late_ui(
    monkeypatch,
):
    started_event = threading.Event()
    started_lock = threading.Lock()
    started_count = 0
    release = threading.Event()
    message_calls = []
    descriptor_updates = []
    compchem_frames = []
    canvas_insertions = []
    template_status = []

    def hold_then_release():
        nonlocal started_count
        with started_lock:
            started_count += 1
            if started_count == 6:
                started_event.set()
        release.wait(5.0)

    class HeldDescriptorWorker(QObject):
        finished = pyqtSignal(int, dict, str)

        def __init__(self, job_id, graph):
            super().__init__()
            self.job_id = job_id

        @pyqtSlot()
        def run(self):
            hold_then_release()
            self.finished.emit(self.job_id, {"logp": 1.0}, "")

    class HeldNameWorker(QObject):
        finished = pyqtSignal(int, object, str)

        def __init__(self, job_id, query):
            super().__init__()
            self.job_id = job_id

        @pyqtSlot()
        def run(self):
            hold_then_release()
            self.finished.emit(self.job_id, None, "late test error")

    class HeldCompChemWorker(QObject):
        frame_ready = pyqtSignal(int, object)
        finished = pyqtSignal(int, object)

        def __init__(self, job_id, graph, spec, coordset=None):
            super().__init__()
            self.job_id = job_id

        @pyqtSlot()
        def run(self):
            hold_then_release()
            self.frame_ready.emit(self.job_id, object())
            self.finished.emit(self.job_id, object())

    class HeldTemplateWorker(QObject):
        finished = pyqtSignal(int, str, str)

        def __init__(self, job_id, graph):
            super().__init__()
            self.job_id = job_id

        @pyqtSlot()
        def run(self):
            hold_then_release()
            self.finished.emit(self.job_id, "C", "")

    class HeldCanvasWorker(QObject):
        finished = pyqtSignal(int, str, str)

        def __init__(self, job_id, graph, mode, name_options):
            super().__init__()
            self.job_id = job_id

        @pyqtSlot()
        def run(self):
            hold_then_release()
            self.finished.emit(self.job_id, "late analysis", "")

    monkeypatch.setattr(main_window_module, "_DescriptorWorker", HeldDescriptorWorker)
    monkeypatch.setattr(main_window_module, "_NameToStructureWorker", HeldNameWorker)
    monkeypatch.setattr(compchem3d_controller, "CompChem3DWorker", HeldCompChemWorker)
    monkeypatch.setattr(template_controller, "_SmilesExportWorker", HeldTemplateWorker)
    monkeypatch.setattr(canvas_structure, "_CanvasAnalysisWorker", HeldCanvasWorker)
    monkeypatch.setattr(
        main_window_module.QMessageBox,
        "critical",
        lambda *args: message_calls.append("critical"),
    )
    monkeypatch.setattr(
        main_window_module.QMessageBox,
        "warning",
        lambda *args: message_calls.append("warning"),
    )
    monkeypatch.setattr(
        main_window_module.QMessageBox,
        "information",
        lambda *args: message_calls.append("information"),
    )

    window = ChemusonWindow()
    window._properties_update_timer.stop()
    _allow_close(window)
    monkeypatch.setattr(
        window.chemical_properties_dock,
        "update_properties",
        lambda rows: descriptor_updates.append(rows),
    )
    monkeypatch.setattr(
        window.compchem_dock,
        "add_frame",
        lambda frame: compchem_frames.append(frame),
    )
    monkeypatch.setattr(
        window.canvas,
        "_insert_analysis_text",
        lambda *args: canvas_insertions.append(args),
    )
    window.show()

    graph = _single_carbon_graph()
    window._start_descriptor_job(graph, [])
    name_job_id = window._start_name_to_structure_job("fake offline query")
    assert name_job_id is not None

    assistant_result = MolecularAssistantResult(
        status=MolecularAssistantStatus.SUCCESS,
        provider_id="test-provider",
        model_id="test-model",
        proposed_smiles="C",
        graph=graph,
        validation_passed=True,
    )

    def generate_without_provider(_description, _config):
        hold_then_release()
        return assistant_result

    assistant = window._molecular_assistant_controller
    assistant._generator = generate_without_provider
    assistant_spy = QSignalSpy(assistant.job_finished)
    assistant.start_job(
        "test description",
        base_url="https://example.invalid/v1",
        model="test-model",
        api_key="test-only-key",
        provider_id="openai",
    )

    compchem = window._compchem3d_controller
    compchem_frames_spy = QSignalSpy(compchem.frame_ready)
    compchem_finished_spy = QSignalSpy(compchem.job_finished)
    compchem.start_job(
        graph,
        CompChemJobSpec(
            operation="optimize",
            backend="rdkit",
            settings=OptimizationSettings(timeout_s=1.0),
        ),
    )

    template_context = TemplateControllerContext(
        parent=window,
        canvas=window.canvas,
        template_library=None,
        show_status=template_status.append,
        refresh_template_views=lambda: None,
        insert_template=lambda *_args: None,
    )
    window._template_controller._start_smiles_export_job(template_context, graph)
    window.canvas._start_analysis_job(graph, None, "name", QPointF(0.0, 0.0))

    try:
        # Every patched worker reaches the gate without invoking a provider,
        # network connector, RDKit subprocess, or other external service.
        assert started_event.wait(5.0)
        active_threads = [
            thread
            for thread in window.findChildren(QThread)
            if thread.isRunning()
        ]
        assert len(active_threads) == 6

        assert not window.close()
        assert window._shutdown_started
        assert not window._async_shutdown_complete
        assert len(window._shutdown_threads) == 6
        assert window.isVisible()

        release.set()
        assert _wait_for(
            lambda: window._async_shutdown_complete and not window.isVisible()
        )
        assert not descriptor_updates
        assert not compchem_frames
        assert not canvas_insertions
        assert template_status == ["Exportando SMILES..."]
        assert not message_calls
        assert len(assistant_spy) == 0
        assert len(compchem_frames_spy) == 0
        assert len(compchem_finished_spy) == 0
        assert assistant.start_job(
            "after shutdown",
            base_url="https://example.invalid/v1",
            model="test-model",
            api_key="test-only-key",
            provider_id="openai",
        ) is None
        assert compchem.start_job(
            graph,
            CompChemJobSpec(
                operation="optimize",
                backend="rdkit",
                settings=OptimizationSettings(timeout_s=1.0),
            ),
        ) is None
    finally:
        release.set()
        if window._shutdown_started and not window._async_shutdown_complete:
            _wait_for(lambda: window._async_shutdown_complete)
        window.close()
        window.deleteLater()
        QApplication.processEvents()
