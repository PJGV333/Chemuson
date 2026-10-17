"""Regresiones del dock 3D / CompChem."""

from __future__ import annotations

import math

import pytest
from PyQt6.QtCore import QEvent, QEventLoop, QPointF, Qt, QTimer
from PyQt6.QtGui import QKeyEvent, QMouseEvent
from PyQt6.QtWidgets import QApplication, QFileDialog, QMessageBox


from chemuson.core.model import MolGraph
from chemuson.geometry3d import (
    Conformer3DResult,
    CoordinateSet3D,
    ForceField,
    OptimizationFrame,
    OptimizationResult,
    OptimizationSettings,
    SceneAtom3D,
    SceneMolecule3D,
)
from chemuson.gui.controllers import CompChem3DController, CompChem3DWorker, CompChemJobSpec
from chemuson.gui.main_window import ChemusonWindow
import chemuson.gui.controllers.compchem3d_controller as compchem_module


@pytest.fixture(scope="module", autouse=True)
def _qapp():
    return QApplication.instance() or QApplication([])


def _chain_graph() -> MolGraph:
    graph = MolGraph()
    a1 = graph.add_atom("C", 0.0, 0.0)
    a2 = graph.add_atom("C", 42.0, 0.0)
    a3 = graph.add_atom("O", 84.0, 0.0)
    graph.add_bond(a1.id, a2.id, order=1)
    graph.add_bond(a2.id, a3.id, order=1)
    return graph


def test_compchem_controller_generates_async_with_fake_backend(monkeypatch) -> None:
    graph = _chain_graph()

    def fake_conformer(received_graph, **_kwargs):
        assert len(received_graph.atoms) == 3
        return Conformer3DResult(
            {1: (0.0, 0.0, 0.0), 2: (1.0, 0.0, 0.0), 3: (2.0, 0.0, 0.0)},
            source="fake",
            cache_key="fake-key",
            method="UFF",
            energy=-1.25,
        )

    monkeypatch.setattr(compchem_module, "conformer_3d_for_graph", fake_conformer)
    controller = CompChem3DController()
    finished: list[object] = []
    loop = QEventLoop()
    controller.job_finished.connect(lambda _job_id, result: (finished.append(result), loop.quit()))

    controller.start_job(
        graph,
        CompChemJobSpec("generate", "rdkit", OptimizationSettings(forcefield=ForceField.UFF)),
    )
    QTimer.singleShot(3000, loop.quit)
    loop.exec()

    assert finished
    result = finished[0]
    assert result.ok
    assert result.coordinates is not None
    assert result.coordinates.source == "fake"


def test_compchem_generate_worker_relays_real_generation_stages(monkeypatch) -> None:
    diagnostics = {
        "error_code": None,
        "stage": "done",
        "graph": {"atom_count": 3, "bond_count": 2},
        "degradation_code": None,
    }

    def fake_conformer(_graph, **kwargs):
        kwargs["on_progress"](
            {
                "event": "progress",
                "stage": "embedding",
                "message": "Buscando una geometría ETKDG reproducible.",
                "graph": diagnostics["graph"],
            }
        )
        return Conformer3DResult(
            {1: (0.0, 0.0, 0.0), 2: (1.5, 0.0, 0.0), 3: (2.8, 0.0, 0.0)},
            source="rdkit",
            cache_key="stage-key",
            method="unoptimized",
            message="missing_forcefield_parameters: no parameters",
            metadata={"converged": False, "diagnostics": diagnostics},
        )

    monkeypatch.setattr(compchem_module, "conformer_3d_for_graph", fake_conformer)
    worker = CompChem3DWorker(
        15,
        _chain_graph(),
        CompChemJobSpec("generate", "rdkit", OptimizationSettings()),
    )
    events = []
    worker.progress_ready.connect(lambda job_id, event: events.append((job_id, event)))
    finished = []
    worker.finished.connect(lambda _job_id, result: finished.append(result))

    worker.run()

    assert len(events) == 1
    assert events[0][0] == 15
    assert events[0][1]["stage"] == "embedding"
    assert finished and finished[0].ok
    assert finished[0].coordinates.method == "unoptimized"
    assert not finished[0].converged
    assert finished[0].metadata["diagnostics"] == diagnostics


def test_compchem_generation_failure_displays_diagnostics_and_retries(monkeypatch) -> None:
    window = ChemusonWindow()
    try:
        for element in ("C", "C", "O"):
            window.canvas.model.add_atom(element, 0.0, 0.0)
        job_id = 91
        window._latest_compchem_job_id = job_id
        window._compchem_job_operations[job_id] = "generate"
        window._compchem_job_backends[job_id] = "rdkit"
        window._compchem_job_snapshots[job_id] = window._current_compchem_snapshot()
        result = OptimizationResult(
            None,
            message="RDKit did not embed the graph.",
            metadata={
                "diagnostics": {
                    "error_code": "embedding_failed",
                    "stage": "embedding",
                    "return_code": -1,
                    "timings_s": {"embedding": 0.25},
                    "stereochemistry": {"specified_atoms": {"1": "S"}},
                    "graph": {"atom_count": 3, "bond_count": 2},
                }
            },
        )
        window._on_compchem_job_finished(job_id, result)

        status = window.compchem_dock.status_label.text()
        assert "embedding_failed" in status
        assert "embedding" in status
        assert "3 átomos" in status
        assert '"return_code": -1' in window.compchem_dock.status_label.toolTip()
        assert '"stereochemistry"' in window.viewport3d_dock.status_label.toolTip()
        assert not window.viewport3d_dock.retry_button.isHidden()
        assert window.compchem_dock.btn_generate.text() == "Reintentar generación 3D"

        starts = []
        monkeypatch.setattr(
            window,
            "_start_compchem_job",
            lambda operation, **kwargs: starts.append((operation, kwargs)),
        )
        window.viewport3d_dock.retry_button.click()
        window.compchem_dock.btn_generate.click()
        assert [operation for operation, _kwargs in starts] == ["generate", "generate"]
    finally:
        window.close()


def test_compchem_discards_stale_result_and_reports_reason() -> None:
    window = ChemusonWindow()
    try:
        canvas = window.canvas
        atom = canvas.model.add_atom("C", 0.0, 0.0)
        original = CoordinateSet3D(
            {atom.id: (0.0, 0.0, 0.0)}, source="test", method="manual"
        )
        canvas.set_geometry3d_coordinates(original)
        job_id = 92
        window._latest_compchem_job_id = job_id
        window._compchem_job_operations[job_id] = "generate"
        window._compchem_job_backends[job_id] = "rdkit"
        window._compchem_job_snapshots[job_id] = window._current_compchem_snapshot()
        canvas.model.add_atom("O", 42.0, 0.0)

        window._on_compchem_job_finished(
            job_id,
            OptimizationResult(
                CoordinateSet3D(
                    {atom.id: (1.0, 0.0, 0.0)}, source="rdkit", method="ETKDGv3"
                ),
                converged=True,
            ),
        )

        assert canvas.coordinate_set_3d == original
        assert "stale_result" in window.compchem_dock.status_label.text()
        assert '"error_code": "stale_result"' in window.compchem_dock.status_label.toolTip()
        assert not window.viewport3d_dock.retry_button.isHidden()
    finally:
        window.close()


def test_scene_build_failure_reports_structured_diagnostic(monkeypatch) -> None:
    window = ChemusonWindow()
    try:
        atom = window.canvas.model.add_atom("C", 0.0, 0.0)
        window._compchem_coordset = CoordinateSet3D(
            {atom.id: (0.0, 0.0, 0.0)}, source="test", method="manual"
        )

        def fail_scene(_graph, _coordinates):
            raise RuntimeError("invalid coordinate mapping")

        monkeypatch.setattr("chemuson.gui.main_window.scene_molecule_from_graph", fail_scene)
        window._refresh_viewport3d_scene()

        assert "scene_failed" in window.viewport3d_dock.status_label.text()
        assert '"error_code": "scene_failed"' in window.viewport3d_dock.status_label.toolTip()
    finally:
        window.close()


def test_viewport_display_failure_reaches_compchem_diagnostics(_qapp):
    window = ChemusonWindow()
    try:
        atom = window.canvas.model.add_atom("C", 0.0, 0.0)
        widget = window.viewport3d_dock.viewport
        widget.set_scene(
            SceneMolecule3D(
                atoms=(SceneAtom3D(atom.id, "C", 0.0, 0.0, 0.0, 1.7, "#404040"),),
                bonds=(),
            )
        )

        class FakeGL:
            def glClearColor(self, *_args):
                pass

            def glClear(self, *_args):
                pass

        class FailedUploadRenderer:
            _gl = FakeGL()

            def upload(self, *_args):
                return False

        widget._renderer = FailedUploadRenderer()
        widget._gl_ready = True
        widget.paintGL()

        assert "display_failed" in window.viewport3d_dock.status_label.text()
        assert "display_failed" in window.compchem_dock.status_label.text()
        assert '"error_code": "display_failed"' in window.compchem_dock.status_label.toolTip()
    finally:
        window.close()


def test_compchem_stage_progress_is_visible_in_both_docks() -> None:
    window = ChemusonWindow()
    try:
        for element in ("C", "C", "O"):
            window.canvas.model.add_atom(element, 0.0, 0.0)
        job_id = 24
        window._latest_compchem_job_id = job_id
        window._compchem_job_operations[job_id] = "generate"
        window._compchem_job_snapshots[job_id] = window._current_compchem_snapshot()

        window._on_compchem_progress(
            job_id,
            {
                "event": "progress",
                "stage": "embedding",
                "message": "Buscando una geometría ETKDG reproducible.",
                "graph": {"atom_count": 3, "bond_count": 2},
            },
        )

        assert "ETKDG" in window.compchem_dock.status_label.text()
        assert "3 átomos" in window.compchem_dock.status_label.text()
        assert "ETKDG" in window.viewport3d_dock.status_label.text()
    finally:
        window._close_approved = True
        window._async_shutdown_complete = True
        window.close()


def test_compchem_worker_handles_backend_failure() -> None:
    graph = _chain_graph()
    worker = CompChem3DWorker(
        1,
        graph,
        CompChemJobSpec("generate", "openbabel", OptimizationSettings(forcefield=ForceField.UFF)),
    )
    finished: list[object] = []
    worker.finished.connect(lambda _job_id, result: finished.append(result))

    worker.run()

    assert finished
    assert not finished[0].ok
    assert "RDKit aislado" in finished[0].message


def test_compchem_worker_relays_live_frames_before_completion(monkeypatch) -> None:
    from chemuson.geometry3d import rdkit_backend

    graph = _chain_graph()
    coordset = CoordinateSet3D(
        {1: (0.0, 0.0, 0.0), 2: (1.0, 0.0, 0.0), 3: (2.0, 0.0, 0.0)},
        source="test",
        method="manual",
    )
    frame = OptimizationFrame(
        step=2,
        coordinates=CoordinateSet3D(
            {1: (0.0, 0.0, 0.0), 2: (1.1, 0.0, 0.0), 3: (2.0, 0.0, 0.0)},
            source="rdkit",
            method="UFF",
            energy=-2.5,
        ),
        energy=-2.5,
    )
    received = {}

    def fake_optimize_stream(_graph, _coordset, _settings, **kwargs):
        received.update(kwargs)
        kwargs["on_frame"](frame)
        return OptimizationResult(
            frame.coordinates,
            converged=True,
            energy=frame.energy,
            method=frame.coordinates.method,
            frames=(frame,),
        )

    monkeypatch.setattr(rdkit_backend, "optimize_stream", fake_optimize_stream)
    worker = CompChem3DWorker(
        9,
        graph,
        CompChemJobSpec(
            "optimize",
            "rdkit",
            OptimizationSettings(forcefield=ForceField.UFF),
            fixed_atom_ids=(1,),
        ),
        coordset,
    )
    events = []
    worker.frame_ready.connect(lambda job_id, value: events.append(("frame", job_id, value)))
    worker.finished.connect(lambda job_id, _result: events.append(("finished", job_id)))

    worker.run()

    assert [event[0] for event in events] == ["frame", "finished"]
    assert events[0][1:] == (9, frame)
    assert received["fixed_atom_ids"] == (1,)
    assert callable(received["on_frame"])


def test_live_mm_coalesces_latest_drag_until_old_worker_finishes() -> None:
    class FakeController:
        def __init__(self):
            self.active = [71]
            self.cancelled = []
            self.started = []

        def active_jobs(self):
            return list(self.active)

        def cancel_job(self, job_id):
            self.cancelled.append(job_id)
            self.active = []
            return True

        def start_job(self, graph, spec, coordinates):
            self.started.append((graph, spec, coordinates))
            self.active = [72]
            return 72

    window = ChemusonWindow()
    try:
        canvas = window.canvas
        atom = canvas.model.add_atom("C", 0.0, 0.0)
        canvas.set_geometry3d_coordinates(
            CoordinateSet3D(
                {atom.id: (0.0, 0.0, 0.0)}, source="test", method="manual"
            )
        )
        fake = FakeController()
        window._compchem3d_controller = fake
        window._compchem_job_backends[71] = "rdkit"
        window._compchem_job_operations[71] = "optimize"
        window._compchem_job_snapshots[71] = (
            canvas.document_token,
            canvas.document_revision - 1,
        )
        window._latest_compchem_job_id = 71
        window.viewport3d_dock.live_mm_checkbox.blockSignals(True)
        window.viewport3d_dock.live_mm_checkbox.setChecked(True)
        window.viewport3d_dock.live_mm_checkbox.blockSignals(False)

        window._start_live_mm_job()
        assert fake.cancelled == [71]
        assert fake.started == []
        assert window._live_mm_pending
        assert window._latest_compchem_job_id == 0

        window._on_compchem_job_finished(
            71,
            OptimizationResult(None, converged=False, method="MMFF94", message="cancelled"),
        )

        assert len(fake.started) == 1
        assert fake.active == [72]
        assert window._latest_compchem_job_id == 72
        assert not window._live_mm_pending
    finally:
        window._live_mm_timer.stop()
        window._close_approved = True
        window._async_shutdown_complete = True
        window.close()


def test_compchem_geometry_is_restored_per_document_on_tab_switch() -> None:
    window = ChemusonWindow()
    try:
        first = window.canvas
        atom = first.model.add_atom("C", 12.0, 24.0)
        first_coordinates = CoordinateSet3D(
            {atom.id: (0.25, 1.5, -0.5)},
            source="test",
            method="manual",
        )
        first.set_geometry3d_coordinates(first_coordinates)
        second = window._create_document_tab(make_current=False)

        window._set_active_canvas(second)
        assert window._compchem_coordset is None
        assert first.coordinate_set_3d == first_coordinates

        window._set_active_canvas(first)
        assert window._compchem_coordset == first_coordinates
        assert first.model.atoms[atom.id].x == 12.0
        assert first.model.atoms[atom.id].y == 24.0
    finally:
        window.close()


def test_2d_topology_changes_mark_3d_stale_and_undo_restores_it() -> None:
    from chemuson.gui.commands import AddAtomCommand, MoveAtomsCommand

    window = ChemusonWindow()
    try:
        canvas = window.canvas
        atom = canvas.model.add_atom("C", 0.0, 0.0)
        coordinates = CoordinateSet3D(
            {atom.id: (0.0, 0.0, 0.0)}, source="test", method="manual"
        )
        canvas.set_geometry3d_coordinates(coordinates)
        canvas.undo_stack.push(
            MoveAtomsCommand(
                canvas.model,
                canvas,
                {atom.id: (atom.x, atom.y)},
                {atom.id: (atom.x + 10.0, atom.y + 5.0)},
            )
        )
        assert not canvas.geometry3d_is_stale
        assert canvas.coordinate_set_3d == coordinates

        canvas.undo_stack.push(
            AddAtomCommand(canvas.model, canvas, "O", 80.0, 0.0, is_explicit=True)
        )
        assert canvas.geometry3d_is_stale
        assert canvas.coordinate_set_3d == coordinates
        assert not window.compchem_dock.btn_export_xyz.isEnabled()

        canvas.undo_stack.undo()
        assert not canvas.geometry3d_is_stale
        assert canvas.coordinate_set_3d == coordinates
        assert window.compchem_dock.btn_export_xyz.isEnabled()
    finally:
        window._close_approved = True
        window._async_shutdown_complete = True
        window.close()


def test_open_viewport_debounces_regeneration_after_topology_changes() -> None:
    class FakeController:
        def __init__(self):
            self.started = []
            self.cancelled = []

        def active_jobs(self):
            return []

        def cancel_job(self, job_id):
            self.cancelled.append(job_id)
            return True

        def start_job(self, graph, spec, coordinates):
            self.started.append((graph, spec, coordinates))
            return len(self.started) + 100

    from PyQt6.QtCore import QEventLoop, QTimer
    from chemuson.gui.commands import AddAtomCommand

    window = ChemusonWindow()
    try:
        canvas = window.canvas
        atom = canvas.model.add_atom("C", 0.0, 0.0)
        canvas.set_geometry3d_coordinates(
            CoordinateSet3D({atom.id: (0.0, 0.0, 0.0)}, "test", "manual")
        )
        fake = FakeController()
        window._compchem3d_controller = fake
        window.viewport3d_dock.isVisible = lambda: True
        timer = window._geometry3d_sync_timer
        timer.setInterval(20)

        for element, offset in (("O", 1), ("N", 2)):
            window.canvas.undo_stack.push(
                AddAtomCommand(
                    canvas.model,
                    canvas,
                    element,
                    float(offset * 42),
                    0.0,
                    is_explicit=True,
                )
            )

        assert canvas.geometry3d_is_stale
        assert window.viewport3d_dock.status_label.text() == "Actualizando geometría 3D…"
        loop = QEventLoop()
        timer.timeout.connect(loop.quit)
        QTimer.singleShot(100, loop.quit)
        loop.exec()

        assert len(fake.started) == 1
        assert len(fake.started[0][0].atoms) == 3
        assert fake.started[0][1].operation == "generate"

        canvas.undo_stack.undo()
        canvas.undo_stack.undo()
        assert not canvas.geometry3d_is_stale
        assert not timer.isActive()

        canvas.undo_stack.redo()
        canvas.undo_stack.redo()
        assert canvas.geometry3d_is_stale
        assert timer.isActive()
        redo_loop = QEventLoop()
        timer.timeout.connect(redo_loop.quit)
        QTimer.singleShot(100, redo_loop.quit)
        redo_loop.exec()
        assert len(fake.started) == 2
        assert len(fake.started[1][0].atoms) == 3
    finally:
        window._geometry3d_sync_timer.stop()
        window._close_approved = True
        window._async_shutdown_complete = True
        window.close()


def test_restoring_canvas_does_not_relabel_stale_geometry_as_current() -> None:
    from chemuson.gui.commands import ChangeAtomCommand

    window = ChemusonWindow()
    try:
        canvas = window.canvas
        atom = canvas.model.add_atom("C", 0.0, 0.0)
        canvas.set_geometry3d_coordinates(
            CoordinateSet3D({atom.id: (0.0, 0.0, 0.0)}, "test", "manual")
        )
        original_signature = canvas.geometry3d_graph_signature
        canvas.undo_stack.push(ChangeAtomCommand(canvas.model, canvas, atom.id, "O"))
        assert canvas.model.get_atom(atom.id).element == "O"
        assert canvas.geometry3d_is_stale

        window._restore_compchem_state_for_canvas(canvas)

        assert canvas.geometry3d_is_stale
        assert canvas.geometry3d_graph_signature == original_signature
    finally:
        window._close_approved = True
        window._async_shutdown_complete = True
        window.close()


def test_layout_edit_does_not_invalidate_running_geometry_job() -> None:
    from chemuson.gui.commands.transform_commands import MoveAtomsCommand

    class FakeController:
        def __init__(self):
            self.cancelled = []

        def cancel_job(self, job_id):
            self.cancelled.append(job_id)
            return True

    window = ChemusonWindow()
    try:
        canvas = window.canvas
        atom = canvas.model.add_atom("C", 0.0, 0.0)
        coordinates = CoordinateSet3D(
            {atom.id: (0.0, 0.0, 0.0)}, source="test", method="manual"
        )
        canvas.set_geometry3d_coordinates(coordinates)
        job_id = 53
        window._compchem3d_controller = FakeController()
        window._latest_compchem_job_id = job_id
        window._compchem_job_operations[job_id] = "generate"
        window._compchem_job_snapshots[job_id] = window._current_compchem_snapshot(canvas)

        canvas.undo_stack.push(
            MoveAtomsCommand(
                canvas.model,
                canvas,
                {atom.id: (0.0, 0.0)},
                {atom.id: (24.0, 12.0)},
            )
        )

        assert window._compchem3d_controller.cancelled == []
        assert window._compchem_job_snapshot_matches(job_id)
    finally:
        window._close_approved = True
        window._async_shutdown_complete = True
        window.close()


def test_compchem_discards_frames_after_document_edit() -> None:
    from chemuson.gui.commands import AddAtomCommand

    window = ChemusonWindow()
    try:
        canvas = window.canvas
        atom = canvas.model.add_atom("C", 12.0, 24.0)
        committed = CoordinateSet3D(
            {atom.id: (0.0, 0.0, 0.0)}, source="test", method="manual"
        )
        canvas.set_geometry3d_coordinates(committed)
        job_id = 41
        window._latest_compchem_job_id = job_id
        window._compchem_job_snapshots[job_id] = (
            canvas.document_token,
            canvas.document_revision,
        )
        window._compchem_job_backends[job_id] = "rdkit"

        canvas.undo_stack.push(
            AddAtomCommand(canvas.model, canvas, "O", 60.0, 60.0, is_explicit=True)
        )
        frame = OptimizationFrame(
            step=1,
            coordinates=CoordinateSet3D(
                {atom.id: (1.0, 0.0, 0.0)}, source="test", method="MMFF94"
            ),
            energy=1.0,
        )
        window._on_compchem_frame_ready(job_id, frame)

        assert window._latest_compchem_job_id == 0
        assert window._compchem_preview_coordset is None
        assert window._compchem_coordset == committed
    finally:
        window._close_approved = True
        window._async_shutdown_complete = True
        window.close()


def test_compchem_frames_preview_without_committing_to_document(monkeypatch) -> None:
    window = ChemusonWindow()
    try:
        canvas = window.canvas
        atom = canvas.model.add_atom("C", 12.0, 24.0)
        committed = CoordinateSet3D(
            {atom.id: (0.0, 0.0, 0.0)}, source="test", method="manual"
        )
        preview = CoordinateSet3D(
            {atom.id: (0.5, 0.0, 0.0)}, source="test", method="MMFF94"
        )
        canvas.set_geometry3d_coordinates(committed)
        job_id = 42
        window._latest_compchem_job_id = job_id
        window._compchem_job_snapshots[job_id] = (
            canvas.document_token,
            canvas.document_revision,
        )
        frame = OptimizationFrame(step=1, coordinates=preview, energy=2.5)
        refreshes = []
        monkeypatch.setattr(window, "_refresh_viewport3d_scene", lambda: refreshes.append(True))
        monkeypatch.setattr(window.viewport3d_dock, "isVisible", lambda: True)

        window._on_compchem_frame_ready(job_id, frame)

        assert window._compchem_preview_coordset == preview
        assert canvas.coordinate_set_3d == committed
        assert refreshes == [True]
    finally:
        window.close()


def test_compchem_final_geometry_is_undoable() -> None:
    window = ChemusonWindow()
    try:
        canvas = window.canvas
        atom = canvas.model.add_atom("C", 12.0, 24.0)
        original = CoordinateSet3D(
            {atom.id: (0.0, 0.0, 0.0)}, source="test", method="manual"
        )
        optimized = CoordinateSet3D(
            {atom.id: (0.8, 0.0, 0.0)}, source="rdkit", method="MMFF94", energy=0.5
        )
        canvas.set_geometry3d_coordinates(original)
        job_id = 43
        window._latest_compchem_job_id = job_id
        window._compchem_job_backends[job_id] = "rdkit"
        window._compchem_job_snapshots[job_id] = (
            canvas.document_token,
            canvas.document_revision,
        )
        result = OptimizationResult(
            coordinates=optimized,
            converged=True,
            energy=0.5,
            method="MMFF94",
            metadata={"energy_unit": "kcal/mol"},
        )

        window._on_compchem_job_finished(job_id, result)

        assert canvas.coordinate_set_3d == optimized
        assert canvas.undo_stack.count() == 1
        canvas.undo_stack.undo()
        assert canvas.coordinate_set_3d == original
        canvas.undo_stack.redo()
        assert canvas.coordinate_set_3d == optimized
    finally:
        window._close_approved = True
        window._async_shutdown_complete = True
        window.close()


def test_3d_builder_adds_atom_and_bond_as_one_undoable_molecular_edit() -> None:
    window = ChemusonWindow()
    try:
        canvas = window.canvas
        parent = canvas.model.add_atom("C", 10.0, 20.0)
        original = CoordinateSet3D(
            {parent.id: (0.0, 0.0, 0.0)}, source="test", method="manual"
        )
        canvas.set_geometry3d_coordinates(original)
        window._pending_3d_build_anchor_id = parent.id

        window._add_3d_atom("C", (2.0, 0.0, 0.0))

        assert len(canvas.model.atoms) == 2
        assert len(canvas.model.bonds) == 1
        built = canvas.coordinate_set_3d
        assert built is not None
        added_id = next(atom_id for atom_id in canvas.model.atoms if atom_id != parent.id)
        assert added_id in built.normalized_positions()
        assert math.dist(
            built.normalized_positions()[parent.id],
            built.normalized_positions()[added_id],
        ) == pytest.approx(1.5)
        assert canvas.undo_stack.count() == 1

        canvas.undo_stack.undo()
        assert list(canvas.model.atoms) == [parent.id]
        assert not canvas.model.bonds
        assert canvas.coordinate_set_3d == original
        canvas.undo_stack.redo()
        assert len(canvas.model.atoms) == 2
        assert canvas.coordinate_set_3d == built
    finally:
        window._close_approved = True
        window._async_shutdown_complete = True
        window.close()


def test_3d_builder_double_click_extends_with_one_single_bond_macro() -> None:
    window = ChemusonWindow()
    try:
        canvas = window.canvas
        parent = canvas.model.add_atom("C", 10.0, 20.0)
        original = CoordinateSet3D(
            {parent.id: (0.0, 0.0, 0.0)}, source="test", method="manual"
        )
        canvas.set_geometry3d_coordinates(original)
        window._pending_3d_build_anchor_id = parent.id
        window.viewport3d_dock.bond_order_combo.setCurrentIndex(2)

        window._on_viewport_build_chain_extension_requested(("C", (3.0, 0.0, 0.0)))

        assert len(canvas.model.atoms) == 2
        assert len(canvas.model.bonds) == 1
        bond = next(iter(canvas.model.bonds.values()))
        assert (bond.order, bond.is_aromatic) == (1, False)
        assert canvas.undo_stack.count() == 1
        built = canvas.coordinate_set_3d
        assert built is not None
        added_id = next(atom_id for atom_id in canvas.model.atoms if atom_id != parent.id)
        assert math.dist(
            built.normalized_positions()[parent.id],
            built.normalized_positions()[added_id],
        ) == pytest.approx(1.5)

        canvas.undo_stack.undo()
        assert list(canvas.model.atoms) == [parent.id]
        assert not canvas.model.bonds
        assert canvas.coordinate_set_3d == original
    finally:
        window._close_approved = True
        window._async_shutdown_complete = True
        window.close()

@pytest.mark.parametrize(("name", "initial_atom_count"), (("ethane", 1), ("butane", 3)))
def test_3d_builder_extends_ethane_and_butane_with_single_bonds(name, initial_atom_count):
    window = ChemusonWindow()
    try:
        canvas = window.canvas
        atom_ids = [
            canvas.model.add_atom("C", 42.0 * index, 0.0).id
            for index in range(initial_atom_count)
        ]
        for first_id, second_id in zip(atom_ids, atom_ids[1:]):
            canvas.model.add_bond(first_id, second_id, order=1)
        original = CoordinateSet3D(
            {atom_id: (1.5 * index, 0.0, 0.0) for index, atom_id in enumerate(atom_ids)},
            source="test",
            method="manual",
        )
        canvas.set_geometry3d_coordinates(original)
        window._pending_3d_build_anchor_id = atom_ids[-1]

        window._on_viewport_build_chain_extension_requested(
            ("C", (1.5 * initial_atom_count, 0.0, 0.0))
        )

        assert len(canvas.model.atoms) == initial_atom_count + 1, name
        assert len(canvas.model.bonds) == initial_atom_count, name
        assert all(bond.order == 1 and not bond.is_aromatic for bond in canvas.model.bonds.values())
        assert canvas.undo_stack.count() == 1
        canvas.undo_stack.undo()
        assert len(canvas.model.atoms) == initial_atom_count
        assert len(canvas.model.bonds) == initial_atom_count - 1
        assert canvas.coordinate_set_3d == original
    finally:
        window._close_approved = True
        window._async_shutdown_complete = True
        window.close()


def test_3d_builder_double_click_without_anchor_is_a_noop() -> None:
    window = ChemusonWindow()
    try:
        canvas = window.canvas
        before = (canvas.undo_stack.count(), canvas.undo_stack.index())

        window._on_viewport_build_chain_extension_requested(("C", (1.5, 0.0, 0.0)))

        assert not canvas.model.atoms
        assert not canvas.model.bonds
        assert (canvas.undo_stack.count(), canvas.undo_stack.index()) == before
        assert "ancla" in window.viewport3d_dock.status_label.text().lower()
    finally:
        window._close_approved = True
        window._async_shutdown_complete = True
        window.close()


def test_build_tool_element_replacement_is_one_undoable_atom_change() -> None:
    window = ChemusonWindow()
    try:
        canvas = window.canvas
        first = canvas.model.add_atom("C", 0.0, 0.0)
        second = canvas.model.add_atom("C", 42.0, 0.0)
        canvas.model.add_bond(first.id, second.id, order=1)
        original = CoordinateSet3D(
            {first.id: (0.0, 0.0, 0.0), second.id: (1.5, 0.0, 0.0)},
            source="test",
            method="manual",
        )
        canvas.set_geometry3d_coordinates(original)

        window._on_viewport_atom_element_change(first.id, "N")

        assert canvas.model.get_atom(first.id).element == "N"
        assert canvas.undo_stack.count() == 1
        canvas.undo_stack.undo()
        assert canvas.model.get_atom(first.id).element == "C"
        assert canvas.coordinate_set_3d == original
    finally:
        window._close_approved = True
        window._async_shutdown_complete = True
        window.close()


def test_3d_builder_links_existing_atoms_as_one_undoable_edit() -> None:
    window = ChemusonWindow()
    try:
        canvas = window.canvas
        first = canvas.model.add_atom("C", 0.0, 0.0)
        second = canvas.model.add_atom("O", 42.0, 0.0)
        canvas.set_geometry3d_coordinates(
            CoordinateSet3D(
                {first.id: (0.0, 0.0, 0.0), second.id: (1.5, 0.0, 0.0)},
                source="test",
                method="manual",
            )
        )
        window.viewport3d_dock.bond_order_combo.setCurrentIndex(0)

        window._on_viewport_build_bond_requested((first.id, second.id))

        assert len(canvas.model.bonds) == 1
        assert canvas.coordinate_set_3d is not None
        assert canvas.geometry3d_is_stale
        assert canvas.undo_stack.count() == 1
        canvas.undo_stack.undo()
        assert not canvas.model.bonds
        assert not canvas.geometry3d_is_stale
    finally:
        window._close_approved = True
        window._async_shutdown_complete = True
        window.close()


def test_3d_builder_updates_existing_bond_order_as_one_undoable_edit() -> None:
    window = ChemusonWindow()
    try:
        canvas = window.canvas
        first = canvas.model.add_atom("C", 0.0, 0.0)
        second = canvas.model.add_atom("C", 42.0, 0.0)
        bond = canvas.model.add_bond(first.id, second.id, order=1)
        canvas.set_geometry3d_coordinates(
            CoordinateSet3D(
                {first.id: (0.0, 0.0, 0.0), second.id: (1.5, 0.0, 0.0)},
                source="test",
                method="manual",
            )
        )
        window.viewport3d_dock.bond_order_combo.setCurrentIndex(1)

        window._on_viewport_build_bond_requested((first.id, second.id))

        assert canvas.model.bonds[bond.id].order == 2
        assert canvas.geometry3d_is_stale
        assert canvas.undo_stack.count() == 1
        canvas.undo_stack.undo()
        assert canvas.model.bonds[bond.id].order == 1
        assert not canvas.geometry3d_is_stale
    finally:
        window._close_approved = True
        window._async_shutdown_complete = True
        window.close()


def test_3d_viewport_representation_control_updates_mesh_style() -> None:
    window = ChemusonWindow()
    try:
        combo = window.viewport3d_dock.representation_combo
        index = combo.findData("space_filling")
        assert index >= 0

        combo.setCurrentIndex(index)

        assert window.viewport3d_dock.viewport.representation == "space_filling"
    finally:
        window._close_approved = True
        window._async_shutdown_complete = True
        window.close()


def test_3d_builder_respects_selected_double_bond_order() -> None:
    window = ChemusonWindow()
    try:
        canvas = window.canvas
        parent = canvas.model.add_atom("C", 0.0, 0.0)
        canvas.set_geometry3d_coordinates(
            CoordinateSet3D(
                {parent.id: (0.0, 0.0, 0.0)}, source="test", method="manual"
            )
        )
        window._pending_3d_build_anchor_id = parent.id
        window.viewport3d_dock.bond_order_combo.setCurrentIndex(1)

        window._add_3d_atom("C", (2.0, 0.0, 0.0))

        bond = next(iter(canvas.model.bonds.values()))
        assert bond.order == 2
        positions = canvas.coordinate_set_3d.normalized_positions()
        added_id = next(atom_id for atom_id in canvas.model.atoms if atom_id != parent.id)
        assert math.dist(positions[parent.id], positions[added_id]) == pytest.approx(1.34)
    finally:
        window._close_approved = True
        window._async_shutdown_complete = True
        window.close()


@pytest.mark.parametrize(
    ("name", "elements"),
    (
        ("ethane", ("C", "C")),
        ("ethanol", ("C", "C", "O")),
        ("butane", ("C", "C", "C", "C")),
    ),
)
def test_simple_molecules_can_be_selected_and_moved_without_manual_atom_ids(
    monkeypatch, name, elements
) -> None:
    from chemuson.geometry3d.model import scene_molecule_from_graph
    from chemuson.gui.viewport3d.picking import PickResult

    window = ChemusonWindow()
    try:
        graph = MolGraph()
        atoms = [
            graph.add_atom(element, index * 42.0, 0.0)
            for index, element in enumerate(elements)
        ]
        for first, second in zip(atoms, atoms[1:]):
            graph.add_bond(first.id, second.id)
        window.canvas.model = graph
        window.canvas._rebuild_items_from_model()
        original_positions = {
            atom.id: (index * 1.5, (index % 2) * 0.2, 0.0)
            for index, atom in enumerate(atoms)
        }
        last = atoms[-1]
        original = CoordinateSet3D(
            original_positions,
            source="test",
            method="manual",
            metadata={
                "generated_atoms": [
                    {
                        "id": f"genH-{name}",
                        "element": "H",
                        "coords": [original_positions[last.id][0], 1.0, 0.0],
                    }
                ],
                "generated_bonds": [
                    {
                        "id": f"genB-{name}",
                        "a1_id": last.id,
                        "a2_id": f"genH-{name}",
                        "order": 1,
                    }
                ],
            },
        )
        window.canvas.set_geometry3d_coordinates(original)
        viewport = window.viewport3d_dock.viewport
        viewport.resize(1000, 800)
        viewport.set_scene(scene_molecule_from_graph(graph, original))
        viewport.set_edit_mode(True)
        viewport.set_edit_tool("move")
        monkeypatch.setattr(
            "chemuson.gui.viewport3d.widget.pick_scene",
            lambda *_args: PickResult("atom", atoms[0].id, 1.0),
        )
        start = QPointF(400.0, 300.0)
        end = QPointF(460.0, 330.0)
        viewport.mousePressEvent(
            QMouseEvent(
                QEvent.Type.MouseButtonPress,
                start,
                Qt.MouseButton.LeftButton,
                Qt.MouseButton.LeftButton,
                Qt.KeyboardModifier.ControlModifier,
            )
        )
        assert set(viewport.selected_atom_ids) == {atom.id for atom in atoms}
        viewport.mouseMoveEvent(
            QMouseEvent(
                QEvent.Type.MouseMove,
                end,
                Qt.MouseButton.NoButton,
                Qt.MouseButton.LeftButton,
                Qt.KeyboardModifier.ControlModifier,
            )
        )
        viewport.mouseReleaseEvent(
            QMouseEvent(
                QEvent.Type.MouseButtonRelease,
                end,
                Qt.MouseButton.LeftButton,
                Qt.MouseButton.NoButton,
                Qt.KeyboardModifier.ControlModifier,
            )
        )

        updated = window.canvas.coordinate_set_3d
        assert updated is not None, name
        deltas = [
            tuple(
                updated.positions[atom.id][axis] - original_positions[atom.id][axis]
                for axis in range(3)
            )
            for atom in atoms
        ]
        assert all(delta == pytest.approx(deltas[0]) for delta in deltas)
        generated_atoms = updated.metadata.get("generated_atoms")
        assert isinstance(generated_atoms, list) and generated_atoms
        generated = generated_atoms[0]
        assert isinstance(generated, dict)
        generated_delta = tuple(
            generated["coords"][axis] - original.metadata["generated_atoms"][0]["coords"][axis]
            for axis in range(3)
        )
        assert generated_delta == pytest.approx(deltas[0])
        assert window.canvas.undo_stack.count() == 1
        window.canvas.undo_stack.undo()
        assert window.canvas.coordinate_set_3d == original
        window.canvas.undo_stack.redo()
        assert window.canvas.coordinate_set_3d == updated
        window.canvas.undo_stack.setClean()
    finally:
        window._close_approved = True
        window._async_shutdown_complete = True
        window.close()


def test_3d_geometry_preview_keeps_and_moves_derived_hydrogens() -> None:
    from chemuson.geometry3d.model import scene_molecule_from_graph

    window = ChemusonWindow()
    try:
        canvas = window.canvas
        carbon = canvas.model.add_atom("C", 0.0, 0.0, is_explicit=True)
        initial = CoordinateSet3D(
            {carbon.id: (0.0, 0.0, 0.0)},
            source="test",
            method="mmff",
            metadata={
                "generated_atoms": [
                    {"id": "genH1", "element": "H", "coords": [1.0, 0.0, 0.0]}
                ],
                "generated_bonds": [
                    {"id": "genB1", "a1_id": carbon.id, "a2_id": "genH1", "order": 1}
                ],
            },
        )
        canvas.set_geometry3d_coordinates(initial)

        window._on_viewport_geometry_preview(
            initial.positions,
            {carbon.id: (0.0, 1.0, 0.0)},
            (carbon.id,),
        )

        preview = window._compchem_preview_coordset
        assert preview is not None
        scene = scene_molecule_from_graph(canvas.model, preview)
        hydrogen = next(atom for atom in scene.atoms if atom.generated)
        assert hydrogen.parent_atom_id == carbon.id
        assert len(scene.atoms) == 2
        assert len(scene.bonds) == 1
        assert (hydrogen.x, hydrogen.y, hydrogen.z) == (1.0, 1.0, 0.0)
    finally:
        window._close_approved = True
        window._async_shutdown_complete = True
        window.close()


def test_dragging_derived_hydrogen_promotes_it_as_one_undoable_edit(monkeypatch) -> None:
    from chemuson.geometry3d.model import scene_molecule_from_graph

    window = ChemusonWindow()
    try:
        canvas = window.canvas
        carbon = canvas.model.add_atom("C", 0.0, 0.0, is_explicit=True)
        original = CoordinateSet3D(
            {carbon.id: (0.0, 0.0, 0.0)},
            source="test",
            method="mmff",
            metadata={
                "generated_atoms": [
                    {"id": "genH1", "element": "H", "coords": [1.0, 0.0, 0.0]}
                ],
                "generated_bonds": [
                    {"id": "genB1", "a1_id": carbon.id, "a2_id": "genH1", "order": 1}
                ],
            },
        )
        canvas.set_geometry3d_coordinates(original)
        window.viewport3d_dock.live_mm_checkbox.setChecked(True)
        window._live_mm_timer.stop()

        window._on_viewport_derived_hydrogen_edit_requested(
            ("genH1", carbon.id, (0.0, 1.0, 0.0))
        )

        assert len(canvas.model.atoms) == 2
        hydrogen = next(atom for atom in canvas.model.atoms.values() if atom.element == "H")
        assert hydrogen.is_explicit
        assert math.hypot(hydrogen.x - carbon.x, hydrogen.y - carbon.y) == pytest.approx(
            canvas.state.bond_length
        )
        bond = next(iter(canvas.model.bonds.values()))
        assert {bond.a1_id, bond.a2_id} == {carbon.id, hydrogen.id}
        coordinates = canvas.coordinate_set_3d
        assert coordinates is not None
        assert coordinates.normalized_positions()[hydrogen.id] == (0.0, 1.0, 0.0)
        assert coordinates.metadata.get("generated_atoms", []) == []
        assert coordinates.metadata.get("generated_bonds", []) == []
        assert not canvas.geometry3d_is_stale
        assert canvas.undo_stack.count() == 1

        canvas.undo_stack.undo()
        assert list(canvas.model.atoms) == [carbon.id]
        assert canvas.coordinate_set_3d == original
        canvas.undo_stack.redo()
        assert len(canvas.model.atoms) == 2
        assert canvas.coordinate_set_3d == coordinates
        assert window._live_mm_fixed_atom_ids == (hydrogen.id,)
        captured = {}
        controller = window._compchem3d_controller
        monkeypatch.setattr(controller, "active_jobs", lambda: [])
        monkeypatch.setattr(
            controller,
            "start_job",
            lambda graph, spec, coordset: captured.update(
                graph=graph, spec=spec, coordset=coordset
            )
            or 9001,
        )
        window._start_live_mm_job()
        window._live_mm_timer.stop()
        assert captured["spec"].fixed_atom_ids == (hydrogen.id,)
        assert captured["coordset"].normalized_positions()[hydrogen.id] == (
            0.0,
            1.0,
            0.0,
        )
        from chemuson.chemio.persistence import PersistenceManager
        from chemuson.gui.canvas import ChemusonCanvas

        payload = PersistenceManager.save_to_dict(canvas)
        restored = ChemusonCanvas()
        PersistenceManager.load_from_dict(payload, restored)
        assert len(restored.model.atoms) == 2
        assert restored.coordinate_set_3d is not None
        restored_hydrogen = next(
            atom for atom in restored.model.atoms.values() if atom.element == "H"
        )
        assert restored.coordinate_set_3d.normalized_positions()[restored_hydrogen.id] == (
            0.0,
            1.0,
            0.0,
        )

        scene = scene_molecule_from_graph(canvas.model, coordinates)
        assert sum(atom.element == "H" for atom in scene.atoms) == 1
    finally:
        window._close_approved = True
        window._async_shutdown_complete = True
        window.close()


def test_promoted_explicit_hydrogen_is_not_duplicated_by_rdkit() -> None:
    from chemuson.chemio.rdkit_io import _rdkit_available, molgraph_to_rdkit
    from chemuson.gui.canvas import ChemusonCanvas

    if not _rdkit_available():
        pytest.skip("RDKit no disponible en este entorno")
    canvas = ChemusonCanvas()
    carbon = canvas.model.add_atom("C", 0.0, 0.0, is_explicit=True)
    hydrogen = canvas.model.add_atom("H", 40.0, 0.0, is_explicit=True)
    canvas.model.add_bond(carbon.id, hydrogen.id, order=1)

    rdkit_mol = molgraph_to_rdkit(canvas.model)

    assert rdkit_mol.GetNumAtoms() == 2


def test_cancelled_3d_drag_restores_derived_hydrogen_coordinates() -> None:
    from chemuson.geometry3d.model import scene_molecule_from_graph

    window = ChemusonWindow()
    try:
        canvas = window.canvas
        carbon = canvas.model.add_atom("C", 0.0, 0.0, is_explicit=True)
        original = CoordinateSet3D(
            {carbon.id: (0.0, 0.0, 0.0)},
            source="test",
            method="mmff",
            metadata={
                "generated_atoms": [
                    {"id": "genH1", "element": "H", "coords": [1.0, 0.0, 0.0]}
                ],
                "generated_bonds": [
                    {"id": "genB1", "a1_id": carbon.id, "a2_id": "genH1", "order": 1}
                ],
            },
        )
        canvas.set_geometry3d_coordinates(original)
        window._on_viewport_geometry_preview(
            original.positions,
            {carbon.id: (0.0, 1.0, 0.0)},
            (carbon.id,),
        )
        preview = window._compchem_preview_coordset
        assert preview is not None
        preview_scene = scene_molecule_from_graph(canvas.model, preview)
        preview_hydrogen = next(atom for atom in preview_scene.atoms if atom.generated)
        assert (preview_hydrogen.x, preview_hydrogen.y, preview_hydrogen.z) == (
            1.0,
            1.0,
            0.0,
        )

        window._on_viewport_geometry_cancelled()

        assert window._compchem_preview_coordset is None
        assert canvas.coordinate_set_3d == original
        restored = canvas.coordinate_set_3d
        assert restored is not None
        restored_scene = scene_molecule_from_graph(canvas.model, restored)
        restored_hydrogen = next(atom for atom in restored_scene.atoms if atom.generated)
        assert (restored_hydrogen.x, restored_hydrogen.y, restored_hydrogen.z) == (
            1.0,
            0.0,
            0.0,
        )
    finally:
        window._close_approved = True
        window._async_shutdown_complete = True
        window.close()


def test_3d_builder_rejects_incompatible_valence_without_partial_edit() -> None:
    window = ChemusonWindow()
    try:
        canvas = window.canvas
        center = canvas.model.add_atom("C", 0.0, 0.0, is_explicit=True)
        hydrogens = [
            canvas.model.add_atom(
                "H", float(index + 1), 0.0, is_explicit=True
            )
            for index in range(4)
        ]
        positions = {center.id: (0.0, 0.0, 0.0)}
        positions.update(
            {atom.id: (float(index + 1), 0.0, 0.0) for index, atom in enumerate(hydrogens)}
        )
        for hydrogen in hydrogens:
            canvas.model.add_bond(center.id, hydrogen.id)
        original_coordinates = CoordinateSet3D(
            positions, source="test", method="manual"
        )
        canvas.set_geometry3d_coordinates(original_coordinates)
        before_atoms = {
            atom_id: (atom.element, atom.implicit_h, atom.has_valence_error)
            for atom_id, atom in canvas.model.atoms.items()
        }
        before_bonds = {
            bond_id: (bond.a1_id, bond.a2_id, bond.order)
            for bond_id, bond in canvas.model.bonds.items()
        }
        before_index = canvas.undo_stack.index()
        before_count = canvas.undo_stack.count()
        window._pending_3d_build_anchor_id = center.id

        window._add_3d_atom("C", (2.0, 0.0, 0.0))

        assert set(canvas.model.atoms) == set(before_atoms)
        assert {
            atom_id: (atom.element, atom.implicit_h, atom.has_valence_error)
            for atom_id, atom in canvas.model.atoms.items()
        } == before_atoms
        assert {
            bond_id: (bond.a1_id, bond.a2_id, bond.order)
            for bond_id, bond in canvas.model.bonds.items()
        } == before_bonds
        assert canvas.coordinate_set_3d == original_coordinates
        assert canvas.undo_stack.index() == before_index
        assert canvas.undo_stack.count() == before_count
        assert "valencia" in window.viewport3d_dock.status_label.text().lower()
    finally:
        window._close_approved = True
        window._async_shutdown_complete = True
        window.close()


def test_3d_bond_command_rejects_missing_endpoint_atomically() -> None:
    from chemuson.gui.commands import AddBondCommand

    window = ChemusonWindow()
    try:
        canvas = window.canvas
        atom = canvas.model.add_atom("C", 0.0, 0.0)
        coordinates = CoordinateSet3D(
            {atom.id: (0.0, 0.0, 0.0)}, source="test", method="manual"
        )
        canvas.set_geometry3d_coordinates(coordinates)
        before_atom_ids = set(canvas.model.atoms)
        before_bond_ids = set(canvas.model.bonds)
        before_history = (canvas.undo_stack.count(), canvas.undo_stack.index())
        command = AddBondCommand(canvas.model, canvas, atom.id, atom.id + 100)

        error = None
        try:
            command.redo()
        except Exception as exc:
            error = exc

        assert set(canvas.model.atoms) == before_atom_ids
        assert set(canvas.model.bonds) == before_bond_ids
        assert canvas.coordinate_set_3d == coordinates
        assert (canvas.undo_stack.count(), canvas.undo_stack.index()) == before_history
        assert isinstance(error, ValueError)
        assert "átomo" in str(error).lower()
    finally:
        window._close_approved = True
        window._async_shutdown_complete = True
        window.close()


def test_3d_bond_command_rejects_duplicate_endpoint_pair_atomically() -> None:
    from chemuson.gui.commands import AddBondCommand

    window = ChemusonWindow()
    try:
        canvas = window.canvas
        first = canvas.model.add_atom("C", 0.0, 0.0)
        second = canvas.model.add_atom("C", 42.0, 0.0)
        existing = canvas.model.add_bond(first.id, second.id, order=1)
        canvas.add_bond_item(existing)
        coordinates = CoordinateSet3D(
            {first.id: (0.0, 0.0, 0.0), second.id: (1.5, 0.0, 0.0)},
            source="test",
            method="manual",
        )
        canvas.set_geometry3d_coordinates(coordinates)
        before_bonds = {
            bond_id: (bond.a1_id, bond.a2_id, bond.order)
            for bond_id, bond in canvas.model.bonds.items()
        }
        before_history = (canvas.undo_stack.count(), canvas.undo_stack.index())
        command = AddBondCommand(canvas.model, canvas, first.id, second.id)

        error = None
        try:
            command.redo()
        except Exception as exc:
            error = exc

        assert {
            bond_id: (bond.a1_id, bond.a2_id, bond.order)
            for bond_id, bond in canvas.model.bonds.items()
        } == before_bonds
        assert canvas.coordinate_set_3d == coordinates
        assert (canvas.undo_stack.count(), canvas.undo_stack.index()) == before_history
        assert isinstance(error, ValueError)
        assert "duplicado" in str(error).lower()
    finally:
        window._close_approved = True
        window._async_shutdown_complete = True
        window.close()


@pytest.mark.parametrize(
    ("operation", "atom_ids", "target", "metric"),
    [
        ("bond_length", (1, 2), 1.8, "distance"),
        ("bond_angle", (1, 2, 3), 100.0, "angle_degrees"),
        ("dihedral", (1, 2, 3, 4), 60.0, "dihedral_degrees"),
    ],
)
def test_geometry_constraint_controls_are_applied_and_undoable(
    operation, atom_ids, target, metric
) -> None:
    from chemuson.geometry3d import editing

    window = ChemusonWindow()
    try:
        graph = MolGraph()
        atoms = [graph.add_atom("C", index * 42.0, 0.0) for index in range(4)]
        for first, second in zip(atoms, atoms[1:]):
            graph.add_bond(first.id, second.id, order=1)
        window.canvas.model = graph
        window.canvas._rebuild_items_from_model()
        original_positions = {
            1: (0.0, 0.0, 0.0),
            2: (1.5, 0.0, 0.0),
            3: (2.1, 1.2, 0.0),
            4: (2.8, 1.4, 1.0),
        }
        generated_parent_id = {
            "bond_length": 2,
            "bond_angle": 3,
            "dihedral": 4,
        }[operation]
        parent_position = original_positions[generated_parent_id]
        original = CoordinateSet3D(
            original_positions,
            source="test",
            method="manual",
            metadata={
                "generated_atoms": [
                    {
                        "id": "genH",
                        "element": "H",
                        "coords": [parent_position[0], parent_position[1] + 1.0, parent_position[2]],
                    }
                ],
                "generated_bonds": [
                    {
                        "id": "genB",
                        "a1_id": generated_parent_id,
                        "a2_id": "genH",
                        "order": 1,
                    }
                ],
            },
        )
        if operation == "bond_length":
            expected = editing.set_bond_length(graph, original, 1, target)
        elif operation == "bond_angle":
            expected = editing.set_bond_angle(graph, original, atom_ids, target)
        else:
            expected = editing.set_dihedral(graph, original, atom_ids, target)
        window.canvas.set_geometry3d_coordinates(original)

        window._on_geometry_constraint_requested(operation, atom_ids, target)

        updated = window.canvas.coordinate_set_3d
        assert updated is not None
        from chemuson.geometry3d.model import scene_molecule_from_graph

        assert len(scene_molecule_from_graph(graph, updated).atoms) == 5
        measured = getattr(editing, metric)(updated, *atom_ids)
        assert measured == pytest.approx(target, abs=1.0e-4)
        actual_generated_atoms = updated.metadata.get("generated_atoms")
        expected_generated_atoms = expected.metadata.get("generated_atoms")
        assert isinstance(actual_generated_atoms, list) and actual_generated_atoms
        assert isinstance(expected_generated_atoms, list) and expected_generated_atoms
        actual_hydrogen = actual_generated_atoms[0]
        expected_hydrogen = expected_generated_atoms[0]
        assert isinstance(actual_hydrogen, dict)
        assert isinstance(expected_hydrogen, dict)
        assert actual_hydrogen["coords"] == pytest.approx(expected_hydrogen["coords"])
        assert window.canvas.undo_stack.count() == 1
        window.canvas.undo_stack.undo()
        assert window.canvas.coordinate_set_3d == original
    finally:
        window._close_approved = True
        window._async_shutdown_complete = True
        window.close()


def test_dihedral_edit_preserves_rigid_generated_hydrogen_transform() -> None:
    from chemuson.geometry3d import editing

    window = ChemusonWindow()
    try:
        graph = MolGraph()
        atoms = [graph.add_atom("C", index * 42.0, 0.0) for index in range(4)]
        for first, second in zip(atoms, atoms[1:]):
            graph.add_bond(first.id, second.id, order=1)
        window.canvas.model = graph
        window.canvas._rebuild_items_from_model()
        original = CoordinateSet3D(
            {
                1: (0.0, 1.0, 0.0),
                2: (1.0, 0.0, 0.0),
                3: (2.0, 0.0, 0.0),
                4: (3.0, 0.0, 1.0),
            },
            source="test",
            method="manual",
            metadata={
                "generated_atoms": [
                    {"id": "genH4", "element": "H", "coords": [3.0, 1.0, 1.0]}
                ],
                "generated_bonds": [
                    {"id": "genB4", "a1_id": 4, "a2_id": "genH4", "order": 1}
                ],
            },
        )
        expected = editing.set_dihedral(graph, original, (1, 2, 3, 4), 60.0)
        window.canvas.set_geometry3d_coordinates(original)

        window._on_geometry_constraint_requested("dihedral", (1, 2, 3, 4), 60.0)

        committed = window.canvas.coordinate_set_3d
        assert committed is not None
        assert committed.positions == expected.positions
        assert committed.metadata["generated_atoms"][0]["coords"] == pytest.approx(
            expected.metadata["generated_atoms"][0]["coords"]
        )
        assert committed.metadata["generated_atoms"][0]["coords"] != pytest.approx(
            original.metadata["generated_atoms"][0]["coords"]
        )
        assert window.canvas.undo_stack.count() == 1

        window.canvas.undo_stack.undo()
        assert window.canvas.coordinate_set_3d == original
        window.canvas.undo_stack.redo()
        assert window.canvas.coordinate_set_3d == committed
    finally:
        window._close_approved = True
        window._async_shutdown_complete = True
        window.close()


@pytest.mark.parametrize(
    ("operation", "atom_ids", "delta_x", "delta_y", "target"),
    (
        ("bond_length", (1, 2), 0.0, -50.0, math.sqrt(2.0) + 0.5),
        ("bond_angle", (1, 2, 3), 0.0, -20.0, 145.0),
        ("dihedral", (1, 2, 3, 4), 40.0, 0.0, 110.0),
    ),
)
def test_geometry_handles_adjust_length_angle_and_dihedral(
    monkeypatch, operation, atom_ids, delta_x, delta_y, target
) -> None:
    from chemuson.geometry3d import editing
    from chemuson.geometry3d.model import scene_molecule_from_graph
    from chemuson.gui.viewport3d.picking import PickResult

    window = ChemusonWindow()
    try:
        graph = MolGraph()
        for index in range(4):
            graph.add_atom("C", index * 42.0, 0.0, atom_id=index + 1)
        for atom_id in range(1, 4):
            graph.add_bond(atom_id, atom_id + 1)
        window.canvas.model = graph
        window.canvas._rebuild_items_from_model()
        original = CoordinateSet3D(
            {
                1: (0.0, 1.0, 0.0),
                2: (1.0, 0.0, 0.0),
                3: (2.0, 0.0, 0.0),
                4: (3.0, 0.0, 1.0),
            },
            source="test",
            method="manual",
            metadata={
                "generated_atoms": [
                    {"id": "genH4", "element": "H", "coords": [3.0, 1.0, 1.0]}
                ],
                "generated_bonds": [
                    {"id": "genB4", "a1_id": 4, "a2_id": "genH4", "order": 1}
                ],
            },
        )
        if operation == "bond_length":
            bond = graph.find_bond_between(*atom_ids)
            assert bond is not None
            expected = editing.set_bond_length(
                graph, original, bond.id, target, moving_endpoint_id=atom_ids[-1]
            )
        elif operation == "bond_angle":
            expected = editing.set_bond_angle(
                graph, original, atom_ids, target, moving_terminal_id=atom_ids[-1]
            )
        else:
            expected = editing.set_dihedral(
                graph, original, atom_ids, target, moving_endpoint_id=atom_ids[2]
            )
        window.canvas.set_geometry3d_coordinates(original)
        viewport = window.viewport3d_dock.viewport
        viewport.resize(1000, 800)
        viewport.set_scene(scene_molecule_from_graph(graph, original))
        viewport.set_edit_mode(True)
        viewport.set_edit_tool("geometry")
        viewport.set_geometry_operation(operation)
        hit_ids = iter(atom_ids)
        monkeypatch.setattr(
            "chemuson.gui.viewport3d.widget.pick_scene",
            lambda *_args: PickResult("atom", next(hit_ids), 1.0),
        )

        for index, atom_id in enumerate(atom_ids):
            point = QPointF(100.0 + index * 30.0, 100.0)
            modifiers = (
                Qt.KeyboardModifier.NoModifier
                if index == 0
                else Qt.KeyboardModifier.ShiftModifier
            )
            viewport.mousePressEvent(
                QMouseEvent(
                    QEvent.Type.MouseButtonPress,
                    point,
                    Qt.MouseButton.LeftButton,
                    Qt.MouseButton.LeftButton,
                    modifiers,
                )
            )
            viewport.mouseReleaseEvent(
                QMouseEvent(
                    QEvent.Type.MouseButtonRelease,
                    point,
                    Qt.MouseButton.LeftButton,
                    Qt.MouseButton.NoButton,
                    modifiers,
                )
            )

        handle = viewport._geometry_handle_screen_position()
        assert handle is not None
        start = QPointF(handle.x(), handle.y())
        end = QPointF(handle.x() + delta_x, handle.y() + delta_y)
        viewport.mousePressEvent(
            QMouseEvent(
                QEvent.Type.MouseButtonPress,
                start,
                Qt.MouseButton.LeftButton,
                Qt.MouseButton.LeftButton,
                Qt.KeyboardModifier.NoModifier,
            )
        )
        viewport.mouseMoveEvent(
            QMouseEvent(
                QEvent.Type.MouseMove,
                end,
                Qt.MouseButton.NoButton,
                Qt.MouseButton.LeftButton,
                Qt.KeyboardModifier.NoModifier,
            )
        )
        expected_unit = "Å" if operation == "bond_length" else "°"
        assert expected_unit in window.viewport3d_dock.status_label.text()
        viewport.mouseReleaseEvent(
            QMouseEvent(
                QEvent.Type.MouseButtonRelease,
                end,
                Qt.MouseButton.LeftButton,
                Qt.MouseButton.NoButton,
                Qt.KeyboardModifier.NoModifier,
            )
        )

        updated = window.canvas.coordinate_set_3d
        assert updated is not None
        measure = {
            "bond_length": editing.distance,
            "bond_angle": editing.angle_degrees,
            "dihedral": editing.dihedral_degrees,
        }[operation]
        assert measure(updated, *atom_ids) == pytest.approx(target, abs=1.0e-4)
        actual_hydrogens = updated.metadata.get("generated_atoms")
        expected_hydrogens = expected.metadata.get("generated_atoms")
        assert isinstance(actual_hydrogens, list) and actual_hydrogens
        assert isinstance(expected_hydrogens, list) and expected_hydrogens
        actual_hydrogen = actual_hydrogens[0]
        expected_hydrogen = expected_hydrogens[0]
        assert isinstance(actual_hydrogen, dict)
        assert isinstance(expected_hydrogen, dict)
        assert actual_hydrogen["coords"] == pytest.approx(expected_hydrogen["coords"])
        assert window.canvas.undo_stack.count() == 1
        window.canvas.undo_stack.setClean()
    finally:
        window._close_approved = True
        window._async_shutdown_complete = True
        window.close()


def test_geometry_handle_rejects_degenerate_bond_without_mutating_coordinates(
    monkeypatch,
) -> None:
    from chemuson.geometry3d.model import scene_molecule_from_graph
    from chemuson.gui.viewport3d.picking import PickResult

    window = ChemusonWindow()
    try:
        graph = MolGraph()
        graph.add_atom("C", 0.0, 0.0, atom_id=1)
        graph.add_atom("C", 42.0, 0.0, atom_id=2)
        graph.add_bond(1, 2)
        window.canvas.model = graph
        window.canvas._rebuild_items_from_model()
        original = CoordinateSet3D(
            {1: (0.0, 0.0, 0.0), 2: (0.0, 0.0, 0.0)},
            source="test",
            method="manual",
        )
        window.canvas.set_geometry3d_coordinates(original)
        viewport = window.viewport3d_dock.viewport
        viewport.resize(1000, 800)
        viewport.set_scene(scene_molecule_from_graph(graph, original))
        viewport.set_edit_mode(True)
        viewport.set_edit_tool("geometry")
        viewport.set_geometry_operation("bond_length")
        hits = iter(
            (
                PickResult("atom", 1, 1.0),
                PickResult("atom", 2, 1.0),
                None,
            )
        )
        monkeypatch.setattr(
            "chemuson.gui.viewport3d.widget.pick_scene",
            lambda *_args: next(hits),
        )
        for index, modifiers in (
            (0, Qt.KeyboardModifier.NoModifier),
            (1, Qt.KeyboardModifier.ShiftModifier),
        ):
            point = QPointF(100.0 + index * 30.0, 100.0)
            viewport.mousePressEvent(
                QMouseEvent(
                    QEvent.Type.MouseButtonPress,
                    point,
                    Qt.MouseButton.LeftButton,
                    Qt.MouseButton.LeftButton,
                    modifiers,
                )
            )
            viewport.mouseReleaseEvent(
                QMouseEvent(
                    QEvent.Type.MouseButtonRelease,
                    point,
                    Qt.MouseButton.LeftButton,
                    Qt.MouseButton.NoButton,
                    modifiers,
                )
            )
        handle = viewport._geometry_handle_screen_position()
        assert handle is not None
        start = QPointF(handle.x(), handle.y())
        end = QPointF(handle.x(), handle.y() - 20.0)
        history = (window.canvas.undo_stack.count(), window.canvas.undo_stack.index())

        viewport.mousePressEvent(
            QMouseEvent(
                QEvent.Type.MouseButtonPress,
                start,
                Qt.MouseButton.LeftButton,
                Qt.MouseButton.LeftButton,
                Qt.KeyboardModifier.NoModifier,
            )
        )
        viewport.mouseMoveEvent(
            QMouseEvent(
                QEvent.Type.MouseMove,
                end,
                Qt.MouseButton.NoButton,
                Qt.MouseButton.LeftButton,
                Qt.KeyboardModifier.NoModifier,
            )
        )
        viewport.mouseReleaseEvent(
            QMouseEvent(
                QEvent.Type.MouseButtonRelease,
                end,
                Qt.MouseButton.LeftButton,
                Qt.MouseButton.NoButton,
                Qt.KeyboardModifier.NoModifier,
            )
        )

        assert not viewport.geometry_drag_active
        assert window.canvas.coordinate_set_3d == original
        assert (window.canvas.undo_stack.count(), window.canvas.undo_stack.index()) == history
        assert "geometría degenerada" in window.viewport3d_dock.status_label.text().lower()
    finally:
        window.canvas.undo_stack.setClean()
        window._close_approved = True
        window._async_shutdown_complete = True
        window.close()


def test_geometry_handle_rejects_cycle_without_mutating_coordinates(monkeypatch) -> None:
    from chemuson.geometry3d.model import scene_molecule_from_graph

    window = ChemusonWindow()
    try:
        graph = MolGraph()
        for atom_id, point in enumerate(((0.0, 0.0, 0.0), (1.5, 0.0, 0.0), (0.75, 1.3, 0.0)), 1):
            graph.add_atom("C", point[0] * 42.0, point[1] * 42.0, atom_id=atom_id)
        graph.add_bond(1, 2)
        graph.add_bond(2, 3)
        graph.add_bond(3, 1)
        window.canvas.model = graph
        window.canvas._rebuild_items_from_model()
        original = CoordinateSet3D(
            {1: (0.0, 0.0, 0.0), 2: (1.5, 0.0, 0.0), 3: (0.75, 1.3, 0.0)},
            source="test",
            method="manual",
        )
        window.canvas.set_geometry3d_coordinates(original)
        viewport = window.viewport3d_dock.viewport
        viewport.resize(1000, 800)
        viewport.set_scene(scene_molecule_from_graph(graph, original))
        viewport.set_edit_mode(True)
        viewport.set_edit_tool("geometry")
        viewport.set_geometry_operation("bond_length")
        from chemuson.gui.viewport3d.picking import PickResult

        hit_ids = iter((1, 2))
        monkeypatch.setattr(
            "chemuson.gui.viewport3d.widget.pick_scene",
            lambda *_args: PickResult("atom", next(hit_ids), 1.0),
        )
        for index, modifiers in (
            (0, Qt.KeyboardModifier.NoModifier),
            (1, Qt.KeyboardModifier.ShiftModifier),
        ):
            point = QPointF(100.0 + index * 30.0, 100.0)
            viewport.mousePressEvent(
                QMouseEvent(
                    QEvent.Type.MouseButtonPress,
                    point,
                    Qt.MouseButton.LeftButton,
                    Qt.MouseButton.LeftButton,
                    modifiers,
                )
            )
            viewport.mouseReleaseEvent(
                QMouseEvent(
                    QEvent.Type.MouseButtonRelease,
                    point,
                    Qt.MouseButton.LeftButton,
                    Qt.MouseButton.NoButton,
                    modifiers,
                )
            )
        handle = viewport._geometry_handle_screen_position()
        assert handle is not None
        start = QPointF(handle.x(), handle.y())
        end = QPointF(handle.x(), handle.y() - 20.0)
        viewport.mousePressEvent(
            QMouseEvent(
                QEvent.Type.MouseButtonPress,
                start,
                Qt.MouseButton.LeftButton,
                Qt.MouseButton.LeftButton,
                Qt.KeyboardModifier.NoModifier,
            )
        )
        viewport.mouseMoveEvent(
            QMouseEvent(
                QEvent.Type.MouseMove,
                end,
                Qt.MouseButton.NoButton,
                Qt.MouseButton.LeftButton,
                Qt.KeyboardModifier.NoModifier,
            )
        )
        viewport.mouseReleaseEvent(
            QMouseEvent(
                QEvent.Type.MouseButtonRelease,
                end,
                Qt.MouseButton.LeftButton,
                Qt.MouseButton.NoButton,
                Qt.KeyboardModifier.NoModifier,
            )
        )

        assert window.canvas.coordinate_set_3d == original
        assert window.canvas.undo_stack.count() == 0
        assert "ciclo" in window.viewport3d_dock.status_label.text().lower()
    finally:
        window.canvas.undo_stack.setClean()
        window._close_approved = True
        window._async_shutdown_complete = True
        window.close()


def test_viewport_geometry_handle_commits_once_and_escape_discards_preview(monkeypatch) -> None:
    from chemuson.geometry3d import editing
    from chemuson.geometry3d.model import scene_molecule_from_graph
    from chemuson.gui.viewport3d.picking import PickResult

    window = ChemusonWindow()
    try:
        graph = _chain_graph()
        window.canvas.model = graph
        window.canvas._rebuild_items_from_model()
        original = CoordinateSet3D(
            {
                1: (0.0, 0.0, 0.0),
                2: (1.5, 0.0, 0.0),
                3: (2.5, 1.0, 0.0),
            },
            source="test",
            method="manual",
        )
        window.canvas.set_geometry3d_coordinates(original)
        viewport = window.viewport3d_dock.viewport
        viewport.resize(1000, 800)
        viewport.set_scene(scene_molecule_from_graph(graph, original))
        viewport.set_edit_mode(True)
        viewport.set_edit_tool("geometry")
        viewport.set_geometry_operation("bond_length")
        hit_ids = iter((1, 2))
        monkeypatch.setattr(
            "chemuson.gui.viewport3d.widget.pick_scene",
            lambda *_args: PickResult("atom", next(hit_ids), 1.0),
        )

        for atom_index, modifiers in (
            (0, Qt.KeyboardModifier.NoModifier),
            (1, Qt.KeyboardModifier.ShiftModifier),
        ):
            point = QPointF(100.0 + atom_index * 100.0, 100.0)
            viewport.mousePressEvent(
                QMouseEvent(
                    QEvent.Type.MouseButtonPress,
                    point,
                    Qt.MouseButton.LeftButton,
                    Qt.MouseButton.LeftButton,
                    modifiers,
                )
            )
            viewport.mouseReleaseEvent(
                QMouseEvent(
                    QEvent.Type.MouseButtonRelease,
                    point,
                    Qt.MouseButton.LeftButton,
                    Qt.MouseButton.NoButton,
                    modifiers,
                )
            )

        handle = viewport._geometry_handle_screen_position()
        assert handle is not None
        start = QPointF(handle.x(), handle.y())
        end = QPointF(handle.x(), handle.y() - 50.0)
        viewport.mousePressEvent(
            QMouseEvent(
                QEvent.Type.MouseButtonPress,
                start,
                Qt.MouseButton.LeftButton,
                Qt.MouseButton.LeftButton,
                Qt.KeyboardModifier.NoModifier,
            )
        )
        viewport.mouseMoveEvent(
            QMouseEvent(
                QEvent.Type.MouseMove,
                end,
                Qt.MouseButton.NoButton,
                Qt.MouseButton.LeftButton,
                Qt.KeyboardModifier.NoModifier,
            )
        )
        viewport.mouseReleaseEvent(
            QMouseEvent(
                QEvent.Type.MouseButtonRelease,
                end,
                Qt.MouseButton.LeftButton,
                Qt.MouseButton.NoButton,
                Qt.KeyboardModifier.NoModifier,
            )
        )

        committed = window.canvas.coordinate_set_3d
        assert committed is not None
        assert editing.distance(committed, 1, 2) == pytest.approx(2.0, abs=1.0e-4)
        assert window.canvas.undo_stack.count() == 1
        window.canvas.undo_stack.undo()
        assert window.canvas.coordinate_set_3d == original
        window.canvas.undo_stack.redo()
        assert window.canvas.coordinate_set_3d == committed

        handle = viewport._geometry_handle_screen_position()
        assert handle is not None
        start = QPointF(handle.x(), handle.y())
        end = QPointF(handle.x(), handle.y() - 20.0)
        viewport.mousePressEvent(
            QMouseEvent(
                QEvent.Type.MouseButtonPress,
                start,
                Qt.MouseButton.LeftButton,
                Qt.MouseButton.LeftButton,
                Qt.KeyboardModifier.NoModifier,
            )
        )
        viewport.mouseMoveEvent(
            QMouseEvent(
                QEvent.Type.MouseMove,
                end,
                Qt.MouseButton.NoButton,
                Qt.MouseButton.LeftButton,
                Qt.KeyboardModifier.NoModifier,
            )
        )
        viewport.keyPressEvent(
            QKeyEvent(
                QEvent.Type.KeyPress,
                Qt.Key.Key_Escape,
                Qt.KeyboardModifier.NoModifier,
            )
        )

        assert window.canvas.coordinate_set_3d == committed
        assert window.canvas.undo_stack.count() == 1
        assert window._compchem_preview_coordset is None
    finally:
        window.canvas.undo_stack.setClean()
        window._close_approved = True
        window._async_shutdown_complete = True
        window.close()


def test_compchem_projection_is_undoable(monkeypatch) -> None:
    window = ChemusonWindow()
    try:
        graph = _chain_graph()
        window.canvas.model = graph
        window.canvas._rebuild_items_from_model()
        before = {atom_id: (atom.x, atom.y) for atom_id, atom in graph.atoms.items()}
        window._compchem_coordset = CoordinateSet3D(
            {
                1: (0.0, 0.0, 0.0),
                2: (1.0, 1.0, 0.0),
                3: (2.0, 0.0, 0.0),
            },
            source="fake",
            method="UFF",
        )
        monkeypatch.setattr(
            QMessageBox,
            "question",
            lambda *_args, **_kwargs: QMessageBox.StandardButton.Yes,
        )

        window._on_compchem_project_to_2d()

        after = {atom_id: (atom.x, atom.y) for atom_id, atom in graph.atoms.items()}
        assert after != before
        assert window.canvas.undo_stack.count() == 1
        window.canvas.undo_stack.undo()
        restored = {atom_id: (atom.x, atom.y) for atom_id, atom in graph.atoms.items()}
        for atom_id, coords in before.items():
            assert restored[atom_id] == pytest.approx(coords)
    finally:
        window.close()


def test_compchem_exports_xyz_and_inputs(monkeypatch, tmp_path) -> None:
    window = ChemusonWindow()
    try:
        graph = _chain_graph()
        window.canvas.model = graph
        window.canvas._rebuild_items_from_model()
        window._compchem_coordset = CoordinateSet3D(
            {1: (0.0, 0.0, 0.0), 2: (1.0, 0.0, 0.0), 3: (2.0, 0.0, 0.0)},
            source="fake",
            method="UFF",
        )
        xyz_path = tmp_path / "mol.xyz"
        monkeypatch.setattr(QFileDialog, "getSaveFileName", lambda *_args, **_kwargs: (str(xyz_path), ""))
        window._on_compchem_export_xyz()
        assert "O 2.00000000" in xyz_path.read_text(encoding="utf-8")

        inp_path = tmp_path / "mol.inp"
        monkeypatch.setattr(QFileDialog, "getSaveFileName", lambda *_args, **_kwargs: (str(inp_path), ""))
        window._on_compchem_export_input("orca")
        assert "* xyz" in inp_path.read_text(encoding="utf-8")
    finally:
        window.close()
