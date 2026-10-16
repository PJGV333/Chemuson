"""Regresiones del dock 3D / CompChem."""

from __future__ import annotations

import math

import pytest
from PyQt6.QtCore import QEventLoop, QTimer
from PyQt6.QtWidgets import QApplication, QFileDialog, QMessageBox


from chemuson.core.model import MolGraph
from chemuson.geometry3d import (
    Conformer3DResult,
    CoordinateSet3D,
    ForceField,
    OptimizationFrame,
    OptimizationResult,
    OptimizationSettings,
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
        original = CoordinateSet3D(
            {
                1: (0.0, 0.0, 0.0),
                2: (1.5, 0.0, 0.0),
                3: (2.1, 1.2, 0.0),
                4: (2.8, 1.4, 1.0),
            },
            source="test",
            method="manual",
        )
        window.canvas.set_geometry3d_coordinates(original)

        window._on_geometry_constraint_requested(operation, atom_ids, target)

        updated = window.canvas.coordinate_set_3d
        assert updated is not None
        measured = getattr(editing, metric)(updated, *atom_ids)
        assert measured == pytest.approx(target, abs=1.0e-4)
        assert window.canvas.undo_stack.count() == 1
        window.canvas.undo_stack.undo()
        assert window.canvas.coordinate_set_3d == original
    finally:
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
