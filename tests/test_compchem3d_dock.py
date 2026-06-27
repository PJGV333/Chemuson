"""Regresiones del dock 3D / CompChem."""

from __future__ import annotations


import pytest
from PyQt6.QtCore import QEventLoop, QTimer
from PyQt6.QtWidgets import QApplication, QFileDialog, QMessageBox


from chemuson.core.model import MolGraph
from chemuson.geometry3d import Conformer3DResult, CoordinateSet3D, ForceField, OptimizationSettings
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
