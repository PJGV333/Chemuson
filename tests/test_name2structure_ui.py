"""Pruebas UI del flujo profesional Name→Structure."""

from __future__ import annotations


import pytest
from PyQt6.QtWidgets import QApplication, QMessageBox


from chemuson.core.model import MolGraph
from chemuson.gui.main_window import ChemusonWindow
from chemuson.name2structure import NameToStructureResult


@pytest.fixture(scope="module", autouse=True)
def _qapp():
    return QApplication.instance() or QApplication([])


@pytest.fixture(autouse=True)
def _isolated_config_home(tmp_path, monkeypatch):
    monkeypatch.setenv("XDG_CONFIG_HOME", str(tmp_path))


def _result() -> NameToStructureResult:
    graph = MolGraph()
    graph.add_atom("O", 0.0, 0.0)
    return NameToStructureResult(
        query="agua",
        graph=graph,
        source="offline-common",
        confidence=0.72,
        smiles="O",
        resolved_name="water",
    )


def test_name_to_structure_result_confirms_and_inserts(monkeypatch) -> None:
    window = ChemusonWindow()
    try:
        monkeypatch.setattr(window, "_schedule_chemical_properties_update", lambda: None)
        inserted = {}
        monkeypatch.setattr(
            window.canvas,
            "_insert_molgraph",
            lambda graph, select_inserted=False: inserted.update(
                {"graph": graph, "select_inserted": select_inserted}
            ),
        )
        prompts: list[str] = []

        def _fake_question(_parent, _title, text, *_args):
            prompts.append(text)
            return QMessageBox.StandardButton.Yes

        monkeypatch.setattr(QMessageBox, "question", _fake_question)

        window._on_name_to_structure_result(1, _result(), "")

        assert inserted["graph"].atoms
        assert inserted["select_inserted"] is True
        assert prompts
        assert "Fuente: offline-common" in prompts[0]
        assert "Confianza: 0.72" in prompts[0]
        assert "SMILES: O" in prompts[0]
        assert window.canvas._last_name_to_structure["query"] == "agua"
        assert window.canvas._last_name_to_structure["resolved_name"] == "water"
    finally:
        window.close()


def test_name_to_structure_result_cancel_does_not_insert(monkeypatch) -> None:
    window = ChemusonWindow()
    try:
        monkeypatch.setattr(window, "_schedule_chemical_properties_update", lambda: None)
        monkeypatch.setattr(
            QMessageBox,
            "question",
            lambda *_args: QMessageBox.StandardButton.No,
        )

        window._on_name_to_structure_result(1, _result(), "")

        assert len(window.canvas.model.atoms) == 0
        assert not hasattr(window.canvas, "_last_name_to_structure")
    finally:
        window.close()


def test_cancelled_name_to_structure_job_ignores_late_result(monkeypatch) -> None:
    window = ChemusonWindow()
    try:
        monkeypatch.setattr(window, "_schedule_chemical_properties_update", lambda: None)
        called = False

        def _fake_question(*_args):
            nonlocal called
            called = True
            return QMessageBox.StandardButton.Yes

        monkeypatch.setattr(QMessageBox, "question", _fake_question)
        window._cancel_name_to_structure_job(9)

        window._on_name_to_structure_result(9, _result(), "")

        assert not called
        assert len(window.canvas.model.atoms) == 0
    finally:
        window.close()
