"""Regresiones para rutas costosas de canvas que deben permanecer incrementales."""

import os
import sys

import pytest
from PyQt6.QtWidgets import QApplication

sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), "..", "src")))

from chemuson.gui.canvas import ChemusonCanvas
from chemuson.gui.commands import AddAtomCommand, MoveAtomsCommand


@pytest.fixture(scope="module", autouse=True)
def _qapp():
    return QApplication.instance() or QApplication([])


def test_add_atom_command_does_not_force_full_scene_resync(monkeypatch) -> None:
    canvas = ChemusonCanvas()
    calls = {"count": 0}

    def _count_sync():
        calls["count"] += 1

    monkeypatch.setattr(canvas, "_sync_scene_with_model", _count_sync)

    canvas.undo_stack.push(AddAtomCommand(canvas.model, canvas, "C", 20.0, 20.0, is_explicit=True))

    assert calls["count"] == 0


def test_move_atoms_command_does_not_force_full_scene_resync(monkeypatch) -> None:
    canvas = ChemusonCanvas()
    a1 = canvas.model.add_atom("C", 20.0, 20.0, is_explicit=True).id
    a2 = canvas.model.add_atom("C", 60.0, 20.0, is_explicit=True).id
    canvas.model.add_bond(a1, a2, order=1)
    canvas._rebuild_items_from_model()

    calls = {"count": 0}

    def _count_sync():
        calls["count"] += 1

    monkeypatch.setattr(canvas, "_sync_scene_with_model", _count_sync)

    before = {a1: (20.0, 20.0), a2: (60.0, 20.0)}
    after = {a1: (40.0, 35.0), a2: (80.0, 35.0)}
    canvas.undo_stack.push(MoveAtomsCommand(canvas.model, canvas, before, after))

    assert calls["count"] == 0
    assert canvas.model.get_atom(a1).x == pytest.approx(40.0)
    assert canvas.model.get_atom(a2).y == pytest.approx(35.0)


def test_paste_selection_payload_batches_structure_validation(monkeypatch) -> None:
    canvas = ChemusonCanvas()
    validations = {"count": 0}

    def _count_validate():
        validations["count"] += 1
        return []

    monkeypatch.setattr(canvas, "validate_structure", _count_validate)

    payload = {
        "atoms": [
            {"id": 1, "element": "C", "x": 0.0, "y": 0.0, "is_explicit": True},
            {"id": 2, "element": "N", "x": 40.0, "y": 0.0, "is_explicit": True},
            {"id": 3, "element": "O", "x": 80.0, "y": 0.0, "is_explicit": True},
        ],
        "bonds": [
            {"a1": 1, "a2": 2, "order": 1},
            {"a1": 2, "a2": 3, "order": 1},
        ],
        "arrows": [],
        "brackets": [],
        "texts": [],
        "images": [],
    }

    canvas._paste_selection_payload(payload)

    assert validations["count"] == 1
