"""Pruebas base de polímeros y Markush."""

from __future__ import annotations

import os
import sys

import pytest
from PyQt6.QtCore import QRectF
from PyQt6.QtWidgets import QApplication

sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), "..", "src")))

from chemuson.chemio.persistence import PersistenceManager
from chemuson.gui.canvas import ChemusonCanvas
from chemuson.gui.commands import ChangeBracketRepeatLabelCommand, ChangeRGroupSubstituentsCommand
from chemuson.markush import sanitize_r_group_substituents, set_r_group_substituents, summarize_markush


@pytest.fixture(scope="module", autouse=True)
def _qapp():
    return QApplication.instance() or QApplication([])


def test_polymer_repeat_label_persists_on_bracket() -> None:
    canvas = ChemusonCanvas()
    bracket = canvas.add_bracket_item(
        QRectF(10.0, 20.0, 80.0, 60.0),
        "[]",
        repeat_label="n",
    )

    data = PersistenceManager.save_to_dict(canvas)
    restored = ChemusonCanvas()
    PersistenceManager.load_from_dict(data, restored)

    assert bracket.repeat_label() == "n"
    assert data["canvas"]["annotations"]["brackets"][0]["repeat_label"] == "n"
    assert len(restored.bracket_items) == 1
    assert restored.bracket_items[0].repeat_label() == "n"


def test_markush_summary_detects_r_groups_and_polymer_repeats() -> None:
    canvas = ChemusonCanvas()
    r1 = canvas.model.add_atom("R1", 0.0, 0.0, is_explicit=True)
    carbon = canvas.model.add_atom("C", 40.0, 0.0)
    canvas.model.add_bond(r1.id, carbon.id, order=1)
    canvas.add_bracket_item(QRectF(0.0, -20.0, 60.0, 40.0), "[]", repeat_label="x")

    summary = summarize_markush(canvas.model, bracket_items=canvas.bracket_items)

    assert summary.has_markush
    assert summary.r_groups[0].atom_id == r1.id
    assert summary.r_groups[0].label == "R1"
    assert summary.polymer_repeats[0].repeat_label == "x"


def test_markush_summary_includes_r_group_substituents_and_dict() -> None:
    canvas = ChemusonCanvas()
    r1 = canvas.model.add_atom("R1", 0.0, 0.0, is_explicit=True, r_group_substituents=("Me", "Et"))

    summary = summarize_markush(canvas.model)
    payload = summary.as_dict()

    assert summary.r_groups[0].allowed_substituents == ("Me", "Et")
    assert payload["r_groups"][0]["allowed_substituents"] == ["Me", "Et"]
    assert payload["r_groups"][0]["atom_id"] == r1.id


def test_r_group_substituents_persist_in_cmsn() -> None:
    canvas = ChemusonCanvas()
    canvas.model.add_atom("R", 0.0, 0.0, is_explicit=True, r_group_substituents=("Me", "Ph"))

    data = PersistenceManager.save_to_dict(canvas)
    restored = ChemusonCanvas()
    PersistenceManager.load_from_dict(data, restored)

    restored_atom = next(iter(restored.model.atoms.values()))
    assert restored_atom.r_group_substituents == ("Me", "Ph")


def test_set_r_group_substituents_sanitizes_and_skips_non_r_atoms() -> None:
    graph = ChemusonCanvas().model
    r = graph.add_atom("R2", 0.0, 0.0, is_explicit=True)
    carbon = graph.add_atom("C", 40.0, 0.0)

    values = set_r_group_substituents(graph, [r.id, carbon.id], "Me, Et; Ph, Me, bad value!")

    assert values == ("Me", "Et", "Ph", "badvalue")
    assert r.r_group_substituents == values
    assert carbon.r_group_substituents == ()


def test_r_group_substituent_change_command_is_undoable() -> None:
    canvas = ChemusonCanvas()
    r = canvas.model.add_atom("R", 0.0, 0.0, is_explicit=True, r_group_substituents=("Me",))

    canvas.undo_stack.push(ChangeRGroupSubstituentsCommand(canvas.model, [r.id], ("Et", "Ph")))

    assert r.r_group_substituents == ("Et", "Ph")
    canvas.undo_stack.undo()
    assert r.r_group_substituents == ("Me",)
    canvas.undo_stack.redo()
    assert r.r_group_substituents == ("Et", "Ph")


def test_repeat_label_is_sanitized_for_persistence() -> None:
    canvas = ChemusonCanvas()
    bracket = canvas.add_bracket_item(QRectF(0.0, 0.0, 40.0, 40.0), "[]", repeat_label="n; rm -rf")

    assert bracket.repeat_label() == "nrm-rf"


def test_sanitize_r_group_substituents_accepts_lists() -> None:
    assert sanitize_r_group_substituents(["Me", "Et", "Me", "Ph "]) == ("Me", "Et", "Ph")


def test_polymer_repeat_change_command_is_undoable() -> None:
    canvas = ChemusonCanvas()
    bracket = canvas.add_bracket_item(QRectF(0.0, 0.0, 40.0, 40.0), "[]")

    canvas.undo_stack.push(ChangeBracketRepeatLabelCommand(canvas, [bracket], "n"))

    assert bracket.repeat_label() == "n"
    canvas.undo_stack.undo()
    assert bracket.repeat_label() == ""
    canvas.undo_stack.redo()
    assert bracket.repeat_label() == "n"
