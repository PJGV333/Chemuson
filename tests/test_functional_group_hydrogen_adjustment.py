"""Regresiones para grupos funcionales simples con H ajustables."""

import os
import re
import sys

import pytest
from PyQt6.QtWidgets import QApplication

sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), "..", "src")))

from chemuson.chemio.persistence import PersistenceManager
from chemuson.gui.canvas import ChemusonCanvas
from chemuson.gui.commands import AddBondCommand, ChangeAtomCommand


@pytest.fixture(scope="module", autouse=True)
def _qapp():
    return QApplication.instance() or QApplication([])


def _label_text(canvas: ChemusonCanvas, atom_id: int) -> str:
    atom = canvas.model.get_atom(atom_id)
    label, _anchor, _offset = canvas._build_atom_label(atom)
    return label


def _label_hydrogen_count(label: str) -> int:
    match = re.search(r"H(\d*)", label)
    if match is None:
        return 0
    digits = match.group(1)
    return int(digits) if digits else 1


def _bonded_placeholder(canvas: ChemusonCanvas) -> tuple[int, int]:
    anchor = canvas.model.add_atom("C", 10.0, 10.0, is_explicit=True).id
    target = canvas.model.add_atom("C", 50.0, 10.0, is_explicit=True).id
    canvas.model.add_bond(anchor, target, order=1)
    canvas._rebuild_items_from_model()
    return anchor, target


def test_nh2_group_shorthand_reduces_hydrogens_when_new_bond_is_added() -> None:
    canvas = ChemusonCanvas()
    _anchor_id, target_id = _bonded_placeholder(canvas)

    canvas.undo_stack.push(ChangeAtomCommand(canvas.model, canvas, target_id, "NH2"))

    atom = canvas.model.get_atom(target_id)
    assert atom.element == "N"
    assert atom.group_h_cap == 2
    assert atom.no_implicit is True
    assert canvas.model.assigned_hydrogen_count(target_id) == 2
    assert _label_hydrogen_count(_label_text(canvas, target_id)) == 2
    assert target_id not in canvas.validate_structure()

    extra = canvas.model.add_atom("C", 90.0, 40.0, is_explicit=True).id
    canvas.add_atom_item(canvas.model.get_atom(extra))
    canvas.undo_stack.push(AddBondCommand(canvas.model, canvas, target_id, extra))

    assert canvas.model.assigned_hydrogen_count(target_id) == 1
    assert _label_hydrogen_count(_label_text(canvas, target_id)) == 1
    assert target_id not in canvas.validate_structure()


@pytest.mark.parametrize(("shorthand", "element"), [("OH", "O"), ("SH", "S")])
def test_oh_and_sh_shorthand_drop_h_after_second_bond(shorthand: str, element: str) -> None:
    canvas = ChemusonCanvas()
    _anchor_id, target_id = _bonded_placeholder(canvas)

    canvas.undo_stack.push(ChangeAtomCommand(canvas.model, canvas, target_id, shorthand))

    atom = canvas.model.get_atom(target_id)
    assert atom.element == element
    assert atom.group_h_cap == 1
    assert canvas.model.assigned_hydrogen_count(target_id) == 1
    assert _label_hydrogen_count(_label_text(canvas, target_id)) == 1

    extra = canvas.model.add_atom("C", 90.0, 40.0, is_explicit=True).id
    canvas.add_atom_item(canvas.model.get_atom(extra))
    canvas.undo_stack.push(AddBondCommand(canvas.model, canvas, target_id, extra))

    assert canvas.model.assigned_hydrogen_count(target_id) == 0
    assert _label_hydrogen_count(_label_text(canvas, target_id)) == 0
    assert _label_text(canvas, target_id) == element
    assert target_id not in canvas.validate_structure()


def test_simple_hydrogen_group_marks_overvalence_after_excess_bonds() -> None:
    canvas = ChemusonCanvas()
    _anchor_id, target_id = _bonded_placeholder(canvas)

    canvas.undo_stack.push(ChangeAtomCommand(canvas.model, canvas, target_id, "NH2"))

    positions = [(90.0, 10.0), (50.0, 50.0), (50.0, -30.0)]
    for x, y in positions:
        extra = canvas.model.add_atom("C", x, y, is_explicit=True).id
        canvas.add_atom_item(canvas.model.get_atom(extra))
        canvas.undo_stack.push(AddBondCommand(canvas.model, canvas, target_id, extra))

    assert canvas.model.assigned_hydrogen_count(target_id) == 0
    errors = canvas.validate_structure()
    assert target_id in errors
    assert canvas.atom_items[target_id]._valence_error is True


def test_simple_hydrogen_group_persists_after_save_load(tmp_path) -> None:
    source = ChemusonCanvas()
    _anchor_id, target_id = _bonded_placeholder(source)
    source.undo_stack.push(ChangeAtomCommand(source.model, source, target_id, "NH2"))

    file_path = tmp_path / "simple_group.cmsn"
    PersistenceManager.save_to_file(str(file_path), source)

    restored = ChemusonCanvas()
    PersistenceManager.load_from_file(str(file_path), restored)

    atom = restored.model.get_atom(target_id)
    assert atom.element == "N"
    assert atom.group_h_cap == 2
    assert atom.no_implicit is True
    assert restored.model.assigned_hydrogen_count(target_id) == 2
    assert _label_hydrogen_count(_label_text(restored, target_id)) == 2


def test_simple_hydrogen_group_survives_internal_copy_paste() -> None:
    source = ChemusonCanvas()
    _anchor_id, target_id = _bonded_placeholder(source)
    source.undo_stack.push(ChangeAtomCommand(source.model, source, target_id, "NH2"))

    source.atom_items[target_id].setSelected(True)
    source._sync_selection_from_scene()
    payload = source._build_selection_payload()
    assert payload is not None

    target = ChemusonCanvas()
    target._last_scene_pos = source.mapToScene(source.viewport().rect().center())
    target._paste_selection_payload(payload)

    inserted = [atom for atom in target.model.atoms.values() if atom.element == "N"]
    assert len(inserted) == 1
    atom = inserted[0]
    assert atom.group_h_cap == 2
    assert atom.no_implicit is True
    assert target.model.assigned_hydrogen_count(atom.id) == 2


@pytest.mark.parametrize(
    ("shorthand", "expected_labels", "expect_error"),
    [
        ("NH2", ["NH", "N", "N"], True),
        ("OH", ["O", "O", "O"], True),
        ("SH", ["S", "S", "S"], False),
    ],
)
def test_legacy_raw_simple_group_aliases_still_update_h_count_and_overvalence(
    shorthand: str,
    expected_labels: list[str],
    expect_error: bool,
) -> None:
    canvas = ChemusonCanvas()
    anchor = canvas.model.add_atom("C", 10.0, 10.0, is_explicit=True).id
    target = canvas.model.add_atom("N", 50.0, 10.0, is_explicit=True).id
    canvas.model.add_bond(anchor, target, order=1)
    canvas._rebuild_items_from_model()

    atom = canvas.model.get_atom(target)
    atom.element = shorthand
    atom.group_h_cap = None
    atom.explicit_h = None
    atom.no_implicit = False
    canvas._refresh_atom_label(target)

    positions = [(90.0, 10.0), (50.0, 50.0), (50.0, -30.0)]
    seen_labels: list[str] = []
    for x, y in positions:
        extra = canvas.model.add_atom("C", x, y, is_explicit=True).id
        canvas.add_atom_item(canvas.model.get_atom(extra))
        canvas.undo_stack.push(AddBondCommand(canvas.model, canvas, target, extra))
        seen_labels.append(_label_text(canvas, target))

    assert seen_labels == expected_labels
    errors = canvas.validate_structure()
    assert (target in errors) is expect_error
    assert canvas.atom_items[target]._valence_error is expect_error


def test_model_add_atom_normalizes_simple_hydrogen_group_aliases() -> None:
    canvas = ChemusonCanvas()
    atom = canvas.model.add_atom("NH2", 40.0, 40.0, is_explicit=True)

    assert atom.element == "N"
    assert atom.group_h_cap == 2
    assert atom.explicit_h is None
    assert atom.no_implicit is True
