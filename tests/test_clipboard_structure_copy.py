"""Regresiones para copiado de estructuras al portapapeles."""


import pytest
from PyQt6.QtGui import QImage
from PyQt6.QtWidgets import QApplication


import chemuson.gui.canvas as canvas_module
from chemuson.gui.canvas import ChemusonCanvas


@pytest.fixture(scope="module", autouse=True)
def _qapp():
    return QApplication.instance() or QApplication([])


def test_selected_structure_copy_dedupes_duplicate_bond_pairs() -> None:
    canvas = ChemusonCanvas()
    a1 = canvas.model.add_atom("C", 20.0, 20.0, is_explicit=True).id
    a2 = canvas.model.add_atom("C", 60.0, 20.0, is_explicit=True).id
    canvas.model.add_bond(a1, a2, order=1)
    canvas.model.add_bond(a1, a2, order=2)
    canvas._rebuild_items_from_model()

    canvas.scene.clearSelection()
    canvas.atom_items[a1].setSelected(True)
    canvas.atom_items[a2].setSelected(True)
    canvas._sync_selection_from_scene()

    atom_ids, bonds = canvas._selected_structure_ids()
    payload = canvas._build_selection_payload()

    assert atom_ids == {a1, a2}
    assert len(bonds) == 1
    assert payload is not None
    assert len(payload["bonds"]) == 1


def test_copy_to_clipboard_large_selection_uses_lightweight_export(monkeypatch) -> None:
    canvas = ChemusonCanvas()
    previous = None
    first = None
    atom_ids = []
    for index in range(22):
        atom = canvas.model.add_atom("C", 20.0 + index * 30.0, 40.0, is_explicit=True)
        atom_ids.append(atom.id)
        if first is None:
            first = atom.id
        if previous is not None:
            canvas.model.add_bond(previous, atom.id, order=1)
        previous = atom.id
    canvas._rebuild_items_from_model()

    canvas.scene.clearSelection()
    for atom_id in atom_ids:
        canvas.atom_items[atom_id].setSelected(True)
    canvas._sync_selection_from_scene()

    calls = {"molfile": 0, "smiles": 0, "svg": 0}

    def _count_molfile(_graph):
        calls["molfile"] += 1
        return "mock-molfile"

    def _count_smiles(_graph):
        calls["smiles"] += 1
        return "mock-smiles"

    def _count_svg(*args, **kwargs):
        calls["svg"] += 1
        return b"<svg/>"

    def _cheap_image(*args, **kwargs):
        image = QImage(16, 16, QImage.Format.Format_ARGB32)
        image.fill(0xFFFFFFFF)
        return image

    monkeypatch.setattr(canvas_module, "molgraph_to_molfile", _count_molfile)
    monkeypatch.setattr(canvas_module, "molgraph_to_smiles", _count_smiles)
    monkeypatch.setattr(canvas, "_render_scene_svg", _count_svg)
    monkeypatch.setattr(canvas, "_render_scene_image", _cheap_image)

    canvas.copy_to_clipboard()

    mime = QApplication.clipboard().mimeData()
    assert calls["molfile"] == 0
    assert calls["smiles"] == 0
    assert calls["svg"] == 0
    assert mime is not None
    assert mime.hasFormat("application/x-chemuson-selection")
