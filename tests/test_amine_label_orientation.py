"""Regresión: aminas ordinarias dibujan NH2/H2N/NH horizontales, sin apilado.

En 0.3.0-beta.1 el H implícito de aminas primarias/secundarias se apilaba
verticalmente (p. ej. H₂ sobre N en un amino y bajo N en el otro de
1,1-diamino-3-fenilpropan-2-ona). La estrategia de cuatro direcciones de
``separated_nh`` queda reservada a contextos conjugados con contrato visual
(amidas, ureas, nitrógenos aromáticos); las aminas ordinarias se restringen
a posiciones horizontales, con la clasificación por conectividad de
``_nitrogen_geometry_context``.
"""

import math

import pytest
from PyQt6.QtCore import QPoint, Qt
from PyQt6.QtTest import QTest
from PyQt6.QtWidgets import QApplication

from chemuson.gui.canvas import ChemusonCanvas
from chemuson.gui.commands import AddBondCommand
from chemuson.gui.theme import apply_theme

BOND_LENGTH = 40.0


@pytest.fixture(scope="module", autouse=True)
def _qapp():
    return QApplication.instance() or QApplication([])


def _pol(x: float, y: float, angle_deg: float, length: float = BOND_LENGTH) -> tuple[float, float]:
    angle = math.radians(angle_deg)
    return (x + length * math.cos(angle), y + length * math.sin(angle))


def _build(atoms, bonds) -> tuple[ChemusonCanvas, dict[str, int]]:
    canvas = ChemusonCanvas()
    ids: dict[str, int] = {}
    for name, element, x, y in atoms:
        ids[name] = canvas.model.add_atom(element, x, y, is_explicit=True).id
    for a, b, order, aromatic in bonds:
        canvas.model.add_bond(ids[a], ids[b], order=order, is_aromatic=aromatic)
    canvas._rebuild_items_from_model()
    return canvas, ids


def _ethylamine():
    atoms = [("c1", "C", 100.0, 100.0), ("n", "N", 140.0, 100.0)]
    bonds = [("c1", "n", 1, False)]
    return atoms, bonds


def _diethylamine():
    atoms = [
        ("c1", "C", 60.0, 120.0),
        ("n", "N", 100.0, 80.0),
        ("c2", "C", 140.0, 120.0),
    ]
    bonds = [("c1", "n", 1, False), ("n", "c2", 1, False)]
    return atoms, bonds


def _aniline():
    center = (300.0, 300.0)
    ring = [
        (center[0] + 40.0 * math.cos(math.radians(60 * i - 90)),
         center[1] + 40.0 * math.sin(math.radians(60 * i - 90)))
        for i in range(6)
    ]
    atoms = [(f"r{i}", "C", x, y) for i, (x, y) in enumerate(ring)]
    bonds = [(f"r{i}", f"r{(i + 1) % 6}", 1, True) for i in range(6)]
    nx, ny = ring[0]
    nx = center[0] + 2.0 * (nx - center[0])
    ny = center[1] + 2.0 * (ny - center[1])
    atoms.append(("n", "N", nx, ny))
    bonds.append(("r0", "n", 1, False))
    return atoms, bonds


def _diamino_ketone():
    x1, y1 = 300.0, 300.0
    n_a = _pol(x1, y1, -75.0)
    n_b = _pol(x1, y1, 75.0)
    c2 = (x1 + BOND_LENGTH, y1)
    o2 = _pol(*c2, -90.0)
    c3 = _pol(*c2, 0.0)
    atoms = [
        ("c1", "C", x1, y1),
        ("na", "N", *n_a),
        ("nb", "N", *n_b),
        ("c2", "C", *c2),
        ("o2", "O", *o2),
        ("c3", "C", *c3),
    ]
    bonds = [
        ("c1", "na", 1, False),
        ("c1", "nb", 1, False),
        ("c1", "c2", 1, False),
        ("c2", "o2", 2, False),
        ("c2", "c3", 1, False),
    ]
    ring = [
        (c3[0] + 40.0 * math.cos(math.radians(60 * i - 90)),
         c3[1] + 40.0 * math.sin(math.radians(60 * i - 90)))
        for i in range(6)
    ]
    atoms += [(f"r{i}", "C", x, y) for i, (x, y) in enumerate(ring)]
    bonds += [("c3", "r0", 1, False)]
    bonds += [(f"r{i}", f"r{(i + 1) % 6}", 1, True) for i in range(6)]
    return atoms, bonds


def _acetamide():
    c2 = (140.0, 500.0)
    atoms = [
        ("c1", "C", 100.0, 500.0),
        ("c2", "C", *c2),
        ("o", "O", *_pol(*c2, -90.0)),
        ("n", "N", *_pol(*c2, 60.0)),
    ]
    bonds = [("c1", "c2", 1, False), ("c2", "o", 2, False), ("c2", "n", 1, False)]
    return atoms, bonds


def _diethylurea():
    cu = (400.0, 500.0)
    n1 = _pol(*cu, 150.0)
    n2 = _pol(*cu, 30.0)
    atoms = [
        ("c0", "C", *_pol(*n1, -60.0)),
        ("n1", "N", *n1),
        ("cu", "C", *cu),
        ("o", "O", *_pol(*cu, -90.0)),
        ("n2", "N", *n2),
        ("c3", "C", *_pol(*n2, -60.0)),
    ]
    bonds = [
        ("c0", "n1", 1, False),
        ("n1", "cu", 1, False),
        ("cu", "o", 2, False),
        ("cu", "n2", 1, False),
        ("n2", "c3", 1, False),
    ]
    return atoms, bonds


def _assert_horizontal(canvas: ChemusonCanvas, atom_id: int) -> None:
    """H implícito a la izquierda/derecha de N, sin apilado vertical y sin
    solapar los enlaces del nitrógeno."""
    item = canvas.atom_items[atom_id]
    position = item.implicit_hydrogen_position
    assert position in {"left", "right"}
    assert item.label.isVisible()
    assert item.implicit_hydrogen_label.isVisible()
    n_center = item.label.mapToScene(item.label.boundingRect().center())
    h_center = item.implicit_hydrogen_label.mapToScene(
        item.implicit_hydrogen_label.boundingRect().center()
    )
    assert abs(h_center.y() - n_center.y()) <= 1.5
    dx = h_center.x() - n_center.x()
    assert (dx > 0) == (position == "right")
    h_shape = item.implicit_hydrogen_label.mapToScene(
        item.implicit_hydrogen_label.shape()
    )
    for bond in canvas.model.bonds.values():
        if atom_id not in (bond.a1_id, bond.a2_id):
            continue
        bond_item = canvas.bond_items[bond.id]
        assert not h_shape.intersects(bond_item.mapToScene(bond_item.shape()))


def test_primary_amine_ethylamine_is_horizontal() -> None:
    canvas, ids = _build(*_ethylamine())
    assert canvas._nitrogen_geometry_context(ids["n"]) == "amine"
    label, _anchor, _offset = canvas._build_atom_label(canvas.model.get_atom(ids["n"]))
    assert label in {"NH2", "H2N"}
    assert canvas.atom_items[ids["n"]]._implicit_hydrogen_text == "H2"
    _assert_horizontal(canvas, ids["n"])


def test_secondary_amine_diethylamine_is_horizontal() -> None:
    canvas, ids = _build(*_diethylamine())
    assert canvas._nitrogen_geometry_context(ids["n"]) == "amine"
    label, _anchor, _offset = canvas._build_atom_label(canvas.model.get_atom(ids["n"]))
    assert label in {"NH", "HN"}
    _assert_horizontal(canvas, ids["n"])


def test_aniline_primary_amine_is_horizontal() -> None:
    canvas, ids = _build(*_aniline())
    assert canvas._nitrogen_geometry_context(ids["n"]) == "amine"
    label, _anchor, _offset = canvas._build_atom_label(canvas.model.get_atom(ids["n"]))
    assert label in {"NH2", "H2N"}
    _assert_horizontal(canvas, ids["n"])


def test_diamino_ketone_both_amines_horizontal_without_stacking() -> None:
    canvas, ids = _build(*_diamino_ketone())
    for name in ("na", "nb"):
        assert canvas._nitrogen_geometry_context(ids[name]) == "amine"
        label, _anchor, _offset = canvas._build_atom_label(canvas.model.get_atom(ids[name]))
        assert label in {"NH2", "H2N"}
        assert canvas.atom_items[ids[name]]._implicit_hydrogen_text == "H2"
        _assert_horizontal(canvas, ids[name])
    # Invariantes: la corrección es tipográfica, no geométrica.
    distances = [
        math.hypot(
            canvas.model.get_atom(bond.a2_id).x - canvas.model.get_atom(bond.a1_id).x,
            canvas.model.get_atom(bond.a2_id).y - canvas.model.get_atom(bond.a1_id).y,
        )
        for bond in canvas.model.bonds.values()
    ]
    assert distances == pytest.approx([BOND_LENGTH] * len(distances), abs=1e-7)


def test_amide_control_keeps_four_direction_strategy() -> None:
    canvas, ids = _build(*_acetamide())
    assert canvas._nitrogen_geometry_context(ids["n"]) == "amide"
    position = canvas.atom_items[ids["n"]].implicit_hydrogen_position
    assert position in {"right", "left", "above", "below"}


def test_urea_control_positions_preserved() -> None:
    canvas, ids = _build(*_diethylurea())
    for name in ("n1", "n2"):
        assert canvas._nitrogen_geometry_context(ids[name]) == "urea"
    positions = {
        name: canvas.atom_items[ids[name]].implicit_hydrogen_position
        for name in ("n1", "n2")
    }
    assert positions == {"n1": "left", "n2": "below"}


def test_charged_ammonium_is_not_split() -> None:
    canvas = ChemusonCanvas()
    c1 = canvas.model.add_atom("C", 600.0, 300.0, is_explicit=True).id
    n = canvas.model.add_atom("N", 640.0, 300.0, is_explicit=True).id
    canvas.model.get_atom(n).charge = 1
    canvas.model.get_atom(n).formal_charge = 1
    canvas.model.add_bond(c1, n, order=1)
    canvas._rebuild_items_from_model()
    assert canvas._nitrogen_geometry_context(n) == "other"
    assert canvas.atom_items[n]._implicit_hydrogen_text == ""


def test_edit_undo_redo_refresh_keeps_amine_horizontal() -> None:
    canvas, ids = _build(*_ethylamine())
    _assert_horizontal(canvas, ids["n"])

    extra = canvas.model.add_atom("C", 120.0, 60.0, is_explicit=True).id
    canvas.add_atom_item(canvas.model.get_atom(extra))
    canvas.undo_stack.push(AddBondCommand(canvas.model, canvas, ids["n"], extra))
    assert canvas.model.implicit_h_count(ids["n"]) == 1
    _assert_horizontal(canvas, ids["n"])

    canvas.undo_stack.undo()
    assert canvas.model.implicit_h_count(ids["n"]) == 2
    _assert_horizontal(canvas, ids["n"])

    canvas.undo_stack.redo()
    assert canvas.model.implicit_h_count(ids["n"]) == 1
    _assert_horizontal(canvas, ids["n"])


def test_drag_move_refresh_keeps_amine_horizontal() -> None:
    canvas, ids = _build(*_diethylamine())
    _assert_horizontal(canvas, ids["n"])

    canvas.resize(900, 700)
    canvas.show()
    QApplication.processEvents()
    item = canvas.atom_items[ids["n"]]
    item.setSelected(True)
    canvas._sync_selection_from_scene()
    start_view = canvas.mapFromScene(item.scenePos())
    end_view = QPoint(start_view.x(), start_view.y() + 30)
    QTest.mousePress(canvas.viewport(), Qt.MouseButton.LeftButton, pos=start_view)
    QTest.mouseMove(canvas.viewport(), end_view)
    QTest.mouseRelease(canvas.viewport(), Qt.MouseButton.LeftButton, pos=end_view)
    QApplication.processEvents()

    assert canvas.model.get_atom(ids["n"]).y == pytest.approx(110.0, abs=1.0)
    _assert_horizontal(canvas, ids["n"])
    canvas.undo_stack.undo()
    QApplication.processEvents()
    _assert_horizontal(canvas, ids["n"])


@pytest.mark.parametrize("theme", ["light", "dark"])
def test_theme_change_keeps_amine_horizontal(theme: str) -> None:
    canvas, ids = _build(*_diethylamine())
    assert apply_theme(canvas, theme) == theme
    _assert_horizontal(canvas, ids["n"])
