"""Regresiones de geometría para clean_2d_fallback en cadenas alifáticas."""

import math

import pytest
from PyQt6.QtCore import QPointF
from PyQt6.QtWidgets import QApplication


from chemuson.gui.canvas import ChemusonCanvas
from chemuson.gui.geom import angle_deg, angle_distance_deg


@pytest.fixture(scope="module", autouse=True)
def _qapp():
    return QApplication.instance() or QApplication([])


def _internal_angle(canvas: ChemusonCanvas, left_id: int, center_id: int, right_id: int) -> float:
    left = canvas.model.get_atom(left_id)
    center = canvas.model.get_atom(center_id)
    right = canvas.model.get_atom(right_id)
    center_pt = QPointF(center.x, center.y)
    a1 = angle_deg(center_pt, QPointF(left.x, left.y))
    a2 = angle_deg(center_pt, QPointF(right.x, right.y))
    return angle_distance_deg(a1, a2)


def test_clean_2d_fallback_linear_alkane_becomes_zigzag() -> None:
    canvas = ChemusonCanvas()
    length = float(canvas.state.bond_length)

    # Esqueleto totalmente recto (caso problemático reportado).
    c1 = canvas.model.add_atom("C", 100.0, 240.0)
    c2 = canvas.model.add_atom("C", 140.0, 240.0)
    c3 = canvas.model.add_atom("C", 180.0, 240.0)
    c4 = canvas.model.add_atom("C", 220.0, 240.0)
    canvas.model.add_bond(c1.id, c2.id, order=1)
    canvas.model.add_bond(c2.id, c3.id, order=1)
    canvas.model.add_bond(c3.id, c4.id, order=1)

    # H explícitos ortogonales al eje de la cadena.
    h_specs = [
        (c1.id, 60.0, 240.0), (c1.id, 100.0, 200.0), (c1.id, 100.0, 280.0),
        (c2.id, 140.0, 200.0), (c2.id, 140.0, 280.0),
        (c3.id, 180.0, 200.0), (c3.id, 180.0, 280.0),
        (c4.id, 260.0, 240.0), (c4.id, 220.0, 200.0), (c4.id, 220.0, 280.0),
    ]
    c4_h_ids: list[int] = []
    for anchor_id, x, y in h_specs:
        h = canvas.model.add_atom("H", x, y, is_explicit=True)
        canvas.model.add_bond(anchor_id, h.id, order=1)
        if anchor_id == c4.id:
            c4_h_ids.append(h.id)

    atom_ids = set(canvas.model.atoms.keys())
    canvas.clean_2d_fallback(atom_ids=atom_ids, iterations=20)

    # En cadena alifática sp3 no debe quedar lineal (180°) en los carbonos internos.
    c2_angle = _internal_angle(canvas, c1.id, c2.id, c3.id)
    c3_angle = _internal_angle(canvas, c2.id, c3.id, c4.id)
    assert c2_angle < 150.0
    assert c3_angle < 150.0

    # Debe conservar longitud base de enlaces C-C.
    for a_id, b_id in ((c1.id, c2.id), (c2.id, c3.id), (c3.id, c4.id)):
        a = canvas.model.get_atom(a_id)
        b = canvas.model.get_atom(b_id)
        dist = math.hypot(b.x - a.x, b.y - a.y)
        assert dist == pytest.approx(length, rel=0.08)

    # CH3 terminal: al menos un H debe quedar opuesto al enlace C-C (hacia afuera).
    c4_atom = canvas.model.get_atom(c4.id)
    c3_atom = canvas.model.get_atom(c3.id)
    incoming = angle_deg(QPointF(c4_atom.x, c4_atom.y), QPointF(c3_atom.x, c3_atom.y))
    outward = (incoming + 180.0) % 360.0
    h_angles = []
    for h_id in c4_h_ids:
        h_atom = canvas.model.get_atom(h_id)
        h_angles.append(angle_deg(QPointF(c4_atom.x, c4_atom.y), QPointF(h_atom.x, h_atom.y)))
    assert min(angle_distance_deg(a, outward) for a in h_angles) < 15.0
