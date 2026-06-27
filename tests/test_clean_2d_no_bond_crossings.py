"""Regresiones de clean_2d_fallback para evitar cruces en árboles acíclicos."""


import pytest
from PyQt6.QtCore import QPointF
from PyQt6.QtWidgets import QApplication


from chemuson.gui.canvas import ChemusonCanvas
from chemuson.gui.geom import angle_deg, angle_distance_deg, segments_intersect


@pytest.fixture(scope="module", autouse=True)
def _qapp():
    return QApplication.instance() or QApplication([])


def _internal_angle(canvas: ChemusonCanvas, left_id: int, center_id: int, right_id: int) -> float:
    left = canvas.model.get_atom(left_id)
    center = canvas.model.get_atom(center_id)
    right = canvas.model.get_atom(right_id)
    a1 = angle_deg(QPointF(center.x, center.y), QPointF(left.x, left.y))
    a2 = angle_deg(QPointF(center.x, center.y), QPointF(right.x, right.y))
    return angle_distance_deg(a1, a2)


def test_clean_2d_fallback_avoids_nonadjacent_bond_crossings():
    """No debe generar cruces entre enlaces no adyacentes en un árbol ramificado."""
    canvas = ChemusonCanvas()

    # Caso determinista que producía cruce/fusión visual en clean_2d_fallback.
    coords = {
        1: (140.0, 200.0),
        2: (300.0, 120.0),
        3: (320.0, 120.0),
        4: (280.0, 320.0),
        5: (100.0, 160.0),
        6: (300.0, 240.0),
        7: (340.0, 300.0),
        8: (300.0, 200.0),
    }
    for atom_id in range(1, 9):
        x, y = coords[atom_id]
        canvas.model.add_atom("C", x, y, atom_id=atom_id)

    edges = [(1, 5), (1, 3), (5, 4), (5, 2), (1, 6), (1, 7), (5, 8)]
    for a1, a2 in edges:
        canvas.model.add_bond(a1, a2, order=1)

    canvas.clean_2d_fallback(atom_ids=set(coords.keys()), iterations=40)

    bonds = list(canvas.model.bonds.values())
    for idx, b1 in enumerate(bonds):
        for b2 in bonds[idx + 1 :]:
            # Ignorar enlaces adyacentes (comparten vértice).
            if {b1.a1_id, b1.a2_id} & {b2.a1_id, b2.a2_id}:
                continue
            a1 = canvas.model.get_atom(b1.a1_id)
            a2 = canvas.model.get_atom(b1.a2_id)
            c1 = canvas.model.get_atom(b2.a1_id)
            c2 = canvas.model.get_atom(b2.a2_id)
            assert not segments_intersect(
                QPointF(a1.x, a1.y),
                QPointF(a2.x, a2.y),
                QPointF(c1.x, c1.y),
                QPointF(c2.x, c2.y),
            )


def test_clean_2d_fallback_branched_sp3_backbone_not_linear():
    """Una cadena sp3 ramificada no debe colapsar a esqueleto lineal 180°."""
    canvas = ChemusonCanvas()

    backbone = []
    for i in range(7):
        atom = canvas.model.add_atom("C", 100.0 + i * 40.0, 220.0)
        backbone.append(atom.id)
    for i in range(6):
        canvas.model.add_bond(backbone[i], backbone[i + 1], order=1)

    # Ramas verticales en cada átomo del backbone.
    for atom_id in backbone:
        atom = canvas.model.get_atom(atom_id)
        up = canvas.model.add_atom("C", atom.x, atom.y - 40.0)
        down = canvas.model.add_atom("C", atom.x, atom.y + 40.0)
        canvas.model.add_bond(atom_id, up.id, order=1)
        canvas.model.add_bond(atom_id, down.id, order=1)

    canvas.clean_2d_fallback(atom_ids=set(canvas.model.atoms.keys()), iterations=40)

    for i in range(1, len(backbone) - 1):
        internal = _internal_angle(canvas, backbone[i - 1], backbone[i], backbone[i + 1])
        assert internal < 155.0
