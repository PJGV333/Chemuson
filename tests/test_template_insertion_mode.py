"""Regresiones para inserción de plantillas por clic."""


import pytest
from PyQt6.QtCore import QPointF
from PyQt6.QtWidgets import QApplication


from chemuson.core.model import BondStyle, MolGraph
from chemuson.gui.canvas import ChemusonCanvas


@pytest.fixture(scope="module", autouse=True)
def _qapp():
    return QApplication.instance() or QApplication([])


def _simple_template_graph() -> MolGraph:
    graph = MolGraph()
    a1 = graph.add_atom("C", 0.0, 0.0, is_explicit=False)
    a2 = graph.add_atom("C", 40.0, 0.0, is_explicit=False)
    graph.add_bond(a1.id, a2.id, order=1, style=BondStyle.PLAIN)
    return graph


def test_insert_template_at_atom_creates_connection_bond():
    """Al colocar sobre un átomo, la plantilla debe quedar unida a la estructura."""
    canvas = ChemusonCanvas()
    anchor = canvas.model.add_atom("C", 200.0, 200.0, is_explicit=False)
    canvas._rebuild_items_from_model()

    before_bonds = len(canvas.model.bonds)
    before_atoms = len(canvas.model.atoms)

    canvas._insert_molgraph_at(
        _simple_template_graph(),
        QPointF(anchor.x, anchor.y),
        attach_to_atom_id=anchor.id,
    )

    assert len(canvas.model.atoms) == before_atoms + 2
    assert len(canvas.model.bonds) == before_bonds + 2
    connected = [
        bond
        for bond in canvas.model.bonds.values()
        if bond.a1_id == anchor.id or bond.a2_id == anchor.id
    ]
    assert len(connected) == 1
