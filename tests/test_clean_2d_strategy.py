"""Estrategia de selección de backend para Limpiar 2D."""



from chemuson.core.model import Bond
from chemuson.gui.controllers.clean2d_controller import Clean2DController


def test_is_acyclic_structure_detects_chain() -> None:
    atom_ids = {1, 2, 3, 4}
    bonds = [
        Bond(id=1, a1_id=1, a2_id=2, order=1),
        Bond(id=2, a1_id=2, a2_id=3, order=1),
        Bond(id=3, a1_id=3, a2_id=4, order=1),
    ]
    assert Clean2DController._is_acyclic_structure(atom_ids, bonds)


def test_is_acyclic_structure_detects_cycle() -> None:
    atom_ids = {1, 2, 3}
    bonds = [
        Bond(id=1, a1_id=1, a2_id=2, order=1),
        Bond(id=2, a1_id=2, a2_id=3, order=1),
        Bond(id=3, a1_id=3, a2_id=1, order=1),
    ]
    assert not Clean2DController._is_acyclic_structure(atom_ids, bonds)
