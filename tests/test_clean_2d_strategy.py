"""Estrategia de selección de backend para Limpiar 2D."""

import os
import sys

sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), "..", "src")))

from chemuson.core.model import Bond
from chemuson.gui.main_window import ChemusonWindow


def test_is_acyclic_structure_detects_chain() -> None:
    atom_ids = {1, 2, 3, 4}
    bonds = [
        Bond(id=1, a1_id=1, a2_id=2, order=1),
        Bond(id=2, a1_id=2, a2_id=3, order=1),
        Bond(id=3, a1_id=3, a2_id=4, order=1),
    ]
    assert ChemusonWindow._is_acyclic_structure(atom_ids, bonds)


def test_is_acyclic_structure_detects_cycle() -> None:
    atom_ids = {1, 2, 3}
    bonds = [
        Bond(id=1, a1_id=1, a2_id=2, order=1),
        Bond(id=2, a1_id=2, a2_id=3, order=1),
        Bond(id=3, a1_id=3, a2_id=1, order=1),
    ]
    assert not ChemusonWindow._is_acyclic_structure(atom_ids, bonds)

