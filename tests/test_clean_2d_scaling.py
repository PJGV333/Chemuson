"""Regresiones de escala para Limpiar 2D con coordenadas de RDKit."""

import os
import sys
import math

import pytest

sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), "..", "src")))

from core.model import Bond
from gui.main_window import ChemusonWindow


def _chain_bonds() -> list[Bond]:
    return [
        Bond(id=1, a1_id=1, a2_id=2, order=1),
        Bond(id=2, a1_id=2, a2_id=3, order=1),
    ]


def _cycle_bonds(size: int) -> list[Bond]:
    bonds: list[Bond] = []
    for idx in range(size):
        a1 = idx + 1
        a2 = ((idx + 1) % size) + 1
        bonds.append(Bond(id=idx + 1, a1_id=a1, a2_id=a2, order=1))
    return bonds


def _regular_polygon(size: int, radius: float, rotation_deg: float = 0.0) -> dict[int, tuple[float, float]]:
    coords: dict[int, tuple[float, float]] = {}
    for idx in range(size):
        angle = math.radians(rotation_deg + (360.0 * idx / float(size)))
        coords[idx + 1] = (radius * math.cos(angle), radius * math.sin(angle))
    return coords


def test_rescale_rdkit_coords_matches_target_bond_length() -> None:
    # Geometría típica de RDKit (~1.5 unidades por enlace).
    rdkit_coords = {
        1: (0.0, 0.0),
        2: (1.5, 0.0),
        3: (3.0, 0.0),
    }
    bonds = _chain_bonds()

    scaled = ChemusonWindow._rescale_coords_to_bond_length(rdkit_coords, bonds, target_length=40.0)
    avg = ChemusonWindow._average_bond_length(scaled, bonds)

    assert avg == pytest.approx(40.0, rel=1e-6)


def test_rescale_rdkit_coords_preserves_center() -> None:
    rdkit_coords = {
        1: (0.0, 0.0),
        2: (1.5, 0.0),
        3: (3.0, 0.0),
    }
    bonds = _chain_bonds()
    before_center = ChemusonWindow._coords_center(rdkit_coords)

    scaled = ChemusonWindow._rescale_coords_to_bond_length(rdkit_coords, bonds, target_length=40.0)
    after_center = ChemusonWindow._coords_center(scaled)

    assert after_center[0] == pytest.approx(before_center[0], rel=1e-12)
    assert after_center[1] == pytest.approx(before_center[1], rel=1e-12)


def test_rescale_rdkit_coords_no_bonds_returns_same_geometry() -> None:
    coords = {
        1: (10.0, 10.0),
        2: (20.0, 15.0),
    }
    scaled = ChemusonWindow._rescale_coords_to_bond_length(coords, [], target_length=40.0)
    assert scaled == coords


def test_partial_blend_of_rotated_cycle_keeps_size_after_rescale() -> None:
    # Simula "before" y "after" con la misma molécula, pero orientaciones distintas.
    # El blend parcial encoge tamaño; la renormalización debe devolverlo al original.
    size = 6
    bonds = _cycle_bonds(size)
    before = _regular_polygon(size=size, radius=40.0, rotation_deg=0.0)
    after_rotated = _regular_polygon(size=size, radius=40.0, rotation_deg=30.0)
    target = ChemusonWindow._average_bond_length(before, bonds)
    assert target > 1e-6

    step_ratio = 0.35
    blended = {}
    for atom_id, (tx, ty) in after_rotated.items():
        bx, by = before[atom_id]
        blended[atom_id] = (
            bx + (tx - bx) * step_ratio,
            by + (ty - by) * step_ratio,
        )

    shrunk = ChemusonWindow._average_bond_length(blended, bonds)
    assert shrunk < target

    restored = ChemusonWindow._rescale_coords_to_bond_length(blended, bonds, target)
    restored_avg = ChemusonWindow._average_bond_length(restored, bonds)
    assert restored_avg == pytest.approx(target, rel=1e-6)
