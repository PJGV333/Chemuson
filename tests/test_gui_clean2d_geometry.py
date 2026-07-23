"""Regresiones sin Qt para la geometría de integración Clean2D."""

from __future__ import annotations

import math

import pytest

from chemuson.core.model import Bond
from chemuson.gui.clean2d_geometry import (
    align_coords_to_reference,
    average_bond_length,
    coords_center,
    project_missing_hydrogen_coords,
    rescale_coords_to_bond_length,
)


def _chain_bonds() -> list[Bond]:
    return [
        Bond(id=1, a1_id=1, a2_id=2, order=1),
        Bond(id=2, a1_id=2, a2_id=3, order=1),
    ]


def _transform(
    coords: dict[int, tuple[float, float]],
    *,
    rotation_deg: float,
    mirror_x: bool = False,
) -> dict[int, tuple[float, float]]:
    angle = math.radians(rotation_deg)
    cos_t = math.cos(angle)
    sin_t = math.sin(angle)
    transformed: dict[int, tuple[float, float]] = {}
    for atom_id, (x, y) in coords.items():
        x = -x if mirror_x else x
        transformed[atom_id] = (
            x * cos_t - y * sin_t + 120.0,
            x * sin_t + y * cos_t - 45.0,
        )
    return transformed


def test_rescale_matches_target_and_preserves_center() -> None:
    coords = {1: (0.0, 0.0), 2: (1.5, 0.0), 3: (3.0, 0.0)}
    bonds = _chain_bonds()

    scaled = rescale_coords_to_bond_length(coords, bonds, 40.0)

    assert average_bond_length(scaled, bonds) == pytest.approx(40.0)
    assert coords_center(scaled) == pytest.approx(coords_center(coords))


def test_rescale_without_bonds_preserves_geometry() -> None:
    coords = {1: (10.0, 10.0), 2: (20.0, 15.0)}

    assert rescale_coords_to_bond_length(coords, [], 40.0) == coords


@pytest.mark.parametrize(
    ("rotation_deg", "mirror_x"),
    [(90.0, False), (35.0, True)],
)
def test_alignment_recovers_rotation_and_reflection(
    rotation_deg: float,
    mirror_x: bool,
) -> None:
    reference = {
        1: (-35.0, -10.0),
        2: (-8.0, 22.0),
        3: (26.0, 18.0),
        4: (34.0, -14.0),
        5: (3.0, -28.0),
    }

    aligned = align_coords_to_reference(
        reference,
        _transform(
            reference,
            rotation_deg=rotation_deg,
            mirror_x=mirror_x,
        ),
    )

    for atom_id, expected in reference.items():
        assert aligned[atom_id] == pytest.approx(expected, abs=1e-6)


def test_single_common_atom_alignment_translates_geometry() -> None:
    reference = {1: (20.0, 30.0)}
    coords = {1: (2.0, 3.0), 2: (7.0, 9.0)}

    assert align_coords_to_reference(reference, coords) == {
        1: (20.0, 30.0),
        2: (25.0, 36.0),
    }


def test_missing_hydrogens_keep_local_vectors() -> None:
    before = {
        1: (0.0, 0.0),
        2: (1.0, 0.0),
        3: (0.0, 1.0),
        4: (0.0, -1.0),
    }
    after = {1: (10.0, 20.0)}
    bonds = [
        Bond(id=1, a1_id=1, a2_id=2, order=1),
        Bond(id=2, a1_id=1, a2_id=3, order=1),
        Bond(id=3, a1_id=1, a2_id=4, order=1),
    ]

    projected = project_missing_hydrogen_coords(
        before,
        after,
        bonds,
        {1: "C", 2: "H", 3: "H", 4: "H"},
    )

    assert projected[2] == pytest.approx((11.0, 20.0))
    assert projected[3] == pytest.approx((10.0, 21.0))
    assert projected[4] == pytest.approx((10.0, 19.0))
