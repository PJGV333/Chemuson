"""Regresiones para proyección de H cuando la depicción 2D usa esqueleto sin H."""

import os
import sys

import pytest

sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), "..", "src")))

from chemuson.core.model import Bond
from chemuson.gui.main_window import ChemusonWindow


def test_project_missing_hydrogen_coords_keeps_local_vectors_for_ch3() -> None:
    before = {
        1: (0.0, 0.0),  # C
        2: (1.0, 0.0),  # H
        3: (0.0, 1.0),  # H
        4: (0.0, -1.0),  # H
    }
    after = {
        1: (10.0, 20.0),
    }
    bonds = [
        Bond(id=1, a1_id=1, a2_id=2, order=1),
        Bond(id=2, a1_id=1, a2_id=3, order=1),
        Bond(id=3, a1_id=1, a2_id=4, order=1),
    ]
    atom_elements = {1: "C", 2: "H", 3: "H", 4: "H"}

    projected = ChemusonWindow._project_missing_hydrogen_coords(before, after, bonds, atom_elements)

    assert projected[2] == pytest.approx((11.0, 20.0), rel=1e-9)
    assert projected[3] == pytest.approx((10.0, 21.0), rel=1e-9)
    assert projected[4] == pytest.approx((10.0, 19.0), rel=1e-9)


def test_project_missing_hydrogen_coords_averages_multiple_anchor_shifts() -> None:
    before = {
        1: (0.0, 0.0),  # C
        2: (4.0, 0.0),  # C
        3: (2.0, 1.0),  # bridging H
    }
    after = {
        1: (1.0, 0.5),  # shift (+1,+0.5)
        2: (5.0, 0.5),  # same shift
    }
    bonds = [
        Bond(id=1, a1_id=1, a2_id=3, order=1),
        Bond(id=2, a1_id=2, a2_id=3, order=1),
    ]
    atom_elements = {1: "C", 2: "C", 3: "H"}

    projected = ChemusonWindow._project_missing_hydrogen_coords(before, after, bonds, atom_elements)

    assert projected[3] == pytest.approx((3.0, 1.5), rel=1e-9)

