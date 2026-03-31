"""Pruebas para el builder de diagramas MO diatómicos."""

from __future__ import annotations

import os
import sys

import pytest

sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), "..", "src")))

from chemuson.gui.energy_diagrams import build_diatomic_mo_diagram


def test_o2_heavy_2p_has_two_unpaired_electrons() -> None:
    diagram = build_diatomic_mo_diagram(
        "O",
        "O",
        12,
        ordering="heavy_2p",
    )

    assert diagram.metadata["ordering"] == "heavy_2p"
    assert diagram.metadata["unpaired_electrons"] == 2


def test_n2_light_2p_has_bond_order_three() -> None:
    diagram = build_diatomic_mo_diagram(
        "N",
        "N",
        10,
        ordering="light_2p",
    )

    assert diagram.metadata["ordering"] == "light_2p"
    assert diagram.metadata["bond_order"] == pytest.approx(3.0)
