"""Pruebas para especies atómicas con excepciones conocidas."""

from __future__ import annotations

import os
import sys

sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), "..", "src")))

from chemuson.gui.energy_diagrams import build_atomic_species_diagram


def test_atomic_species_builder_oxygen_neutral_has_eight_electrons() -> None:
    diagram = build_atomic_species_diagram("O")

    assert diagram.metadata["electron_count"] == 8
    assert sum(sum(level.occupancies) for level in diagram.levels) == 8


def test_atomic_species_builder_na_plus_has_ten_electrons() -> None:
    diagram = build_atomic_species_diagram("Na", charge=1)

    assert diagram.metadata["electron_count"] == 10
    assert diagram.metadata["configuration_string"] == "1s2 2s2 2p6"


def test_atomic_species_builder_cr_uses_known_exception_when_enabled() -> None:
    diagram = build_atomic_species_diagram("Cr", use_known_exceptions=True)

    assert diagram.metadata["configuration_string"].endswith("3d5 4s1")
    assert diagram.metadata["builder"]["name"] == "build_atomic_species_diagram"


def test_atomic_species_builder_cu_uses_known_exception_when_enabled() -> None:
    diagram = build_atomic_species_diagram("Cu", use_known_exceptions=True)

    assert diagram.metadata["configuration_string"].endswith("3d10 4s1")
    assert diagram.metadata["builder"]["params"]["use_known_exceptions"] is True
