"""Pruebas para el builder de diagramas de campo ligando."""

from __future__ import annotations



from chemuson.gui.energy_diagrams import build_ligand_field_diagram


def test_octahedral_d6_high_and_low_spin_differ_in_unpaired_electrons() -> None:
    high_spin = build_ligand_field_diagram(6, geometry="octahedral", spin_mode="high")
    low_spin = build_ligand_field_diagram(6, geometry="octahedral", spin_mode="low")

    assert high_spin.metadata["unpaired_electrons"] == 4
    assert low_spin.metadata["unpaired_electrons"] == 0
    assert high_spin.metadata["unpaired_electrons"] > low_spin.metadata["unpaired_electrons"]


def test_tetrahedral_d5_high_spin_maximizes_unpaired_electrons() -> None:
    diagram = build_ligand_field_diagram(5, geometry="tetrahedral", spin_mode="high")

    assert diagram.metadata["unpaired_electrons"] == 5
