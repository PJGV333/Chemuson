"""Pruebas para presets semánticos de diagramas electrónicos."""

from __future__ import annotations



from chemuson.gui.diagram_presets import (
    build_semantic_diagram_from_preset,
    list_atomic_presets,
    list_ligand_field_presets,
    list_molecular_orbital_presets,
)


def test_preset_registry_contains_expected_builtins() -> None:
    atomic_names = {preset.display_name for preset in list_atomic_presets()}
    mo_names = {preset.display_name for preset in list_molecular_orbital_presets()}
    ligand_names = {preset.display_name for preset in list_ligand_field_presets()}

    assert {"H", "C", "N", "O", "Na", "Cl", "Fe", "Co", "Ni", "Cr", "Cu"} <= atomic_names
    assert {"B2", "C2", "N2", "O2", "F2", "O2-", "O2^2-"} <= mo_names
    assert {
        "d4 octahedral high spin",
        "d4 octahedral low spin",
        "d5 octahedral high spin",
        "d6 octahedral high spin",
        "d6 octahedral low spin",
        "d8 square planar",
        "d5 tetrahedral high spin",
    } <= ligand_names


def test_o2_preset_yields_two_unpaired_electrons() -> None:
    diagram = build_semantic_diagram_from_preset("mo:O2")

    assert diagram.metadata["unpaired_electrons"] == 2


def test_n2_preset_yields_bond_order_three() -> None:
    diagram = build_semantic_diagram_from_preset("mo:N2")

    assert diagram.metadata["bond_order"] == 3.0


def test_d6_octahedral_low_spin_has_fewer_unpaired_than_high_spin() -> None:
    high_spin = build_semantic_diagram_from_preset("ligand_field:d6_octahedral_high")
    low_spin = build_semantic_diagram_from_preset("ligand_field:d6_octahedral_low")

    assert low_spin.metadata["unpaired_electrons"] < high_spin.metadata["unpaired_electrons"]
