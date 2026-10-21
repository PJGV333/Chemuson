"""Pruebas para presets semánticos de diagramas electrónicos."""

from __future__ import annotations

import hashlib
import json

from chemuson.gui.diagram_presets import (
    build_semantic_diagram_from_preset,
    diagram_preset_menu_label,
    list_atomic_presets,
    list_ligand_field_presets,
    list_molecular_orbital_presets,
)


def test_preset_registry_contains_expected_builtins() -> None:
    atomic_presets = list_atomic_presets()
    mo_presets = list_molecular_orbital_presets()
    ligand_presets = list_ligand_field_presets()
    atomic_names = {preset.display_name for preset in atomic_presets}
    mo_names = {preset.display_name for preset in mo_presets}
    ligand_names = {preset.display_name for preset in ligand_presets}

    assert atomic_names == {"H", "C", "N", "O", "Na", "Cl", "Fe", "Co", "Ni", "Cr", "Cu"}
    assert mo_names == {"B2", "C2", "N2", "O2", "F2", "O2-", "O2^2-"}
    assert ligand_names == {
        "d4 octahedral high spin",
        "d4 octahedral low spin",
        "d5 octahedral high spin",
        "d6 octahedral high spin",
        "d6 octahedral low spin",
        "d8 square planar",
        "d5 tetrahedral high spin",
    }
    assert {preset.name for preset in atomic_presets} == {
        "atomic:H", "atomic:C", "atomic:N", "atomic:O", "atomic:Na", "atomic:Cl",
        "atomic:Fe", "atomic:Co", "atomic:Ni", "atomic:Cr", "atomic:Cu",
    }
    assert {preset.name for preset in mo_presets} == {
        "mo:B2", "mo:C2", "mo:N2", "mo:O2", "mo:F2", "mo:O2-", "mo:O2^2-",
    }
    assert {preset.name for preset in ligand_presets} == {
        "ligand_field:d4_octahedral_high",
        "ligand_field:d4_octahedral_low",
        "ligand_field:d5_octahedral_high",
        "ligand_field:d6_octahedral_high",
        "ligand_field:d6_octahedral_low",
        "ligand_field:d8_square_planar",
        "ligand_field:d5_tetrahedral_high",
    }

    expected_menu_labels = {
        "mo:B2": "B₂",
        "mo:C2": "C₂",
        "mo:N2": "N₂",
        "mo:O2": "O₂",
        "mo:F2": "F₂",
        "mo:O2-": "O₂⁻",
        "mo:O2^2-": "O₂²⁻",
        "ligand_field:d4_octahedral_high": "d⁴ octahedral high spin",
        "ligand_field:d4_octahedral_low": "d⁴ octahedral low spin",
        "ligand_field:d5_octahedral_high": "d⁵ octahedral high spin",
        "ligand_field:d6_octahedral_high": "d⁶ octahedral high spin",
        "ligand_field:d6_octahedral_low": "d⁶ octahedral low spin",
        "ligand_field:d8_square_planar": "d⁸ square planar",
        "ligand_field:d5_tetrahedral_high": "d⁵ tetrahedral high spin",
    }
    assert {
        name: diagram_preset_menu_label(name) for name in expected_menu_labels
    } == expected_menu_labels
    assert all(
        diagram_preset_menu_label(preset.name) == preset.display_name
        for preset in atomic_presets
    )

    # Independent baseline hashes captured from the beta base before the UI-only format.
    expected_diagram_fingerprints = {
        "mo:B2": "82d690663ba625639016ee4e87ec82e9c6854ef75a5998abf65ba700205d7d7d",
        "mo:C2": "9f2d77ff652a11e98137d9eadb511c91ca8e85d722a3bfe8d3b6d3325a51b17b",
        "mo:N2": "b922eeda182c547e18b6280fef85e829558969032b361cf8383ed2a24cc792ee",
        "mo:O2": "de792ae73e0cb316c4e265b420786bacebdbee5f29caf7c2c19a34907c5ec94a",
        "mo:F2": "d4df73bde47af683ceffe4e16c91abb1c6e12186c69d9a4aea1ea9f3496662e5",
        "mo:O2-": "c8b1a57cb32640aa60c0b96a309abc2a256a494b115617f950abc07a35c0c55d",
        "mo:O2^2-": "402a47a1e4889cf39420835bc658a7acd467358f5a6efd618a8fc4523784b04b",
        "ligand_field:d4_octahedral_high": "bdd3037f51633e83af99956c76f5f3ea16a8b905a085b42e08200edab1497cdf",
        "ligand_field:d4_octahedral_low": "b31864815321f93bd03017b4a89ed9973a50ae48f815179086e7f02c2aa9d43f",
        "ligand_field:d5_octahedral_high": "a7e99ee9884d50f011538557da4fcfb5d50735cd4b1718dc6d45a3d149c1b985",
        "ligand_field:d6_octahedral_high": "2af2cc0e51ca4d593d02f606e77a6d3d84452595f302a05802e06849c3a319e9",
        "ligand_field:d6_octahedral_low": "393d625dfd35fa71543d13dd30fc4710bfe23faaa5ebda5756daa2a9fd012693",
        "ligand_field:d8_square_planar": "b083cd09643ad9711e9ff64576ca6f52d8c42b6c19857927de3eea6c3dbccb83",
        "ligand_field:d5_tetrahedral_high": "280a3b5a7bca2873ed688472289959e7d69e38344c0d8ab18a67200ce07753d8",
    }
    actual_diagram_fingerprints = {}
    for name in expected_diagram_fingerprints:
        diagram = build_semantic_diagram_from_preset(name).to_json_dict()
        canonical = json.dumps(
            diagram,
            sort_keys=True,
            separators=(",", ":"),
            ensure_ascii=False,
        ).encode("utf-8")
        actual_diagram_fingerprints[name] = hashlib.sha256(canonical).hexdigest()
    assert actual_diagram_fingerprints == expected_diagram_fingerprints


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
