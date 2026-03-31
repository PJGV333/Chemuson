"""Presets de inserción rápida para diagramas electrónicos semánticos."""

from __future__ import annotations

from dataclasses import dataclass
from typing import Literal

from chemuson.gui.diagram_models import SemanticDiagram
from chemuson.gui.energy_diagrams import (
    build_atomic_species_diagram,
    build_diatomic_mo_diagram,
    build_ligand_field_diagram,
)


PresetFamily = Literal["atomic", "molecular_orbital", "ligand_field"]


@dataclass(frozen=True)
class AtomicDiagramPreset:
    name: str
    display_name: str
    symbol: str
    charge: int = 0
    expanded_subshells: bool = True
    title: str | None = None
    use_known_exceptions: bool = True
    family: PresetFamily = "atomic"

    def build(self) -> SemanticDiagram:
        return build_atomic_species_diagram(
            symbol=self.symbol,
            charge=self.charge,
            expanded_subshells=self.expanded_subshells,
            title=self.title,
            use_known_exceptions=self.use_known_exceptions,
        )


@dataclass(frozen=True)
class MolecularOrbitalDiagramPreset:
    name: str
    display_name: str
    left_label: str
    right_label: str
    total_electrons: int
    ordering: Literal["light_2p", "heavy_2p"] = "heavy_2p"
    include_core_1s: bool = False
    title: str | None = None
    family: PresetFamily = "molecular_orbital"

    def build(self) -> SemanticDiagram:
        return build_diatomic_mo_diagram(
            left_label=self.left_label,
            right_label=self.right_label,
            total_electrons=self.total_electrons,
            ordering=self.ordering,
            include_core_1s=self.include_core_1s,
            title=self.title,
        )


@dataclass(frozen=True)
class LigandFieldDiagramPreset:
    name: str
    display_name: str
    d_electrons: int
    geometry: Literal["octahedral", "tetrahedral", "square_planar"] = "octahedral"
    spin_mode: Literal["high", "low"] = "high"
    title: str | None = None
    family: PresetFamily = "ligand_field"

    def build(self) -> SemanticDiagram:
        return build_ligand_field_diagram(
            d_electrons=self.d_electrons,
            geometry=self.geometry,
            spin_mode=self.spin_mode,
            title=self.title,
        )


ATOMIC_PRESETS: tuple[AtomicDiagramPreset, ...] = (
    AtomicDiagramPreset("atomic:H", "H", "H"),
    AtomicDiagramPreset("atomic:C", "C", "C"),
    AtomicDiagramPreset("atomic:N", "N", "N"),
    AtomicDiagramPreset("atomic:O", "O", "O"),
    AtomicDiagramPreset("atomic:Na", "Na", "Na"),
    AtomicDiagramPreset("atomic:Cl", "Cl", "Cl"),
    AtomicDiagramPreset("atomic:Fe", "Fe", "Fe"),
    AtomicDiagramPreset("atomic:Co", "Co", "Co"),
    AtomicDiagramPreset("atomic:Ni", "Ni", "Ni"),
    AtomicDiagramPreset("atomic:Cr", "Cr", "Cr"),
    AtomicDiagramPreset("atomic:Cu", "Cu", "Cu"),
)

MOLECULAR_ORBITAL_PRESETS: tuple[MolecularOrbitalDiagramPreset, ...] = (
    MolecularOrbitalDiagramPreset("mo:B2", "B2", "B", "B", 6, "light_2p", False, "B2 MO Diagram"),
    MolecularOrbitalDiagramPreset("mo:C2", "C2", "C", "C", 8, "light_2p", False, "C2 MO Diagram"),
    MolecularOrbitalDiagramPreset("mo:N2", "N2", "N", "N", 10, "light_2p", False, "N2 MO Diagram"),
    MolecularOrbitalDiagramPreset("mo:O2", "O2", "O", "O", 12, "heavy_2p", False, "O2 MO Diagram"),
    MolecularOrbitalDiagramPreset("mo:F2", "F2", "F", "F", 14, "heavy_2p", False, "F2 MO Diagram"),
    MolecularOrbitalDiagramPreset("mo:O2-", "O2-", "O", "O", 13, "heavy_2p", False, "O2- MO Diagram"),
    MolecularOrbitalDiagramPreset("mo:O2^2-", "O2^2-", "O", "O", 14, "heavy_2p", False, "O2^2- MO Diagram"),
)

LIGAND_FIELD_PRESETS: tuple[LigandFieldDiagramPreset, ...] = (
    LigandFieldDiagramPreset("ligand_field:d4_octahedral_high", "d4 octahedral high spin", 4, "octahedral", "high"),
    LigandFieldDiagramPreset("ligand_field:d4_octahedral_low", "d4 octahedral low spin", 4, "octahedral", "low"),
    LigandFieldDiagramPreset("ligand_field:d5_octahedral_high", "d5 octahedral high spin", 5, "octahedral", "high"),
    LigandFieldDiagramPreset("ligand_field:d6_octahedral_high", "d6 octahedral high spin", 6, "octahedral", "high"),
    LigandFieldDiagramPreset("ligand_field:d6_octahedral_low", "d6 octahedral low spin", 6, "octahedral", "low"),
    LigandFieldDiagramPreset("ligand_field:d8_square_planar", "d8 square planar", 8, "square_planar", "low"),
    LigandFieldDiagramPreset("ligand_field:d5_tetrahedral_high", "d5 tetrahedral high spin", 5, "tetrahedral", "high"),
)


DIAGRAM_PRESETS = {
    **{preset.name: preset for preset in ATOMIC_PRESETS},
    **{preset.name: preset for preset in MOLECULAR_ORBITAL_PRESETS},
    **{preset.name: preset for preset in LIGAND_FIELD_PRESETS},
}


def get_diagram_preset(name: str):
    """Devuelve un descriptor de preset por nombre estable."""
    return DIAGRAM_PRESETS.get(str(name or "").strip())


def list_atomic_presets() -> list[AtomicDiagramPreset]:
    return list(ATOMIC_PRESETS)


def list_molecular_orbital_presets() -> list[MolecularOrbitalDiagramPreset]:
    return list(MOLECULAR_ORBITAL_PRESETS)


def list_ligand_field_presets() -> list[LigandFieldDiagramPreset]:
    return list(LIGAND_FIELD_PRESETS)


def build_semantic_diagram_from_preset(name: str) -> SemanticDiagram:
    """Construye un diagrama semántico desde un nombre de preset."""
    preset = get_diagram_preset(name)
    if preset is None:
        raise KeyError(f"Unknown semantic diagram preset: {name}")
    return preset.build()
