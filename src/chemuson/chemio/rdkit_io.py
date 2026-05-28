"""Conversión entre el grafo interno y formatos químicos externos.

Este módulo integra RDKit (si está disponible) para exportar/importar
SMILES, MOL y SVG. Incluye rutas de respaldo para exportar sin RDKit.
"""

from __future__ import annotations

import os
import re
from dataclasses import dataclass
from typing import Dict, Iterable, Optional, Tuple

from chemuson.core.model import ATOMIC_NUMBERS, BondStyle, BondStereo, MolGraph, bond_is_structural

Chem = None
AllChem = None
rdMolDraw2D = None
rdMolStandardize = None
_RDKIT_IMPORT_ATTEMPTED = False


def _direct_rdkit_enabled() -> bool:
    """Indica si se permite importar RDKit en el proceso principal."""
    return str(os.environ.get("CHEMUSON_ENABLE_DIRECT_RDKIT", "")).strip().lower() in {
        "1",
        "true",
        "yes",
        "on",
    }


def _load_rdkit_modules() -> bool:
    """Carga RDKit de forma perezosa solo cuando se habilita explícitamente.

    RDKit puede abortar el proceso si sus extensiones nativas no coinciden con
    el intérprete activo. La ruta normal de Chemuson usa workers aislados; esta
    carga directa queda reservada para entornos controlados.
    """
    global Chem, AllChem, rdMolDraw2D, rdMolStandardize, _RDKIT_IMPORT_ATTEMPTED
    if Chem is not None and AllChem is not None and rdMolDraw2D is not None:
        return True
    if _RDKIT_IMPORT_ATTEMPTED or not _direct_rdkit_enabled():
        return False
    _RDKIT_IMPORT_ATTEMPTED = True
    try:
        from rdkit import Chem as _Chem
        from rdkit.Chem import AllChem as _AllChem
        from rdkit.Chem.Draw import rdMolDraw2D as _rdMolDraw2D

        try:
            from rdkit.Chem.MolStandardize import rdMolStandardize as _rdMolStandardize
        except Exception:  # pragma: no cover - optional dependency at runtime
            _rdMolStandardize = None
    except Exception:  # pragma: no cover - optional dependency at runtime
        Chem = None
        AllChem = None
        rdMolDraw2D = None
        rdMolStandardize = None
        return False
    Chem = _Chem
    AllChem = _AllChem
    rdMolDraw2D = _rdMolDraw2D
    rdMolStandardize = _rdMolStandardize
    return True


def _rdkit_available() -> bool:
    """Indica si RDKit está disponible en tiempo de ejecución."""
    return _load_rdkit_modules()


def _bond_priority_for_export(bond) -> int:
    """Prioriza el enlace más informativo cuando hay duplicados por par de átomos."""
    if getattr(bond, "is_aromatic", False):
        return 50
    style = getattr(bond, "style", BondStyle.PLAIN)
    style_bonus = 5 if style == BondStyle.COORDINATION else 0
    display_order = getattr(bond, "display_order", None)
    display_bonus = int(display_order or 0)
    return int(getattr(bond, "order", 1) or 1) * 10 + style_bonus + display_bonus


def _unique_bonds_for_export(molgraph: MolGraph):
    """Devuelve enlaces únicos por par de átomos para exportación robusta."""
    unique: Dict[tuple[int, int], object] = {}
    for bond in sorted(molgraph.bonds.values(), key=lambda item: item.id):
        if not bond_is_structural(bond):
            continue
        pair = (min(int(bond.a1_id), int(bond.a2_id)), max(int(bond.a1_id), int(bond.a2_id)))
        existing = unique.get(pair)
        if existing is None or _bond_priority_for_export(bond) > _bond_priority_for_export(existing):
            unique[pair] = bond
    return list(unique.values())


def _graph_has_duplicate_bond_pairs(molgraph: MolGraph) -> bool:
    """Indica si el grafo contiene dos enlaces para el mismo par de átomos."""
    seen: set[tuple[int, int]] = set()
    for bond in molgraph.bonds.values():
        if not bond_is_structural(bond):
            continue
        pair = (min(int(bond.a1_id), int(bond.a2_id)), max(int(bond.a1_id), int(bond.a2_id)))
        if pair in seen:
            return True
        seen.add(pair)
    return False


def _require_rdkit():
    """Lanza un error si RDKit no está instalado."""
    if not _rdkit_available():
        raise RuntimeError("RDKit no disponible")


def _atomic_number_for_symbol(symbol: str) -> int:
    """Resuelve número atómico usando tabla interna (sin invocar RDKit)."""
    text = str(symbol or "").strip()
    if not text:
        return 0
    if text in ATOMIC_NUMBERS:
        return int(ATOMIC_NUMBERS[text])
    if len(text) == 1:
        normalized = text.upper()
    else:
        normalized = text[0].upper() + text[1:].lower()
    return int(ATOMIC_NUMBERS.get(normalized, 0) or 0)


_SIMPLE_GROUP_EXPORT_LABELS: dict[tuple[str, int], str] = {
    ("N", 2): "NH2",
    ("N", 1): "NH",
    ("O", 1): "OH",
    ("S", 1): "SH",
}

_RDKIT_QUERY_SYMBOLS = {"*", "A", "Q", "L", "LP", "Du", "R", "R#"}


@dataclass(frozen=True)
class _ExportAtomSpec:
    element: str
    charge: int = 0
    explicit_h: Optional[int] = None
    group_h_cap: Optional[int] = None
    no_implicit: bool = False


@dataclass(frozen=True)
class _ExportBondSpec:
    a1: int
    a2: int
    order: int = 1
    is_aromatic: bool = False


@dataclass(frozen=True)
class _ExportAbbreviationSpec:
    atoms: tuple[_ExportAtomSpec, ...]
    bonds: tuple[_ExportBondSpec, ...]
    anchor_index: int = 0


def _atom_spec(element: str, charge: int = 0) -> _ExportAtomSpec:
    """Atajo legible para especificaciones de abreviaturas exportables."""
    return _ExportAtomSpec(element=element, charge=int(charge))


_EXPORT_ABBREVIATIONS: dict[str, _ExportAbbreviationSpec] = {
    "Me": _ExportAbbreviationSpec((_atom_spec("C"),), ()),
    "CH3": _ExportAbbreviationSpec((_atom_spec("C"),), ()),
    "Et": _ExportAbbreviationSpec(
        (_atom_spec("C"), _atom_spec("C")),
        (_ExportBondSpec(0, 1),),
    ),
    "C2H5": _ExportAbbreviationSpec(
        (_atom_spec("C"), _atom_spec("C")),
        (_ExportBondSpec(0, 1),),
    ),
    "iPr": _ExportAbbreviationSpec(
        (_atom_spec("C"), _atom_spec("C"), _atom_spec("C")),
        (_ExportBondSpec(0, 1), _ExportBondSpec(0, 2)),
    ),
    "tBu": _ExportAbbreviationSpec(
        (_atom_spec("C"), _atom_spec("C"), _atom_spec("C"), _atom_spec("C")),
        (_ExportBondSpec(0, 1), _ExportBondSpec(0, 2), _ExportBondSpec(0, 3)),
    ),
    "OMe": _ExportAbbreviationSpec(
        (_atom_spec("O"), _atom_spec("C")),
        (_ExportBondSpec(0, 1),),
    ),
    "MeO": _ExportAbbreviationSpec(
        (_atom_spec("O"), _atom_spec("C")),
        (_ExportBondSpec(0, 1),),
    ),
    "OCH3": _ExportAbbreviationSpec(
        (_atom_spec("O"), _atom_spec("C")),
        (_ExportBondSpec(0, 1),),
    ),
    "CH3O": _ExportAbbreviationSpec(
        (_atom_spec("O"), _atom_spec("C")),
        (_ExportBondSpec(0, 1),),
    ),
    "OEt": _ExportAbbreviationSpec(
        (_atom_spec("O"), _atom_spec("C"), _atom_spec("C")),
        (_ExportBondSpec(0, 1), _ExportBondSpec(1, 2)),
    ),
    "EtO": _ExportAbbreviationSpec(
        (_atom_spec("O"), _atom_spec("C"), _atom_spec("C")),
        (_ExportBondSpec(0, 1), _ExportBondSpec(1, 2)),
    ),
    "Ph": _ExportAbbreviationSpec(
        (
            _atom_spec("C"),
            _atom_spec("C"),
            _atom_spec("C"),
            _atom_spec("C"),
            _atom_spec("C"),
            _atom_spec("C"),
        ),
        (
            _ExportBondSpec(0, 1, is_aromatic=True),
            _ExportBondSpec(1, 2, is_aromatic=True),
            _ExportBondSpec(2, 3, is_aromatic=True),
            _ExportBondSpec(3, 4, is_aromatic=True),
            _ExportBondSpec(4, 5, is_aromatic=True),
            _ExportBondSpec(5, 0, is_aromatic=True),
        ),
    ),
    "CN": _ExportAbbreviationSpec(
        (_atom_spec("C"), _atom_spec("N")),
        (_ExportBondSpec(0, 1, order=3),),
    ),
    "CHO": _ExportAbbreviationSpec(
        (_atom_spec("C"), _atom_spec("O")),
        (_ExportBondSpec(0, 1, order=2),),
    ),
    "COOH": _ExportAbbreviationSpec(
        (_atom_spec("C"), _atom_spec("O"), _atom_spec("O")),
        (_ExportBondSpec(0, 1, order=2), _ExportBondSpec(0, 2)),
    ),
    "CO2H": _ExportAbbreviationSpec(
        (_atom_spec("C"), _atom_spec("O"), _atom_spec("O")),
        (_ExportBondSpec(0, 1, order=2), _ExportBondSpec(0, 2)),
    ),
    "NO2": _ExportAbbreviationSpec(
        (_atom_spec("N", charge=1), _atom_spec("O"), _atom_spec("O", charge=-1)),
        (_ExportBondSpec(0, 1, order=2), _ExportBondSpec(0, 2)),
    ),
    "SO3H": _ExportAbbreviationSpec(
        (_atom_spec("S"), _atom_spec("O"), _atom_spec("O"), _atom_spec("O")),
        (
            _ExportBondSpec(0, 1, order=2),
            _ExportBondSpec(0, 2, order=2),
            _ExportBondSpec(0, 3),
        ),
    ),
    "Ac": _ExportAbbreviationSpec(
        (_atom_spec("C"), _atom_spec("O"), _atom_spec("C")),
        (_ExportBondSpec(0, 1, order=2), _ExportBondSpec(0, 2)),
    ),
    "Boc": _ExportAbbreviationSpec(
        (
            _atom_spec("C"),
            _atom_spec("O"),
            _atom_spec("O"),
            _atom_spec("C"),
            _atom_spec("C"),
            _atom_spec("C"),
            _atom_spec("C"),
        ),
        (
            _ExportBondSpec(0, 1, order=2),
            _ExportBondSpec(0, 2),
            _ExportBondSpec(2, 3),
            _ExportBondSpec(3, 4),
            _ExportBondSpec(3, 5),
            _ExportBondSpec(3, 6),
        ),
    ),
    "Ts": _ExportAbbreviationSpec(
        (
            _atom_spec("S"),
            _atom_spec("O"),
            _atom_spec("O"),
            _atom_spec("C"),
            _atom_spec("C"),
            _atom_spec("C"),
            _atom_spec("C"),
            _atom_spec("C"),
            _atom_spec("C"),
            _atom_spec("C"),
        ),
        (
            _ExportBondSpec(0, 1, order=2),
            _ExportBondSpec(0, 2, order=2),
            _ExportBondSpec(0, 3),
            _ExportBondSpec(3, 4, is_aromatic=True),
            _ExportBondSpec(4, 5, is_aromatic=True),
            _ExportBondSpec(5, 6, is_aromatic=True),
            _ExportBondSpec(6, 7, is_aromatic=True),
            _ExportBondSpec(7, 8, is_aromatic=True),
            _ExportBondSpec(8, 3, is_aromatic=True),
            _ExportBondSpec(6, 9),
        ),
    ),
}


_EXPORT_ABBREVIATION_ALIASES: dict[str, str] = {
    label.lower(): label for label in _EXPORT_ABBREVIATIONS
}


def _structural_degree_map(molgraph: MolGraph) -> dict[int, int]:
    """Cuenta enlaces estructurales por átomo."""
    degree: dict[int, int] = {int(atom_id): 0 for atom_id in molgraph.atoms}
    for bond in molgraph.bonds.values():
        if not bond_is_structural(bond):
            continue
        degree[int(bond.a1_id)] = degree.get(int(bond.a1_id), 0) + 1
        degree[int(bond.a2_id)] = degree.get(int(bond.a2_id), 0) + 1
    return degree


def _median_structural_bond_length(molgraph: MolGraph) -> float:
    """Estima una longitud de enlace robusta a partir del grafo."""
    lengths: list[float] = []
    for bond in molgraph.bonds.values():
        if not bond_is_structural(bond):
            continue
        a1 = molgraph.atoms.get(int(bond.a1_id))
        a2 = molgraph.atoms.get(int(bond.a2_id))
        if a1 is None or a2 is None:
            continue
        dx = float(a2.x) - float(a1.x)
        dy = float(a2.y) - float(a1.y)
        length = (dx * dx + dy * dy) ** 0.5
        if length > 1e-6:
            lengths.append(length)
    if not lengths:
        return 40.0
    lengths.sort()
    mid = len(lengths) // 2
    if len(lengths) % 2:
        return float(lengths[mid])
    return (float(lengths[mid - 1]) + float(lengths[mid])) * 0.5


def _atoms_are_export_merge_compatible(a, b) -> bool:
    """Determina si dos átomos casi superpuestos pueden representar el mismo nodo."""
    if str(getattr(a, "element", "")) != str(getattr(b, "element", "")):
        return False
    if str(getattr(a, "element", "")) == "H":
        return False
    if int(getattr(a, "charge", 0) or 0) != int(getattr(b, "charge", 0) or 0):
        return False
    if getattr(a, "isotope", None) != getattr(b, "isotope", None):
        return False
    if int(getattr(a, "radical_electrons", 0) or 0) != int(
        getattr(b, "radical_electrons", 0) or 0
    ):
        return False
    if bool(getattr(a, "is_coordination_center", False)) or bool(
        getattr(b, "is_coordination_center", False)
    ):
        return False
    return True


def _atom_export_merge_priority(atom, degree: int) -> tuple[int, int, int, int, int]:
    """Prioriza qué átomo conservar al colapsar duplicados de exportación."""
    explicit_h = getattr(atom, "explicit_h", None)
    return (
        int(degree),
        1 if bool(getattr(atom, "is_explicit", False)) else 0,
        1 if explicit_h is not None else 0,
        1 if getattr(atom, "group_h_cap", None) is not None else 0,
        -int(getattr(atom, "id", 0)),
    )


def _find(parent: dict[int, int], atom_id: int) -> int:
    """Busca la raíz union-find con compresión de caminos."""
    root = int(atom_id)
    while parent[root] != root:
        root = parent[root]
    while parent[int(atom_id)] != int(atom_id):
        next_id = parent[int(atom_id)]
        parent[int(atom_id)] = root
        atom_id = next_id
    return root


def _merge_near_duplicate_atoms_for_export(molgraph: MolGraph) -> MolGraph:
    """Colapsa átomos casi superpuestos solo en una copia para SMILES.

    Esto cubre el caso de edicion donde un cierre de anillo se ve conectado
    porque dos átomos caen sobre el mismo punto, pero internamente son nodos
    distintos. El umbral es pequeño para no inferir enlaces por proximidad.
    """
    if len(molgraph.atoms) < 2:
        return molgraph
    median_length = _median_structural_bond_length(molgraph)
    threshold = min(6.0, max(1.5, median_length * 0.15))
    threshold_sq = threshold * threshold

    bonded_pairs = {
        (min(int(bond.a1_id), int(bond.a2_id)), max(int(bond.a1_id), int(bond.a2_id)))
        for bond in molgraph.bonds.values()
        if bond_is_structural(bond)
    }
    degree = _structural_degree_map(molgraph)
    parent: dict[int, int] = {int(atom_id): int(atom_id) for atom_id in molgraph.atoms}
    atoms = sorted(molgraph.atoms.values(), key=lambda item: int(item.id))

    for idx, atom_a in enumerate(atoms):
        for atom_b in atoms[idx + 1 :]:
            pair = (
                min(int(atom_a.id), int(atom_b.id)),
                max(int(atom_a.id), int(atom_b.id)),
            )
            if pair in bonded_pairs:
                continue
            if not _atoms_are_export_merge_compatible(atom_a, atom_b):
                continue
            dx = float(atom_a.x) - float(atom_b.x)
            dy = float(atom_a.y) - float(atom_b.y)
            if dx * dx + dy * dy > threshold_sq:
                continue
            root_a = _find(parent, int(atom_a.id))
            root_b = _find(parent, int(atom_b.id))
            if root_a == root_b:
                continue
            keep = root_a
            drop = root_b
            keep_atom = molgraph.atoms[keep]
            drop_atom = molgraph.atoms[drop]
            drop_priority = _atom_export_merge_priority(drop_atom, degree.get(drop, 0))
            keep_priority = _atom_export_merge_priority(keep_atom, degree.get(keep, 0))
            if drop_priority > keep_priority:
                keep, drop = drop, keep
            parent[drop] = keep

    groups: dict[int, list[int]] = {}
    for atom_id in molgraph.atoms:
        root = _find(parent, int(atom_id))
        groups.setdefault(root, []).append(int(atom_id))
    if all(len(ids) == 1 for ids in groups.values()):
        return molgraph

    merged = MolGraph()
    id_map: dict[int, int] = {}
    for root, _atom_ids in sorted(groups.items()):
        keep_atom = molgraph.atoms[root]
        new_id = _copy_atom_for_export(merged, keep_atom, atom_id=int(root))
        for atom_id in groups[root]:
            id_map[int(atom_id)] = new_id

    unique_bonds: dict[tuple[int, int], object] = {}
    for bond in _unique_bonds_for_export(molgraph):
        a1_id = id_map.get(int(bond.a1_id), int(bond.a1_id))
        a2_id = id_map.get(int(bond.a2_id), int(bond.a2_id))
        if a1_id == a2_id:
            continue
        pair = (min(a1_id, a2_id), max(a1_id, a2_id))
        existing = unique_bonds.get(pair)
        if existing is None or _bond_priority_for_export(bond) > _bond_priority_for_export(existing):
            unique_bonds[pair] = bond

    for bond in sorted(unique_bonds.values(), key=lambda item: int(getattr(item, "id", 0))):
        a1_id = id_map.get(int(bond.a1_id), int(bond.a1_id))
        a2_id = id_map.get(int(bond.a2_id), int(bond.a2_id))
        if a1_id != a2_id:
            _copy_bond_for_export(merged, bond, a1_id, a2_id)

    return merged


def _pseudoatom_label_for_export(atom) -> str | None:
    """Devuelve etiqueta pseudoatómica cuando conviene exportar un dummy atom."""
    element = str(getattr(atom, "element", "") or "").strip()
    if not element:
        return None
    group_h_cap = getattr(atom, "group_h_cap", None)
    if group_h_cap is not None:
        label = _SIMPLE_GROUP_EXPORT_LABELS.get((element, int(group_h_cap)))
        if label:
            return label
    if element in _RDKIT_QUERY_SYMBOLS:
        return None
    if _atomic_number_for_symbol(element) <= 0:
        return element
    return None


def _export_abbreviation_for_label(label: str) -> _ExportAbbreviationSpec | None:
    """Devuelve la expansión química de una abreviatura conocida."""
    text = str(label or "").strip()
    if not text:
        return None
    canonical = _EXPORT_ABBREVIATION_ALIASES.get(text.lower())
    if canonical is None:
        return None
    return _EXPORT_ABBREVIATIONS.get(canonical)


def _copy_atom_for_export(
    target: MolGraph,
    source_atom,
    *,
    atom_id: int | None = None,
    atom_spec: _ExportAtomSpec | None = None,
) -> int:
    """Copia un átomo al grafo de exportación, opcionalmente reemplazando su elemento."""
    element = source_atom.element if atom_spec is None else atom_spec.element
    charge = int(getattr(source_atom, "charge", 0) or 0)
    explicit_h = getattr(source_atom, "explicit_h", None)
    group_h_cap = getattr(source_atom, "group_h_cap", None)
    no_implicit = bool(getattr(source_atom, "no_implicit", False))
    if atom_spec is not None:
        charge += int(atom_spec.charge)
        explicit_h = atom_spec.explicit_h
        group_h_cap = atom_spec.group_h_cap
        no_implicit = bool(atom_spec.no_implicit)
    new_atom = target.add_atom(
        element,
        float(getattr(source_atom, "x", 0.0)),
        float(getattr(source_atom, "y", 0.0)),
        atom_id=atom_id,
        charge=charge,
        isotope=getattr(source_atom, "isotope", None) if atom_spec is None else None,
        radical_electrons=(
            int(getattr(source_atom, "radical_electrons", 0) or 0)
            if atom_spec is None
            else 0
        ),
        oxidation_state=(
            getattr(source_atom, "oxidation_state", None) if atom_spec is None else None
        ),
        stereo_cip=getattr(source_atom, "stereo_cip", None) if atom_spec is None else None,
        stereo_axial=getattr(source_atom, "stereo_axial", None) if atom_spec is None else None,
        stereo_helical=getattr(source_atom, "stereo_helical", None) if atom_spec is None else None,
        stereo_si_re=getattr(source_atom, "stereo_si_re", None) if atom_spec is None else None,
        explicit_h=explicit_h,
        group_h_cap=group_h_cap,
        mapping=getattr(source_atom, "mapping", None) if atom_spec is None else None,
        is_query=bool(getattr(source_atom, "is_query", False)) if atom_spec is None else False,
        is_explicit=(
            bool(getattr(source_atom, "is_explicit", False))
            if atom_spec is None
            else element != "C"
        ),
        no_implicit=no_implicit,
        label_scale=getattr(source_atom, "label_scale", None) if atom_spec is None else None,
        is_coordination_center=(
            bool(getattr(source_atom, "is_coordination_center", False))
            if atom_spec is None
            else False
        ),
        sphere_radius=getattr(source_atom, "sphere_radius", None) if atom_spec is None else None,
        sphere_color=getattr(source_atom, "sphere_color", None) if atom_spec is None else None,
        sphere_filled=bool(getattr(source_atom, "sphere_filled", True)),
        sphere_transparent=bool(getattr(source_atom, "sphere_transparent", False)),
        opacity=getattr(source_atom, "opacity", None) if atom_spec is None else None,
    )
    return int(new_atom.id)


def _copy_bond_for_export(target: MolGraph, source_bond, a1_id: int, a2_id: int) -> None:
    """Copia un enlace original al grafo de exportación."""
    target.add_bond(
        a1_id,
        a2_id,
        order=int(getattr(source_bond, "order", 1) or 1),
        bond_id=int(getattr(source_bond, "id")),
        style=getattr(source_bond, "style", BondStyle.PLAIN),
        stereo=getattr(source_bond, "stereo", BondStereo.NONE),
        stereo_ez=getattr(source_bond, "stereo_ez", None),
        stereo_axial=getattr(source_bond, "stereo_axial", None),
        stereo_endo_exo=getattr(source_bond, "stereo_endo_exo", None),
        stereo_helical=getattr(source_bond, "stereo_helical", None),
        is_aromatic=bool(getattr(source_bond, "is_aromatic", False)),
        display_order=getattr(source_bond, "display_order", None),
        is_query=bool(getattr(source_bond, "is_query", False)),
        ring_id=getattr(source_bond, "ring_id", None),
        length_px=getattr(source_bond, "length_px", None),
        stroke_px=getattr(source_bond, "stroke_px", None),
        color=getattr(source_bond, "color", None),
        donor_atom_id=getattr(source_bond, "donor_atom_id", None),
        flex_curve_1=getattr(source_bond, "flex_curve_1", None),
        flex_curve_2=getattr(source_bond, "flex_curve_2", None),
        pi_offset_sign=getattr(source_bond, "pi_offset_sign", None),
        opacity=getattr(source_bond, "opacity", None),
    )


def _expand_export_abbreviations(molgraph: MolGraph) -> MolGraph:
    """Expande abreviaturas comunes a átomos reales para exportar SMILES.

    SMILES no conserva aliases de dibujo como `OMe`: RDKit y ChemDraw los
    interpretan como átomos dummy (`*`). La expansión mantiene el grafo visible
    intacto y solo materializa grupos conocidos en una copia para exportación.
    """
    expansions: dict[int, _ExportAbbreviationSpec] = {}
    for atom in molgraph.atoms.values():
        spec = _export_abbreviation_for_label(getattr(atom, "element", ""))
        if spec is not None:
            expansions[int(atom.id)] = spec
    if not expansions:
        return molgraph

    expanded = MolGraph()
    id_map: dict[int, int] = {}
    pending: list[tuple[int, _ExportAbbreviationSpec, list[int]]] = []

    for atom in sorted(molgraph.atoms.values(), key=lambda item: item.id):
        spec = expansions.get(int(atom.id))
        if spec is None:
            id_map[int(atom.id)] = _copy_atom_for_export(expanded, atom, atom_id=int(atom.id))
            continue
        anchor_spec = spec.atoms[int(spec.anchor_index)]
        anchor_id = _copy_atom_for_export(
            expanded,
            atom,
            atom_id=int(atom.id),
            atom_spec=anchor_spec,
        )
        id_map[int(atom.id)] = anchor_id
        local_ids: list[int] = [0 for _ in spec.atoms]
        local_ids[int(spec.anchor_index)] = anchor_id
        pending.append((int(atom.id), spec, local_ids))

    for original_id, spec, local_ids in pending:
        source_atom = molgraph.atoms[original_id]
        for idx, atom_spec in enumerate(spec.atoms):
            if idx == int(spec.anchor_index):
                continue
            local_ids[idx] = _copy_atom_for_export(
                expanded,
                source_atom,
                atom_spec=atom_spec,
            )

    for bond in _unique_bonds_for_export(molgraph):
        a1_id = id_map.get(int(bond.a1_id))
        a2_id = id_map.get(int(bond.a2_id))
        if a1_id is None or a2_id is None or a1_id == a2_id:
            continue
        _copy_bond_for_export(expanded, bond, a1_id, a2_id)

    for _original_id, spec, local_ids in pending:
        for bond_spec in spec.bonds:
            a1_id = local_ids[int(bond_spec.a1)]
            a2_id = local_ids[int(bond_spec.a2)]
            if a1_id == a2_id:
                continue
            expanded.add_bond(
                a1_id,
                a2_id,
                order=int(bond_spec.order),
                is_aromatic=bool(bond_spec.is_aromatic),
            )

    return expanded


def expand_abbreviations_for_calculation(molgraph: MolGraph) -> MolGraph:
    """Expande superátomos conocidos en una copia para cálculo/exportación."""
    return _expand_export_abbreviations(molgraph)


@dataclass
class StrictValidationResult:
    """Resultado de validación/normalización estricta basada en RDKit."""

    is_valid: bool
    errors: list[str]
    normalized_smiles: Optional[str] = None
    normalized_graph: Optional[MolGraph] = None


def molgraph_to_rdkit_with_map(molgraph: MolGraph):
    """Convierte un `MolGraph` a RDKit `Mol` y mapea IDs internos.

    Args:
        molgraph: Grafo molecular de Chemuson.

    Returns:
        Tupla `(mol, id_map)` donde `mol` es un RDKit Mol y `id_map` mapea
        IDs de átomos del grafo a índices RDKit.

    Raises:
        RuntimeError: Si RDKit no está disponible.

    Side Effects:
        No modifica el grafo; construye objetos RDKit en memoria.
    """
    _require_rdkit()
    rw = Chem.RWMol()
    id_map: Dict[int, int] = {}

    def _bond_priority(bond_type) -> int:
        if bond_type == Chem.BondType.AROMATIC:
            return 5
        if bond_type == Chem.BondType.TRIPLE:
            return 4
        if bond_type == Chem.BondType.DOUBLE:
            return 3
        if bond_type == Chem.BondType.SINGLE:
            return 2
        if bond_type == Chem.BondType.DATIVE:
            return 1
        return 0

    for atom in sorted(molgraph.atoms.values(), key=lambda a: a.id):
        element = atom.element
        pseudo_label = _pseudoatom_label_for_export(atom)
        if pseudo_label is not None:
            rd_atom = Chem.Atom(0)
            rd_atom.SetProp("atomLabel", pseudo_label)
            rd_atom.SetProp("dummyLabel", pseudo_label)
        else:
            atomic_number = _atomic_number_for_symbol(element)
            rd_atom = Chem.Atom(int(atomic_number))
        rd_atom.SetFormalCharge(atom.charge)
        if atom.isotope is not None:
            rd_atom.SetIsotope(atom.isotope)
        radical_electrons = int(getattr(atom, "radical_electrons", 0) or 0)
        if radical_electrons > 0:
            rd_atom.SetNumRadicalElectrons(radical_electrons)
        if getattr(atom, "stereo_axial", None):
            rd_atom.SetProp("_ChemusonStereoAxial", str(atom.stereo_axial))
        if getattr(atom, "stereo_helical", None):
            rd_atom.SetProp("_ChemusonStereoHelical", str(atom.stereo_helical))
        if getattr(atom, "stereo_si_re", None):
            rd_atom.SetProp("_ChemusonStereoSiRe", str(atom.stereo_si_re))
        rd_idx = rw.AddAtom(rd_atom)
        id_map[atom.id] = rd_idx

    for bond in _unique_bonds_for_export(molgraph):
        # RDKit distingue enlaces aromáticos y órdenes discretos.
        begin_id = bond.a1_id
        end_id = bond.a2_id
        if bond.style == BondStyle.COORDINATION:
            bond_type = Chem.BondType.DATIVE
            donor = bond.donor_atom_id
            if donor == bond.a2_id:
                begin_id, end_id = bond.a2_id, bond.a1_id
            elif donor == bond.a1_id:
                begin_id, end_id = bond.a1_id, bond.a2_id
        elif bond.is_aromatic:
            rw.GetAtomWithIdx(id_map[bond.a1_id]).SetIsAromatic(True)
            rw.GetAtomWithIdx(id_map[bond.a2_id]).SetIsAromatic(True)
            bond_type = Chem.BondType.AROMATIC
        elif bond.order == 2:
            bond_type = Chem.BondType.DOUBLE
        elif bond.order == 3:
            bond_type = Chem.BondType.TRIPLE
        else:
            bond_type = Chem.BondType.SINGLE
        begin_idx = id_map[begin_id]
        end_idx = id_map[end_id]
        if begin_idx == end_idx:
            continue

        rd_bond = rw.GetBondBetweenAtoms(begin_idx, end_idx)
        if rd_bond is None:
            try:
                rw.AddBond(begin_idx, end_idx, bond_type)
                rd_bond = rw.GetBondBetweenAtoms(begin_idx, end_idx)
            except Exception:
                # Defensive skip
                continue
        elif _bond_priority(bond_type) > _bond_priority(rd_bond.GetBondType()):
            rd_bond.SetBondType(bond_type)
        if bond_type == Chem.BondType.AROMATIC:
            rw.GetAtomWithIdx(id_map[bond.a1_id]).SetIsAromatic(True)
            rw.GetAtomWithIdx(id_map[bond.a2_id]).SetIsAromatic(True)
        if rd_bond is not None:
            if getattr(bond, "stereo_axial", None):
                rd_bond.SetProp("_ChemusonStereoAxial", str(bond.stereo_axial))
            if getattr(bond, "stereo_endo_exo", None):
                rd_bond.SetProp("_ChemusonStereoEndoExo", str(bond.stereo_endo_exo))
            if getattr(bond, "stereo_helical", None):
                rd_bond.SetProp("_ChemusonStereoHelical", str(bond.stereo_helical))

    mol = rw.GetMol()
    # Preservar coordenadas 2D del editor.
    conf = Chem.Conformer(mol.GetNumAtoms())
    for atom_id, idx in id_map.items():
        atom = molgraph.atoms[atom_id]
        conf.SetAtomPosition(idx, (atom.x, atom.y, 0.0))
    mol.AddConformer(conf, assignId=True)
    return mol, id_map


def molgraph_to_rdkit(molgraph: MolGraph):
    """Convierte un `MolGraph` a RDKit `Mol` ignorando el mapeo de IDs.

    Args:
        molgraph: Grafo molecular de Chemuson.

    Returns:
        Mol de RDKit con conformación 2D.
    """
    mol, _ = molgraph_to_rdkit_with_map(molgraph)
    return mol


def _prepare_molgraph_for_smiles_export(molgraph: MolGraph) -> MolGraph:
    """Aplica normalizaciones seguras sobre una copia conceptual para SMILES."""
    merged = _merge_near_duplicate_atoms_for_export(molgraph)
    return _expand_export_abbreviations(merged)


def molgraph_to_smiles(molgraph: MolGraph) -> str:
    """Genera SMILES desde un `MolGraph`.

    Usa RDKit en un worker aislado si está disponible; si falla, utiliza un
    escritor interno. El worker evita que fallos nativos o bucles de RDKit
    congelen la interfaz gráfica.

    Args:
        molgraph: Grafo molecular de Chemuson.

    Returns:
        SMILES canónico o aproximado en el camino de respaldo.
    """
    export_graph = _prepare_molgraph_for_smiles_export(molgraph)
    if _graph_requires_rdkit_fallback(export_graph):
        return _molgraph_to_smiles_fallback(export_graph)
    if _rdkit_available():
        try:
            from chemuson.chemio.rdkit_safe import molgraph_to_smiles_isolated

            smiles, error = molgraph_to_smiles_isolated(export_graph)
            if not error and smiles:
                return smiles
        except Exception:
            # RDKit can reject some hypervalent depictions (e.g., interhalogens, noble gases).
            # Fall back to the internal writer so the editor can still export something useful.
            return _molgraph_to_smiles_fallback(export_graph)
    return _molgraph_to_smiles_fallback(export_graph)


def molgraph_to_smiles_isolated_or_error(
    molgraph: MolGraph,
    timeout_s: float = 8.0,
) -> str:
    """Genera SMILES con RDKit aislado y falla rápido si no hay resultado.

    Esta ruta está pensada para acciones de UI: evita caer al escritor interno,
    que puede ser costoso en estructuras cíclicas grandes.
    """
    from chemuson.chemio.rdkit_safe import molgraph_to_smiles_isolated

    export_graph = _prepare_molgraph_for_smiles_export(molgraph)
    smiles, error = molgraph_to_smiles_isolated(export_graph, timeout_s=timeout_s)
    if error or not smiles:
        raise RuntimeError(error or "empty_smiles")
    return smiles


def molgraph_to_molfile(molgraph: MolGraph) -> str:
    """Genera un bloque MOL (V2000) desde un `MolGraph`.

    Args:
        molgraph: Grafo molecular de Chemuson.

    Returns:
        Cadena MOL. Si RDKit falla, usa un exportador interno básico.
    """
    if _graph_requires_rdkit_fallback(molgraph):
        return _molgraph_to_molfile_fallback(molgraph)
    if _rdkit_available():
        try:
            mol = molgraph_to_rdkit(molgraph)
            return Chem.MolToMolBlock(mol)
        except Exception:
            return _molgraph_to_molfile_fallback(molgraph)
    return _molgraph_to_molfile_fallback(molgraph)


def molgraph_to_svg(molgraph: MolGraph, size: Tuple[int, int] = (300, 200)) -> str:
    """Renderiza un `MolGraph` como SVG con RDKit.

    Args:
        molgraph: Grafo molecular de Chemuson.
        size: Tamaño del lienzo (ancho, alto) en píxeles.

    Returns:
        SVG como cadena.
    """
    _require_rdkit()
    mol = molgraph_to_rdkit(molgraph)
    drawer = rdMolDraw2D.MolDraw2DSVG(size[0], size[1])
    rdMolDraw2D.PrepareAndDrawMolecule(drawer, mol)
    drawer.FinishDrawing()
    return drawer.GetDrawingText()


def _looks_like_counts_line(line: str) -> bool:
    """Heurística para identificar la línea de conteo de un MOL V2000/V3000."""
    text = str(line).rstrip()
    if not text:
        return False
    if "V2000" in text or "V3000" in text:
        return True
    padded = f"{text:<6}"[:6]
    left = padded[:3].strip()
    right = padded[3:6].strip()
    return left.isdigit() and right.isdigit()


def _parse_counts(line: str) -> Optional[tuple[int, int]]:
    """Extrae número de átomos/enlaces desde la línea de conteo."""
    text = str(line or "")
    if not text:
        return None
    try:
        return int(text[0:3]), int(text[3:6])
    except Exception:
        pass
    ints = re.findall(r"-?\d+", text)
    if len(ints) < 2:
        return None
    try:
        return int(ints[0]), int(ints[1])
    except Exception:
        return None


def _find_counts_line_index(lines: list[str]) -> Optional[int]:
    """Localiza la línea de conteo CTAB en las primeras líneas del bloque."""
    limit = min(len(lines), 8)
    for idx in range(limit):
        if _looks_like_counts_line(lines[idx]):
            if _parse_counts(lines[idx]) is not None:
                return idx
    return None


def _is_float_token(token: str) -> bool:
    """Verifica si un token representa un número real."""
    try:
        float(token)
        return True
    except Exception:
        return False


def _parse_atom_line(line: str) -> tuple[float, float, str]:
    """Parsea coordenadas XY y símbolo desde una línea de átomo MOL."""
    tokens = line.split()
    if (
        len(tokens) >= 4
        and _is_float_token(tokens[0])
        and _is_float_token(tokens[1])
        and _is_float_token(tokens[2])
    ):
        x = float(tokens[0])
        y = float(tokens[1])
        symbol = tokens[3].strip() or "C"
        return x, y, symbol
    try:
        x = float(line[0:10])
        y = float(line[10:20])
    except Exception as exc:
        raise ValueError(f"Línea de átomo inválida: {line!r}") from exc
    symbol = line[31:34].strip() if len(line) >= 34 else ""
    if not symbol:
        symbol = "C"
    return x, y, symbol


def _parse_bond_line(line: str) -> tuple[int, int, int]:
    """Parsea índices de átomos y tipo de enlace desde una línea MOL."""
    try:
        return int(line[0:3]), int(line[3:6]), int(line[6:9])
    except Exception:
        tokens = line.split()
        if len(tokens) < 3:
            raise ValueError(f"Línea de enlace inválida: {line!r}")
        return int(tokens[0]), int(tokens[1]), int(tokens[2])


def _fallback_bond_order_and_aromatic(bond_type: int) -> tuple[int, bool]:
    """Normaliza tipo de enlace MOL a orden/aromaticidad para parser interno."""
    order = 1
    aromatic = False
    if bond_type == 2:
        order = 2
    elif bond_type == 3:
        order = 3
    elif bond_type == 4:
        aromatic = True
        order = 1
    return order, aromatic


def _fallback_bond_rank(order: int, aromatic: bool) -> int:
    """Ranking para resolver duplicados por par de átomos en parser MOL."""
    if aromatic:
        return 100
    return int(order)


def _mol_symbol_requires_fallback(symbol: str) -> bool:
    """Indica si un símbolo no es elemento estándar y conviene parseo local."""
    clean = str(symbol or "").strip()
    if not clean:
        return True
    if clean in ATOMIC_NUMBERS:
        return False
    if clean in {"*", "A", "Q", "L", "LP", "Du", "R", "R#"}:
        return False
    return True


def _graph_requires_rdkit_fallback(molgraph: MolGraph) -> bool:
    """Determina si conviene evitar RDKit por símbolos no estándar."""
    for atom in molgraph.atoms.values():
        if bool(getattr(atom, "is_coordination_center", False)):
            return True
        if _mol_symbol_requires_fallback(atom.element):
            return True
    return _graph_has_duplicate_bond_pairs(molgraph)


def _should_use_molfile_fallback(molfile: str) -> bool:
    """Detecta MOL con pseudoátomos/estructura no compatible con RDKit parser."""
    text = str(molfile or "")
    if not text:
        return True
    lines = text.splitlines()
    counts_idx = _find_counts_line_index(lines)
    if counts_idx is None:
        return True
    counts = _parse_counts(lines[counts_idx])
    if counts is None:
        return True
    atom_count, _ = counts
    atom_start = counts_idx + 1
    if atom_count < 0 or atom_start + atom_count > len(lines):
        return True
    for i in range(atom_count):
        _, _, symbol = _parse_atom_line(lines[atom_start + i])
        if _mol_symbol_requires_fallback(symbol):
            return True
    bond_start = atom_start + atom_count
    if bond_start > len(lines):
        return True
    seen_pairs: set[tuple[int, int]] = set()
    idx = bond_start
    while idx < len(lines):
        line = lines[idx]
        if line.startswith("M  END"):
            break
        if line.startswith("M  "):
            idx += 1
            continue
        try:
            a1_idx, a2_idx, _bond_type = _parse_bond_line(line)
        except Exception:
            idx += 1
            continue
        if a1_idx == a2_idx:
            return True
        pair = (min(a1_idx, a2_idx), max(a1_idx, a2_idx))
        if pair in seen_pairs:
            # RDKit puede disparar "bond already exists" al parsear este caso.
            return True
        seen_pairs.add(pair)
        idx += 1
    # Molblocks con alias `A  n` pueden venir de pseudoátomos serializados
    # como `*` + etiqueta. RDKit tiende a devolver `*` y se pierde la etiqueta;
    # forzamos parser interno para preservar el alias original.
    for line in lines:
        if line.startswith("A  "):
            return True
    return False


def _molfile_to_molgraph_fallback(molfile: str) -> MolGraph:
    """Importa MOL mediante parser interno tolerante (sin RDKit)."""
    text = str(molfile or "")
    if not text.strip():
        raise ValueError("Mol inválido")
    lines = text.splitlines()
    counts_idx = _find_counts_line_index(lines)
    if counts_idx is None:
        raise ValueError("Mol inválido")
    counts = _parse_counts(lines[counts_idx])
    if counts is None:
        raise ValueError("Mol inválido")
    atom_count, bond_count = counts
    if atom_count < 0 or bond_count < 0:
        raise ValueError("Mol inválido")

    atom_start = counts_idx + 1
    bond_start = atom_start + atom_count
    if bond_start + bond_count > len(lines):
        raise ValueError("Mol inválido")

    graph = MolGraph()
    atom_ids: list[int] = []
    for offset in range(atom_count):
        x, y, symbol = _parse_atom_line(lines[atom_start + offset])
        atom = graph.add_atom(symbol, x, y, is_explicit=(symbol != "C"))
        atom_ids.append(atom.id)

    unique_bonds: dict[tuple[int, int], tuple[int, bool]] = {}
    for offset in range(bond_count):
        a1_idx, a2_idx, bond_type = _parse_bond_line(lines[bond_start + offset])
        if not (1 <= a1_idx <= len(atom_ids) and 1 <= a2_idx <= len(atom_ids)):
            continue
        if a1_idx == a2_idx:
            continue
        order, aromatic = _fallback_bond_order_and_aromatic(int(bond_type))
        pair = (min(a1_idx, a2_idx), max(a1_idx, a2_idx))
        previous = unique_bonds.get(pair)
        if previous is None or _fallback_bond_rank(order, aromatic) > _fallback_bond_rank(
            previous[0], previous[1]
        ):
            unique_bonds[pair] = (order, aromatic)

    for (a1_idx, a2_idx), (order, aromatic) in sorted(unique_bonds.items()):
        graph.add_bond(
            atom_ids[a1_idx - 1],
            atom_ids[a2_idx - 1],
            order=order,
            style=BondStyle.PLAIN,
            stereo=BondStereo.NONE,
            is_aromatic=aromatic,
        )

    aliases: dict[int, str] = {}
    idx = bond_start + bond_count
    while idx < len(lines):
        line = lines[idx]
        if line.startswith("M  END"):
            break
        if line.startswith("A  "):
            try:
                atom_idx = int(line[3:].strip())
            except Exception:
                parts = line.split()
                atom_idx = int(parts[1]) if len(parts) > 1 else -1
            if idx + 1 < len(lines):
                aliases[atom_idx] = lines[idx + 1].strip()
                idx += 2
                continue
        if line.startswith("M  CHG"):
            tokens = line.split()
            if len(tokens) >= 3:
                try:
                    count = int(tokens[2])
                except Exception:
                    count = 0
                cursor = 3
                for _ in range(count):
                    if cursor + 1 >= len(tokens):
                        break
                    try:
                        atom_idx = int(tokens[cursor])
                        charge = int(tokens[cursor + 1])
                    except Exception:
                        break
                    cursor += 2
                    if 1 <= atom_idx <= len(atom_ids):
                        graph.get_atom(atom_ids[atom_idx - 1]).charge = charge
        if line.startswith("M  ISO"):
            tokens = line.split()
            if len(tokens) >= 3:
                try:
                    count = int(tokens[2])
                except Exception:
                    count = 0
                cursor = 3
                for _ in range(count):
                    if cursor + 1 >= len(tokens):
                        break
                    try:
                        atom_idx = int(tokens[cursor])
                        mass = int(tokens[cursor + 1])
                    except Exception:
                        break
                    cursor += 2
                    if 1 <= atom_idx <= len(atom_ids):
                        graph.get_atom(atom_ids[atom_idx - 1]).isotope = int(mass)
        if line.startswith("M  RAD"):
            tokens = line.split()
            if len(tokens) >= 3:
                try:
                    count = int(tokens[2])
                except Exception:
                    count = 0
                cursor = 3
                for _ in range(count):
                    if cursor + 1 >= len(tokens):
                        break
                    try:
                        atom_idx = int(tokens[cursor])
                        rad_code = int(tokens[cursor + 1])
                    except Exception:
                        break
                    cursor += 2
                    if 1 <= atom_idx <= len(atom_ids):
                        electrons = {1: 2, 2: 1, 3: 2}.get(rad_code, 1)
                        graph.get_atom(atom_ids[atom_idx - 1]).radical_electrons = int(electrons)
        idx += 1

    for atom_idx, label in aliases.items():
        if 1 <= atom_idx <= len(atom_ids):
            atom = graph.get_atom(atom_ids[atom_idx - 1])
            atom.element = label
            atom.is_explicit = label != "C"

    _scale_to_default(graph)
    return graph


def normalize_molblock_header(molfile: str) -> str:
    """Normaliza cabecera CTAB para tolerar MOL con línea de comentario omitida.

    Algunos bloques guardados en bibliotecas pueden perder una línea de
    encabezado al aplicar `strip()`, desplazando la línea de conteo de CTAB.
    Esta función inserta una línea vacía cuando detecta ese patrón.
    """
    text = str(molfile or "")
    if not text:
        return ""
    normalized = text.replace("\r\n", "\n").replace("\r", "\n")
    lines = normalized.split("\n")
    if len(lines) >= 3 and _looks_like_counts_line(lines[2]):
        line4 = lines[3] if len(lines) > 3 else ""
        if not _looks_like_counts_line(line4):
            lines.insert(2, "")
    return "\n".join(lines)


def molfile_to_molgraph(molfile: str) -> MolGraph:
    """Importa un bloque MOL (V2000) a `MolGraph`.

    Args:
        molfile: Texto MOL.

    Returns:
        Grafo molecular equivalente.
    """
    normalized = normalize_molblock_header(molfile)
    try:
        # El parser interno preserva pseudoátomos, isótopos y enlaces tal como
        # vienen en el CTAB, y evita divergencias de sanitización de RDKit.
        return _molfile_to_molgraph_fallback(normalized)
    except Exception as fallback_error:
        if not _rdkit_available():
            raise fallback_error
        try:
            # Intentamos con RDKit sin sanitización para evitar crashes por aromaticidad.
            mol = Chem.MolFromMolBlock(normalized, sanitize=False)
            if mol is not None:
                try:
                    # Sanitización defensiva
                    Chem.SanitizeMol(mol, Chem.SanitizeFlags.SANITIZE_ALL ^ Chem.SanitizeFlags.SANITIZE_KEKULIZE)
                except Exception:
                    pass
                return rdkit_to_molgraph(mol)
        except Exception:
            pass
        raise fallback_error


def smiles_to_molgraph(smiles: str) -> MolGraph:
    """Importa un SMILES a `MolGraph` usando RDKit.

    Args:
        smiles: Cadena SMILES.

    Returns:
        Grafo molecular equivalente.
    """
    _require_rdkit()
    mol = Chem.MolFromSmiles(smiles)
    return rdkit_to_molgraph(mol)


def rdkit_to_molgraph(mol) -> MolGraph:
    """Convierte un RDKit Mol a `MolGraph`.

    Args:
        mol: Instancia de RDKit Mol.

    Returns:
        Grafo molecular con coordenadas 2D.

    Raises:
        ValueError: Si el Mol es `None`.
    """
    _require_rdkit()
    if mol is None:
        raise ValueError("Mol inválido")
    try:
        Chem.AssignStereochemistry(mol, cleanIt=True, force=True)
        Chem.FindPotentialStereoBonds(mol)
    except Exception:
        pass
    if mol.GetNumConformers() == 0:
        AllChem.Compute2DCoords(mol)
    conf = mol.GetConformer()

    graph = MolGraph()
    idx_map: Dict[int, int] = {}

    for atom in mol.GetAtoms():
        idx = atom.GetIdx()
        pos = conf.GetAtomPosition(idx)
        symbol = atom.GetSymbol()
        new_atom = graph.add_atom(
            symbol,
            pos.x,
            pos.y,
            is_explicit=symbol != "C",
        )
        new_atom.charge = atom.GetFormalCharge()
        if atom.GetIsotope():
            new_atom.isotope = atom.GetIsotope()
        try:
            new_atom.radical_electrons = int(atom.GetNumRadicalElectrons() or 0)
        except Exception:
            new_atom.radical_electrons = 0
        if atom.HasProp("_CIPCode"):
            try:
                new_atom.stereo_cip = atom.GetProp("_CIPCode")
            except Exception:
                pass
        for prop_name, attr_name in (
            ("_ChemusonStereoAxial", "stereo_axial"),
            ("_ChemusonStereoHelical", "stereo_helical"),
            ("_ChemusonStereoSiRe", "stereo_si_re"),
        ):
            if atom.HasProp(prop_name):
                try:
                    setattr(new_atom, attr_name, atom.GetProp(prop_name))
                except Exception:
                    pass
        idx_map[idx] = new_atom.id

    for bond in mol.GetBonds():
        order = 1
        style = BondStyle.PLAIN
        donor_atom_id = None
        bond_type = bond.GetBondType()
        if bond_type == Chem.BondType.DOUBLE:
            order = 2
        elif bond_type == Chem.BondType.TRIPLE:
            order = 3
        elif bond_type == Chem.BondType.DATIVE:
            style = BondStyle.COORDINATION
            donor_atom_id = idx_map[bond.GetBeginAtomIdx()]
        new_bond = graph.add_bond(
            idx_map[bond.GetBeginAtomIdx()],
            idx_map[bond.GetEndAtomIdx()],
            order,
            style=style,
            stereo=BondStereo.NONE,
            is_aromatic=bond.GetIsAromatic(),
            donor_atom_id=donor_atom_id,
        )
        stereo = bond.GetStereo()
        if stereo == Chem.BondStereo.STEREOE:
            new_bond.stereo_ez = "E"
        elif stereo == Chem.BondStereo.STEREOZ:
            new_bond.stereo_ez = "Z"
        for prop_name, attr_name in (
            ("_ChemusonStereoAxial", "stereo_axial"),
            ("_ChemusonStereoEndoExo", "stereo_endo_exo"),
            ("_ChemusonStereoHelical", "stereo_helical"),
        ):
            if bond.HasProp(prop_name):
                try:
                    setattr(new_bond, attr_name, bond.GetProp(prop_name))
                except Exception:
                    pass

    _scale_to_default(graph)
    return graph


def strict_validate_and_normalize(molgraph: MolGraph) -> StrictValidationResult:
    """Valida y normaliza de forma estricta un `MolGraph` usando RDKit.

    La validación estricta ejecuta sanitización de RDKit y, cuando está
    disponible, aplica un pipeline de normalización (`Normalize`/`Reionize`).
    Está orientada a chequeos avanzados de especies iónicas y complejos.
    """
    _require_rdkit()
    mol, _ = molgraph_to_rdkit_with_map(molgraph)

    try:
        Chem.SanitizeMol(mol)
    except Exception as exc:
        errors = [f"SanitizeMol: {exc}"]
        try:
            problems = Chem.DetectChemistryProblems(mol)
        except Exception:
            problems = []
        for problem in problems:
            message = ""
            if hasattr(problem, "Message"):
                try:
                    message = problem.Message()
                except Exception:
                    message = str(problem)
            else:
                message = str(problem)
            errors.append(message)
        return StrictValidationResult(is_valid=False, errors=errors)

    normalized = Chem.Mol(mol)
    if rdMolStandardize is not None:
        normalized = rdMolStandardize.Cleanup(normalized)
        normalized = rdMolStandardize.Normalize(normalized)
        normalized = rdMolStandardize.Reionize(normalized)
        Chem.SanitizeMol(normalized)

    smiles = Chem.MolToSmiles(normalized, canonical=True)
    normalized_graph = rdkit_to_molgraph(Chem.Mol(normalized))
    return StrictValidationResult(
        is_valid=True,
        errors=[],
        normalized_smiles=smiles,
        normalized_graph=normalized_graph,
    )


def kekulize_display_orders(
    molgraph: MolGraph, seed_atoms: Optional[Iterable[int]] = None
) -> Optional[Dict[int, int]]:
    """Calcula órdenes de dibujo para enlaces aromáticos.

    Usa la kekulización de RDKit y, opcionalmente, restringe el cálculo
    a los átomos conectados a una lista de semillas.

    Args:
        molgraph: Grafo molecular de Chemuson.
        seed_atoms: IDs de átomos desde los cuales expandir el cálculo.

    Returns:
        Diccionario de `bond_id -> display_order` o `None` si falla.
    """
    if not _rdkit_available():
        return None
    # Con pseudoátomos (p. ej. OH/CH2OH) evitamos la ruta RDKit para no
    # disparar errores nativos de tabla periódica.
    if any(_mol_symbol_requires_fallback(atom.element) for atom in molgraph.atoms.values()):
        return None
    try:
        mol, id_map = molgraph_to_rdkit_with_map(molgraph)
        Chem.Kekulize(mol, clearAromaticFlags=True)
    except Exception:
        return None

    aromatic_atoms: Optional[set[int]] = None
    if seed_atoms is not None:
        # Expandimos las semillas a todo el componente aromático conectado.
        seeds = set(seed_atoms)
        if seeds:
            adjacency: dict[int, list[int]] = {}
            for bond in molgraph.bonds.values():
                if bond.is_aromatic:
                    adjacency.setdefault(bond.a1_id, []).append(bond.a2_id)
                    adjacency.setdefault(bond.a2_id, []).append(bond.a1_id)
            aromatic_atoms = set()
            stack = [atom_id for atom_id in seeds if atom_id in adjacency]
            while stack:
                node = stack.pop()
                if node in aromatic_atoms:
                    continue
                aromatic_atoms.add(node)
                for neighbor in adjacency.get(node, []):
                    if neighbor not in aromatic_atoms:
                        stack.append(neighbor)

    display_orders: Dict[int, int] = {}
    for bond in molgraph.bonds.values():
        if not bond.is_aromatic:
            continue
        if aromatic_atoms is not None and bond.a1_id not in aromatic_atoms:
            continue
        rd_bond = mol.GetBondBetweenAtoms(id_map[bond.a1_id], id_map[bond.a2_id])
        if rd_bond is None:
            continue
        bond_type = rd_bond.GetBondType()
        display_orders[bond.id] = 2 if bond_type == Chem.BondType.DOUBLE else 1
    return display_orders


def _scale_to_default(graph: MolGraph, target: float = 40.0) -> None:
    """Escala el grafo para aproximar una longitud de enlace objetivo.

    Args:
        graph: Grafo molecular a escalar.
        target: Longitud promedio deseada (px).

    Side Effects:
        Modifica las coordenadas de los átomos del grafo.
    """
    if not graph.bonds:
        return
    lengths = []
    for bond in graph.bonds.values():
        a1 = graph.atoms[bond.a1_id]
        a2 = graph.atoms[bond.a2_id]
        dx = a2.x - a1.x
        dy = a2.y - a1.y
        lengths.append((dx * dx + dy * dy) ** 0.5)
    avg = sum(lengths) / len(lengths)
    if avg <= 0:
        return
    scale = target / avg
    for atom in graph.atoms.values():
        atom.x *= scale
        atom.y *= scale


@dataclass
class _SmilesNode:
    """Nodo auxiliar para emitir SMILES en el escritor de respaldo."""
    atom_id: int
    symbol: str
    children: list[tuple[str, "_SmilesNode"]]
    ring_closures: list[tuple[str, int]]


def _molgraph_to_smiles_fallback(molgraph: MolGraph) -> str:
    """Emite un SMILES básico sin RDKit.

    Este exportador es deliberadamente simple: maneja aromaticidad básica,
    ramificaciones y cierres de anillo con numeración incremental.

    Args:
        molgraph: Grafo molecular de Chemuson.

    Returns:
        Cadena SMILES aproximada.
    """
    if not molgraph.atoms:
        return ""
    bonds = _unique_bonds_for_export(molgraph)
    adjacency: dict[int, list[tuple[int, Bond]]] = {}
    for bond in bonds:
        adjacency.setdefault(bond.a1_id, []).append((bond.a2_id, bond))
        adjacency.setdefault(bond.a2_id, []).append((bond.a1_id, bond))
    for neighbors in adjacency.values():
        neighbors.sort(key=lambda item: item[0])

    aromatic_atoms: set[int] = set()
    for bond in bonds:
        if bond.is_aromatic:
            aromatic_atoms.add(bond.a1_id)
            aromatic_atoms.add(bond.a2_id)

    nodes_by_id: dict[int, _SmilesNode] = {}
    edge_handled: set[tuple[int, int]] = set()
    ring_counter = 1

    organic_subset = {"B", "C", "N", "O", "P", "S", "F", "Cl", "Br", "I", "b", "c", "n", "o", "p", "s"}

    assigned_h_counter = getattr(molgraph, "assigned_hydrogen_count", None)

    def assigned_hydrogens(atom_id: int) -> int:
        """Devuelve H asignados en el átomo, excluyendo H explícitos como nodos."""
        if callable(assigned_h_counter):
            try:
                return int(max(0, assigned_h_counter(atom_id)))
            except Exception:
                pass
        atom = molgraph.atoms[atom_id]
        return int(getattr(atom, "explicit_h", 0) or 0)

    def simple_aromatic_atom(atom_id: int) -> bool:
        """Indica si el átomo puede emitirse como aromático sin corchetes."""
        atom = molgraph.atoms[atom_id]
        return (
            atom_id in aromatic_atoms
            and atom.element in {"B", "C", "N", "O", "P", "S"}
            and atom.isotope is None
            and not atom.charge
            and assigned_hydrogens(atom_id) <= 0
            and atom.mapping is None
        )

    def bond_symbol(bond: Bond) -> str:
        """Convierte un enlace en su símbolo SMILES."""
        if bond.is_aromatic:
            if simple_aromatic_atom(bond.a1_id) and simple_aromatic_atom(bond.a2_id):
                return ""
            return ":"
        if bond.order == 2:
            return "="
        if bond.order == 3:
            return "#"
        return ""

    def atom_symbol(atom_id: int) -> str:
        """Construye el símbolo SMILES del átomo con carga/isótopo."""
        atom = molgraph.atoms[atom_id]
        element = atom.element
        aromatic = atom_id in aromatic_atoms and element in {"B", "C", "N", "O", "P", "S"}
        symbol = element.lower() if aromatic else element
        charge = ""
        if atom.charge > 0:
            charge = f"+{atom.charge}" if atom.charge > 1 else "+"
        elif atom.charge < 0:
            magnitude = abs(atom.charge)
            charge = f"-{magnitude}" if magnitude > 1 else "-"
        isotope = f"{atom.isotope}" if atom.isotope is not None else ""
        explicit_h = ""
        assigned_h = assigned_hydrogens(atom_id)
        if assigned_h > 0:
            explicit_h = "H" if assigned_h == 1 else f"H{assigned_h}"
        if atom.mapping is not None:
            return f"[{isotope}{symbol}{explicit_h}{charge}:{atom.mapping}]"
        if isotope or explicit_h or charge or symbol not in organic_subset:
            return f"[{isotope}{symbol}{explicit_h}{charge}]"
        return symbol

    def build(node_id: int, parent_id: Optional[int]) -> _SmilesNode:
        """Construye el árbol DFS y registra cierres de anillo."""
        nonlocal ring_counter
        node = nodes_by_id.get(node_id)
        if node is None:
            node = _SmilesNode(
                atom_id=node_id,
                symbol=atom_symbol(node_id),
                children=[],
                ring_closures=[],
            )
            nodes_by_id[node_id] = node
        for neighbor_id, bond in adjacency.get(node_id, []):
            if neighbor_id == parent_id:
                continue
            edge_key = (min(node_id, neighbor_id), max(node_id, neighbor_id))
            if neighbor_id not in nodes_by_id:
                child = build(neighbor_id, node_id)
                node.children.append((bond_symbol(bond), child))
            else:
                # Enlace ya visto: representamos un cierre de anillo.
                if edge_key in edge_handled:
                    continue
                edge_handled.add(edge_key)
                ring_id = ring_counter
                ring_counter += 1
                node.ring_closures.append((bond_symbol(bond), ring_id))
                other = nodes_by_id[neighbor_id]
                other.ring_closures.append((bond_symbol(bond), ring_id))
        return node

    def ring_token(symbol: str, ring_id: int) -> str:
        """Convierte ID de anillo en token SMILES, usando % para >=10."""
        if ring_id >= 10:
            return f"{symbol}%{ring_id}"
        return f"{symbol}{ring_id}"

    def emit(node: _SmilesNode) -> str:
        """Emite SMILES recorriendo el árbol con ramificaciones."""
        subtree_sizes: dict[int, int] = {}

        def subtree_size(current: _SmilesNode) -> int:
            cached = subtree_sizes.get(current.atom_id)
            if cached is not None:
                return cached
            size = 1
            for _symbol, child in current.children:
                size += subtree_size(child)
            subtree_sizes[current.atom_id] = size
            return size

        def main_child_index(children: list[tuple[str, _SmilesNode]]) -> int:
            def bond_priority(symbol: str) -> int:
                if symbol == "":
                    return 3
                if symbol == ":":
                    return 2
                if symbol == "=":
                    return 1
                if symbol == "#":
                    return 0
                return -1

            return max(
                range(len(children)),
                key=lambda idx: (
                    subtree_size(children[idx][1]),
                    bond_priority(children[idx][0]),
                    -children[idx][1].atom_id,
                ),
            )

        text = node.symbol
        for symbol, ring_id in node.ring_closures:
            text += ring_token(symbol, ring_id)
        if not node.children:
            return text
        main_idx = main_child_index(node.children)
        for idx, (symbol, child) in enumerate(node.children):
            branch = f"{symbol}{emit(child)}"
            if idx == main_idx:
                continue
            text += f"({branch})"
        main_symbol, main_child = node.children[main_idx]
        text += f"{main_symbol}{emit(main_child)}"
        return text

    seen: set[int] = set()
    parts: list[str] = []
    for atom_id in sorted(molgraph.atoms.keys()):
        if atom_id in seen:
            continue
        node = build(atom_id, None)
        for node_id in nodes_by_id:
            seen.add(node_id)
        parts.append(emit(node))
    return ".".join(parts)


def _molgraph_to_molfile_fallback(molgraph: MolGraph) -> str:
    """Serializa un `MolGraph` a un MOL V2000 mínimo sin RDKit.

    Args:
        molgraph: Grafo molecular de Chemuson.

    Returns:
        Cadena con el bloque MOL (V2000) básico.
    """
    atoms = [molgraph.atoms[atom_id] for atom_id in sorted(molgraph.atoms.keys())]
    bonds = _unique_bonds_for_export(molgraph)
    atom_index = {atom.id: idx + 1 for idx, atom in enumerate(atoms)}
    aliases: list[tuple[int, str]] = []

    lines = [
        "Chemuson",
        "Chemuson",
        "",
        f"{len(atoms):>3}{len(bonds):>3}  0  0  0  0  0  0  0  0999 V2000",
    ]

    for idx, atom in enumerate(atoms, start=1):
        symbol = atom.element.strip() or "C"
        symbol_for_line = symbol
        if _mol_symbol_requires_fallback(symbol) or len(symbol) > 3:
            symbol_for_line = "*"
            aliases.append((idx, symbol))
        lines.append(
            f"{atom.x:>10.4f}{atom.y:>10.4f}{0.0:>10.4f} {symbol_for_line:<3} 0  0  0  0  0  0  0  0  0  0  0  0"
        )

    for bond in bonds:
        bond_type = 4 if bond.is_aromatic else max(1, min(3, bond.order))
        lines.append(
            f"{atom_index[bond.a1_id]:>3}{atom_index[bond.a2_id]:>3}{bond_type:>3}  0  0  0  0"
        )

    charges: list[tuple[int, int]] = []
    isotopes: list[tuple[int, int]] = []
    radicals: list[tuple[int, int]] = []
    for atom in atoms:
        if atom.charge:
            charges.append((atom_index[atom.id], atom.charge))
        if atom.isotope is not None:
            isotopes.append((atom_index[atom.id], int(atom.isotope)))
        radical_electrons = int(getattr(atom, "radical_electrons", 0) or 0)
        if radical_electrons > 0:
            # Molfile RAD code: 1=singlet, 2=doublet, 3=triplet.
            rad_code = {1: 2, 2: 1}.get(radical_electrons, 3)
            radicals.append((atom_index[atom.id], rad_code))
    for i in range(0, len(charges), 8):
        chunk = charges[i : i + 8]
        parts = [f"{len(chunk):>3}"]
        for idx, charge in chunk:
            parts.append(f"{idx:>3}{charge:>3}")
        lines.append(f"M  CHG{''.join(parts)}")
    for i in range(0, len(isotopes), 8):
        chunk = isotopes[i : i + 8]
        parts = [f"{len(chunk):>3}"]
        for idx, mass in chunk:
            parts.append(f"{idx:>3}{mass:>3}")
        lines.append(f"M  ISO{''.join(parts)}")
    for i in range(0, len(radicals), 8):
        chunk = radicals[i : i + 8]
        parts = [f"{len(chunk):>3}"]
        for idx, rad_code in chunk:
            parts.append(f"{idx:>3}{rad_code:>3}")
        lines.append(f"M  RAD{''.join(parts)}")

    for idx, label in aliases:
        lines.append(f"A  {idx:>3}")
        lines.append(label)

    lines.append("M  END")
    return "\n".join(lines)
