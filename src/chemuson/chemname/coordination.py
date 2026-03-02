"""Nomenclatura mínima para complejos de coordinación (experimental)."""

from __future__ import annotations

from collections import defaultdict
from dataclasses import dataclass
from math import acos, sqrt

from chemuson.chemcalc.valence import implicit_h_count

from .molview import MolView

_METAL_NAMES = {
    "Fe": "iron",
    "Ni": "nickel",
    "Co": "cobalt",
    "Cu": "copper",
    "Zn": "zinc",
    "Pt": "platinum",
    "Pd": "palladium",
    "Ag": "silver",
    "Au": "gold",
}

_METALS = set(_METAL_NAMES.keys())

_SIMPLE_MULTIPLIER = {
    2: "di",
    3: "tri",
    4: "tetra",
    5: "penta",
    6: "hexa",
}

_COMPLEX_MULTIPLIER = {
    2: "bis",
    3: "tris",
    4: "tetrakis",
    5: "pentakis",
    6: "hexakis",
}

_HALO_LIGANDS = {
    "F": "fluoro",
    "Cl": "chloro",
    "Br": "bromo",
    "I": "iodo",
}


@dataclass(frozen=True)
class _Ligand:
    """Representa un ligando identificado alrededor del metal."""

    name: str
    charge: int
    atoms: frozenset[int]
    anchor: int


def detect_coordination_name(view: MolView) -> str | None:
    """Nombra complejos de coordinación simples; retorna `None` si no aplica."""
    metals = [atom_id for atom_id in view.atoms() if view.element(atom_id) in _METALS]
    if len(metals) != 1:
        return None
    metal = metals[0]
    metal_neighbors = [nbr for nbr in view.neighbors(metal) if view.element(nbr) != "H"]
    if not metal_neighbors:
        return None

    ligands: list[_Ligand] = []
    consumed_atoms: set[int] = {metal}
    consumed_neighbors: set[int] = set()

    for cp_ligand in _detect_eta5_cyclopentadienyl(view, metal):
        ligands.append(cp_ligand)
        consumed_atoms |= set(cp_ligand.atoms)
        consumed_neighbors.add(cp_ligand.anchor)

    for nbr in metal_neighbors:
        if nbr in consumed_neighbors:
            continue
        component = _component_without_metal(view, start=nbr, metal=metal)
        if component & consumed_atoms:
            continue
        ligand = _classify_ligand_component(view, metal, nbr, component)
        if ligand is None:
            return None
        ligands.append(ligand)
        consumed_atoms |= set(component)

    if not ligands:
        return None

    ligand_text = _render_ligands(ligands)
    prefix = _cis_trans_prefix(view, metal, ligands)
    metal_name = _METAL_NAMES.get(view.element(metal), view.element(metal).lower())
    oxidation = _oxidation_state(view, ligands)
    if oxidation is None:
        return f"{prefix}{ligand_text}{metal_name}"
    return f"{prefix}{ligand_text}{metal_name}({_roman(oxidation)})"


def _component_without_metal(view: MolView, start: int, metal: int) -> set[int]:
    """Obtiene el componente conectado al ligando excluyendo el centro metálico."""
    stack = [start]
    visited: set[int] = set()
    while stack:
        atom_id = stack.pop()
        if atom_id in visited or atom_id == metal:
            continue
        visited.add(atom_id)
        for nbr in view.neighbors(atom_id):
            if nbr == metal or view.element(nbr) == "H" or nbr in visited:
                continue
            stack.append(nbr)
    return visited


def _classify_ligand_component(
    view: MolView,
    metal: int,
    anchor: int,
    component: set[int],
) -> _Ligand | None:
    """Clasifica ligandos puntuales comunes (CO, NH3, H2O, haluros, CN)."""
    if not component:
        return None
    if len(component) == 1:
        atom_id = anchor
        element = view.element(atom_id)
        if element in _HALO_LIGANDS:
            return _Ligand(_HALO_LIGANDS[element], -1, frozenset(component), anchor)
        if element == "N":
            h_total = implicit_h_count(view, atom_id) + view.explicit_h(atom_id)
            heavy_neighbors = [nbr for nbr in view.neighbors(atom_id) if view.element(nbr) != "H"]
            if set(heavy_neighbors) == {metal} and h_total >= 3:
                return _Ligand("ammine", 0, frozenset(component), anchor)
        if element == "O":
            h_total = implicit_h_count(view, atom_id) + view.explicit_h(atom_id)
            heavy_neighbors = [nbr for nbr in view.neighbors(atom_id) if view.element(nbr) != "H"]
            if set(heavy_neighbors) == {metal} and h_total >= 2:
                return _Ligand("aqua", 0, frozenset(component), anchor)

    if len(component) == 2 and view.element(anchor) == "C":
        other = next(iter(component - {anchor}))
        if view.element(other) == "O":
            order = view.bond_order_between(anchor, other)
            if order in {2, 3}:
                return _Ligand("carbonyl", 0, frozenset(component), anchor)
        if view.element(other) == "N":
            order = view.bond_order_between(anchor, other)
            if order == 3:
                return _Ligand("cyano", -1, frozenset(component), anchor)

    return None


def _detect_eta5_cyclopentadienyl(view: MolView, metal: int) -> list[_Ligand]:
    """Detecta anillos Cp enlazados al metal con hapticidad η5 (MVP)."""
    ligands: list[_Ligand] = []
    candidate_atoms = {
        atom_id
        for atom_id in view.neighbors(metal)
        if view.element(atom_id) == "C"
    }
    if len(candidate_atoms) < 5:
        return ligands

    adjacency: dict[int, set[int]] = {atom_id: set() for atom_id in candidate_atoms}
    for atom_id in candidate_atoms:
        for nbr in view.neighbors(atom_id):
            if nbr not in candidate_atoms:
                continue
            if not view.bond_is_aromatic(atom_id, nbr):
                continue
            adjacency[atom_id].add(nbr)

    visited: set[int] = set()
    for atom_id in sorted(candidate_atoms):
        if atom_id in visited:
            continue
        stack = [atom_id]
        component: set[int] = set()
        while stack:
            current = stack.pop()
            if current in visited:
                continue
            visited.add(current)
            component.add(current)
            for nbr in adjacency.get(current, set()):
                if nbr not in visited:
                    stack.append(nbr)
        if len(component) != 5:
            continue
        if any(len(adjacency.get(node, set()) & component) != 2 for node in component):
            continue
        attached = [node for node in component if metal in view.neighbors(node)]
        if len(attached) < 3:
            continue
        anchor = min(attached)
        ligands.append(_Ligand("η5-cyclopentadienyl", -1, frozenset(component), anchor))
    return ligands


def _render_ligands(ligands: list[_Ligand]) -> str:
    """Renderiza lista de ligandos ordenada alfabéticamente con multiplicadores."""
    grouped: dict[str, list[_Ligand]] = defaultdict(list)
    for lig in ligands:
        grouped[lig.name].append(lig)

    blocks: list[str] = []
    for name in sorted(grouped.keys()):
        count = len(grouped[name])
        if count == 1:
            blocks.append(name)
            continue
        if _is_complex_ligand_name(name):
            multiplier = _COMPLEX_MULTIPLIER.get(count)
            if multiplier is None:
                return ""
            blocks.append(f"{multiplier}({name})")
            continue
        multiplier = _SIMPLE_MULTIPLIER.get(count)
        if multiplier is None:
            return ""
        blocks.append(f"{multiplier}{name}")
    return "".join(blocks)


def _is_complex_ligand_name(name: str) -> bool:
    """Define si un ligando requiere multiplicador tipo bis/tris."""
    return ("η" in name) or ("-" in name) or ("(" in name) or (")" in name)


def _cis_trans_prefix(view: MolView, metal: int, ligands: list[_Ligand]) -> str:
    """Infere descriptor cis/trans en complejos tetracoordinados simples."""
    if len(ligands) != 4:
        return ""
    by_name: dict[str, list[_Ligand]] = defaultdict(list)
    for lig in ligands:
        by_name[lig.name].append(lig)
    pair = next((lst for lst in by_name.values() if len(lst) == 2), None)
    if pair is None:
        return ""

    center = view._get_atom(metal)  # noqa: SLF001 - lectura controlada para coordenadas.
    if center is None:
        return ""
    cx = float(getattr(center, "x", 0.0))
    cy = float(getattr(center, "y", 0.0))
    v1 = _vector_to_atom(view, pair[0].anchor, cx, cy)
    v2 = _vector_to_atom(view, pair[1].anchor, cx, cy)
    if v1 is None or v2 is None:
        return ""
    angle = _angle_deg(v1, v2)
    if angle >= 150.0:
        return "trans-"
    return "cis-"


def _vector_to_atom(view: MolView, atom_id: int, cx: float, cy: float) -> tuple[float, float] | None:
    """Devuelve vector desde el centro metálico hacia un átomo ancla."""
    atom = view._get_atom(atom_id)  # noqa: SLF001 - lectura controlada para coordenadas.
    if atom is None:
        return None
    x = float(getattr(atom, "x", 0.0))
    y = float(getattr(atom, "y", 0.0))
    return (x - cx, y - cy)


def _angle_deg(v1: tuple[float, float], v2: tuple[float, float]) -> float:
    """Calcula ángulo entre dos vectores en grados."""
    n1 = sqrt(v1[0] * v1[0] + v1[1] * v1[1])
    n2 = sqrt(v2[0] * v2[0] + v2[1] * v2[1])
    if n1 <= 1e-9 or n2 <= 1e-9:
        return 0.0
    dot = (v1[0] * v2[0] + v1[1] * v2[1]) / (n1 * n2)
    dot = max(-1.0, min(1.0, dot))
    return acos(dot) * 180.0 / 3.141592653589793


def _oxidation_state(view: MolView, ligands: list[_Ligand]) -> int | None:
    """Estima estado de oxidación con carga total del complejo y ligandos."""
    metal_id = next((atom_id for atom_id in view.atoms() if view.element(atom_id) in _METALS), None)
    if metal_id is None:
        return None
    metal_formal = view.formal_charge(metal_id)
    lig_charge = sum(lig.charge for lig in ligands)
    return int(metal_formal - lig_charge)


def _roman(value: int) -> str:
    """Convierte enteros pequeños a número romano (incluye 0)."""
    if value == 0:
        return "0"
    negative = value < 0
    n = abs(int(value))
    table = [
        (1000, "M"),
        (900, "CM"),
        (500, "D"),
        (400, "CD"),
        (100, "C"),
        (90, "XC"),
        (50, "L"),
        (40, "XL"),
        (10, "X"),
        (9, "IX"),
        (5, "V"),
        (4, "IV"),
        (1, "I"),
    ]
    out = []
    for arabic, roman in table:
        while n >= arabic:
            out.append(roman)
            n -= arabic
    rendered = "".join(out) if out else "0"
    return f"-{rendered}" if negative else rendered
