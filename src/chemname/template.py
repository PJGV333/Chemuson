from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
import json
from typing import Dict, Iterable, List

from .errors import ChemNameNotSupported


@dataclass(frozen=True)
class TemplateAtom:
    element: str
    aromatic: bool
    degree: int


@dataclass(frozen=True)
class TemplateBond:
    a: int
    b: int
    order: int
    aromatic: bool


@dataclass(frozen=True)
class TemplateMol:
    atoms: List[TemplateAtom]
    bonds: List[TemplateBond]
    locant_by_atom_idx: Dict[int, int]


def load_template(path: str | Path) -> TemplateMol:
    path = Path(path)
    if path.suffix.lower() in {".json"}:
        return _load_template_json(path)
    if path.suffix.lower() in {".mol"}:
        return _load_template_mol(path)
    raise ChemNameNotSupported("Unsupported template format")


def _load_template_json(path: Path) -> TemplateMol:
    data = json.loads(path.read_text(encoding="utf-8"))
    atoms_data = data.get("atoms", [])
    bonds_data = data.get("bonds", [])
    locants_data = data.get("locants", {})
    if not atoms_data or not bonds_data:
        raise ChemNameNotSupported("Invalid template JSON")

    bonds: List[TemplateBond] = []
    for bond in bonds_data:
        a = int(bond["a"])
        b = int(bond["b"])
        order = int(bond.get("order", 1))
        aromatic = bool(bond.get("arom", False))
        bonds.append(TemplateBond(a=a, b=b, order=order, aromatic=aromatic))

    degree = [0 for _ in atoms_data]
    aromatic_atoms = [False for _ in atoms_data]
    for bond in bonds:
        degree[bond.a] += 1
        degree[bond.b] += 1
        if bond.aromatic:
            aromatic_atoms[bond.a] = True
            aromatic_atoms[bond.b] = True

    atoms: List[TemplateAtom] = []
    for idx, atom in enumerate(atoms_data):
        element = str(atom.get("el", "")).strip()
        if not element:
            raise ChemNameNotSupported("Invalid template atom")
        elem = element[0].upper() + element[1:].lower()
        aromatic = bool(atom.get("arom", False)) or aromatic_atoms[idx]
        atoms.append(TemplateAtom(element=elem, aromatic=aromatic, degree=degree[idx]))

    locant_by_atom_idx: Dict[int, int] = {
        int(k): int(v) for k, v in locants_data.items()
    }
    return TemplateMol(atoms=atoms, bonds=bonds, locant_by_atom_idx=locant_by_atom_idx)


def _load_template_mol(path: Path) -> TemplateMol:
    lines = path.read_text(encoding="utf-8", errors="replace").splitlines()
    if len(lines) < 4:
        raise ChemNameNotSupported("Invalid mol template")

    counts_idx = None
    for idx, line in enumerate(lines[:10]):
        if "V2000" in line:
            counts_idx = idx
            break
    if counts_idx is None:
        counts_idx = 3

    counts_line = lines[counts_idx]
    try:
        atom_count = int(counts_line[0:3])
        bond_count = int(counts_line[3:6])
    except Exception:
        parts = counts_line.split()
        if len(parts) < 2:
            raise ChemNameNotSupported("Invalid mol counts line")
        atom_count = int(parts[0])
        bond_count = int(parts[1])

    atom_lines = lines[counts_idx + 1 : counts_idx + 1 + atom_count]
    bond_lines = lines[counts_idx + 1 + atom_count : counts_idx + 1 + atom_count + bond_count]
    if len(atom_lines) != atom_count or len(bond_lines) != bond_count:
        raise ChemNameNotSupported("Invalid mol template size")

    atoms_element: List[str] = []
    atom_map: Dict[int, int] = {}
    for idx, line in enumerate(atom_lines):
        parts = line.split()
        if len(parts) < 4:
            raise ChemNameNotSupported("Invalid mol atom line")
        elem = parts[3].strip()
        elem = elem[0].upper() + elem[1:].lower()
        atoms_element.append(elem)

        map_num = None
        if len(line) >= 63:
            field = line[60:63].strip()
            if field.isdigit():
                map_num = int(field)
        if (map_num is None or map_num == 0) and len(parts) >= 17 and parts[-1].lstrip("-").isdigit():
            candidate = int(parts[-1])
            if candidate != 0:
                map_num = candidate
        if map_num and map_num > 0:
            atom_map[idx] = map_num

    bonds: List[TemplateBond] = []
    degree = [0 for _ in atoms_element]
    aromatic_atoms = [False for _ in atoms_element]
    for line in bond_lines:
        parts = line.split()
        if len(parts) < 3:
            raise ChemNameNotSupported("Invalid mol bond line")
        a = int(parts[0]) - 1
        b = int(parts[1]) - 1
        order = int(parts[2])
        aromatic = False
        if order == 4:
            aromatic = True
            order = 1
        bonds.append(TemplateBond(a=a, b=b, order=order, aromatic=aromatic))
        degree[a] += 1
        degree[b] += 1
        if aromatic:
            aromatic_atoms[a] = True
            aromatic_atoms[b] = True

    atoms: List[TemplateAtom] = []
    for idx, elem in enumerate(atoms_element):
        atoms.append(TemplateAtom(element=elem, aromatic=aromatic_atoms[idx], degree=degree[idx]))

    return TemplateMol(atoms=atoms, bonds=bonds, locant_by_atom_idx=atom_map)
