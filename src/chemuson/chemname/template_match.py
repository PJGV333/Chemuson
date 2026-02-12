"""Emparejamiento de plantillas moleculares con grafos dibujados."""

from __future__ import annotations

from typing import Dict, Iterable, List

from .errors import ChemNameNotSupported
from .molview import MolView
from .rings import find_rings_simple, _ring_aromatic_basic
from .template import TemplateBond, TemplateMol


def match_template_exact(
    template: TemplateMol,
    view: MolView,
    atom_ids: Iterable[int] | None = None,
) -> List[Dict[int, int]]:
    """Busca mapeos exactos entre una plantilla y el grafo.

    Args:
        template: Plantilla con átomos, enlaces y locantes.
        view: Vista del grafo molecular.
        atom_ids: Subconjunto de átomos candidatos (opcional).

    Returns:
        Lista de mapeos `template_idx -> mol_atom_id` válidos.
    """
    atom_ids = list(atom_ids) if atom_ids is not None else view.atoms()
    atom_set = set(atom_ids)
    if len(atom_ids) != len(template.atoms):
        return []

    mol_element = {atom_id: view.element(atom_id) for atom_id in atom_ids}
    mol_degree = {
        atom_id: len(
            [
                nbr
                for nbr in view.neighbors(atom_id)
                if nbr in atom_set and view.element(nbr) != "H"
            ]
        )
        for atom_id in atom_ids
    }
    mol_aromatic = _aromatic_atoms(view, atom_set)

    mol_bonds: Dict[tuple[int, int], tuple[int, bool]] = {}
    for a in atom_ids:
        for b in view.neighbors(a):
            if b not in atom_set:
                continue
            key = (a, b) if a < b else (b, a)
            if key in mol_bonds:
                continue
            order = view.bond_order_between(a, b)
            aromatic = view.bond_is_aromatic(a, b)
            if not aromatic and mol_aromatic.get(a) and mol_aromatic.get(b) and order in {1, 2}:
                aromatic = True
            mol_bonds[key] = (order, aromatic)

    template_adj: Dict[int, List[TemplateBond]] = {i: [] for i in range(len(template.atoms))}
    bonded_pairs = set()
    for bond in template.bonds:
        template_adj[bond.a].append(bond)
        template_adj[bond.b].append(TemplateBond(bond.b, bond.a, bond.order, bond.aromatic))
        pair = (bond.a, bond.b) if bond.a < bond.b else (bond.b, bond.a)
        bonded_pairs.add(pair)

    # Ordenamos por grado/aromaticidad para podar el backtracking.
    order = sorted(
        range(len(template.atoms)),
        key=lambda idx: (
            -template.atoms[idx].degree,
            0 if template.atoms[idx].aromatic else 1,
            template.atoms[idx].element,
            idx,
        ),
    )

    candidates: Dict[int, List[int]] = {}
    for t_idx in range(len(template.atoms)):
        t_atom = template.atoms[t_idx]
        options = []
        for m_id in atom_ids:
            if mol_element[m_id] != t_atom.element:
                continue
            if mol_degree[m_id] != t_atom.degree:
                continue
            if t_atom.aromatic != bool(mol_aromatic.get(m_id)):
                continue
            options.append(m_id)
        if not options:
            return []
        candidates[t_idx] = options

    results: List[Dict[int, int]] = []
    mapping: Dict[int, int] = {}
    used = set()

    def bond_compatible(tbond: TemplateBond, mol_key: tuple[int, int]) -> bool:
        """Comprueba compatibilidad de un enlace plantilla con el grafo."""
        if mol_key not in mol_bonds:
            return False
        order, aromatic = mol_bonds[mol_key]
        if tbond.aromatic:
            return aromatic or order in {1, 2}
        return order == tbond.order

    def backtrack(pos: int) -> None:
        """Búsqueda en profundidad con poda por enlaces incompatibles."""
        if pos == len(order):
            if _mapping_is_exact(mapping, bonded_pairs, mol_bonds):
                results.append(dict(mapping))
            return
        t_idx = order[pos]
        for m_id in candidates[t_idx]:
            if m_id in used:
                continue
            ok = True
            for tbond in template_adj[t_idx]:
                if tbond.b not in mapping:
                    continue
                m_other = mapping[tbond.b]
                key = (m_id, m_other) if m_id < m_other else (m_other, m_id)
                if not bond_compatible(tbond, key):
                    ok = False
                    break
            if not ok:
                continue
            # Evitar falsos positivos: no se permiten enlaces extra fuera
            # de los definidos en la plantilla.
            for t_other, m_other in mapping.items():
                if t_other == t_idx:
                    continue
                pair = (t_idx, t_other) if t_idx < t_other else (t_other, t_idx)
                if pair in bonded_pairs:
                    continue
                key = (m_id, m_other) if m_id < m_other else (m_other, m_id)
                if key in mol_bonds:
                    ok = False
                    break
            if not ok:
                continue
            mapping[t_idx] = m_id
            used.add(m_id)
            backtrack(pos + 1)
            used.remove(m_id)
            mapping.pop(t_idx, None)

    backtrack(0)
    return results


def select_template_mapping(
    template: TemplateMol,
    view: MolView,
    mappings: Iterable[Dict[int, int]],
) -> Dict[int, int] | None:
    """Selecciona el mejor mapeo según locantes y orden estable.

    Args:
        template: Plantilla de referencia.
        view: Vista del grafo molecular.
        mappings: Iterable de mapeos candidatos.

    Returns:
        Mapeo elegido o `None` si no hay candidatos.
    """
    best = None
    best_key = None
    for mapping in mappings:
        locants = _mapping_substituent_locants(template, view, mapping)
        has_locants = bool(locants)
        loc_key = tuple(sorted(locants))
        map_key = _mapping_key(mapping)
        key = (0, loc_key, map_key) if has_locants else (1, (), map_key)
        if best_key is None or key < best_key:
            best_key = key
            best = mapping
    return best


def _mapping_substituent_locants(
    template: TemplateMol, view: MolView, mapping: Dict[int, int]
) -> List[int]:
    """Extrae los locantes de sustitución dados un mapeo."""
    inverse = {mol_id: t_idx for t_idx, mol_id in mapping.items()}
    ring_set = set(mapping.values())
    locants: List[int] = []
    for mol_id in ring_set:
        for nbr in view.neighbors(mol_id):
            if nbr in ring_set:
                continue
            if view.element(nbr) == "H":
                continue
            t_idx = inverse[mol_id]
            locant = template.locant_by_atom_idx.get(t_idx)
            if locant is None:
                continue
            locants.append(locant)
    return locants


def _mapping_key(mapping: Dict[int, int]) -> tuple[int, ...]:
    """Clave determinista de un mapeo para desempates."""
    return tuple(mapping[idx] for idx in sorted(mapping.keys()))


def _mapping_is_exact(
    mapping: Dict[int, int],
    bonded_pairs: set[tuple[int, int]],
    mol_bonds: Dict[tuple[int, int], tuple[int, bool]],
) -> bool:
    """Verifica que el mapeo respete exactamente enlaces y no-enlaces."""
    t_indices = list(mapping.keys())
    for i in range(len(t_indices)):
        for j in range(i + 1, len(t_indices)):
            t_a = t_indices[i]
            t_b = t_indices[j]
            pair = (t_a, t_b) if t_a < t_b else (t_b, t_a)
            m_a = mapping[t_a]
            m_b = mapping[t_b]
            m_pair = (m_a, m_b) if m_a < m_b else (m_b, m_a)
            if pair in bonded_pairs:
                if m_pair not in mol_bonds:
                    return False
            else:
                if m_pair in mol_bonds:
                    return False
    return True


def _aromatic_atoms(view: MolView, atom_set: set[int]) -> Dict[int, bool]:
    """Marca átomos aromáticos dentro del subconjunto dado."""
    aromatic: Dict[int, bool] = {atom_id: False for atom_id in atom_set}
    rings = find_rings_simple(view)
    for ring in rings:
        if not ring.issubset(atom_set):
            continue
        if _ring_aromatic_basic(view, ring):
            for atom_id in ring:
                aromatic[atom_id] = True
    return aromatic
