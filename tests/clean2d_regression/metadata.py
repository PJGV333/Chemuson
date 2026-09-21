from __future__ import annotations

from collections.abc import Mapping
from typing import TYPE_CHECKING, Any

from chemuson.core.model import MolGraph

if TYPE_CHECKING:
    from .cases import Clean2DRegressionCase


SIZE_CLASSES = frozenset({"simple", "medium", "large", "complex-scale"})


def validate_case_taxonomy(case: Clean2DRegressionCase) -> None:
    """Validate stable, observational taxonomy declared by a corpus case."""

    if not str(case.name).strip():
        raise AssertionError("regression case name is required")
    if case.size_class not in SIZE_CLASSES:
        raise AssertionError(f"{case.name}: invalid size class {case.size_class!r}")
    if not str(case.family).strip():
        raise AssertionError(f"{case.name}: family is required")
    if not case.tags or any(not str(tag).strip() for tag in case.tags):
        raise AssertionError(f"{case.name}: non-empty classification tags are required")


def derive_topology_metadata(graph: MolGraph) -> dict[str, Any]:
    """Return deterministic topology facts without changing or routing a graph.

    Block and connector decomposition is intentionally left unavailable until a
    later campaign owns that contract; `None` is more honest than a guessed
    value in a benchmark report.
    """

    atom_ids = tuple(sorted(int(atom_id) for atom_id in graph.atoms))
    adjacency: dict[int, set[int]] = {atom_id: set() for atom_id in atom_ids}
    for bond in graph.bonds.values():
        if bond.a1_id not in adjacency or bond.a2_id not in adjacency:
            continue
        adjacency[bond.a1_id].add(bond.a2_id)
        adjacency[bond.a2_id].add(bond.a1_id)

    connected_components = _count_components(adjacency)
    atom_count = len(atom_ids)
    bond_count = sum(1 for bond in graph.bonds.values() if bond.a1_id in adjacency and bond.a2_id in adjacency)
    ring_count = max(0, bond_count - atom_count + connected_components)
    heavy_atom_count = sum(
        1 for atom_id in atom_ids if str(getattr(graph.atoms[atom_id], "element", "")).upper() != "H"
    )
    branch_point_count = sum(1 for neighbors in adjacency.values() if len(neighbors) >= 3)

    return {
        "atom_count": atom_count,
        "heavy_atom_count": heavy_atom_count,
        "bond_count": bond_count,
        "ring_count": ring_count,
        "connected_components": connected_components,
        "branch_point_count": branch_point_count,
        "rigid_block_count": None,
        "rotatable_connector_count": None,
        "macrocycle_count": None,
    }


def _count_components(adjacency: Mapping[int, set[int]]) -> int:
    remaining = set(adjacency)
    components = 0
    while remaining:
        components += 1
        pending = [remaining.pop()]
        while pending:
            atom_id = pending.pop()
            for neighbor in adjacency[atom_id]:
                if neighbor in remaining:
                    remaining.remove(neighbor)
                    pending.append(neighbor)
    return components
