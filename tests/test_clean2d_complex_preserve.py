from __future__ import annotations

from dataclasses import asdict
import math
import os
import sys

sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), "..", "src")))

from chemuson.clean2d import (
    assert_clean2d_invariants,
    classify_clean2d_complexity,
    run_clean2d_engine,
    stereo_layout_signature,
)
from chemuson.core.model import BondStereo, BondStyle, MolGraph


def _add_ring(
    graph: MolGraph,
    cx: float,
    cy: float,
    *,
    aromatic: bool = True,
    radius: float = 40.0,
    oxygen_index: int | None = None,
) -> list[int]:
    atoms: list[int] = []
    for idx in range(6):
        angle = math.radians(idx * 60.0)
        element = "O" if oxygen_index == idx else "C"
        atoms.append(graph.add_atom(element, cx + radius * math.cos(angle), cy + radius * math.sin(angle)).id)
    for idx in range(6):
        graph.add_bond(atoms[idx], atoms[(idx + 1) % 6], order=1, is_aromatic=aromatic)
    return atoms


def _connect_with_zigzag(graph: MolGraph, left_atom: int, right_atom: int, y_offset: float = 26.46) -> list[int]:
    x0, y0 = graph.atoms[left_atom].x, graph.atoms[left_atom].y
    x4, y4 = graph.atoms[right_atom].x, graph.atoms[right_atom].y
    dx = (x4 - x0) / 4.0
    atoms = [
        graph.add_atom("C", x0 + dx, y0 + y_offset).id,
        graph.add_atom("C", x0 + 2.0 * dx, (y0 + y4) / 2.0).id,
        graph.add_atom("C", x0 + 3.0 * dx, y4 + y_offset).id,
    ]
    graph.add_bond(left_atom, atoms[0], order=1)
    graph.add_bond(atoms[0], atoms[1], order=1)
    graph.add_bond(atoms[1], atoms[2], order=1)
    graph.add_bond(atoms[2], right_atom, order=1)
    return atoms


def _build_tetrandrine_like_graph() -> MolGraph:
    graph = MolGraph()
    ring_a = _add_ring(graph, 0.0, 0.0)
    ring_b = _add_ring(graph, 200.0, 0.0)
    linker = _connect_with_zigzag(graph, ring_a[0], ring_b[3], 26.46)
    graph.atoms[linker[0]].stereo_cip = "R"
    graph.atoms[linker[2]].stereo_cip = "S"
    graph.add_bond(linker[0], graph.add_atom("O", graph.atoms[linker[0]].x, graph.atoms[linker[0]].y + 40.0).id, style=BondStyle.WEDGE, stereo=BondStereo.UP)
    graph.add_bond(linker[2], graph.add_atom("N", graph.atoms[linker[2]].x, graph.atoms[linker[2]].y + 40.0).id, style=BondStyle.HASHED, stereo=BondStereo.DOWN)
    return graph


def _build_vancomycin_like_graph() -> MolGraph:
    graph = MolGraph()
    rings = [_add_ring(graph, idx * 200.0, 0.0) for idx in range(3)]
    sugar = _add_ring(graph, 600.0, 0.0, aromatic=False, oxygen_index=0)
    for base, element, cx, cy, style, stereo in (
        (rings[0][0], "N", 0.0, 0.0, BondStyle.WEDGE, BondStereo.UP),
        (rings[1][1], "O", 200.0, 0.0, BondStyle.PLAIN, BondStereo.NONE),
        (rings[2][0], "O", 400.0, 0.0, BondStyle.HASHED, BondStereo.DOWN),
        (rings[2][5], "N", 400.0, 0.0, BondStyle.PLAIN, BondStereo.NONE),
        (sugar[0], "O", 600.0, 0.0, BondStyle.PLAIN, BondStereo.NONE),
    ):
        bx, by = graph.atoms[base].x, graph.atoms[base].y
        length = math.hypot(bx - cx, by - cy) or 1.0
        atom = graph.add_atom(element, bx + (bx - cx) / length * 40.0, by + (by - cy) / length * 40.0).id
        graph.add_bond(base, atom, order=1, style=style, stereo=stereo)
    graph.atoms[rings[0][0]].stereo_cip = "R"
    graph.atoms[rings[2][0]].stereo_cip = "S"
    return graph


def _build_simple_graph() -> MolGraph:
    graph = MolGraph()
    ring = _add_ring(graph, 0.0, 0.0)
    chain = [graph.add_atom("C", 80.0 + idx * 40.0, 0.0).id for idx in range(3)]
    graph.add_bond(ring[0], chain[0], order=1)
    graph.add_bond(chain[0], chain[1], order=1)
    graph.add_bond(chain[1], chain[2], order=1)
    return graph


def _coords(graph: MolGraph) -> dict[int, tuple[float, float]]:
    return {atom_id: (atom.x, atom.y) for atom_id, atom in graph.atoms.items()}


def _mean_displacement(before: dict[int, tuple[float, float]], after: dict[int, tuple[float, float]]) -> float:
    common = set(before) & set(after)
    return sum(math.dist(before[atom_id], after[atom_id]) for atom_id in common) / len(common)


def _chemical_signature(graph: MolGraph) -> tuple[tuple[tuple[int, dict], ...], tuple[tuple[int, dict], ...]]:
    atoms = []
    for atom_id, atom in sorted(graph.atoms.items()):
        data = asdict(atom)
        data.pop("x")
        data.pop("y")
        atoms.append((atom_id, data))
    bonds = [(bond_id, asdict(bond)) for bond_id, bond in sorted(graph.bonds.items())]
    return tuple(atoms), tuple(bonds)


def test_complexity_profile_detects_preserve_only_for_tetrandrine_like_graph() -> None:
    graph = _build_tetrandrine_like_graph()

    profile = classify_clean2d_complexity(graph)

    assert profile.preserve_only is True
    assert profile.reason
    assert profile.aromatic_ring_count >= 2
    assert profile.stereo_center_count >= 2
    assert profile.linker_count >= 1
    assert profile.block_count >= profile.aromatic_ring_count + profile.stereo_center_count


def test_complexity_profile_detects_preserve_only_for_vancomycin_like_graph() -> None:
    graph = _build_vancomycin_like_graph()

    profile = classify_clean2d_complexity(graph)

    assert profile.preserve_only is True
    assert profile.aromatic_ring_count >= 3 or profile.has_hierarchical_blocks
    assert profile.stereo_center_count >= 2


def test_quick_clean_preserves_good_complex_layout() -> None:
    graph = _build_vancomycin_like_graph()
    before = _coords(graph)
    signature = _chemical_signature(graph)

    result = run_clean2d_engine(graph, set(graph.atoms), mode="quick", target_bond_length=40.0)

    assert result.selected is not None
    assert result.selected.source == "current"
    assert "compleja preservada" in result.message or "no se redibujó" in result.message
    assert _mean_displacement(before, result.selected.coords) < 1e-9
    assert _chemical_signature(graph) == signature
    assert_clean2d_invariants(graph, graph, before, result.selected.coords, atom_ids=set(graph.atoms))


def test_quick_clean_complex_never_uses_global_redraw_sources() -> None:
    graph = _build_tetrandrine_like_graph()
    graph.atoms[7].x += 8.0
    graph.atoms[8].y -= 6.0

    result = run_clean2d_engine(graph, set(graph.atoms), mode="quick", target_bond_length=40.0)
    sources = {candidate.source for candidate in (*result.candidates, *result.rejected)}

    assert not sources & {"rdkit_isolated", "rdkit_direct", "internal_templates", "v2", "clean2d_v2", "local_graph"}
    assert sources <= {"current", "complex_preserve"}


def test_complex_preserve_does_not_change_stereo_signature() -> None:
    graph = _build_tetrandrine_like_graph()
    before = _coords(graph)
    before_signature = stereo_layout_signature(graph, before, set(graph.atoms))

    result = run_clean2d_engine(graph, set(graph.atoms), mode="quick", target_bond_length=40.0)
    after_coords = result.selected.coords if result.selected is not None else before

    assert stereo_layout_signature(graph, after_coords, set(graph.atoms)) == before_signature


def test_simple_clean2d_regression_still_uses_existing_flow() -> None:
    graph = _build_simple_graph()

    result = run_clean2d_engine(graph, set(graph.atoms), mode="quick", target_bond_length=40.0)

    assert result.selected is not None
    assert result.selected.metadata.get("complex_preserve_only") is not True
