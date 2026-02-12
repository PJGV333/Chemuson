"""Pruebas de numeración visual de átomos y estructuras."""

from __future__ import annotations

from chemuson.core.model import MolGraph
from chemuson.gui.numbering import compute_atom_numbers, compute_structure_numbers


def _build_simple_ethanol_like() -> MolGraph:
    """Crea una cadena simple C-C-O para pruebas de orden."""
    graph = MolGraph()
    a1 = graph.add_atom("C", 10.0, 40.0)
    a2 = graph.add_atom("C", 30.0, 38.0)
    a3 = graph.add_atom("O", 50.0, 39.0)
    graph.add_bond(a1.id, a2.id, order=1)
    graph.add_bond(a2.id, a3.id, order=1)
    return graph


def test_compute_atom_numbers_xy_order_is_stable_with_small_moves():
    """La numeración se mantiene si no cambia el orden relativo."""
    graph = _build_simple_ethanol_like()
    numbers_before = compute_atom_numbers(graph)

    # Pequeño movimiento que no altera orden x asc / y desc.
    atom = graph.get_atom(2)
    atom.x += 0.2
    atom.y += 0.1
    numbers_after = compute_atom_numbers(graph)

    assert numbers_before == numbers_after
    assert numbers_after == {1: 1, 2: 2, 3: 3}


def test_compute_atom_numbers_xy_tie_breaks_by_y_desc():
    """Con mismo x, el átomo más alto (mayor y) recibe número menor."""
    graph = MolGraph()
    a_low = graph.add_atom("C", 20.0, 10.0)
    a_high = graph.add_atom("O", 20.0, 30.0)
    graph.add_bond(a_low.id, a_high.id, order=1)

    numbers = compute_atom_numbers(graph)
    assert numbers[a_high.id] == 1
    assert numbers[a_low.id] == 2


def test_compute_structure_numbers_for_disconnected_components():
    """Numera estructuras por centroide (x asc, y desc)."""
    graph = MolGraph()

    # Estructura izquierda (debe ser #1)
    l1 = graph.add_atom("C", 10.0, 20.0)
    l2 = graph.add_atom("O", 20.0, 20.0)
    graph.add_bond(l1.id, l2.id, order=1)

    # Estructura derecha (debe ser #2)
    r1 = graph.add_atom("C", 110.0, 20.0)
    r2 = graph.add_atom("N", 120.0, 20.0)
    graph.add_bond(r1.id, r2.id, order=1)

    structures = compute_structure_numbers(graph)
    assert len(structures) == 2
    assert structures[0].number == 1
    assert structures[1].number == 2
    assert set(structures[0].atom_ids) == {l1.id, l2.id}
    assert set(structures[1].atom_ids) == {r1.id, r2.id}


def test_compute_atom_numbers_resets_per_structure():
    """La numeración atómica se reinicia en cada estructura desconectada."""
    graph = MolGraph()
    # Componente A
    a1 = graph.add_atom("C", 10.0, 20.0)
    a2 = graph.add_atom("C", 20.0, 20.0)
    graph.add_bond(a1.id, a2.id, order=1)
    # Componente B
    b1 = graph.add_atom("N", 110.0, 20.0)
    b2 = graph.add_atom("O", 120.0, 20.0)
    graph.add_bond(b1.id, b2.id, order=1)

    numbers = compute_atom_numbers(graph)
    assert numbers[a1.id] == 1
    assert numbers[a2.id] == 2
    assert numbers[b1.id] == 1
    assert numbers[b2.id] == 2


def test_compute_structure_numbers_includes_isolated_atoms():
    """Un átomo sin enlaces cuenta como estructura propia."""
    graph = MolGraph()
    iso = graph.add_atom("Cl", 200.0, 200.0)

    structures = compute_structure_numbers(graph)
    assert len(structures) == 1
    assert structures[0].number == 1
    assert structures[0].atom_ids == (iso.id,)
