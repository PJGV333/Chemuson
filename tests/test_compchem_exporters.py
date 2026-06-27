"""Pruebas de exportadores de química computacional."""

from __future__ import annotations



from chemuson.compchem.exporters import export_gaussian_input, export_nwchem_input, export_orca_input
from chemuson.core.model import MolGraph
from chemuson.geometry3d.model import CoordinateSet3D


def _water() -> tuple[MolGraph, CoordinateSet3D]:
    graph = MolGraph()
    o = graph.add_atom("O", 0.0, 0.0)
    h1 = graph.add_atom("H", 1.0, 0.0, is_explicit=True)
    h2 = graph.add_atom("H", -1.0, 0.0, is_explicit=True)
    graph.add_bond(o.id, h1.id)
    graph.add_bond(o.id, h2.id)
    coords = CoordinateSet3D(
        {o.id: (0.0, 0.0, 0.0), h1.id: (0.96, 0.0, 0.0), h2.id: (-0.24, 0.93, 0.0)},
        source="test",
        method="manual",
    )
    return graph, coords


def test_export_orca_input() -> None:
    graph, coords = _water()

    text = export_orca_input(graph, coords, charge=0, multiplicity=1, method="B3LYP", basis="def2-SVP", cores=4, calculation="opt")

    assert "! B3LYP def2-SVP opt" in text
    assert "%pal nprocs 4 end" in text
    assert "* xyz 0 1" in text
    assert "O    0.00000000" in text


def test_export_gaussian_input() -> None:
    graph, coords = _water()

    text = export_gaussian_input(graph, coords, charge=1, multiplicity=2, memory_mb=4000, cores=2, calculation="freq")

    assert "%mem=4000MB" in text
    assert "%nprocshared=2" in text
    assert "# B3LYP/6-31G(d) freq" in text
    assert "1 2" in text


def test_export_nwchem_input() -> None:
    graph, coords = _water()

    text = export_nwchem_input(graph, coords, charge=-1, multiplicity=2, calculation="sp")

    assert "memory 2000 mb" in text
    assert "charge -1" in text
    assert "mult 2" in text
    assert "task dft energy" in text
