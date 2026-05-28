from __future__ import annotations

"""Exportador de entradas ORCA."""

from chemuson.compchem.exporters.common import calculation_keyword, geometry_lines
from chemuson.core.model import MolGraph
from chemuson.geometry3d.model import CoordinateSet3D


def export_orca_input(
    graph: MolGraph,
    coordinates: CoordinateSet3D,
    *,
    charge: int = 0,
    multiplicity: int = 1,
    method: str = "B3LYP",
    basis: str = "def2-SVP",
    memory_mb: int = 2000,
    cores: int = 1,
    calculation: str = "sp",
) -> str:
    calc = calculation_keyword(calculation)
    lines = [
        f"! {method} {basis} {calc}",
        f"%pal nprocs {max(1, int(cores))} end",
        f"%maxcore {max(1, int(memory_mb))}",
        "",
        f"* xyz {int(charge)} {int(multiplicity)}",
    ]
    lines.extend(geometry_lines(graph, coordinates))
    lines.append("*")
    return "\n".join(lines) + "\n"
