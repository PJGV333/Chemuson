from __future__ import annotations

"""Exportador de entradas Gaussian."""

from chemuson.compchem.exporters.common import calculation_keyword, geometry_lines
from chemuson.core.model import MolGraph
from chemuson.geometry3d.model import CoordinateSet3D


def export_gaussian_input(
    graph: MolGraph,
    coordinates: CoordinateSet3D,
    *,
    charge: int = 0,
    multiplicity: int = 1,
    method: str = "B3LYP",
    basis: str = "6-31G(d)",
    memory_mb: int = 2000,
    cores: int = 1,
    calculation: str = "sp",
) -> str:
    calc = calculation_keyword(calculation)
    route_calc = "" if calc == "sp" else f" {calc}"
    lines = [
        f"%mem={max(1, int(memory_mb))}MB",
        f"%nprocshared={max(1, int(cores))}",
        f"# {method}/{basis}{route_calc}",
        "",
        "ChemUSON generated input",
        "",
        f"{int(charge)} {int(multiplicity)}",
    ]
    lines.extend(geometry_lines(graph, coordinates))
    lines.extend(["", ""])
    return "\n".join(lines)
