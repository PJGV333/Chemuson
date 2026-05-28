from __future__ import annotations

"""Exportador de entradas NWChem."""

from chemuson.compchem.exporters.common import calculation_keyword, geometry_lines
from chemuson.core.model import MolGraph
from chemuson.geometry3d.model import CoordinateSet3D


def export_nwchem_input(
    graph: MolGraph,
    coordinates: CoordinateSet3D,
    *,
    charge: int = 0,
    multiplicity: int = 1,
    method: str = "b3lyp",
    basis: str = "6-31g*",
    memory_mb: int = 2000,
    cores: int = 1,
    calculation: str = "sp",
) -> str:
    calc = calculation_keyword(calculation).split()[0]
    task_calc = "energy" if calc == "sp" else calc
    lines = [
        "start chemuson",
        f"memory {max(1, int(memory_mb))} mb",
        f"charge {int(charge)}",
        "geometry units angstroms noautosym",
    ]
    lines.extend(geometry_lines(graph, coordinates))
    lines.extend(
        [
            "end",
            "basis",
            f"  * library {basis}",
            "end",
            "dft",
            f"  xc {method}",
            f"  mult {int(multiplicity)}",
            "end",
            f"set nproc {max(1, int(cores))}",
            f"task dft {task_calc}",
        ]
    )
    return "\n".join(lines) + "\n"
