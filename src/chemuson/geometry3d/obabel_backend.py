from __future__ import annotations

"""Backend opcional Open Babel por subproceso."""

import shutil
import subprocess
import tempfile
from pathlib import Path

from chemuson.core.model import MolGraph
from chemuson.geometry3d.coordinates import planar_coordinate_set
from chemuson.geometry3d.export_xyz import molgraph_to_xyz, xyz_to_positions
from chemuson.geometry3d.model import CoordinateSet3D, ForceField, OptimizationResult, OptimizationSettings


SUPPORTED_FORCEFIELDS = (ForceField.GAFF, ForceField.GHEMICAL, ForceField.MMFF94, ForceField.MMFF94S, ForceField.UFF)


def obabel_path() -> str | None:
    """Devuelve ruta de `obabel` si está instalado."""
    return shutil.which("obabel")


def is_available() -> bool:
    return obabel_path() is not None


def supported_forcefields() -> tuple[ForceField, ...]:
    """Lista declarativa de force fields que Open Babel puede aceptar."""
    return SUPPORTED_FORCEFIELDS if is_available() else ()


def optimize(
    graph: MolGraph,
    coordset: CoordinateSet3D | None = None,
    settings: OptimizationSettings | None = None,
) -> OptimizationResult:
    """Optimiza con `obabel --minimize`; falla de forma controlada si no existe."""
    exe = obabel_path()
    if exe is None:
        return OptimizationResult(None, method="openbabel", message="Open Babel no disponible: ejecutable `obabel` no encontrado")
    settings = settings or OptimizationSettings(forcefield=ForceField.UFF)
    if settings.forcefield not in SUPPORTED_FORCEFIELDS:
        return OptimizationResult(None, method="openbabel", message=f"Force field no soportado por Open Babel: {settings.forcefield.value}")
    coordset = coordset or planar_coordinate_set(graph)
    atom_ids = [int(atom.id) for atom in sorted(graph.atoms.values(), key=lambda item: item.id)]
    with tempfile.TemporaryDirectory(prefix="chemuson-obabel-") as tmpdir:
        input_path = Path(tmpdir) / "input.xyz"
        output_path = Path(tmpdir) / "output.xyz"
        input_path.write_text(molgraph_to_xyz(graph, coordset), encoding="utf-8")
        cmd = [
            exe,
            str(input_path),
            "-O",
            str(output_path),
            "--minimize",
            "--ff",
            settings.forcefield.value,
            "--steps",
            str(max(1, int(settings.max_iters))),
        ]
        try:
            proc = subprocess.run(cmd, text=True, capture_output=True, timeout=max(1.0, settings.timeout_s), check=False)
        except subprocess.TimeoutExpired:
            return OptimizationResult(None, method="openbabel", message="Open Babel timeout")
        except Exception as exc:
            return OptimizationResult(None, method="openbabel", message=str(exc))
        if proc.returncode != 0:
            return OptimizationResult(None, method="openbabel", message=(proc.stderr or proc.stdout or "Open Babel failed").strip())
        positions = xyz_to_positions(output_path.read_text(encoding="utf-8"), atom_ids)
    if not positions:
        return OptimizationResult(None, method="openbabel", message="Open Babel no devolvió coordenadas XYZ válidas")
    result_coords = CoordinateSet3D(
        positions=positions,
        source="openbabel",
        method=settings.forcefield.value,
        metadata={"steps": int(settings.max_iters)},
    )
    return OptimizationResult(result_coords, converged=True, method=settings.forcefield.value)
