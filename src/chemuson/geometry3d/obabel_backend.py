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
        positions, metadata, error = _run_obabel_mol(exe, graph, coordset, atom_ids, settings, Path(tmpdir))
        if error:
            fallback_positions, fallback_metadata, fallback_error = _run_obabel_xyz(exe, graph, coordset, atom_ids, settings, Path(tmpdir))
            metadata = {"mol_error": error, "mol_metadata": metadata, "xyz_metadata": fallback_metadata, "format": "xyz-fallback"}
            positions = fallback_positions
            error = fallback_error
        else:
            metadata["format"] = "mol"
    if not positions:
        return OptimizationResult(
            None,
            method="openbabel",
            message=error or "Open Babel no devolvió coordenadas válidas",
            metadata=metadata,
        )
    if len(positions) != len(atom_ids):
        return OptimizationResult(
            None,
            method="openbabel",
            message=f"Open Babel devolvió {len(positions)} átomos; se esperaban {len(atom_ids)}",
            metadata=metadata,
        )
    result_coords = CoordinateSet3D(
        positions=positions,
        source="openbabel",
        method=settings.forcefield.value,
        metadata={"steps": int(settings.max_iters), **metadata},
    )
    return OptimizationResult(result_coords, converged=True, method=settings.forcefield.value)


def _run_obabel_mol(
    exe: str,
    graph: MolGraph,
    coordset: CoordinateSet3D,
    atom_ids: list[int],
    settings: OptimizationSettings,
    tmpdir: Path,
) -> tuple[dict[int, tuple[float, float, float]], dict[str, str], str | None]:
    input_path = tmpdir / "input.mol"
    output_path = tmpdir / "output.mol"
    input_path.write_text(_molblock_with_coordinates(graph, coordset), encoding="utf-8")
    proc, error = _run_obabel_command(exe, input_path, output_path, settings)
    metadata = _process_metadata(proc)
    if error:
        return {}, metadata, error
    return _positions_from_molblock(output_path.read_text(encoding="utf-8"), atom_ids), metadata, None


def _run_obabel_xyz(
    exe: str,
    graph: MolGraph,
    coordset: CoordinateSet3D,
    atom_ids: list[int],
    settings: OptimizationSettings,
    tmpdir: Path,
) -> tuple[dict[int, tuple[float, float, float]], dict[str, str], str | None]:
    input_path = tmpdir / "input.xyz"
    output_path = tmpdir / "output.xyz"
    input_path.write_text(molgraph_to_xyz(graph, coordset), encoding="utf-8")
    proc, error = _run_obabel_command(exe, input_path, output_path, settings)
    metadata = _process_metadata(proc)
    if error:
        return {}, metadata, error
    return xyz_to_positions(output_path.read_text(encoding="utf-8"), atom_ids), metadata, None


def _run_obabel_command(
    exe: str,
    input_path: Path,
    output_path: Path,
    settings: OptimizationSettings,
) -> tuple[subprocess.CompletedProcess[str] | None, str | None]:
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
        return None, "Open Babel timeout"
    except Exception as exc:
        return None, str(exc)
    if proc.returncode != 0:
        return proc, (proc.stderr or proc.stdout or "Open Babel failed").strip()
    if not output_path.exists():
        return proc, "Open Babel no creó archivo de salida"
    return proc, None


def _process_metadata(proc: subprocess.CompletedProcess[str] | None) -> dict[str, str]:
    if proc is None:
        return {}
    return {"stdout": (proc.stdout or "").strip(), "stderr": (proc.stderr or "").strip()}


def _molblock_with_coordinates(graph: MolGraph, coordset: CoordinateSet3D) -> str:
    positions = coordset.normalized_positions()
    atoms = [atom for atom in sorted(graph.atoms.values(), key=lambda item: item.id) if int(atom.id) in positions]
    atom_index = {int(atom.id): idx + 1 for idx, atom in enumerate(atoms)}
    bonds = [
        bond
        for bond in sorted(graph.bonds.values(), key=lambda item: item.id)
        if int(bond.a1_id) in atom_index and int(bond.a2_id) in atom_index
    ]
    lines = ["ChemUSON", "  ChemUSON 3D", "", f"{len(atoms):>3}{len(bonds):>3}  0  0  0  0            999 V2000"]
    charge_records: list[tuple[int, int]] = []
    for atom in atoms:
        x, y, z = positions[int(atom.id)]
        lines.append(f"{x:>10.4f}{y:>10.4f}{z:>10.4f} {atom.element:<3} 0  0  0  0  0  0  0  0  0  0  0  0")
        charge = int(getattr(atom, "charge", 0) or 0)
        if charge:
            charge_records.append((atom_index[int(atom.id)], charge))
    for bond in bonds:
        order = 4 if bool(getattr(bond, "is_aromatic", False)) else max(1, min(3, int(getattr(bond, "order", 1) or 1)))
        lines.append(f"{atom_index[int(bond.a1_id)]:>3}{atom_index[int(bond.a2_id)]:>3}{order:>3}  0  0  0  0")
    for chunk_start in range(0, len(charge_records), 8):
        chunk = charge_records[chunk_start : chunk_start + 8]
        line = f"M  CHG{len(chunk):>3}"
        for atom_idx, charge in chunk:
            line += f"{atom_idx:>4}{charge:>4}"
        lines.append(line)
    lines.append("M  END")
    return "\n".join(lines) + "\n"


def _positions_from_molblock(molblock: str, atom_ids: list[int]) -> dict[int, tuple[float, float, float]]:
    lines = str(molblock or "").splitlines()
    if len(lines) < 4:
        return {}
    try:
        atom_count = int(lines[3][0:3])
    except Exception:
        return {}
    positions: dict[int, tuple[float, float, float]] = {}
    for atom_id, line in zip(atom_ids, lines[4 : 4 + atom_count]):
        try:
            positions[int(atom_id)] = (float(line[0:10]), float(line[10:20]), float(line[20:30]))
        except Exception:
            continue
    return positions
