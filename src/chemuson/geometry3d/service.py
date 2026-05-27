from __future__ import annotations

"""Pipeline híbrido para conformaciones 3D reales y proyección 2D."""

from dataclasses import dataclass
from concurrent.futures import Future, ThreadPoolExecutor
import hashlib
import math
from typing import Iterable

from chemuson.core.model import MolGraph, bond_affects_valence


@dataclass(frozen=True)
class Conformer3DResult:
    """Resultado trazable de una generación de conformación 3D."""

    atom_positions: dict[int, tuple[float, float, float]]
    source: str
    cache_key: str
    method: str = "rdkit"
    energy: float | None = None
    message: str = ""
    from_cache: bool = False

    @property
    def ok(self) -> bool:
        return bool(self.atom_positions)


@dataclass(frozen=True)
class Rotation3D:
    """Rotación trackball en radianes, sin límite artificial de tilt."""

    pitch: float = 0.0
    yaw: float = 0.0
    roll: float = 0.0


@dataclass(frozen=True)
class DepthCue:
    """Estilo derivado de profundidad para renderizado posterior."""

    z: float
    normalized_z: float
    opacity: float
    stroke_scale: float


@dataclass(frozen=True)
class ProjectedAtom3D:
    """Átomo 3D proyectado a coordenadas 2D con depth cueing."""

    atom_id: int
    x: float
    y: float
    cue: DepthCue


_CONFORMER_CACHE: dict[str, Conformer3DResult] = {}
_MAX_CACHE_SIZE = 128
_EXECUTOR = ThreadPoolExecutor(max_workers=2, thread_name_prefix="chemuson-3d")


def conformer_3d_for_graph(
    graph: MolGraph,
    *,
    timeout_s: float = 8.0,
    force: bool = False,
) -> Conformer3DResult:
    """Genera coordenadas 3D reales mediante RDKit con cache por hash químico."""
    cache_key = chemical_graph_hash(graph)
    if not force:
        cached = _CONFORMER_CACHE.get(cache_key)
        if cached is not None:
            return Conformer3DResult(
                atom_positions=dict(cached.atom_positions),
                source=cached.source,
                cache_key=cached.cache_key,
                method=cached.method,
                energy=cached.energy,
                message=cached.message,
                from_cache=True,
            )

    try:
        from chemuson.chemio.rdkit_safe import conformer_3d_isolated

        positions, metadata, error = conformer_3d_isolated(graph, timeout_s=timeout_s)
    except Exception as exc:
        return Conformer3DResult({}, "rdkit", cache_key, message=str(exc))

    if error or not positions:
        return Conformer3DResult({}, "rdkit", cache_key, message=str(error or "empty_conformer"))

    result = Conformer3DResult(
        atom_positions=dict(positions),
        source=str(metadata.get("source", "rdkit")),
        cache_key=cache_key,
        method=str(metadata.get("method", "rdkit")),
        energy=_float_or_none(metadata.get("energy")),
    )
    _store_cache(cache_key, result)
    return result


def conformer_3d_for_graph_async(
    graph: MolGraph,
    *,
    timeout_s: float = 8.0,
    force: bool = False,
) -> Future[Conformer3DResult]:
    """Programa generación 3D fuera del hilo de UI y devuelve un ``Future``."""
    return _EXECUTOR.submit(
        conformer_3d_for_graph,
        graph,
        timeout_s=timeout_s,
        force=force,
    )


def project_conformer_to_2d(
    atom_positions: dict[int, tuple[float, float, float]],
    *,
    rotation: Rotation3D | None = None,
    center: tuple[float, float] = (0.0, 0.0),
    scale: float = 42.0,
    depth_opacity_range: tuple[float, float] = (0.48, 1.0),
    depth_stroke_range: tuple[float, float] = (0.72, 1.22),
) -> list[ProjectedAtom3D]:
    """Proyecta una conformación 3D a 2D con opacidad/ancho dependientes de Z."""
    rotation = rotation or Rotation3D()
    if not atom_positions:
        return []

    centroid = _centroid3(atom_positions.values())
    rotated: list[tuple[int, float, float, float]] = []
    for atom_id, point in sorted(atom_positions.items()):
        x, y, z = _rotate_point3(
            point[0] - centroid[0],
            point[1] - centroid[1],
            point[2] - centroid[2],
            rotation,
        )
        rotated.append((int(atom_id), x, y, z))

    z_values = [z for _atom_id, _x, _y, z in rotated]
    z_min = min(z_values)
    z_span = max(z_values) - z_min
    opacity_min, opacity_max = depth_opacity_range
    stroke_min, stroke_max = depth_stroke_range
    cx, cy = center

    projected: list[ProjectedAtom3D] = []
    for atom_id, x, y, z in rotated:
        zn = 0.5 if z_span <= 1e-9 else (z - z_min) / z_span
        cue = DepthCue(
            z=z,
            normalized_z=zn,
            opacity=_lerp(opacity_min, opacity_max, zn),
            stroke_scale=_lerp(stroke_min, stroke_max, zn),
        )
        projected.append(ProjectedAtom3D(atom_id, cx + x * scale, cy - y * scale, cue))
    return projected


def chemical_graph_hash(graph: MolGraph) -> str:
    """Hash estable de conectividad y atributos químicos principales."""
    parts: list[str] = []
    for atom in sorted(graph.atoms.values(), key=lambda item: item.id):
        parts.append(
            "A"
            f"{int(atom.id)}:{atom.element}:{int(getattr(atom, 'charge', 0) or 0)}:"
            f"{getattr(atom, 'isotope', None) or ''}:{int(getattr(atom, 'radical_electrons', 0) or 0)}"
        )
    for bond in sorted(graph.bonds.values(), key=lambda item: item.id):
        if not bond_affects_valence(bond):
            continue
        a1, a2 = sorted((int(bond.a1_id), int(bond.a2_id)))
        parts.append(
            "B"
            f"{a1}-{a2}:{int(getattr(bond, 'order', 1) or 1)}:"
            f"{int(bool(getattr(bond, 'is_aromatic', False)))}"
        )
    return hashlib.sha256("|".join(parts).encode("utf-8")).hexdigest()


def _store_cache(cache_key: str, result: Conformer3DResult) -> None:
    if len(_CONFORMER_CACHE) >= _MAX_CACHE_SIZE:
        oldest = next(iter(_CONFORMER_CACHE), None)
        if oldest is not None:
            _CONFORMER_CACHE.pop(oldest, None)
    _CONFORMER_CACHE[cache_key] = result


def _rotate_point3(x: float, y: float, z: float, rotation: Rotation3D) -> tuple[float, float, float]:
    cx = math.cos(rotation.pitch)
    sx = math.sin(rotation.pitch)
    cy = math.cos(rotation.yaw)
    sy = math.sin(rotation.yaw)
    cz = math.cos(rotation.roll)
    sz = math.sin(rotation.roll)

    y, z = y * cx - z * sx, y * sx + z * cx
    x, z = x * cy + z * sy, -x * sy + z * cy
    x, y = x * cz - y * sz, x * sz + y * cz
    return x, y, z


def _centroid3(points: Iterable[tuple[float, float, float]]) -> tuple[float, float, float]:
    values = list(points)
    count = len(values) or 1
    return (
        sum(x for x, _y, _z in values) / count,
        sum(y for _x, y, _z in values) / count,
        sum(z for _x, _y, z in values) / count,
    )


def _lerp(a: float, b: float, t: float) -> float:
    return float(a) + (float(b) - float(a)) * max(0.0, min(1.0, float(t)))


def _float_or_none(value: object) -> float | None:
    try:
        return float(value)
    except Exception:
        return None
