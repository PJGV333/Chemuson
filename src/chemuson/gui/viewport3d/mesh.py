"""Generación pura de mallas de esferas y varillas."""

from __future__ import annotations

from dataclasses import dataclass
import math

from .geometry import ATOM_RADIUS_SCALE, BOND_RADIUS, Vec3, cross, length, normalize, sub


@dataclass(frozen=True, slots=True)
class MeshData:
    vertices: tuple[float, ...]

    @property
    def vertex_count(self) -> int:
        return len(self.vertices) // 9


def _color(value: str) -> tuple[float, float, float]:
    text = str(value or "#B0B0B0").lstrip("#")
    if len(text) != 6:
        return (0.69, 0.69, 0.69)
    try:
        return tuple(int(text[index : index + 2], 16) / 255.0 for index in (0, 2, 4))  # type: ignore[return-value]
    except ValueError:
        return (0.69, 0.69, 0.69)


def _selected_color() -> tuple[float, float, float]:
    return (1.0, 0.78, 0.12)


def _append_vertex(output: list[float], point: Vec3, normal: Vec3, color: tuple[float, float, float]) -> None:
    output.extend((*point, *normal, *color))


def _append_triangle(
    output: list[float],
    p1: Vec3,
    n1: Vec3,
    p2: Vec3,
    n2: Vec3,
    p3: Vec3,
    n3: Vec3,
    color: tuple[float, float, float],
) -> None:
    _append_vertex(output, p1, n1, color)
    _append_vertex(output, p2, n2, color)
    _append_vertex(output, p3, n3, color)


def _append_sphere(
    output: list[float],
    center: Vec3,
    radius: float,
    color: tuple[float, float, float],
    segments: int = 20,
    rings: int = 10,
) -> None:
    if radius <= 1.0e-9 or not math.isfinite(radius):
        return
    for ring in range(rings):
        theta_1 = math.pi * ring / rings
        theta_2 = math.pi * (ring + 1) / rings
        for segment in range(segments):
            phi_1 = 2.0 * math.pi * segment / segments
            phi_2 = 2.0 * math.pi * (segment + 1) / segments
            normals = (
                (math.sin(theta_1) * math.cos(phi_1), math.cos(theta_1), math.sin(theta_1) * math.sin(phi_1)),
                (math.sin(theta_2) * math.cos(phi_1), math.cos(theta_2), math.sin(theta_2) * math.sin(phi_1)),
                (math.sin(theta_2) * math.cos(phi_2), math.cos(theta_2), math.sin(theta_2) * math.sin(phi_2)),
                (math.sin(theta_1) * math.cos(phi_2), math.cos(theta_1), math.sin(theta_1) * math.sin(phi_2)),
            )
            points: tuple[Vec3, Vec3, Vec3, Vec3] = tuple(
                tuple(center[index] + radius * normal[index] for index in range(3))  # type: ignore[misc]
                for normal in normals
            )  # type: ignore[assignment]
            _append_triangle(output, points[0], normals[0], points[1], normals[1], points[2], normals[2], color)
            _append_triangle(output, points[0], normals[0], points[2], normals[2], points[3], normals[3], color)


def _append_cylinder(
    output: list[float],
    start: Vec3,
    end: Vec3,
    radius: float,
    color: tuple[float, float, float],
    segments: int = 12,
) -> None:
    axis = sub(end, start)
    axis_length = length(axis)
    if axis_length <= 1.0e-9 or not math.isfinite(axis_length):
        return
    direction = normalize(axis)
    reference = (0.0, 1.0, 0.0) if abs(direction[1]) < 0.9 else (1.0, 0.0, 0.0)
    side = normalize(cross(direction, reference), (1.0, 0.0, 0.0))
    other = normalize(cross(direction, side), (0.0, 0.0, 1.0))
    for segment in range(segments):
        angle_1 = 2.0 * math.pi * segment / segments
        angle_2 = 2.0 * math.pi * (segment + 1) / segments
        radial_1: Vec3 = normalize((
            math.cos(angle_1) * side[0] + math.sin(angle_1) * other[0],
            math.cos(angle_1) * side[1] + math.sin(angle_1) * other[1],
            math.cos(angle_1) * side[2] + math.sin(angle_1) * other[2],
        ))
        radial_2: Vec3 = normalize((
            math.cos(angle_2) * side[0] + math.sin(angle_2) * other[0],
            math.cos(angle_2) * side[1] + math.sin(angle_2) * other[1],
            math.cos(angle_2) * side[2] + math.sin(angle_2) * other[2],
        ))
        start_1: Vec3 = (
            start[0] + radius * radial_1[0],
            start[1] + radius * radial_1[1],
            start[2] + radius * radial_1[2],
        )
        start_2: Vec3 = (
            start[0] + radius * radial_2[0],
            start[1] + radius * radial_2[1],
            start[2] + radius * radial_2[2],
        )
        end_1: Vec3 = (
            end[0] + radius * radial_1[0],
            end[1] + radius * radial_1[1],
            end[2] + radius * radial_1[2],
        )
        end_2: Vec3 = (
            end[0] + radius * radial_2[0],
            end[1] + radius * radial_2[1],
            end[2] + radius * radial_2[2],
        )
        _append_triangle(output, start_1, radial_1, end_1, radial_1, end_2, radial_2, color)
        _append_triangle(output, start_1, radial_1, end_2, radial_2, start_2, radial_2, color)


def build_scene_mesh(scene, selected_atom_id: int | str | None = None, selected_bond_id: int | str | None = None) -> MeshData:
    output: list[float] = []
    atoms = {atom.id: atom for atom in getattr(scene, "atoms", ())}
    for atom in getattr(scene, "atoms", ()):
        center = (float(atom.x), float(atom.y), float(atom.z))
        radius = float(getattr(atom, "radius", 0.0)) * ATOM_RADIUS_SCALE
        color = _selected_color() if atom.id == selected_atom_id else _color(atom.color)
        _append_sphere(output, center, radius, color)
    for bond in getattr(scene, "bonds", ()):
        start_atom = atoms.get(bond.a1_id)
        end_atom = atoms.get(bond.a2_id)
        if start_atom is None or end_atom is None:
            continue
        start = (float(start_atom.x), float(start_atom.y), float(start_atom.z))
        end = (float(end_atom.x), float(end_atom.y), float(end_atom.z))
        color = _selected_color() if bond.id == selected_bond_id else (0.72, 0.72, 0.76)
        _append_cylinder(output, start, end, BOND_RADIUS, color)
    return MeshData(tuple(output))
