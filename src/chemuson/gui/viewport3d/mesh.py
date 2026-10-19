"""Generación pura de mallas de esferas y varillas."""

from __future__ import annotations

from dataclasses import dataclass
import math
from typing import Any, Sequence

from .geometry import (
    ATOM_RADIUS_SCALE,
    BOND_RADIUS,
    Vec3,
    add,
    cross,
    dot,
    length,
    normalize,
    scale,
    sub,
)


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


def bond_lateral_direction(
    start: Vec3,
    end: Vec3,
    start_id: int | str,
    end_id: int | str,
    atoms: dict[int | str, Any],
    adjacency: dict[int | str, list[int | str]],
) -> Vec3:
    axis = normalize(sub(end, start))
    for anchor_id in (start_id, end_id):
        anchor = atoms[anchor_id]
        anchor_position = (float(anchor.x), float(anchor.y), float(anchor.z))
        for neighbor_id in adjacency.get(anchor_id, ()):
            if neighbor_id in (start_id, end_id):
                continue
            neighbor = atoms.get(neighbor_id)
            if neighbor is None:
                continue
            neighbor_position = (float(neighbor.x), float(neighbor.y), float(neighbor.z))
            neighbor_vector = sub(neighbor_position, anchor_position)
            projected = sub(neighbor_vector, scale(axis, dot(neighbor_vector, axis)))
            if length(projected) > 1.0e-9:
                return normalize(projected)
    reference = (0.0, 1.0, 0.0) if abs(axis[1]) < 0.9 else (1.0, 0.0, 0.0)
    return normalize(cross(axis, reference), (1.0, 0.0, 0.0))


def normalize_representation(value: str) -> str:
    representation = str(value).strip().lower()
    if representation not in {"ball_and_stick", "wireframe", "space_filling"}:
        raise ValueError(f"Representación 3D no soportada: {value}")
    return representation


def bond_lane_layout(
    order: int,
    is_aromatic: bool = False,
    representation: str = "ball_and_stick",
) -> tuple[tuple[float, ...], float]:
    """Devuelve carriles coherentes con el estilo usado por malla y picking."""
    style = normalize_representation(representation)
    if style == "space_filling":
        return (0.0,), 0.045
    if style == "wireframe":
        offsets, _radius = bond_lane_layout(order, is_aromatic)
        return offsets, 0.035
    try:
        normalized_order = max(1, min(3, int(order or 1)))
    except (TypeError, ValueError):
        normalized_order = 1
    if is_aromatic:
        return (-BOND_RADIUS * 0.75, BOND_RADIUS * 0.75), BOND_RADIUS * 0.58
    if normalized_order == 1:
        return (0.0,), BOND_RADIUS
    if normalized_order == 2:
        return (-BOND_RADIUS * 1.6, BOND_RADIUS * 1.6), BOND_RADIUS
    return (-BOND_RADIUS * 2.0, 0.0, BOND_RADIUS * 2.0), BOND_RADIUS * 0.82


def _bond_surface_trim(atom_radius: float, lane_offset: float, bond_radius: float) -> float:
    """Trim a bond lane to where its outer edge meets the atom sphere."""
    radial_extent = abs(lane_offset) + bond_radius
    if (
        atom_radius <= radial_extent
        or atom_radius <= 0.0
        or not math.isfinite(atom_radius)
        or not math.isfinite(radial_extent)
    ):
        return 0.0
    return math.sqrt(max(0.0, atom_radius * atom_radius - radial_extent * radial_extent))


def _append_bond(
    output: list[float],
    start: Vec3,
    end: Vec3,
    lateral: Vec3,
    order: int,
    is_aromatic: bool,
    color: tuple[float, float, float],
    representation: str = "ball_and_stick",
    start_atom_radius: float = 0.0,
    end_atom_radius: float = 0.0,
) -> None:
    offsets, radius = bond_lane_layout(order, is_aromatic, representation)
    axis = sub(end, start)
    axis_length = length(axis)
    if axis_length <= 1.0e-9 or not math.isfinite(axis_length):
        return
    direction = scale(axis, 1.0 / axis_length)
    for offset in offsets:
        displacement = scale(lateral, offset)
        start_trim = 0.0
        end_trim = 0.0
        if representation == "ball_and_stick":
            start_trim = _bond_surface_trim(start_atom_radius, offset, radius)
            end_trim = _bond_surface_trim(end_atom_radius, offset, radius)
        if start_trim + end_trim >= axis_length:
            continue
        lane_start = add(add(start, scale(direction, start_trim)), displacement)
        lane_end = add(sub(end, scale(direction, end_trim)), displacement)
        _append_cylinder(output, lane_start, lane_end, radius, color)


def build_scene_mesh(
    scene,
    selected_atom_id: int | str | None = None,
    selected_bond_id: int | str | None = None,
    representation: str = "ball_and_stick",
    selected_atom_ids: Sequence[int] = (),
) -> MeshData:
    style = normalize_representation(representation)
    selected_ids = set(selected_atom_ids)
    output: list[float] = []
    atoms = {atom.id: atom for atom in getattr(scene, "atoms", ())}
    bonds = tuple(getattr(scene, "bonds", ()))
    adjacency: dict[int | str, list[int | str]] = {atom_id: [] for atom_id in atoms}
    for bond in bonds:
        if bond.a1_id in adjacency and bond.a2_id in adjacency:
            adjacency[bond.a1_id].append(bond.a2_id)
            adjacency[bond.a2_id].append(bond.a1_id)
    for atom in getattr(scene, "atoms", ()):
        center = (float(atom.x), float(atom.y), float(atom.z))
        if style == "wireframe":
            radius = 0.06
        elif style == "space_filling":
            radius = float(getattr(atom, "radius", 0.0))
        else:
            radius = float(getattr(atom, "radius", 0.0)) * ATOM_RADIUS_SCALE
        is_selected = (
            atom.id == selected_atom_id
            or atom.id in selected_ids
            or (
                bool(getattr(atom, "generated", False))
                and getattr(atom, "parent_atom_id", None) in selected_ids
            )
        )
        color = _selected_color() if is_selected else _color(atom.color)
        _append_sphere(output, center, radius, color)
    for bond in bonds:
        start_atom = atoms.get(bond.a1_id)
        end_atom = atoms.get(bond.a2_id)
        if start_atom is None or end_atom is None:
            continue
        start = (float(start_atom.x), float(start_atom.y), float(start_atom.z))
        end = (float(end_atom.x), float(end_atom.y), float(end_atom.z))
        color = _selected_color() if bond.id == selected_bond_id else (0.72, 0.72, 0.76)
        lateral = bond_lateral_direction(start, end, bond.a1_id, bond.a2_id, atoms, adjacency)
        try:
            order = max(1, min(3, int(getattr(bond, "order", 1) or 1)))
        except (TypeError, ValueError):
            order = 1
        _append_bond(
            output,
            start,
            end,
            lateral,
            order,
            bool(getattr(bond, "is_aromatic", False)),
            color,
            style,
            float(getattr(start_atom, "radius", 0.0)) * ATOM_RADIUS_SCALE
            if style == "ball_and_stick"
            else 0.0,
            float(getattr(end_atom, "radius", 0.0)) * ATOM_RADIUS_SCALE
            if style == "ball_and_stick"
            else 0.0,
        )
    return MeshData(tuple(output))
