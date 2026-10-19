from __future__ import annotations

import math
from collections.abc import Iterable, Sequence

from chemuson.core.model import MolGraph, bond_is_structural
from chemuson.geometry3d.model import (
    CoordinateSet3D,
    transformed_generated_hydrogen_metadata,
)


Vector3 = tuple[float, float, float]


def translate_atoms(
    coordinates: CoordinateSet3D,
    atom_ids: Iterable[int],
    delta: Vector3,
) -> CoordinateSet3D:
    """Traslada un conjunto de átomos sin mutar el snapshot de entrada."""
    positions = coordinates.normalized_positions()
    selected = {int(atom_id) for atom_id in atom_ids}
    offset = _finite_vector(delta, "delta")
    missing = selected.difference(positions)
    if missing:
        raise ValueError(f"missing_coordinates:{','.join(map(str, sorted(missing)))}")
    if not selected:
        return coordinates

    dx, dy, dz = offset
    for atom_id in selected:
        x, y, z = positions[atom_id]
        positions[atom_id] = (x + dx, y + dy, z + dz)
    metadata = dict(coordinates.metadata)
    metadata.update(
        transformed_generated_hydrogen_metadata(
            coordinates,
            selected,
            lambda point: (point[0] + dx, point[1] + dy, point[2] + dz),
        )
    )
    return CoordinateSet3D(
        positions=positions,
        source=coordinates.source,
        method=coordinates.method,
        energy=None,
        metadata=metadata,
    )


def _finite_vector(value: Iterable[float], name: str) -> Vector3:
    try:
        result = tuple(float(component) for component in value)
    except (TypeError, ValueError) as exc:
        raise ValueError(f"invalid_{name}") from exc
    if len(result) != 3 or not all(math.isfinite(component) for component in result):
        raise ValueError(f"invalid_{name}")
    return result[0], result[1], result[2]


def distance(coordinates: CoordinateSet3D, first_id: int, second_id: int) -> float:
    """Distancia euclídea entre dos átomos, en ångströms."""
    first, second = _points(coordinates, (first_id, second_id))
    return _norm(_subtract(second, first))


def angle_degrees(
    coordinates: CoordinateSet3D,
    first_id: int,
    center_id: int,
    last_id: int,
) -> float:
    """Ángulo `first-center-last` entre 0 y 180 grados."""
    first, center, last = _points(coordinates, (first_id, center_id, last_id))
    first_vector = _subtract(first, center)
    last_vector = _subtract(last, center)
    denominator = _norm(first_vector) * _norm(last_vector)
    if denominator <= 1e-12:
        raise ValueError("undefined_angle")
    cosine = max(-1.0, min(1.0, _dot(first_vector, last_vector) / denominator))
    return math.degrees(math.acos(cosine))


def dihedral_degrees(
    coordinates: CoordinateSet3D,
    first_id: int,
    second_id: int,
    third_id: int,
    fourth_id: int,
) -> float:
    """Dihedro firmado `first-second-third-fourth`, normalizado a [-180, 180]."""
    first, second, third, fourth = _points(
        coordinates, (first_id, second_id, third_id, fourth_id)
    )
    axis = _normalize(_subtract(third, second), "undefined_dihedral")
    first_vector = _subtract(first, second)
    fourth_vector = _subtract(fourth, third)
    first_projected = _subtract(first_vector, _scale(axis, _dot(first_vector, axis)))
    fourth_projected = _subtract(fourth_vector, _scale(axis, _dot(fourth_vector, axis)))
    first_projected = _normalize(first_projected, "undefined_dihedral")
    fourth_projected = _normalize(fourth_projected, "undefined_dihedral")
    x = _dot(first_projected, fourth_projected)
    y = _dot(_cross(axis, first_projected), fourth_projected)
    return _normalize_degrees(math.degrees(math.atan2(y, x)))


def connected_component(
    graph: MolGraph,
    start_atom_id: int,
    excluded_bond_id: int | None = None,
) -> frozenset[int]:
    """Componente estructural alcanzable al omitir opcionalmente un enlace."""
    start = int(start_atom_id)
    if start not in graph.atoms:
        raise ValueError("missing_atom")
    adjacency = {atom_id: set() for atom_id in graph.atoms}
    for bond in graph.bonds.values():
        if bond.id == excluded_bond_id or not bond_is_structural(bond):
            continue
        if bond.a1_id not in adjacency or bond.a2_id not in adjacency:
            continue
        adjacency[bond.a1_id].add(bond.a2_id)
        adjacency[bond.a2_id].add(bond.a1_id)

    visited = {start}
    pending = [start]
    while pending:
        current = pending.pop()
        for neighbor in adjacency[current].difference(visited):
            visited.add(neighbor)
            pending.append(neighbor)
    return frozenset(visited)


def bond_is_bridge(graph: MolGraph, bond_id: int) -> bool:
    """True si eliminar el enlace separa sus extremos en el grafo estructural."""
    bond = _bond(graph, bond_id)
    return bond.a2_id not in connected_component(graph, bond.a1_id, bond.id)


def set_bond_length(
    graph: MolGraph,
    coordinates: CoordinateSet3D,
    bond_id: int,
    target_length: float,
    moving_endpoint_id: int | None = None,
) -> CoordinateSet3D:
    """Ajusta una distancia trasladando rígidamente el lado móvil de un puente."""
    bond = _editable_bridge(graph, bond_id)
    target = _positive_finite(target_length, "invalid_bond_length")
    moving = _moving_endpoint(bond.a1_id, bond.a2_id, moving_endpoint_id)
    fixed = bond.a1_id if moving == bond.a2_id else bond.a2_id
    moving_point, fixed_point = _points(coordinates, (moving, fixed))
    axis = _normalize(_subtract(moving_point, fixed_point), "undefined_bond_axis")
    current = _norm(_subtract(moving_point, fixed_point))
    shift = _scale(axis, target - current)
    fragment = connected_component(graph, moving, bond.id)
    return translate_atoms(coordinates, fragment, shift)


def set_bond_angle(
    graph: MolGraph,
    coordinates: CoordinateSet3D,
    atom_ids: Sequence[int],
    target_degrees: float,
    moving_terminal_id: int | None = None,
) -> CoordinateSet3D:
    """Ajusta un ángulo moviendo rígidamente el lado de un enlace terminal."""
    if len(atom_ids) != 3:
        raise ValueError("angle_requires_three_atoms")
    first_id, center_id, last_id = (int(atom_id) for atom_id in atom_ids)
    target = _bounded_finite(target_degrees, 0.0, 180.0, "invalid_angle")
    moving = _moving_endpoint(first_id, last_id, moving_terminal_id or last_id)
    if moving == center_id:
        raise ValueError("invalid_moving_endpoint")
    _require_bond(graph, first_id, center_id)
    _require_bond(graph, center_id, last_id)
    moving_bond = _require_editable_single_bond(graph, center_id, moving)
    if not bond_is_bridge(graph, moving_bond.id):
        raise ValueError("bond_not_bridge")

    first, center, last = _points(coordinates, (first_id, center_id, last_id))
    current = angle_degrees(coordinates, first_id, center_id, last_id)
    reference = first if moving == last_id else last
    moving_point = last if moving == last_id else first
    normal = _cross(_subtract(reference, center), _subtract(moving_point, center))
    if _norm(normal) <= 1e-12:
        normal = _orthogonal_axis(_subtract(moving_point, center))
    fragment = connected_component(graph, moving, moving_bond.id)
    delta = target - current
    plus = rotate_atoms_about_axis(coordinates, fragment, center, _add(center, normal), delta)
    minus = rotate_atoms_about_axis(coordinates, fragment, center, _add(center, normal), -delta)
    plus_error = abs(angle_degrees(plus, first_id, center_id, last_id) - target)
    minus_error = abs(angle_degrees(minus, first_id, center_id, last_id) - target)
    return plus if plus_error <= minus_error else minus


def set_dihedral(
    graph: MolGraph,
    coordinates: CoordinateSet3D,
    atom_ids: Sequence[int],
    target_degrees: float,
    moving_endpoint_id: int | None = None,
) -> CoordinateSet3D:
    """Ajusta el dihedro rotando rígidamente un lado de un enlace puente simple."""
    if len(atom_ids) != 4:
        raise ValueError("dihedral_requires_four_atoms")
    first_id, second_id, third_id, fourth_id = (int(atom_id) for atom_id in atom_ids)
    target = _normalize_degrees(_finite_scalar(target_degrees, "invalid_dihedral"))
    central_bond = _require_editable_single_bond(graph, second_id, third_id)
    _require_bond(graph, first_id, second_id)
    _require_bond(graph, third_id, fourth_id)
    if not bond_is_bridge(graph, central_bond.id):
        raise ValueError("bond_not_bridge")
    moving = _moving_endpoint(second_id, third_id, moving_endpoint_id or third_id)
    fragment = connected_component(graph, moving, central_bond.id)
    current = dihedral_degrees(coordinates, first_id, second_id, third_id, fourth_id)
    delta = _normalize_degrees(target - current)
    second, third = _points(coordinates, (second_id, third_id))
    axis_start, axis_end = (second, third) if moving == third_id else (third, second)
    plus = rotate_atoms_about_axis(coordinates, fragment, axis_start, axis_end, delta)
    minus = rotate_atoms_about_axis(coordinates, fragment, axis_start, axis_end, -delta)
    plus_error = abs(_normalize_degrees(
        dihedral_degrees(plus, first_id, second_id, third_id, fourth_id) - target
    ))
    minus_error = abs(_normalize_degrees(
        dihedral_degrees(minus, first_id, second_id, third_id, fourth_id) - target
    ))
    return plus if plus_error <= minus_error else minus


def rotate_atoms_about_axis(
    coordinates: CoordinateSet3D,
    atom_ids: Iterable[int],
    axis_start: Vector3,
    axis_end: Vector3,
    angle_degrees_value: float,
) -> CoordinateSet3D:
    """Rota átomos alrededor de un eje mediante la fórmula de Rodrigues."""
    origin = _finite_vector(axis_start, "axis_start")
    endpoint = _finite_vector(axis_end, "axis_end")
    axis = _normalize(_subtract(endpoint, origin), "undefined_rotation_axis")
    angle = math.radians(_finite_scalar(angle_degrees_value, "invalid_rotation"))
    selected = {int(atom_id) for atom_id in atom_ids}
    positions = coordinates.normalized_positions()
    missing = selected.difference(positions)
    if missing:
        raise ValueError(f"missing_coordinates:{','.join(map(str, sorted(missing)))}")
    cosine = math.cos(angle)
    sine = math.sin(angle)

    def rotate_point(point: Vector3) -> Vector3:
        relative = _subtract(point, origin)
        rotated = _add(
            _add(_scale(relative, cosine), _scale(_cross(axis, relative), sine)),
            _scale(axis, _dot(axis, relative) * (1.0 - cosine)),
        )
        return _add(origin, rotated)

    for atom_id in selected:
        positions[atom_id] = rotate_point(positions[atom_id])
    metadata = dict(coordinates.metadata)
    metadata.update(
        transformed_generated_hydrogen_metadata(coordinates, selected, rotate_point)
    )
    return CoordinateSet3D(
        positions=positions,
        source=coordinates.source,
        method=coordinates.method,
        energy=None,
        metadata=metadata,
    )


def _editable_bridge(graph: MolGraph, bond_id: int):
    bond = _require_editable_single_bond_by_id(graph, bond_id)
    if not bond_is_bridge(graph, bond.id):
        raise ValueError("bond_not_bridge")
    return bond


def _require_editable_single_bond_by_id(graph: MolGraph, bond_id: int):
    bond = _bond(graph, bond_id)
    if (
        bond.order != 1
        or bond.is_aromatic
        or bond.is_query
        or not bond_is_structural(bond)
    ):
        raise ValueError("bond_not_editable")
    return bond


def _require_editable_single_bond(graph: MolGraph, first_id: int, second_id: int):
    bond = _require_bond(graph, first_id, second_id)
    return _require_editable_single_bond_by_id(graph, bond.id)


def _require_bond(graph: MolGraph, first_id: int, second_id: int):
    for bond in graph.bonds.values():
        if {bond.a1_id, bond.a2_id} == {int(first_id), int(second_id)}:
            if bond_is_structural(bond):
                return bond
    raise ValueError("missing_bond")


def _bond(graph: MolGraph, bond_id: int):
    try:
        return graph.bonds[int(bond_id)]
    except (KeyError, TypeError, ValueError) as exc:
        raise ValueError("missing_bond") from exc


def _moving_endpoint(first_id: int, second_id: int, requested: int | None) -> int:
    if requested is None:
        return int(second_id)
    moving = int(requested)
    if moving not in {int(first_id), int(second_id)}:
        raise ValueError("invalid_moving_endpoint")
    return moving


def _points(coordinates: CoordinateSet3D, atom_ids: Iterable[int]) -> tuple[Vector3, ...]:
    positions = coordinates.normalized_positions()
    points: list[Vector3] = []
    for atom_id in atom_ids:
        try:
            points.append(_finite_vector(positions[int(atom_id)], "coordinate"))
        except KeyError as exc:
            raise ValueError(f"missing_coordinate:{atom_id}") from exc
    return tuple(points)


def _finite_scalar(value: float, name: str) -> float:
    try:
        result = float(value)
    except (TypeError, ValueError) as exc:
        raise ValueError(name) from exc
    if not math.isfinite(result):
        raise ValueError(name)
    return result


def _positive_finite(value: float, name: str) -> float:
    result = _finite_scalar(value, name)
    if result <= 0.0:
        raise ValueError(name)
    return result


def _bounded_finite(value: float, lower: float, upper: float, name: str) -> float:
    result = _finite_scalar(value, name)
    if not lower < result < upper:
        raise ValueError(name)
    return result


def _normalize_degrees(value: float) -> float:
    wrapped = (float(value) + 180.0) % 360.0 - 180.0
    return 180.0 if wrapped == -180.0 and value > 0 else wrapped


def _add(first: Vector3, second: Vector3) -> Vector3:
    return first[0] + second[0], first[1] + second[1], first[2] + second[2]


def _subtract(first: Vector3, second: Vector3) -> Vector3:
    return first[0] - second[0], first[1] - second[1], first[2] - second[2]


def _scale(vector: Vector3, factor: float) -> Vector3:
    return vector[0] * factor, vector[1] * factor, vector[2] * factor


def _dot(first: Vector3, second: Vector3) -> float:
    return first[0] * second[0] + first[1] * second[1] + first[2] * second[2]


def _cross(first: Vector3, second: Vector3) -> Vector3:
    return (
        first[1] * second[2] - first[2] * second[1],
        first[2] * second[0] - first[0] * second[2],
        first[0] * second[1] - first[1] * second[0],
    )


def _norm(vector: Vector3) -> float:
    return math.sqrt(_dot(vector, vector))


def _normalize(vector: Vector3, error: str) -> Vector3:
    length = _norm(vector)
    if length <= 1e-12:
        raise ValueError(error)
    return _scale(vector, 1.0 / length)


def _orthogonal_axis(vector: Vector3) -> Vector3:
    unit = _normalize(vector, "undefined_rotation_axis")
    candidate = (1.0, 0.0, 0.0) if abs(unit[0]) < 0.8 else (0.0, 1.0, 0.0)
    return _normalize(_cross(unit, candidate), "undefined_rotation_axis")
