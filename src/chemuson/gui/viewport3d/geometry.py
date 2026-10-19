"""Geometría vectorial pura usada por cámara, mallas y picking."""

from __future__ import annotations

from dataclasses import dataclass
import math
from typing import Iterable

Vec3 = tuple[float, float, float]
EPSILON = 1.0e-9
ATOM_RADIUS_SCALE = 0.30
BOND_RADIUS = 0.12


def as_vec3(value: Iterable[float]) -> Vec3:
    values = tuple(float(component) for component in value)
    if len(values) != 3 or not all(math.isfinite(component) for component in values):
        raise ValueError("Se esperaba un vector 3D finito")
    return values  # type: ignore[return-value]


def add(a: Vec3, b: Vec3) -> Vec3:
    return (a[0] + b[0], a[1] + b[1], a[2] + b[2])


def sub(a: Vec3, b: Vec3) -> Vec3:
    return (a[0] - b[0], a[1] - b[1], a[2] - b[2])


def scale(vector: Vec3, factor: float) -> Vec3:
    return (vector[0] * factor, vector[1] * factor, vector[2] * factor)


def dot(a: Vec3, b: Vec3) -> float:
    return a[0] * b[0] + a[1] * b[1] + a[2] * b[2]


def cross(a: Vec3, b: Vec3) -> Vec3:
    return (
        a[1] * b[2] - a[2] * b[1],
        a[2] * b[0] - a[0] * b[2],
        a[0] * b[1] - a[1] * b[0],
    )


def length(vector: Vec3) -> float:
    return math.sqrt(dot(vector, vector))


def normalize(vector: Vec3, fallback: Vec3 = (0.0, 0.0, 1.0)) -> Vec3:
    magnitude = length(vector)
    if magnitude <= EPSILON or not math.isfinite(magnitude):
        return fallback
    return scale(vector, 1.0 / magnitude)


def is_finite_vec(vector: Vec3) -> bool:
    return all(math.isfinite(component) for component in vector)


@dataclass(frozen=True, slots=True)
class Bounds3D:
    minimum: Vec3
    maximum: Vec3

    @property
    def center(self) -> Vec3:
        return scale(add(self.minimum, self.maximum), 0.5)

    @property
    def radius(self) -> float:
        return 0.5 * length(sub(self.maximum, self.minimum))


@dataclass(frozen=True, slots=True)
class Ray3D:
    origin: Vec3
    direction: Vec3

    def __post_init__(self) -> None:
        if not is_finite_vec(self.origin) or not is_finite_vec(self.direction):
            raise ValueError("El rayo debe contener valores finitos")
        object.__setattr__(self, "direction", normalize(self.direction))


def scene_bounds(scene) -> Bounds3D | None:
    points: list[Vec3] = []
    for atom in getattr(scene, "atoms", ()):
        point = (float(atom.x), float(atom.y), float(atom.z))
        radius = max(0.0, float(getattr(atom, "radius", 0.0))) * ATOM_RADIUS_SCALE
        if not is_finite_vec(point) or not math.isfinite(radius):
            continue
        points.extend(
            (
                (point[0] - radius, point[1] - radius, point[2] - radius),
                (point[0] + radius, point[1] + radius, point[2] + radius),
            )
        )
    if not points:
        return None
    minimum = tuple(min(point[index] for point in points) for index in range(3))
    maximum = tuple(max(point[index] for point in points) for index in range(3))
    return Bounds3D(minimum=minimum, maximum=maximum)  # type: ignore[arg-type]


def ray_sphere_intersection(ray: Ray3D, center: Vec3, radius: float) -> float | None:
    roots = _ray_sphere_roots(ray, center, radius)
    return roots[0] if roots else None


def _ray_sphere_roots(ray: Ray3D, center: Vec3, radius: float) -> tuple[float, ...]:
    if radius <= EPSILON or not math.isfinite(radius) or not is_finite_vec(center):
        return ()
    offset = sub(ray.origin, center)
    half_b = dot(offset, ray.direction)
    discriminant = half_b * half_b - (dot(offset, offset) - radius * radius)
    if discriminant < 0.0:
        return ()
    root = math.sqrt(max(0.0, discriminant))
    near = -half_b - root
    far = -half_b + root
    return tuple(value for value in (near, far) if value >= 0.0)


def ray_capsule_intersection(ray: Ray3D, start: Vec3, end: Vec3, radius: float) -> float | None:
    """Devuelve la primera distancia positiva hasta la superficie de cápsula."""
    if radius <= EPSILON or not math.isfinite(radius) or not is_finite_vec(start) or not is_finite_vec(end):
        return None
    segment = sub(end, start)
    segment_length = length(segment)
    if segment_length <= EPSILON:
        return ray_sphere_intersection(ray, start, radius)

    axis = scale(segment, 1.0 / segment_length)
    relative = sub(ray.origin, start)
    axial_origin = dot(relative, axis)
    axial_direction = dot(ray.direction, axis)
    radial_origin = sub(relative, scale(axis, axial_origin))
    radial_direction = sub(ray.direction, scale(axis, axial_direction))

    candidates: list[float] = []
    quadratic_a = dot(radial_direction, radial_direction)
    quadratic_b = 2.0 * dot(radial_origin, radial_direction)
    quadratic_c = dot(radial_origin, radial_origin) - radius * radius
    if quadratic_a > 1.0e-14:
        discriminant = quadratic_b * quadratic_b - 4.0 * quadratic_a * quadratic_c
        tolerance = 1.0e-12 * max(1.0, quadratic_b * quadratic_b, abs(4.0 * quadratic_a * quadratic_c))
        if discriminant >= -tolerance:
            root = math.sqrt(max(0.0, discriminant))
            q = -0.5 * (quadratic_b + math.copysign(root, quadratic_b))
            if abs(q) > EPSILON:
                roots = (q / quadratic_a, quadratic_c / q)
            else:
                roots = (-quadratic_b / (2.0 * quadratic_a),)
            axial_tolerance = EPSILON * max(1.0, segment_length, radius)
            for distance in roots:
                if distance < 0.0 or not math.isfinite(distance):
                    continue
                axial = axial_origin + distance * axial_direction
                if -axial_tolerance <= axial <= segment_length + axial_tolerance:
                    candidates.append(distance)

    axial_tolerance = EPSILON * max(1.0, segment_length, radius)
    for endpoint, is_start in ((start, True), (end, False)):
        for distance in _ray_sphere_roots(ray, endpoint, radius):
            hit = add(ray.origin, scale(ray.direction, distance))
            axial = dot(sub(hit, start), axis)
            if (is_start and axial <= axial_tolerance) or (not is_start and axial >= segment_length - axial_tolerance):
                candidates.append(distance)

    return min(candidates) if candidates else None
