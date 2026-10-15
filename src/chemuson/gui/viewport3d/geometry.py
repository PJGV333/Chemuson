"""Geometría vectorial pura usada por cámara, mallas y picking."""

from __future__ import annotations

from dataclasses import dataclass
import math
from typing import Iterable

Vec3 = tuple[float, float, float]
EPSILON = 1.0e-9
ATOM_RADIUS_SCALE = 0.38
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
    if radius <= EPSILON or not math.isfinite(radius) or not is_finite_vec(center):
        return None
    offset = sub(ray.origin, center)
    half_b = dot(offset, ray.direction)
    discriminant = half_b * half_b - (dot(offset, offset) - radius * radius)
    if discriminant < 0.0:
        return None
    root = math.sqrt(max(0.0, discriminant))
    near = -half_b - root
    far = -half_b + root
    if near >= 0.0:
        return near
    if far >= 0.0:
        return far
    return None


def ray_capsule_intersection(ray: Ray3D, start: Vec3, end: Vec3, radius: float) -> float | None:
    """Devuelve la distancia del punto de aproximación rayo-segmento.

    Para picking el orden de impactos es más importante que resolver una
    ecuación de cilindro completa. La distancia mínima a la cápsula finita es
    estable también para enlaces degenerados y rayos casi paralelos.
    """
    if radius <= EPSILON or not math.isfinite(radius):
        return None
    segment = sub(end, start)
    segment_length_sq = dot(segment, segment)
    if segment_length_sq <= EPSILON:
        return ray_sphere_intersection(ray, start, radius)

    relative = sub(ray.origin, start)
    ray_segment_dot = dot(ray.direction, segment)
    ray_relative_dot = dot(ray.direction, relative)
    segment_relative_dot = dot(segment, relative)
    denominator = segment_length_sq - ray_segment_dot * ray_segment_dot

    if abs(denominator) > EPSILON:
        ray_t = (ray_segment_dot * segment_relative_dot - segment_length_sq * ray_relative_dot) / denominator
        segment_u = (segment_relative_dot - ray_segment_dot * ray_relative_dot) / denominator
    else:
        ray_t = 0.0
        segment_u = segment_relative_dot / segment_length_sq

    if ray_t < 0.0:
        ray_t = 0.0
        segment_u = segment_relative_dot / segment_length_sq
    segment_u = max(0.0, min(1.0, segment_u))
    if ray_t > 0.0 and (segment_u <= 0.0 or segment_u >= 1.0):
        endpoint = start if segment_u <= 0.0 else end
        ray_t = max(0.0, dot(sub(endpoint, ray.origin), ray.direction))

    ray_point = add(ray.origin, scale(ray.direction, ray_t))
    segment_point = add(start, scale(segment, segment_u))
    if length(sub(ray_point, segment_point)) <= radius + EPSILON:
        return ray_t
    return None
