"""Picking CPU coherente con la cámara del visor."""

from __future__ import annotations

from dataclasses import dataclass
import math

from .camera import CameraProjection, CameraState
from .geometry import ATOM_RADIUS_SCALE, Ray3D, add, ray_capsule_intersection, ray_sphere_intersection, scale
from .mesh import bond_lane_layout, bond_lateral_direction


@dataclass(frozen=True, slots=True)
class PickResult:
    kind: str
    primitive_id: int | str
    distance: float


def ray_from_screen(
    camera: CameraState,
    viewport_width: float,
    viewport_height: float,
    pixel_x: float,
    pixel_y: float,
) -> Ray3D:
    if (
        viewport_width <= 0.0
        or viewport_height <= 0.0
        or not all(
            math.isfinite(value)
            for value in (viewport_width, viewport_height, pixel_x, pixel_y)
        )
    ):
        raise ValueError("Viewport o coordenada de cursor no válidos")
    eye, forward, right, up = camera.ray_parameters(viewport_width, viewport_height)
    ndc_x = (2.0 * float(pixel_x) / float(viewport_width)) - 1.0
    ndc_y = 1.0 - (2.0 * float(pixel_y) / float(viewport_height))
    aspect = float(viewport_width) / float(viewport_height)
    if camera.projection == CameraProjection.ORTHOGRAPHIC:
        half_height = camera.ortho_scale * 0.5
        origin = add(
            add(eye, scale(right, ndc_x * half_height * aspect)),
            scale(up, ndc_y * half_height),
        )
        return Ray3D(origin, forward)
    half_fov = math.tan(math.radians(camera.fov_y_degrees) * 0.5)
    direction = add(
        forward,
        add(scale(right, ndc_x * aspect * half_fov), scale(up, ndc_y * half_fov)),
    )
    return Ray3D(eye, direction)


def pick_scene(
    scene,
    camera: CameraState,
    viewport_width: float,
    viewport_height: float,
    pixel_x: float,
    pixel_y: float,
) -> PickResult | None:
    ray = ray_from_screen(camera, viewport_width, viewport_height, pixel_x, pixel_y)
    atoms = {atom.id: atom for atom in getattr(scene, "atoms", ())}
    bonds = tuple(getattr(scene, "bonds", ()))
    adjacency: dict[int | str, list[int | str]] = {atom_id: [] for atom_id in atoms}
    for bond in bonds:
        if bond.a1_id in adjacency and bond.a2_id in adjacency:
            adjacency[bond.a1_id].append(bond.a2_id)
            adjacency[bond.a2_id].append(bond.a1_id)
    best: PickResult | None = None

    def consider(kind: str, primitive_id: int | str, distance: float | None) -> None:
        nonlocal best
        if distance is None or not math.isfinite(distance) or distance < 0.0:
            return
        candidate = PickResult(kind, primitive_id, float(distance))
        if best is None or candidate.distance < best.distance - 1.0e-8:
            best = candidate
            return
        if abs(candidate.distance - best.distance) <= 1.0e-8:
            priority = {"atom": 0, "bond": 1}
            if priority[candidate.kind] < priority[best.kind]:
                best = candidate

    for atom in getattr(scene, "atoms", ()):
        center = (float(atom.x), float(atom.y), float(atom.z))
        radius = max(0.0, float(getattr(atom, "radius", 0.0))) * ATOM_RADIUS_SCALE
        consider("atom", atom.id, ray_sphere_intersection(ray, center, radius))

    for bond in bonds:
        start_atom = atoms.get(bond.a1_id)
        end_atom = atoms.get(bond.a2_id)
        if start_atom is None or end_atom is None:
            continue
        start = (float(start_atom.x), float(start_atom.y), float(start_atom.z))
        end = (float(end_atom.x), float(end_atom.y), float(end_atom.z))
        lateral = bond_lateral_direction(start, end, bond.a1_id, bond.a2_id, atoms, adjacency)
        offsets, radius = bond_lane_layout(
            getattr(bond, "order", 1),
            bool(getattr(bond, "is_aromatic", False)),
        )
        bond_distance = None
        for offset in offsets:
            displacement = scale(lateral, offset)
            lane_distance = ray_capsule_intersection(
                ray,
                add(start, displacement),
                add(end, displacement),
                radius,
            )
            if lane_distance is not None and (
                bond_distance is None or lane_distance < bond_distance
            ):
                bond_distance = lane_distance
        consider("bond", bond.id, bond_distance)
    return best
