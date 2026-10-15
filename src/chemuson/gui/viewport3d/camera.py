"""Cámara orbital y proyecciones del visor 3D."""

from __future__ import annotations

from dataclasses import dataclass
from enum import Enum
import math

from .geometry import Vec3, add, cross, normalize, scale, scene_bounds, sub


class CameraProjection(str, Enum):
    PERSPECTIVE = "perspective"
    ORTHOGRAPHIC = "orthographic"


@dataclass(slots=True)
class CameraState:
    target: Vec3 = (0.0, 0.0, 0.0)
    distance: float = 10.0
    yaw: float = 0.0
    pitch: float = 0.25
    fov_y_degrees: float = 45.0
    ortho_scale: float = 10.0
    projection: CameraProjection = CameraProjection.PERSPECTIVE

    MIN_DISTANCE = 0.05
    MAX_DISTANCE = 1.0e5
    MIN_ORTHO_SCALE = 0.05
    MAX_ORTHO_SCALE = 1.0e5
    MIN_PITCH = -math.pi / 2.0 + 0.02
    MAX_PITCH = math.pi / 2.0 - 0.02

    def __post_init__(self) -> None:
        self.distance = self._clamp_distance(self.distance)
        self.ortho_scale = self._clamp_ortho_scale(self.ortho_scale)
        self.pitch = max(self.MIN_PITCH, min(self.MAX_PITCH, float(self.pitch)))
        self.fov_y_degrees = max(5.0, min(120.0, float(self.fov_y_degrees)))

    @staticmethod
    def _clamp_distance(value: float) -> float:
        return max(CameraState.MIN_DISTANCE, min(CameraState.MAX_DISTANCE, float(value)))

    @staticmethod
    def _clamp_ortho_scale(value: float) -> float:
        return max(CameraState.MIN_ORTHO_SCALE, min(CameraState.MAX_ORTHO_SCALE, float(value)))

    @property
    def eye(self) -> Vec3:
        cos_pitch = math.cos(self.pitch)
        return add(
            self.target,
            (
                self.distance * cos_pitch * math.sin(self.yaw),
                self.distance * math.sin(self.pitch),
                self.distance * cos_pitch * math.cos(self.yaw),
            ),
        )

    @property
    def basis(self) -> tuple[Vec3, Vec3, Vec3]:
        forward = normalize(sub(self.target, self.eye), (0.0, 0.0, -1.0))
        right = normalize(cross(forward, (0.0, 1.0, 0.0)), (1.0, 0.0, 0.0))
        up = normalize(cross(right, forward), (0.0, 1.0, 0.0))
        return forward, right, up

    def orbit(self, delta_x: float, delta_y: float, viewport_height: float) -> None:
        if viewport_height <= 0.0 or not math.isfinite(viewport_height):
            return
        sensitivity = 2.0 * math.pi / max(240.0, float(viewport_height))
        self.yaw -= float(delta_x) * sensitivity
        self.pitch = max(
            self.MIN_PITCH,
            min(self.MAX_PITCH, self.pitch - float(delta_y) * sensitivity),
        )

    def pan(self, delta_x: float, delta_y: float, viewport_height: float) -> None:
        if viewport_height <= 0.0 or not math.isfinite(viewport_height):
            return
        _forward, right, up = self.basis
        scale_per_pixel = (
            self.ortho_scale if self.projection == CameraProjection.ORTHOGRAPHIC else self.distance
        ) / float(viewport_height)
        self.target = add(
            self.target,
            add(scale(right, -float(delta_x) * scale_per_pixel), scale(up, float(delta_y) * scale_per_pixel)),
        )

    def zoom(self, wheel_delta: float) -> None:
        if not math.isfinite(wheel_delta):
            return
        factor = math.exp(-float(wheel_delta) * 0.001)
        if self.projection == CameraProjection.ORTHOGRAPHIC:
            self.ortho_scale = self._clamp_ortho_scale(self.ortho_scale * factor)
        else:
            self.distance = self._clamp_distance(self.distance * factor)

    def set_projection(self, projection: CameraProjection) -> None:
        self.projection = CameraProjection(projection)

    def fit_scene(self, scene, padding: float = 1.35) -> bool:
        bounds = scene_bounds(scene)
        if bounds is None:
            return False
        padding = max(1.0, float(padding))
        self.target = bounds.center
        radius = max(0.05, bounds.radius) * padding
        if self.projection == CameraProjection.ORTHOGRAPHIC:
            self.ortho_scale = self._clamp_ortho_scale(2.0 * radius)
        else:
            half_fov = math.radians(self.fov_y_degrees) * 0.5
            self.distance = self._clamp_distance(radius / max(0.05, math.tan(half_fov)))
        return True

    def ray_parameters(self, viewport_width: float, viewport_height: float) -> tuple[Vec3, Vec3, Vec3, Vec3]:
        if viewport_width <= 0.0 or viewport_height <= 0.0:
            raise ValueError("El viewport debe tener dimensiones positivas")
        forward, right, up = self.basis
        return self.eye, forward, right, up
