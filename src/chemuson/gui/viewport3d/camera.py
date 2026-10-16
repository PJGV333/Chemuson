"""Cámara orbital y proyecciones del visor 3D."""

from __future__ import annotations

from dataclasses import dataclass
from enum import Enum
import math

from .geometry import Vec3, add, cross, dot, normalize, scale, scene_bounds, sub


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
    fit_radius: float = 0.0

    MIN_DISTANCE = 0.05
    MAX_DISTANCE = 1.0e12
    MIN_ORTHO_SCALE = 0.05
    MAX_ORTHO_SCALE = 1.0e12
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

    def fit_scene(
        self,
        scene,
        padding: float = 1.35,
        viewport_width: float = 1.0,
        viewport_height: float = 1.0,
    ) -> bool:
        bounds = scene_bounds(scene)
        if (
            bounds is None
            or viewport_width <= 0.0
            or viewport_height <= 0.0
            or not math.isfinite(viewport_width)
            or not math.isfinite(viewport_height)
            or not math.isfinite(padding)
        ):
            return False
        padding = max(1.0, float(padding))
        radius = max(0.05, bounds.radius) * padding
        if not math.isfinite(radius):
            return False
        aspect = float(viewport_width) / float(viewport_height)
        if not math.isfinite(aspect) or aspect <= 0.0:
            return False
        if self.projection == CameraProjection.ORTHOGRAPHIC:
            requested_scale = 2.0 * radius * max(1.0, 1.0 / aspect)
            requested_distance = max(self.distance, radius * 1.1)
            if (
                not math.isfinite(requested_scale)
                or requested_scale > self.MAX_ORTHO_SCALE
                or requested_distance > self.MAX_DISTANCE
            ):
                return False
            self.ortho_scale = self._clamp_ortho_scale(requested_scale)
            self.distance = self._clamp_distance(requested_distance)
        else:
            vertical_half_angle = math.radians(self.fov_y_degrees) * 0.5
            horizontal_half_angle = math.atan(math.tan(vertical_half_angle) * aspect)
            limiting_half_angle = min(vertical_half_angle, horizontal_half_angle)
            requested_distance = radius / max(1.0e-9, math.sin(limiting_half_angle))
            if not math.isfinite(requested_distance) or requested_distance > self.MAX_DISTANCE:
                return False
            self.distance = self._clamp_distance(requested_distance)
        self.target = bounds.center
        self.fit_radius = radius
        return True

    def clip_planes(self) -> tuple[float, float]:
        """Return depth planes that contain the currently fitted scene sphere."""
        margin = max(0.05, self.fit_radius * 1.1)
        near = max(0.001, self.distance - margin)
        far = max(near + 0.01, self.distance + margin)
        return near, far

    def ray_parameters(self, viewport_width: float, viewport_height: float) -> tuple[Vec3, Vec3, Vec3, Vec3]:
        if viewport_width <= 0.0 or viewport_height <= 0.0:
            raise ValueError("El viewport debe tener dimensiones positivas")
        forward, right, up = self.basis
        return self.eye, forward, right, up

    def screen_to_world(
        self,
        pixel_x: float,
        pixel_y: float,
        viewport_width: float,
        viewport_height: float,
        *,
        plane_point: Vec3 | None = None,
    ) -> Vec3:
        """Unprojects a pixel sobre un plano perpendicular a la cámara."""
        values = (pixel_x, pixel_y, viewport_width, viewport_height)
        if (
            viewport_width <= 0.0
            or viewport_height <= 0.0
            or not all(math.isfinite(float(value)) for value in values)
        ):
            raise ValueError("Viewport o coordenada de cursor no válidos")
        eye, forward, right, up = self.ray_parameters(viewport_width, viewport_height)
        if plane_point is None:
            plane = self.target
        else:
            plane = (
                float(plane_point[0]),
                float(plane_point[1]),
                float(plane_point[2]),
            )
        if len(plane) != 3 or not all(math.isfinite(value) for value in plane):
            raise ValueError("Punto del plano no válido")
        ndc_x = 2.0 * float(pixel_x) / float(viewport_width) - 1.0
        ndc_y = 1.0 - 2.0 * float(pixel_y) / float(viewport_height)
        aspect = float(viewport_width) / float(viewport_height)
        if self.projection == CameraProjection.ORTHOGRAPHIC:
            half_height = self.ortho_scale * 0.5
            origin = add(
                add(eye, scale(right, ndc_x * half_height * aspect)),
                scale(up, ndc_y * half_height),
            )
            return add(origin, scale(forward, dot(sub(plane, origin), forward)))

        half_fov = math.tan(math.radians(self.fov_y_degrees) * 0.5)
        direction = normalize(
            add(
                forward,
                add(
                    scale(right, ndc_x * aspect * half_fov),
                    scale(up, ndc_y * half_fov),
                ),
            )
        )
        denominator = dot(direction, forward)
        if denominator <= 1.0e-12:
            raise ValueError("El plano no intersecta el rayo de la cámara")
        distance = dot(sub(plane, eye), forward) / denominator
        if distance < 0.0 or not math.isfinite(distance):
            raise ValueError("El plano no está delante de la cámara")
        return add(eye, scale(direction, distance))

    def screen_delta_to_world(
        self,
        delta_x: float,
        delta_y: float,
        viewport_width: float,
        viewport_height: float,
        *,
        position: Vec3,
    ) -> Vec3:
        """Convierte desplazamiento de cursor a unidades del plano de pantalla."""
        values = (delta_x, delta_y, viewport_width, viewport_height, *position)
        if (
            viewport_width <= 0.0
            or viewport_height <= 0.0
            or not all(math.isfinite(float(value)) for value in values)
        ):
            raise ValueError("Viewport o desplazamiento de cursor no válidos")
        forward, right, up = self.basis
        if self.projection == CameraProjection.ORTHOGRAPHIC:
            units_per_pixel = self.ortho_scale / float(viewport_height)
        else:
            position3d = (
                float(position[0]),
                float(position[1]),
                float(position[2]),
            )
            depth = dot(sub(position3d, self.eye), forward)
            if depth <= 0.0:
                raise ValueError("El átomo no está delante de la cámara")
            units_per_pixel = (
                2.0 * depth * math.tan(math.radians(self.fov_y_degrees) * 0.5)
                / float(viewport_height)
            )
        return add(
            scale(right, float(delta_x) * units_per_pixel),
            scale(up, -float(delta_y) * units_per_pixel),
        )
