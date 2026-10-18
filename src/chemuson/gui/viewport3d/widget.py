"""Widget Qt6 del visor: eventos, cámara y selección local."""

from __future__ import annotations

from dataclasses import dataclass, replace
import math
import time
from typing import Any, cast

from PyQt6.QtCore import QPoint, Qt, QTimer, pyqtSignal
from PyQt6.QtGui import (
    QHideEvent,
    QBrush,
    QColor,
    QMouseEvent,
    QOpenGLContext,
    QPainter,
    QPen,
    QResizeEvent,
    QShowEvent,
    QSurfaceFormat,
    QWheelEvent,
)
from PyQt6.QtOpenGLWidgets import QOpenGLWidget
from PyQt6.QtWidgets import QApplication

from chemuson.geometry3d import SceneAtom3D, SceneBond3D, SceneMolecule3D
from chemuson.geometry3d.editing import angle_degrees, dihedral_degrees, distance
from chemuson.geometry3d.model import (
    CoordinateSet3D,
    VisualFramePlayback,
    scene_with_coordinates,
)

from .camera import CameraProjection, CameraState
from .geometry import add, dot
from .renderer import (
    GL_COLOR_BUFFER_BIT,
    GL_DEPTH_BUFFER_BIT,
    GL_DEPTH_TEST,
    GL_LESS,
    MoleculeRenderer,
)
from .mesh import normalize_representation
from .picking import PickResult, pick_scene


def project_camera_axes(camera: CameraState) -> tuple[tuple[str, float, float, float], ...]:
    """Project positive world axes into the current camera's screen basis."""
    forward, right, up = camera.basis
    world_axes = (
        ("X", (1.0, 0.0, 0.0)),
        ("Y", (0.0, 1.0, 0.0)),
        ("Z", (0.0, 0.0, 1.0)),
    )
    return tuple(
        (label, dot(axis, right), -dot(axis, up), dot(axis, forward))
        for label, axis in world_axes
    )


@dataclass
class _GeometryDragState:
    before: dict[int, tuple[float, float, float]]
    latest: dict[int, tuple[float, float, float]]
    selected_atom_ids: tuple[int, ...]
    anchor_atom_id: int
    start_position: QPoint
    moved: bool = False
    depth_mode: bool = False
    depth_button: Qt.MouseButton = Qt.MouseButton.MiddleButton
    geometry_operation: str | None = None
    geometry_atom_ids: tuple[int, ...] = ()
    geometry_start_value: float | None = None


@dataclass
class _BuildDragState:
    anchor_atom_id: int
    anchor_position: tuple[float, float, float]
    start_position: QPoint
    current_position: QPoint
    moved: bool = False


@dataclass(frozen=True)
class _BuildPreviewState:
    anchor_atom_id: int
    anchor_position: tuple[float, float, float]
    end_position: tuple[float, float, float]
    element: str
    target_atom_id: int | None


@dataclass(frozen=True)
class _PendingBuildClick:
    element: str
    position: tuple[float, float, float]


@dataclass
class _DerivedHydrogenDragState:
    atom_id: str
    parent_atom_id: int
    before: tuple[float, float, float]
    latest: tuple[float, float, float]
    start_position: QPoint
    moved: bool = False
    depth_mode: bool = False
    depth_button: Qt.MouseButton = Qt.MouseButton.MiddleButton


class MolecularViewportWidget(QOpenGLWidget):
    """Viewport 3D con estado de escena y selección estrictamente local."""

    primitiveSelected = pyqtSignal(str, object)
    selectionCleared = pyqtSignal()
    statusChanged = pyqtSignal(str)
    editModeChanged = pyqtSignal(bool)
    editToolChanged = pyqtSignal(str)
    geometryPreviewChanged = pyqtSignal(object, object, object)
    geometryEditCommitted = pyqtSignal(object, object, object)
    geometryEditCancelled = pyqtSignal()
    geometryEditStarted = pyqtSignal()
    geometryHandleStarted = pyqtSignal(str, object)
    geometryHandlePreviewRequested = pyqtSignal(str, object, float)
    buildCancelled = pyqtSignal()
    buildPointRequested = pyqtSignal(object)
    buildChainExtensionRequested = pyqtSignal(object)
    buildAnchorSelected = pyqtSignal(object)
    buildBondRequested = pyqtSignal(object)
    buildElementChangeRequested = pyqtSignal(int, str)
    buildPreviewChanged = pyqtSignal(object)
    liveMMPlaybackDrained = pyqtSignal()
    derivedHydrogenPreviewChanged = pyqtSignal(object, object, object)
    derivedHydrogenEditRequested = pyqtSignal(object)

    def __init__(self, parent=None) -> None:
        fmt = QSurfaceFormat()
        fmt.setVersion(4, 1)
        fmt.setProfile(QSurfaceFormat.OpenGLContextProfile.CoreProfile)
        fmt.setDepthBufferSize(24)
        fmt.setSamples(4)
        super().__init__(parent)
        self.setFormat(fmt)
        self.setFocusPolicy(Qt.FocusPolicy.StrongFocus)
        self.setMouseTracking(True)
        self.setToolTip(
            "Navegación: clic selecciona, arrastre izquierdo orbita, "
            "arrastre derecho desplaza y rueda amplía/reduce. "
            "Construir: arrastre desde un átomo para formar un enlace 3D; "
            "doble clic en un átomo sustituye su elemento y doble clic en espacio "
            "extiende el ancla con un enlace sencillo. Escape cancela. "
            "Mover: arrastre izquierdo sobre átomo para editarlo; arrastre vacío "
            "orbita la cámara; botón central mueve en profundidad."
        )
        self.camera = CameraState()
        self.scene = SceneMolecule3D(atoms=(), bonds=(), source="")
        self.selected_atom_id: int | str | None = None
        self.selected_atom_ids: list[int] = []
        self.selected_bond_id: int | str | None = None
        self.edit_mode = False
        self.edit_tool = "move"
        self.geometry_operation = "bond_length"
        self.build_element = "C"
        self.build_bond_order = (1, False)
        self.representation = "ball_and_stick"
        self._measure_atom_ids: list[int] = []
        self._renderer = MoleculeRenderer()
        self._last_mouse_position: QPoint | None = None
        self._dragging = False
        self._camera_drag_origin_empty = False
        self._orbit_guide_active = False
        self._geometry_edit_notified = False
        self._geometry_drag_state: _GeometryDragState | None = None
        self._build_drag_state: _BuildDragState | None = None
        self._derived_hydrogen_drag_state: _DerivedHydrogenDragState | None = None
        self._build_preview: _BuildPreviewState | None = None
        self._pending_build_click: _PendingBuildClick | None = None
        self._suppress_next_build_release = False
        self._build_click_timer = QTimer(self)
        self._build_click_timer.setSingleShot(True)
        self._build_click_timer.setInterval(max(1, QApplication.doubleClickInterval()))
        self._build_click_timer.timeout.connect(self._commit_pending_build_click)
        self._live_mm_playback = VisualFramePlayback()
        self._live_mm_playback_enabled = False
        self._live_mm_telemetry: dict[str, object] | None = None
        self._live_mm_drained_emitted = False
        self._live_mm_tick_count = 0
        self._live_mm_tick_total_s = 0.0
        self._live_mm_animation_timer = QTimer(self)
        self._live_mm_animation_timer.setTimerType(Qt.TimerType.PreciseTimer)
        self._live_mm_animation_timer.setInterval(16)
        self._live_mm_animation_timer.timeout.connect(self._advance_live_mm_playback)
        self._press_hit: PickResult | None = None
        self._gl_ready = False
        self._mesh_dirty = True
        self._scene_diagnostics: dict[str, object] = {
            "stage": "scene_empty",
            "atom_count": 0,
            "bond_count": 0,
            "gl_context_ready": False,
            "upload_status": "not_requested",
            "render_status": "not_requested",
            "error_code": None,
            "status_message": "",
        }
        self._auto_fit = True
        self._active_context: QOpenGLContext | None = None
        self._connected_context: QOpenGLContext | None = None
        self.statusChanged.emit("Esperando contexto OpenGL…")

    @property
    def gl_ready(self) -> bool:
        return self._gl_ready

    @property
    def scene_diagnostics(self) -> dict[str, object]:
        return dict(self._scene_diagnostics)

    def _report_scene_display_status(self, message: str) -> None:
        if self._scene_diagnostics.get("status_message") == message:
            return
        self._scene_diagnostics["status_message"] = message
        self.statusChanged.emit(message)

    @property
    def live_mm_playback_enabled(self) -> bool:
        return self._live_mm_playback_enabled

    @property
    def live_mm_playback_active(self) -> bool:
        return self._live_mm_playback.active

    @property
    def live_mm_playback_stats(self) -> dict[str, float | int]:
        return {
            "ticks": self._live_mm_tick_count,
            "mean_tick_ms": (
                1000.0 * self._live_mm_tick_total_s / self._live_mm_tick_count
                if self._live_mm_tick_count
                else 0.0
            ),
            "dropped_frames": self._live_mm_playback.dropped_frames,
            "pending_frames": self._live_mm_playback.pending_count,
        }

    @property
    def live_mm_telemetry(self) -> dict[str, object] | None:
        return dict(self._live_mm_telemetry) if self._live_mm_telemetry is not None else None

    def begin_live_mm_playback(self, initial: CoordinateSet3D | None) -> None:
        self._live_mm_animation_timer.stop()
        self._live_mm_playback.begin(initial)
        self._live_mm_playback_enabled = True
        self._live_mm_drained_emitted = False
        self._live_mm_tick_count = 0
        self._live_mm_tick_total_s = 0.0

    def enqueue_live_mm_frame(self, frame: object) -> bool:
        coordinates = getattr(frame, "coordinates", None)
        if not self._live_mm_playback_enabled or not isinstance(coordinates, CoordinateSet3D):
            return False
        accepted = self._live_mm_playback.enqueue(coordinates)
        if not accepted:
            self._apply_live_mm_coordinates(coordinates)
            return False
        if not self._live_mm_animation_timer.isActive():
            self._live_mm_animation_timer.start()
        return True

    def finish_live_mm_playback(self, final: CoordinateSet3D | None) -> None:
        if not self._live_mm_playback_enabled:
            self.begin_live_mm_playback(final)
        self._live_mm_playback.finish(final)
        if self._live_mm_playback.drained:
            self._complete_live_mm_playback()
        elif not self._live_mm_animation_timer.isActive():
            self._live_mm_animation_timer.start()

    def cancel_live_mm_playback(
        self,
        restore: CoordinateSet3D | None = None,
        *,
        clear_telemetry: bool = True,
    ) -> None:
        self._live_mm_animation_timer.stop()
        self._live_mm_playback.cancel()
        self._live_mm_playback_enabled = False
        self._live_mm_drained_emitted = False
        if clear_telemetry:
            self._live_mm_telemetry = None
        if restore is not None:
            self._apply_live_mm_coordinates(restore)
        else:
            self.update()

    def set_live_mm_telemetry(self, payload: object | None) -> None:
        self._live_mm_telemetry = dict(payload) if isinstance(payload, dict) else None
        self.update()

    def _apply_live_mm_coordinates(self, coordinates: CoordinateSet3D) -> None:
        self.scene = scene_with_coordinates(self.scene, coordinates)
        self._mesh_dirty = True
        self.update()

    def _advance_live_mm_playback(self) -> None:
        started = time.perf_counter()
        coordinates = self._live_mm_playback.sample(time.monotonic())
        if coordinates is not None:
            self._apply_live_mm_coordinates(coordinates)
        if self._live_mm_telemetry is not None:
            self._live_mm_telemetry["visual_interpolating"] = self._live_mm_playback.active
            self._live_mm_telemetry["dropped_frames"] = self._live_mm_playback.dropped_frames
        self._live_mm_tick_count += 1
        self._live_mm_tick_total_s += time.perf_counter() - started
        if self._live_mm_playback.drained:
            self._complete_live_mm_playback()

    def _complete_live_mm_playback(self) -> None:
        self._live_mm_animation_timer.stop()
        self._live_mm_playback_enabled = False
        if self._live_mm_telemetry is not None:
            self._live_mm_telemetry["visual_interpolating"] = False
        self.update()
        if not self._live_mm_drained_emitted:
            self._live_mm_drained_emitted = True
            self.liveMMPlaybackDrained.emit()

    def _notify_geometry_edit_started(self) -> None:
        if self._geometry_edit_notified:
            return
        self._geometry_edit_notified = True
        self.geometryEditStarted.emit()
        build_state = self._build_drag_state
        if build_state is not None:
            anchor = next(
                (atom for atom in self.scene.atoms if atom.id == build_state.anchor_atom_id),
                None,
            )
            if anchor is not None:
                build_state.anchor_position = (
                    float(anchor.x),
                    float(anchor.y),
                    float(anchor.z),
                )
        state = self._geometry_drag_state
        if state is not None:
            state.before = {
                int(atom.id): (float(atom.x), float(atom.y), float(atom.z))
                for atom in self.scene.atoms
                if isinstance(atom.id, int)
            }
            state.latest = dict(state.before)
            if state.geometry_operation is not None:
                try:
                    state.geometry_start_value = self._geometry_measurement(
                        state.geometry_operation,
                        state.geometry_atom_ids,
                    )
                except (KeyError, ValueError):
                    state.geometry_start_value = None
        derived_state = self._derived_hydrogen_drag_state
        if derived_state is not None:
            atom = next(
                (item for item in self.scene.atoms if item.id == derived_state.atom_id),
                None,
            )
            if atom is not None:
                point = (float(atom.x), float(atom.y), float(atom.z))
                derived_state.before = point
                derived_state.latest = point

    @property
    def geometry_drag_active(self) -> bool:
        return self._geometry_drag_state is not None

    def initializeGL(self) -> None:
        context = self.context()
        if context is not None and context != self._connected_context:
            # PyQt's runtime accepts the connection type; its stubs omit it.
            cast(Any, context.aboutToBeDestroyed).connect(
                self._cleanup_gl,
                Qt.ConnectionType.DirectConnection,
            )
            self._connected_context = context
        self._active_context = context
        ready, message = self._renderer.initialize()
        self._gl_ready = ready
        self._scene_diagnostics["gl_context_ready"] = ready
        self.statusChanged.emit(message)
        if ready:
            self._mesh_dirty = True
            if self.scene.atoms:
                self._scene_diagnostics["upload_status"] = "pending"
                self._scene_diagnostics["render_status"] = "pending"
        elif self.scene.atoms:
            self._scene_diagnostics["upload_status"] = "unavailable"
            self._scene_diagnostics["error_code"] = "display_failed"
            self._report_scene_display_status(
                f"display_failed: no se pudo inicializar OpenGL ({message})"
            )

    def showEvent(self, a0: QShowEvent | None) -> None:
        if a0 is None:
            return
        super().showEvent(a0)
        if self._gl_ready or self.context() is None:
            return
        try:
            self.makeCurrent()
            self.initializeGL()
        except Exception as exc:
            self.statusChanged.emit(f"No se pudo restaurar el contexto OpenGL: {exc}")
        finally:
            try:
                self.doneCurrent()
            except Exception as exc:
                self.statusChanged.emit(f"No se pudo liberar el contexto OpenGL: {exc}")

    def hideEvent(self, a0: QHideEvent | None) -> None:
        if a0 is None:
            return
        if self._active_context is not None:
            self._cleanup_gl()
        super().hideEvent(a0)

    def _cleanup_gl(self) -> None:
        context = self._active_context
        if context is not None:
            try:
                self.makeCurrent()
                self._renderer.cleanup()
            except Exception as exc:
                self.statusChanged.emit(f"No se pudieron liberar los recursos OpenGL: {exc}")
            try:
                self.doneCurrent()
            except Exception as exc:
                self.statusChanged.emit(f"No se pudo liberar el contexto OpenGL: {exc}")
        if context == self._connected_context:
            self._connected_context = None
        self._active_context = None
        self._gl_ready = False
        self._mesh_dirty = True
        self.update()

    def paintGL(self) -> None:
        if not self._gl_ready:
            if self.scene.atoms:
                self._scene_diagnostics["upload_status"] = "unavailable"
                self._scene_diagnostics["error_code"] = "display_failed"
                self._report_scene_display_status(
                    "display_failed: escena preparada, pero el contexto OpenGL no está listo."
                )
            return
        gl = self._renderer._gl
        if gl is None:
            if self.scene.atoms:
                self._scene_diagnostics["upload_status"] = "unavailable"
                self._scene_diagnostics["error_code"] = "display_failed"
                self._report_scene_display_status(
                    "display_failed: renderer OpenGL no disponible para cargar la escena."
                )
            return
        gl.glEnable(GL_DEPTH_TEST)
        gl.glDepthMask(True)
        gl.glDepthFunc(GL_LESS)
        gl.glClearDepth(1.0)
        gl.glClearColor(0.055, 0.065, 0.08, 1.0)
        gl.glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT)
        if self._mesh_dirty:
            try:
                uploaded = self._renderer.upload(
                    self._scene_with_build_preview(),
                    self.selected_atom_id,
                    self.selected_bond_id,
                    self.representation,
                    self.selected_atom_ids,
                )
            except Exception as exc:
                uploaded = False
                self._scene_diagnostics["upload_error"] = str(exc)
            self._mesh_dirty = not uploaded
            if not uploaded:
                self._scene_diagnostics["upload_status"] = "failed"
                self._scene_diagnostics["error_code"] = "display_failed"
                self._report_scene_display_status(
                    "display_failed: falló la carga de la malla al renderer OpenGL."
                )
                return
            self._scene_diagnostics["upload_status"] = "uploaded"
        try:
            self._renderer.render(self.camera, self.width(), self.height())
        except Exception as exc:
            self._scene_diagnostics["render_status"] = "failed"
            self._scene_diagnostics["render_error"] = str(exc)
            self._scene_diagnostics["error_code"] = "display_failed"
            self._report_scene_display_status(
                "display_failed: el renderer no pudo solicitar el dibujo OpenGL."
            )
            return
        self._scene_diagnostics["render_status"] = "requested"
        self._report_scene_display_status(
            "render_requested: escena enviada al renderer OpenGL; aceptación visual pendiente."
        )
        self._paint_geometry_handle()
        self._paint_live_mm_hud()
        self._paint_orientation_axes()
        self._paint_orbit_guide()

    def _paint_orientation_axes(self) -> None:
        origin = QPoint(48, max(48, self.height() - 48))
        painter = QPainter()
        if not painter.begin(self):
            return
        painter.setRenderHint(QPainter.RenderHint.Antialiasing)
        painter.setPen(QPen(QColor(170, 185, 198, 180), 1))
        painter.setBrush(QBrush(QColor(12, 18, 26, 190)))
        painter.drawEllipse(origin, 41, 41)
        colors = {"X": QColor(245, 91, 91), "Y": QColor(87, 216, 156), "Z": QColor(93, 169, 246)}
        length_pixels = 32.0
        for label, screen_x, screen_y, _depth in project_camera_axes(self.camera):
            end = QPoint(
                origin.x() + round(screen_x * length_pixels),
                origin.y() + round(screen_y * length_pixels),
            )
            color = colors[label]
            painter.setPen(QPen(color, 2))
            painter.drawLine(origin, end)
            if end == origin:
                painter.drawEllipse(origin, 3, 3)
                label_position = origin + QPoint(7, -7)
            else:
                label_position = end + QPoint(4, -3)
            painter.drawText(label_position, label)
        painter.end()

    def _paint_orbit_guide(self) -> None:
        if not self._orbit_guide_active:
            return
        center_x, center_y = self.width() // 2, self.height() // 2
        painter = QPainter()
        if not painter.begin(self):
            return
        painter.setRenderHint(QPainter.RenderHint.Antialiasing)
        painter.setPen(QPen(QColor(95, 205, 245, 210), 3))
        painter.drawArc(center_x - 42, center_y - 42, 84, 84, 20 * 16, 125 * 16)
        painter.drawArc(center_x - 42, center_y - 42, 84, 84, 200 * 16, 125 * 16)
        painter.setPen(QPen(QColor(220, 240, 250, 230), 1))
        painter.drawText(center_x - 48, center_y + 58, "Cámara · orbitar")
        painter.end()

    def _live_mm_hud_lines(self) -> tuple[str, ...] | None:
        telemetry = self._live_mm_telemetry
        if telemetry is None:
            return None
        forcefield = str(telemetry.get("forcefield", "RDKit"))
        unit = str(telemetry.get("energy_unit", "kcal/mol"))
        energy = telemetry.get("energy")
        delta_energy = telemetry.get("delta_energy")
        energy_text = (
            f"{float(energy):.6f} {unit}"
            if isinstance(energy, (int, float)) and math.isfinite(float(energy))
            else "—"
        )
        delta_text = (
            f"{float(delta_energy):+.6f} {unit}"
            if isinstance(delta_energy, (int, float)) and math.isfinite(float(delta_energy))
            else "—"
        )

        def telemetry_int(key: str) -> int:
            value = telemetry.get(key, 0)
            if isinstance(value, (int, float, str)):
                try:
                    return max(0, int(value))
                except (TypeError, ValueError, OverflowError):
                    return 0
            return 0

        step = telemetry_int("step")
        max_steps = telemetry_int("max_steps")
        fixed_count = telemetry_int("fixed_atom_count")
        constraints = telemetry.get("constraints", {})
        coordinate_unit = (
            str(constraints.get("coordinate_unit", "Å"))
            if isinstance(constraints, dict)
            else "Å"
        )
        restriction_text = (
            f"{fixed_count} átomos fijos · posiciones {coordinate_unit}"
            if fixed_count
            else "Restricciones: ninguna"
        )
        interpolation = bool(telemetry.get("visual_interpolating", False))
        return (
            f"Live MM · {forcefield}",
            f"E {energy_text} · ΔE (muestras) {delta_text}",
            f"Paso {step}/{max_steps} · {restriction_text}",
            "Interpolación visual · gradiente no calculado"
            if interpolation
            else "ΔE entre muestras · gradiente no calculado",
        )

    def _paint_live_mm_hud(self) -> None:
        lines = self._live_mm_hud_lines()
        if lines is None:
            return
        width = min(max(300, self.width() - 24), 440)
        painter = QPainter()
        if not painter.begin(self):
            return
        painter.setRenderHint(QPainter.RenderHint.Antialiasing)
        painter.setPen(QPen(QColor(220, 230, 238, 235), 1))
        painter.setBrush(QBrush(QColor(12, 18, 26, 220)))
        painter.drawRoundedRect(12, 12, width, 84, 6, 6)
        painter.setPen(QPen(QColor(240, 245, 250, 245)))
        for index, line in enumerate(lines):
            painter.drawText(22, 31 + index * 16, line)
        painter.end()

    def _paint_geometry_handle(self) -> None:
        handle = self._geometry_handle_screen_position()
        if handle is None:
            return
        atom_ids = tuple(self.selected_atom_ids)
        atom = next((item for item in self.scene.atoms if item.id == atom_ids[-1]), None)
        if atom is None:
            return
        anchor_screen = self.camera.world_to_screen(
            (float(atom.x), float(atom.y), float(atom.z)),
            max(1, self.width()),
            max(1, self.height()),
        )
        if anchor_screen is None:
            return
        try:
            value = self._geometry_measurement(self.geometry_operation, atom_ids)
        except (KeyError, ValueError):
            return
        label = {
            "bond_length": f"{value:.2f} Å",
            "bond_angle": f"{value:.1f}°",
            "dihedral": f"{value:.1f}°",
        }[self.geometry_operation]
        anchor = QPoint(round(anchor_screen[0]), round(anchor_screen[1]))
        painter = QPainter()
        if painter.begin(self):
            painter.setRenderHint(QPainter.RenderHint.Antialiasing)
            painter.setPen(QPen(QColor(70, 205, 245, 235), 2))
            painter.setBrush(QBrush(QColor(26, 123, 170, 230)))
            painter.drawLine(anchor, handle)
            painter.drawEllipse(handle, 7, 7)
            painter.drawText(handle + QPoint(10, -8), label)
            painter.end()

    def set_scene(
        self,
        scene: SceneMolecule3D | None,
        *,
        preserve_interaction: bool = False,
    ) -> None:
        self.cancel_live_mm_playback(clear_telemetry=scene is None)
        self._clear_build_preview()
        self.scene = scene or SceneMolecule3D(atoms=(), bonds=(), source="")
        self._scene_diagnostics = {
            "stage": "scene_ready" if self.scene.atoms else "scene_empty",
            "atom_count": len(self.scene.atoms),
            "bond_count": len(self.scene.bonds),
            "gl_context_ready": self._gl_ready,
            "upload_status": "pending" if self.scene.atoms else "not_requested",
            "render_status": "pending" if self.scene.atoms else "not_requested",
            "error_code": None,
            "status_message": "",
        }
        if preserve_interaction and self._geometry_drag_state is not None:
            self._geometry_drag_state.latest = {
                int(atom.id): (float(atom.x), float(atom.y), float(atom.z))
                for atom in self.scene.atoms
                if isinstance(atom.id, int)
            }
        else:
            self.selected_atom_id = None
            self.selected_atom_ids.clear()
            self.selected_bond_id = None
            self._geometry_drag_state = None
            self._derived_hydrogen_drag_state = None
            self._measure_atom_ids.clear()
        if self._auto_fit:
            self._fit_scene_to_viewport()
        self._mesh_dirty = True
        self.update()

    def _clear_build_preview(self) -> None:
        self._build_drag_state = None
        if self._build_preview is None:
            return
        self._build_preview = None
        self._mesh_dirty = True
        self.buildPreviewChanged.emit(None)
        self.update()

    def _scene_with_build_preview(self) -> SceneMolecule3D:
        preview = self._build_preview
        if preview is None:
            return self.scene
        atoms = list(self.scene.atoms)
        bonds = list(self.scene.bonds)
        if preview.target_atom_id is not None:
            target = next(
                (atom for atom in self.scene.atoms if atom.id == preview.target_atom_id),
                None,
            )
            if target is None:
                return self.scene
            preview_atom_id = f"preview-snap:{target.id}"
            atoms.append(
                SceneAtom3D(
                    preview_atom_id,
                    target.element,
                    target.x,
                    target.y,
                    target.z,
                    target.radius + 0.08,
                    "#5FCFF2",
                    generated=True,
                    parent_atom_id=(target.id if isinstance(target.id, int) else None),
                )
            )
            has_bond = any(
                {bond.a1_id, bond.a2_id}
                == {preview.anchor_atom_id, preview.target_atom_id}
                for bond in self.scene.bonds
            )
            if not has_bond:
                bonds.append(
                    SceneBond3D(
                        "preview-build-bond",
                        preview.anchor_atom_id,
                        preview_atom_id,
                        order=self.build_bond_order[0],
                        is_aromatic=self.build_bond_order[1],
                    )
                )
        else:
            preview_atom_id = "preview-build-end"
            atoms.append(
                SceneAtom3D(
                    preview_atom_id,
                    preview.element,
                    *preview.end_position,
                    0.62,
                    "#5FCFF2",
                )
            )
            bonds.append(
                SceneBond3D(
                    "preview-build-bond",
                    preview.anchor_atom_id,
                    preview_atom_id,
                    order=self.build_bond_order[0],
                    is_aromatic=self.build_bond_order[1],
                )
            )
        return SceneMolecule3D(
            atoms=tuple(atoms),
            bonds=tuple(bonds),
            source=self.scene.source,
        )

    def _commit_pending_build_click(self) -> None:
        self._build_click_timer.stop()
        pending = self._pending_build_click
        self._pending_build_click = None
        if pending is not None:
            self.buildPointRequested.emit((pending.element, pending.position))

    def _cancel_pending_build_click(self) -> None:
        self._build_click_timer.stop()
        self._pending_build_click = None

    def _set_pending_build_click(self, element: str, point: tuple[float, float, float]) -> None:
        self._pending_build_click = _PendingBuildClick(element, point)
        self._build_click_timer.start()

    def _update_build_preview(self, position: QPoint) -> None:
        state = self._build_drag_state
        if state is None:
            return
        hit = self._pick_at(position)
        target_id: int | None = None
        end_position: tuple[float, float, float] | None = None
        if (
            hit is not None
            and hit.kind == "atom"
            and isinstance(hit.primitive_id, int)
            and hit.primitive_id != state.anchor_atom_id
        ):
            target = next(
                (atom for atom in self.scene.atoms if atom.id == hit.primitive_id),
                None,
            )
            if target is not None:
                target_id = hit.primitive_id
                end_position = (float(target.x), float(target.y), float(target.z))
        if end_position is None:
            try:
                end_position = self.camera.screen_to_world(
                    position.x(),
                    position.y(),
                    max(1, self.width()),
                    max(1, self.height()),
                    plane_point=state.anchor_position,
                )
            except ValueError as exc:
                self.statusChanged.emit(str(exc))
                return
        state.current_position = position
        preview = _BuildPreviewState(
            state.anchor_atom_id,
            state.anchor_position,
            end_position,
            self.build_element,
            target_id,
        )
        if preview == self._build_preview:
            return
        self._build_preview = preview
        self._mesh_dirty = True
        self.buildPreviewChanged.emit(preview)
        self.update()

    def set_build_bond_order(self, order: int, is_aromatic: bool = False) -> None:
        normalized_order = int(order)
        if normalized_order not in {1, 2, 3}:
            raise ValueError(f"Orden de enlace 3D no soportado: {order}")
        value = (normalized_order, bool(is_aromatic))
        if value == self.build_bond_order:
            return
        self.build_bond_order = value
        if self._build_preview is not None:
            self._mesh_dirty = True
            self.update()

    def set_edit_mode(self, enabled: bool) -> None:
        self._cancel_pending_build_click()
        self.edit_mode = bool(enabled)
        self.editModeChanged.emit(self.edit_mode)
        self._geometry_drag_state = None
        self._derived_hydrogen_drag_state = None
        self._clear_build_preview()
        self._dragging = False

    def set_edit_tool(self, tool: str) -> None:
        self._cancel_pending_build_click()
        value = str(tool).strip().lower()
        if value not in {"select", "move", "build", "bond", "measure", "geometry"}:
            raise ValueError(f"Herramienta 3D no soportada: {tool}")
        self.edit_tool = value
        self._measure_atom_ids.clear()
        self._geometry_drag_state = None
        self._derived_hydrogen_drag_state = None
        self._clear_build_preview()
        self._dragging = False
        self.editToolChanged.emit(value)
        if value == "geometry":
            self.statusChanged.emit(
                "Selección geométrica: elija 2/3/4 átomos en orden con Shift y arrastre el handle."
            )

    def set_geometry_operation(self, operation: str) -> None:
        value = str(operation).strip().lower()
        if value not in {"bond_length", "bond_angle", "dihedral"}:
            raise ValueError(f"Magnitud geométrica no soportada: {operation}")
        self.geometry_operation = value
        self.update()

    def _record_measurement_pick(self, hit: PickResult | None) -> None:
        if hit is None or hit.kind != "atom" or not isinstance(hit.primitive_id, int):
            self.statusChanged.emit("Seleccione átomos documentales para medir.")
            return
        atom_id = hit.primitive_id
        if self._measure_atom_ids and self._measure_atom_ids[-1] == atom_id:
            return
        self._measure_atom_ids.append(atom_id)
        if len(self._measure_atom_ids) > 4:
            self._measure_atom_ids = [atom_id]
        points = {
            int(atom.id): (float(atom.x), float(atom.y), float(atom.z))
            for atom in self.scene.atoms
            if isinstance(atom.id, int)
        }
        coordinates = CoordinateSet3D(points, source="viewport", method="measure")
        ids = self._measure_atom_ids
        try:
            if len(ids) == 1:
                self.statusChanged.emit(f"Medición: átomo {atom_id} seleccionado.")
            elif len(ids) == 2:
                value = distance(coordinates, ids[0], ids[1])
                self.statusChanged.emit(
                    f"Distancia {ids[0]}–{ids[1]}: {value:.2f} Å"
                )
            elif len(ids) == 3:
                value = angle_degrees(coordinates, ids[0], ids[1], ids[2])
                self.statusChanged.emit(
                    f"Ángulo {ids[0]}–{ids[1]}–{ids[2]}: {value:.1f}°"
                )
            else:
                value = dihedral_degrees(
                    coordinates, ids[0], ids[1], ids[2], ids[3]
                )
                self.statusChanged.emit(
                    f"Dihedro {ids[0]}–{ids[1]}–{ids[2]}–{ids[3]}: {value:.1f}°"
                )
                self._measure_atom_ids.clear()
        except (KeyError, ValueError) as exc:
            self._measure_atom_ids.clear()
            self.statusChanged.emit(f"No se pudo medir la geometría: {exc}")

    def set_build_element(self, element: str) -> None:
        self.build_element = str(element).strip().capitalize() or "C"

    def set_representation(self, representation: str) -> None:
        value = normalize_representation(representation)
        if value == self.representation:
            return
        self.representation = value
        self._mesh_dirty = True
        self.update()

    def _apply_geometry_positions(
        self, positions: dict[int, tuple[float, float, float]]
    ) -> None:
        old_positions = {
            int(atom.id): (float(atom.x), float(atom.y), float(atom.z))
            for atom in self.scene.atoms
            if isinstance(atom.id, int)
        }
        parent_deltas = {
            atom_id: tuple(
                positions[atom_id][index] - old_positions[atom_id][index]
                for index in range(3)
            )
            for atom_id in positions.keys() & old_positions.keys()
        }
        atoms = tuple(
            replace(
                atom,
                x=positions[int(atom.id)][0],
                y=positions[int(atom.id)][1],
                z=positions[int(atom.id)][2],
            )
            if isinstance(atom.id, int) and int(atom.id) in positions
            else replace(
                atom,
                x=float(atom.x) + parent_deltas[atom.parent_atom_id][0],
                y=float(atom.y) + parent_deltas[atom.parent_atom_id][1],
                z=float(atom.z) + parent_deltas[atom.parent_atom_id][2],
            )
            if atom.generated and atom.parent_atom_id in parent_deltas
            else atom
            for atom in self.scene.atoms
        )
        self.scene = replace(self.scene, atoms=atoms)
        self._mesh_dirty = True
        self.update()

    def _apply_derived_hydrogen_position(
        self, atom_id: str, position: tuple[float, float, float]
    ) -> None:
        self.scene = replace(
            self.scene,
            atoms=tuple(
                replace(atom, x=position[0], y=position[1], z=position[2])
                if atom.id == atom_id and atom.generated
                else atom
                for atom in self.scene.atoms
            ),
        )
        self._mesh_dirty = True
        self.update()

    def clear_scene(self) -> None:
        self.set_scene(None)

    def set_projection(self, projection: CameraProjection) -> None:
        self.camera.set_projection(projection)
        if self._auto_fit:
            self._fit_scene_to_viewport()
        self.update()

    def fit_scene(self) -> None:
        self._auto_fit = True
        self._fit_scene_to_viewport()
        self.update()

    def _fit_scene_to_viewport(self) -> None:
        if self.width() <= 0 or self.height() <= 0:
            return
        self.camera.fit_scene(
            self.scene,
            viewport_width=self.width(),
            viewport_height=self.height(),
        )

    def resizeEvent(self, e: QResizeEvent | None) -> None:
        if e is None:
            return
        super().resizeEvent(e)
        if self._auto_fit:
            self._fit_scene_to_viewport()

    def _mark_manual_navigation(self) -> None:
        self._auto_fit = False

    def _set_selection(
        self,
        hit: PickResult | None,
        *,
        extend: bool = False,
    ) -> None:
        if hit is None:
            self.selected_atom_id = None
            self.selected_atom_ids.clear()
            self.selected_bond_id = None
            self.selectionCleared.emit()
        elif hit.kind == "atom":
            atom_id = hit.primitive_id
            if extend and isinstance(atom_id, int):
                if atom_id in self.selected_atom_ids:
                    self.selected_atom_ids.remove(atom_id)
                else:
                    self.selected_atom_ids.append(atom_id)
                self.selected_atom_id = (
                    atom_id
                    if atom_id in self.selected_atom_ids
                    else (self.selected_atom_ids[-1] if self.selected_atom_ids else None)
                )
            else:
                self.selected_atom_ids = [atom_id] if isinstance(atom_id, int) else []
                self.selected_atom_id = atom_id
            self.selected_bond_id = None
            self.primitiveSelected.emit(hit.kind, hit.primitive_id)
        else:
            self.selected_atom_id = None
            self.selected_atom_ids.clear()
            self.selected_bond_id = hit.primitive_id
            self.primitiveSelected.emit(hit.kind, hit.primitive_id)
        self._mesh_dirty = True
        self.update()

    def _geometry_handle_screen_position(self) -> QPoint | None:
        expected = {"bond_length": 2, "bond_angle": 3, "dihedral": 4}.get(
            self.geometry_operation
        )
        if (
            self.edit_tool != "geometry"
            or expected is None
            or len(self.selected_atom_ids) != expected
        ):
            return None
        atom = next(
            (atom for atom in self.scene.atoms if atom.id == self.selected_atom_ids[-1]),
            None,
        )
        if atom is None:
            return None
        screen = self.camera.world_to_screen(
            (float(atom.x), float(atom.y), float(atom.z)),
            max(1, self.width()),
            max(1, self.height()),
        )
        if screen is None:
            return None
        return QPoint(round(screen[0] + 28.0), round(screen[1] - 24.0))

    def _geometry_measurement(self, operation: str, atom_ids: tuple[int, ...]) -> float:
        positions = {
            int(atom.id): (float(atom.x), float(atom.y), float(atom.z))
            for atom in self.scene.atoms
            if isinstance(atom.id, int)
        }
        coordinates = CoordinateSet3D(
            positions=positions,
            source="viewport",
            method="geometry-handle",
        )
        if operation == "bond_length":
            return distance(coordinates, *atom_ids)
        if operation == "bond_angle":
            return angle_degrees(coordinates, *atom_ids)
        return dihedral_degrees(coordinates, *atom_ids)

    def _begin_geometry_handle_drag(self, position: QPoint) -> bool:
        handle = self._geometry_handle_screen_position()
        if handle is None or (handle - position).manhattanLength() > 14:
            return False
        operation = self.geometry_operation
        atom_ids = tuple(self.selected_atom_ids)
        before = {
            int(atom.id): (float(atom.x), float(atom.y), float(atom.z))
            for atom in self.scene.atoms
            if isinstance(atom.id, int)
        }
        try:
            initial_value = self._geometry_measurement(operation, atom_ids)
        except (KeyError, ValueError) as exc:
            self.statusChanged.emit(f"Handle geométrico no disponible: {exc}")
            return False
        if not atom_ids or atom_ids[-1] not in before:
            return False
        self._geometry_drag_state = _GeometryDragState(
            before=dict(before),
            latest=dict(before),
            selected_atom_ids=atom_ids,
            anchor_atom_id=atom_ids[-1],
            start_position=position,
            geometry_operation=operation,
            geometry_atom_ids=atom_ids,
            geometry_start_value=initial_value,
        )
        self.geometryHandleStarted.emit(operation, atom_ids)
        self._mark_manual_navigation()
        return True

    def _connected_component(self, atom_id: int) -> set[int]:
        atom_ids = {
            int(atom.id) for atom in self.scene.atoms if isinstance(atom.id, int)
        }
        if atom_id not in atom_ids:
            return set()
        neighbors: dict[int, set[int]] = {value: set() for value in atom_ids}
        for bond in self.scene.bonds:
            if (
                isinstance(bond.a1_id, int)
                and isinstance(bond.a2_id, int)
                and bond.a1_id in atom_ids
                and bond.a2_id in atom_ids
            ):
                neighbors[bond.a1_id].add(bond.a2_id)
                neighbors[bond.a2_id].add(bond.a1_id)
        component = {atom_id}
        pending = [atom_id]
        while pending:
            current = pending.pop()
            for neighbor in neighbors[current] - component:
                component.add(neighbor)
                pending.append(neighbor)
        return component

    def _pick_at(self, position: QPoint) -> PickResult | None:
        try:
            return pick_scene(
                self.scene,
                self.camera,
                self.width(),
                self.height(),
                position.x(),
                position.y(),
                self.representation,
            )
        except ValueError as exc:
            self.statusChanged.emit(str(exc))
            return None

    def _pick(self, position: QPoint) -> None:
        self._set_selection(self._pick_at(position))

    def mousePressEvent(self, a0: QMouseEvent | None) -> None:
        if a0 is None:
            return
        event = a0
        if self._pending_build_click is not None:
            self._commit_pending_build_click()
        self._last_mouse_position = event.position().toPoint()
        self._dragging = False
        self._camera_drag_origin_empty = False
        self._orbit_guide_active = False
        self._geometry_drag_state = None
        self._derived_hydrogen_drag_state = None
        self._geometry_edit_notified = False
        self._clear_build_preview()
        self._press_hit = None
        left_edit = self.edit_mode and event.button() == Qt.MouseButton.LeftButton
        depth_edit = (
            self.edit_mode
            and self.edit_tool == "move"
            and (
                event.button() == Qt.MouseButton.MiddleButton
                or (
                    event.button() == Qt.MouseButton.LeftButton
                    and bool(event.modifiers() & Qt.KeyboardModifier.AltModifier)
                )
            )
        )
        if left_edit or depth_edit:
            position = event.position().toPoint()
            geometry_handle_started = (
                left_edit
                and self.edit_tool == "geometry"
                and self._begin_geometry_handle_drag(position)
            )
            if geometry_handle_started:
                hit = None
            else:
                hit = self._pick_at(position)
                self._press_hit = hit
                self._camera_drag_origin_empty = hit is None and event.button() == Qt.MouseButton.LeftButton
            if left_edit and not geometry_handle_started:
                extend = bool(event.modifiers() & Qt.KeyboardModifier.ShiftModifier)
                if self.edit_tool == "geometry" and not extend:
                    if hit is None:
                        self._set_selection(None)
                    elif hit.primitive_id not in self.selected_atom_ids:
                        self._set_selection(hit)
                    else:
                        self.selected_atom_id = hit.primitive_id
                else:
                    self._set_selection(hit, extend=extend)
            elif depth_edit and hit is not None and hit.kind == "atom":
                if (
                    not isinstance(hit.primitive_id, int)
                    or hit.primitive_id not in self.selected_atom_ids
                ):
                    self._set_selection(hit)
            if (
                left_edit
                and self.edit_tool in {"build", "bond"}
                and hit is not None
                and hit.kind == "atom"
            ):
                anchor_id = hit.primitive_id
                anchor = next(
                    (atom for atom in self.scene.atoms if atom.id == anchor_id),
                    None,
                )
                if isinstance(anchor_id, int) and anchor is not None:
                    self._build_drag_state = _BuildDragState(
                        anchor_atom_id=anchor_id,
                        anchor_position=(float(anchor.x), float(anchor.y), float(anchor.z)),
                        start_position=position,
                        current_position=position,
                    )
            if (
                left_edit
                and hit is not None
                and hit.kind == "atom"
                and isinstance(hit.primitive_id, int)
                and event.modifiers() & Qt.KeyboardModifier.ControlModifier
            ):
                component = self._connected_component(hit.primitive_id)
                self.selected_atom_ids = [
                    int(atom.id)
                    for atom in self.scene.atoms
                    if isinstance(atom.id, int) and atom.id in component
                ]
                self._mesh_dirty = True
                self.update()
            if self.edit_tool == "move" and hit is not None and hit.kind == "atom":
                anchor_id = hit.primitive_id
                positions = {
                    int(atom.id): (float(atom.x), float(atom.y), float(atom.z))
                    for atom in self.scene.atoms
                    if isinstance(atom.id, int)
                }
                if isinstance(anchor_id, int) and anchor_id in positions:
                    selected_ids = tuple(
                        atom_id
                        for atom_id in self.selected_atom_ids
                        if atom_id in positions
                    )
                    if anchor_id in selected_ids:
                        self._geometry_drag_state = _GeometryDragState(
                            before=dict(positions),
                            latest=dict(positions),
                            selected_atom_ids=selected_ids,
                            anchor_atom_id=anchor_id,
                            start_position=position,
                            depth_mode=depth_edit,
                            depth_button=event.button(),
                        )
                        self._mark_manual_navigation()
            if self.edit_tool == "move" and hit is not None and hit.kind == "atom":
                derived_atom = next(
                    (
                        atom
                        for atom in self.scene.atoms
                        if atom.id == hit.primitive_id and atom.generated
                    ),
                    None,
                )
                if (
                    derived_atom is not None
                    and isinstance(derived_atom.id, str)
                    and derived_atom.parent_atom_id is not None
                ):
                    position3d = (
                        float(derived_atom.x),
                        float(derived_atom.y),
                        float(derived_atom.z),
                    )
                    self._derived_hydrogen_drag_state = _DerivedHydrogenDragState(
                        atom_id=derived_atom.id,
                        parent_atom_id=derived_atom.parent_atom_id,
                        before=position3d,
                        latest=position3d,
                        start_position=position,
                        depth_mode=depth_edit,
                        depth_button=event.button(),
                    )
                    self._mark_manual_navigation()
        elif event.button() == Qt.MouseButton.LeftButton:
            position = event.position().toPoint()
            hit = self._pick_at(position)
            self._press_hit = hit
            self._camera_drag_origin_empty = hit is None
        self.setFocus()
        event.accept()

    def mouseDoubleClickEvent(self, a0: QMouseEvent | None) -> None:
        if a0 is None:
            return
        event = a0
        if (
            self.edit_mode
            and self.edit_tool in {"build", "bond"}
            and event.button() == Qt.MouseButton.LeftButton
        ):
            position = event.position().toPoint()
            pending = self._pending_build_click
            self._cancel_pending_build_click()
            hit = self._pick_at(position)
            if hit is not None and hit.kind == "atom" and isinstance(hit.primitive_id, int):
                self.buildElementChangeRequested.emit(hit.primitive_id, self.build_element)
            elif self.edit_tool == "build":
                try:
                    world_position = (
                        pending.position
                        if pending is not None
                        else self.camera.screen_to_world(
                            position.x(),
                            position.y(),
                            max(1, self.width()),
                            max(1, self.height()),
                        )
                    )
                except ValueError as exc:
                    self.statusChanged.emit(str(exc))
                else:
                    element = pending.element if pending is not None else self.build_element
                    self.buildChainExtensionRequested.emit((element, world_position))
            self._suppress_next_build_release = True
            self._press_hit = None
            self._last_mouse_position = None
            event.accept()
            return
        super().mouseDoubleClickEvent(event)

    def _orbit_empty_space(self, delta: QPoint) -> None:
        self.camera.orbit(delta.x(), delta.y(), max(1, self.height()))
        self._mark_manual_navigation()
        self._dragging = True
        if not self._orbit_guide_active:
            self.statusChanged.emit("Cámara · rotación desde espacio vacío; los átomos no se mueven.")
        self._orbit_guide_active = True
        self.update()

    def mouseMoveEvent(self, a0: QMouseEvent | None) -> None:
        if a0 is None:
            return
        event = a0
        if self._last_mouse_position is None:
            return
        current = event.position().toPoint()
        delta = current - self._last_mouse_position
        if delta.isNull():
            return
        buttons = event.buttons()
        if (
            buttons & Qt.MouseButton.LeftButton
            and self.edit_mode
            and self.edit_tool in {"build", "bond"}
        ):
            build_state = self._build_drag_state
            if build_state is not None:
                total_delta = current - build_state.start_position
                if abs(total_delta.x()) + abs(total_delta.y()) >= 4:
                    starting_build_drag = not build_state.moved
                    build_state.current_position = current
                    build_state.moved = True
                    self._notify_geometry_edit_started()
                    if starting_build_drag:
                        self.statusChanged.emit("Construcción 3D · previsualización esfera/enlace.")
                    self._update_build_preview(current)
                    self._dragging = True
            elif abs(delta.x()) + abs(delta.y()) >= 4:
                if self._camera_drag_origin_empty:
                    self._orbit_empty_space(delta)
                else:
                    self._dragging = True
            self._last_mouse_position = current
            event.accept()
            return
        if (
            buttons & Qt.MouseButton.LeftButton
            and self.edit_mode
            and self.edit_tool in {"select", "measure"}
        ):
            if abs(delta.x()) + abs(delta.y()) >= 4:
                if self._camera_drag_origin_empty:
                    self._orbit_empty_space(delta)
                else:
                    self._dragging = True
            self._last_mouse_position = current
            event.accept()
            return
        derived_state = self._derived_hydrogen_drag_state
        derived_button = (
            derived_state.depth_button
            if derived_state is not None and derived_state.depth_mode
            else Qt.MouseButton.LeftButton
        )
        if (
            buttons & derived_button
            and self.edit_mode
            and self.edit_tool == "move"
            and derived_state is not None
        ):
            total_delta = current - derived_state.start_position
            if abs(total_delta.x()) + abs(total_delta.y()) >= 4:
                try:
                    if derived_state.depth_mode:
                        world_delta = self.camera.screen_depth_delta_to_world(
                            total_delta.y(),
                            max(1, self.height()),
                            position=derived_state.before,
                        )
                    else:
                        world_delta = self.camera.screen_delta_to_world(
                            total_delta.x(),
                            total_delta.y(),
                            max(1, self.width()),
                            max(1, self.height()),
                            position=derived_state.before,
                        )
                except ValueError as exc:
                    self.statusChanged.emit(str(exc))
                    return
                starting_atom_drag = not derived_state.moved
                if starting_atom_drag:
                    self._notify_geometry_edit_started()
                    self.statusChanged.emit("Geometría · moviendo H derivado con su átomo padre.")
                latest = add(derived_state.before, world_delta)
                derived_state.latest = latest
                derived_state.moved = True
                self._apply_derived_hydrogen_position(derived_state.atom_id, latest)
                self.derivedHydrogenPreviewChanged.emit(
                    derived_state.atom_id,
                    derived_state.parent_atom_id,
                    latest,
                )
                self._dragging = True
            self._last_mouse_position = current
            event.accept()
            return
        state = self._geometry_drag_state
        if (
            buttons & Qt.MouseButton.LeftButton
            and self.edit_mode
            and self.edit_tool == "geometry"
            and state is not None
            and state.geometry_operation is not None
        ):
            total_delta = current - state.start_position
            if abs(total_delta.x()) + abs(total_delta.y()) >= 4:
                if not state.moved:
                    self._notify_geometry_edit_started()
                start_value = state.geometry_start_value
                if start_value is None:
                    self.statusChanged.emit("Handle geométrico sin valor inicial.")
                    return
                if state.geometry_operation == "bond_length":
                    target_value = max(0.1, start_value - total_delta.y() * 0.01)
                    label = f"Distancia: {target_value:.2f} Å"
                elif state.geometry_operation == "bond_angle":
                    target_value = max(0.1, min(179.9, start_value - total_delta.y() * 0.5))
                    label = f"Ángulo: {target_value:.1f}°"
                else:
                    target_value = (start_value + total_delta.x() * 0.5 + 180.0) % 360.0 - 180.0
                    label = f"Dihedro: {target_value:.1f}°"
                state.moved = True
                self.statusChanged.emit(label)
                self.geometryHandlePreviewRequested.emit(
                    state.geometry_operation,
                    state.geometry_atom_ids,
                    target_value,
                )
                self.update()
            self._last_mouse_position = current
            self._dragging = state.moved
            event.accept()
            return
        state = self._geometry_drag_state
        drag_button = (
            state.depth_button
            if state is not None and state.depth_mode
            else Qt.MouseButton.LeftButton
        )
        if (
            buttons & drag_button
            and self.edit_mode
            and self.edit_tool == "move"
            and state is not None
        ):
            total_delta = current - state.start_position
            if total_delta.isNull():
                return
            anchor = state.before[state.anchor_atom_id]
            try:
                if state.depth_mode:
                    world_delta = self.camera.screen_depth_delta_to_world(
                        total_delta.y(),
                        max(1, self.height()),
                        position=anchor,
                    )
                else:
                    world_delta = self.camera.screen_delta_to_world(
                        total_delta.x(),
                        total_delta.y(),
                        max(1, self.width()),
                        max(1, self.height()),
                        position=anchor,
                    )
            except ValueError as exc:
                self.statusChanged.emit(str(exc))
                return
            starting_atom_drag = not state.moved
            if starting_atom_drag:
                self._notify_geometry_edit_started()
                self.statusChanged.emit("Geometría · moviendo átomo/fragmento en 3D.")
            positions = dict(state.before)
            for atom_id in state.selected_atom_ids:
                origin = state.before[atom_id]
                positions[atom_id] = add(origin, world_delta)
            state.latest = positions
            state.moved = True
            self._apply_geometry_positions(positions)
            self.geometryPreviewChanged.emit(
                dict(state.before), dict(positions), state.selected_atom_ids
            )
            self._last_mouse_position = current
            self._dragging = True
            event.accept()
            return
        self._last_mouse_position = current
        if buttons & Qt.MouseButton.LeftButton:
            if self._camera_drag_origin_empty:
                self._orbit_empty_space(delta)
            else:
                self.camera.orbit(delta.x(), delta.y(), max(1, self.height()))
                self._mark_manual_navigation()
                self._dragging = True
                self.statusChanged.emit("Cámara · orbitar; los átomos no se editan.")
        elif buttons & Qt.MouseButton.RightButton:
            self._orbit_guide_active = False
            self.camera.pan(delta.x(), delta.y(), max(1, self.height()))
            self._mark_manual_navigation()
            self._dragging = True
            self.statusChanged.emit("Cámara · desplazar.")
        if self._dragging:
            self.update()
        event.accept()

    def mouseReleaseEvent(self, a0: QMouseEvent | None) -> None:
        if a0 is None:
            return
        event = a0
        if event.button() == Qt.MouseButton.LeftButton and self._suppress_next_build_release:
            self._suppress_next_build_release = False
            had_orbit_guide = self._orbit_guide_active
            self._last_mouse_position = None
            self._dragging = False
            self._camera_drag_origin_empty = False
            self._orbit_guide_active = False
            self._geometry_edit_notified = False
            self._geometry_drag_state = None
            self._derived_hydrogen_drag_state = None
            self._clear_build_preview()
            self._press_hit = None
            if had_orbit_guide:
                self.update()
            event.accept()
            return
        position = event.position().toPoint()
        if event.button() in {Qt.MouseButton.LeftButton, Qt.MouseButton.MiddleButton} and (
            event.button() == Qt.MouseButton.LeftButton
            or (self.edit_mode and self.edit_tool == "move")
        ):
            state = self._geometry_drag_state
            derived_state = self._derived_hydrogen_drag_state
            if self.edit_mode and self.edit_tool in {"build", "bond"}:
                build_state = self._build_drag_state
                if build_state is not None and build_state.moved:
                    build_edit_requested = False
                    hit = self._pick_at(position)
                    if hit is not None and hit.kind == "atom" and isinstance(
                        hit.primitive_id, int
                    ):
                        if hit.primitive_id != build_state.anchor_atom_id:
                            self.buildBondRequested.emit(
                                (build_state.anchor_atom_id, hit.primitive_id)
                            )
                            build_edit_requested = True
                    elif self.edit_tool == "build":
                        try:
                            world_point = self.camera.screen_to_world(
                                position.x(),
                                position.y(),
                                max(1, self.width()),
                                max(1, self.height()),
                                plane_point=build_state.anchor_position,
                            )
                        except ValueError as exc:
                            self.statusChanged.emit(str(exc))
                        else:
                            self.buildAnchorSelected.emit(build_state.anchor_atom_id)
                            self.buildPointRequested.emit(
                                (self.build_element, world_point)
                            )
                            build_edit_requested = True
                    else:
                        self.statusChanged.emit("La herramienta Enlazar requiere otro átomo.")
                    if not build_edit_requested:
                        self.geometryEditCancelled.emit()
                elif not self._dragging:
                    hit = self._press_hit
                    if hit is not None and hit.kind == "atom" and isinstance(
                        hit.primitive_id, int
                    ):
                        self.buildAnchorSelected.emit(hit.primitive_id)
                    elif hit is None and self.edit_tool == "build":
                        try:
                            world_point = self.camera.screen_to_world(
                                position.x(),
                                position.y(),
                                max(1, self.width()),
                                max(1, self.height()),
                            )
                        except ValueError as exc:
                            self.statusChanged.emit(str(exc))
                        else:
                            self._set_pending_build_click(self.build_element, world_point)
            elif self.edit_mode and self.edit_tool == "geometry" and (
                state is not None
                and state.geometry_operation is not None
                and state.moved
            ):
                self.geometryEditCommitted.emit(
                    dict(state.before),
                    dict(state.latest),
                    state.selected_atom_ids,
                )
            elif self.edit_mode and self.edit_tool == "measure" and not self._dragging:
                self._record_measurement_pick(self._press_hit)
            elif (
                self.edit_mode
                and self.edit_tool == "move"
                and derived_state is not None
                and derived_state.moved
            ):
                self.derivedHydrogenEditRequested.emit(
                    (
                        derived_state.atom_id,
                        derived_state.parent_atom_id,
                        derived_state.latest,
                    )
                )
            elif (
                self.edit_mode
                and self.edit_tool == "move"
                and state is not None
                and state.moved
            ):
                self.geometryEditCommitted.emit(
                    dict(state.before),
                    dict(state.latest),
                    state.selected_atom_ids,
                )
            elif not self.edit_mode and not self._dragging:
                self._pick(position)
        self._last_mouse_position = None
        self._dragging = False
        had_orbit_guide = self._orbit_guide_active
        self._camera_drag_origin_empty = False
        self._orbit_guide_active = False
        self._geometry_edit_notified = False
        self._geometry_drag_state = None
        self._derived_hydrogen_drag_state = None
        self._clear_build_preview()
        self._press_hit = None
        if had_orbit_guide:
            self.update()
        event.accept()

    def wheelEvent(self, a0: QWheelEvent | None) -> None:
        if a0 is None:
            return
        event = a0
        self.camera.zoom(event.angleDelta().y())
        self._mark_manual_navigation()
        self.update()
        event.accept()

    def keyPressEvent(self, a0) -> None:
        if a0 is None:
            return
        event = a0
        if event.key() == Qt.Key.Key_Escape:
            derived_state = self._derived_hydrogen_drag_state
            if self.edit_mode and derived_state is not None:
                self._apply_derived_hydrogen_position(
                    derived_state.atom_id,
                    derived_state.before,
                )
                self._derived_hydrogen_drag_state = None
                self._dragging = False
                self._last_mouse_position = None
                self._press_hit = None
                self.geometryEditCancelled.emit()
                event.accept()
                return
            state = self._geometry_drag_state
            if self.edit_mode and state is not None:
                self._apply_geometry_positions(state.before)
                self._geometry_drag_state = None
                self._dragging = False
                self._last_mouse_position = None
                self._press_hit = None
                self.geometryEditCancelled.emit()
                event.accept()
                return
            if self.edit_mode and self.edit_tool == "measure":
                self._measure_atom_ids.clear()
                self.statusChanged.emit("Medición cancelada.")
                event.accept()
                return
            if self.edit_mode and self.edit_tool in {"build", "bond"}:
                self._cancel_pending_build_click()
                self._press_hit = None
                self._clear_build_preview()
                self.buildCancelled.emit()
                event.accept()
                return
        if event.key() in {Qt.Key.Key_Delete, Qt.Key.Key_Backspace}:
            selected = next(
                (
                    atom
                    for atom in self.scene.atoms
                    if atom.id == self.selected_atom_id and atom.generated
                ),
                None,
            )
            if selected is not None:
                self.statusChanged.emit(
                    "No se puede borrar un H derivado; promuévalo como átomo explícito primero."
                )
                event.accept()
                return
        if event.modifiers() == Qt.KeyboardModifier.NoModifier:
            shortcut_tools = {
                Qt.Key.Key_S: "select",
                Qt.Key.Key_D: "build",
                Qt.Key.Key_M: "move",
                Qt.Key.Key_B: "bond",
                Qt.Key.Key_R: "measure",
                Qt.Key.Key_G: "geometry",
            }
            tool = shortcut_tools.get(cast(Qt.Key, event.key()))
            if tool is not None:
                self.set_edit_mode(True)
                self.set_edit_tool(tool)
                event.accept()
                return
        if event.key() == Qt.Key.Key_F:
            self.fit_scene()
            event.accept()
            return
        if event.key() == Qt.Key.Key_P:
            projection = (
                CameraProjection.ORTHOGRAPHIC
                if self.camera.projection == CameraProjection.PERSPECTIVE
                else CameraProjection.PERSPECTIVE
            )
            self.set_projection(projection)
            event.accept()
            return
        super().keyPressEvent(event)
