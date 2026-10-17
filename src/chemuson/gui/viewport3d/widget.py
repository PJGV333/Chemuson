"""Widget Qt6 del visor: eventos, cámara y selección local."""

from __future__ import annotations

from dataclasses import dataclass, replace
from typing import Any, cast

from PyQt6.QtCore import QPoint, Qt, pyqtSignal
from PyQt6.QtGui import (
    QHideEvent,
    QMouseEvent,
    QOpenGLContext,
    QResizeEvent,
    QShowEvent,
    QSurfaceFormat,
    QWheelEvent,
)
from PyQt6.QtOpenGLWidgets import QOpenGLWidget

from chemuson.geometry3d import SceneMolecule3D

from .camera import CameraProjection, CameraState
from .geometry import add
from .renderer import GL_COLOR_BUFFER_BIT, GL_DEPTH_BUFFER_BIT, MoleculeRenderer
from .picking import PickResult, pick_scene


@dataclass
class _GeometryDragState:
    before: dict[int, tuple[float, float, float]]
    latest: dict[int, tuple[float, float, float]]
    selected_atom_ids: tuple[int, ...]
    anchor_atom_id: int
    start_position: QPoint
    moved: bool = False


class MolecularViewportWidget(QOpenGLWidget):
    """Viewport 3D con estado de escena y selección estrictamente local."""

    primitiveSelected = pyqtSignal(str, object)
    selectionCleared = pyqtSignal()
    statusChanged = pyqtSignal(str)
    geometryPreviewChanged = pyqtSignal(object, object, object)
    geometryEditCommitted = pyqtSignal(object, object, object)
    geometryEditCancelled = pyqtSignal()
    buildCancelled = pyqtSignal()
    buildPointRequested = pyqtSignal(object)
    buildAnchorSelected = pyqtSignal(object)

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
        self.camera = CameraState()
        self.scene = SceneMolecule3D(atoms=(), bonds=(), source="")
        self.selected_atom_id: int | str | None = None
        self.selected_atom_ids: set[int] = set()
        self.selected_bond_id: int | str | None = None
        self.edit_mode = False
        self.edit_tool = "move"
        self.build_element = "C"
        self._renderer = MoleculeRenderer()
        self._last_mouse_position: QPoint | None = None
        self._dragging = False
        self._geometry_drag_state: _GeometryDragState | None = None
        self._press_hit: PickResult | None = None
        self._gl_ready = False
        self._mesh_dirty = True
        self._auto_fit = True
        self._active_context: QOpenGLContext | None = None
        self._connected_context: QOpenGLContext | None = None
        self.statusChanged.emit("Esperando contexto OpenGL…")

    @property
    def gl_ready(self) -> bool:
        return self._gl_ready

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
        self.statusChanged.emit(message)
        if ready:
            self._mesh_dirty = True

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
            return
        gl = self._renderer._gl
        if gl is None:
            return
        gl.glClearColor(0.055, 0.065, 0.08, 1.0)
        gl.glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT)
        if self._mesh_dirty:
            self._mesh_dirty = not self._renderer.upload(
                self.scene,
                self.selected_atom_id,
                self.selected_bond_id,
            )
        self._renderer.render(self.camera, self.width(), self.height())

    def set_scene(
        self,
        scene: SceneMolecule3D | None,
        *,
        preserve_interaction: bool = False,
    ) -> None:
        self.scene = scene or SceneMolecule3D(atoms=(), bonds=(), source="")
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
        if self._auto_fit:
            self._fit_scene_to_viewport()
        self._mesh_dirty = True
        self.update()

    def set_edit_mode(self, enabled: bool) -> None:
        self.edit_mode = bool(enabled)
        self._geometry_drag_state = None
        self._dragging = False

    def set_edit_tool(self, tool: str) -> None:
        value = str(tool).strip().lower()
        if value not in {"move", "build"}:
            raise ValueError(f"Herramienta 3D no soportada: {tool}")
        self.edit_tool = value
        self._geometry_drag_state = None
        self._dragging = False

    def set_build_element(self, element: str) -> None:
        self.build_element = str(element).strip().capitalize() or "C"

    def _apply_geometry_positions(
        self, positions: dict[int, tuple[float, float, float]]
    ) -> None:
        atoms = tuple(
            replace(
                atom,
                x=positions[int(atom.id)][0],
                y=positions[int(atom.id)][1],
                z=positions[int(atom.id)][2],
            )
            if isinstance(atom.id, int) and int(atom.id) in positions
            else atom
            for atom in self.scene.atoms
        )
        self.scene = replace(self.scene, atoms=atoms)
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
                    self.selected_atom_ids.add(atom_id)
                self.selected_atom_id = (
                    atom_id
                    if atom_id in self.selected_atom_ids
                    else next(iter(self.selected_atom_ids), None)
                )
            else:
                self.selected_atom_ids = (
                    {atom_id} if isinstance(atom_id, int) else set()
                )
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
        self._last_mouse_position = event.position().toPoint()
        self._dragging = False
        self._geometry_drag_state = None
        self._press_hit = None
        if self.edit_mode and event.button() == Qt.MouseButton.LeftButton:
            position = event.position().toPoint()
            hit = self._pick_at(position)
            self._press_hit = hit
            extend = bool(event.modifiers() & Qt.KeyboardModifier.ShiftModifier)
            self._set_selection(hit, extend=extend)
            if (
                hit is not None
                and hit.kind == "atom"
                and isinstance(hit.primitive_id, int)
                and event.modifiers() & Qt.KeyboardModifier.ControlModifier
            ):
                self.selected_atom_ids = self._connected_component(hit.primitive_id)
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
                        sorted(self.selected_atom_ids.intersection(positions))
                    )
                    if anchor_id in selected_ids:
                        self._geometry_drag_state = _GeometryDragState(
                            before=dict(positions),
                            latest=dict(positions),
                            selected_atom_ids=selected_ids,
                            anchor_atom_id=anchor_id,
                            start_position=position,
                        )
                        self._mark_manual_navigation()
        self.setFocus()
        event.accept()

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
        state = self._geometry_drag_state
        if (
            buttons & Qt.MouseButton.LeftButton
            and self.edit_mode
            and self.edit_tool == "move"
            and state is not None
        ):
            total_delta = current - state.start_position
            if total_delta.isNull():
                return
            anchor = state.before[state.anchor_atom_id]
            try:
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
            positions = dict(state.latest)
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
            self.camera.orbit(delta.x(), delta.y(), max(1, self.height()))
            self._mark_manual_navigation()
            self._dragging = True
        elif buttons & Qt.MouseButton.RightButton:
            self.camera.pan(delta.x(), delta.y(), max(1, self.height()))
            self._mark_manual_navigation()
            self._dragging = True
        if self._dragging:
            self.update()
        event.accept()

    def mouseReleaseEvent(self, a0: QMouseEvent | None) -> None:
        if a0 is None:
            return
        event = a0
        position = event.position().toPoint()
        if event.button() == Qt.MouseButton.LeftButton:
            state = self._geometry_drag_state
            if self.edit_mode and self.edit_tool == "build" and not self._dragging:
                hit = self._press_hit
                if hit is not None and hit.kind == "atom" and isinstance(
                    hit.primitive_id, int
                ):
                    self.buildAnchorSelected.emit(hit.primitive_id)
                elif hit is None:
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
                        self.buildPointRequested.emit(
                            (self.build_element, world_point)
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
        self._geometry_drag_state = None
        self._press_hit = None
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
            if self.edit_mode and self.edit_tool == "build":
                self._press_hit = None
                self.buildCancelled.emit()
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
