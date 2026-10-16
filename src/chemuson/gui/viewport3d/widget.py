"""Widget Qt6 del visor: eventos, cámara y selección local."""

from __future__ import annotations

from PyQt6.QtCore import QPoint, Qt, pyqtSignal
from PyQt6.QtGui import QMouseEvent, QOpenGLContext, QResizeEvent, QSurfaceFormat, QWheelEvent
from PyQt6.QtOpenGLWidgets import QOpenGLWidget

from chemuson.geometry3d import SceneMolecule3D

from .camera import CameraProjection, CameraState
from .renderer import GL_COLOR_BUFFER_BIT, GL_DEPTH_BUFFER_BIT, MoleculeRenderer
from .picking import PickResult, pick_scene


class MolecularViewportWidget(QOpenGLWidget):
    """Viewport 3D con estado de escena y selección estrictamente local."""

    primitiveSelected = pyqtSignal(str, object)
    selectionCleared = pyqtSignal()
    statusChanged = pyqtSignal(str)

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
        self.selected_bond_id: int | str | None = None
        self._renderer = MoleculeRenderer()
        self._last_mouse_position: QPoint | None = None
        self._dragging = False
        self._gl_ready = False
        self._mesh_dirty = True
        self._auto_fit = True
        self._active_context: QOpenGLContext | None = None
        self._connected_context: QOpenGLContext | None = None
        self.statusChanged.emit("Esperando contexto OpenGL…")

    @property
    def gl_ready(self) -> bool:
        return self._gl_ready

    def initializeGL(self) -> None:
        context = self.context()
        if context is not None and context != self._connected_context:
            context.aboutToBeDestroyed.connect(self._cleanup_gl)
            self._connected_context = context
        self._active_context = context
        ready, message = self._renderer.initialize()
        self._gl_ready = ready
        self.statusChanged.emit(message)
        if ready:
            self._mesh_dirty = True

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

    def set_scene(self, scene: SceneMolecule3D | None) -> None:
        self.scene = scene or SceneMolecule3D(atoms=(), bonds=(), source="")
        self.selected_atom_id = None
        self.selected_bond_id = None
        if self._auto_fit:
            self._fit_scene_to_viewport()
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

    def _set_selection(self, hit: PickResult | None) -> None:
        if hit is None:
            self.selected_atom_id = None
            self.selected_bond_id = None
            self.selectionCleared.emit()
        elif hit.kind == "atom":
            self.selected_atom_id = hit.primitive_id
            self.selected_bond_id = None
            self.primitiveSelected.emit(hit.kind, hit.primitive_id)
        else:
            self.selected_atom_id = None
            self.selected_bond_id = hit.primitive_id
            self.primitiveSelected.emit(hit.kind, hit.primitive_id)
        self._mesh_dirty = True
        self.update()

    def _pick(self, position: QPoint) -> None:
        try:
            hit = pick_scene(
                self.scene,
                self.camera,
                self.width(),
                self.height(),
                position.x(),
                position.y(),
            )
        except ValueError as exc:
            self.statusChanged.emit(str(exc))
            return
        self._set_selection(hit)

    def mousePressEvent(self, a0: QMouseEvent | None) -> None:
        if a0 is None:
            return
        event = a0
        self._last_mouse_position = event.position().toPoint()
        self._dragging = False
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
        self._last_mouse_position = current
        buttons = event.buttons()
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
        if event.button() == Qt.MouseButton.LeftButton and not self._dragging:
            self._pick(position)
        self._last_mouse_position = None
        self._dragging = False
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
