from __future__ import annotations

import pytest
from PyQt6.QtCore import QObject, QPoint, QPointF, QEvent, Qt, pyqtSignal
from PyQt6.QtGui import QMouseEvent, QResizeEvent
from PyQt6.QtWidgets import QApplication

from chemuson.geometry3d import SceneAtom3D, SceneMolecule3D
from chemuson.gui.main_window import ChemusonWindow
from chemuson.gui.viewport3d import MolecularViewportDock, MolecularViewportWidget
from chemuson.gui.viewport3d.picking import PickResult


@pytest.fixture(scope="module", autouse=True)
def _qapp():
    return QApplication.instance() or QApplication([])


def test_viewport_widget_keeps_selection_local() -> None:
    widget = MolecularViewportWidget()
    scene = SceneMolecule3D(
        atoms=(SceneAtom3D(1, "C", 0.0, 0.0, 0.0, 1.7, "#404040"),),
        bonds=(),
    )
    widget.set_scene(scene)
    assert widget.scene is scene
    assert widget.selected_atom_id is None
    assert widget.selected_bond_id is None
    assert not widget.gl_ready


def test_viewport_dock_is_hidden_until_explicitly_opened() -> None:
    dock = MolecularViewportDock()
    assert not dock.isVisible()
    assert dock.viewport.scene.atoms == ()
    dock.close()


def test_main_window_exposes_one_viewport_action_without_side_panel_duplication() -> None:
    window = ChemusonWindow()
    try:
        assert not window.viewport3d_dock.isVisible()
        assert window.action_toggle_viewport3d.text() == "Visor molecular 3D"
        assert window.action_toggle_viewport3d in window.view_menu.actions()
        assert window.side_panel.PAGE_KEYS == (
            "inspector",
            "validation",
            "properties",
            "templates",
            "appearance",
            "spectroscopy",
            "compchem",
        )
        window.action_toggle_viewport3d.trigger()
        assert not window.viewport3d_dock.isHidden()
        window.action_toggle_viewport3d.trigger()
        assert window.viewport3d_dock.isHidden()
    finally:
        window.close()


class _ContextSignals(QObject):
    aboutToBeDestroyed = pyqtSignal()


def test_scene_and_selection_defer_gpu_upload_until_paint(monkeypatch) -> None:
    widget = MolecularViewportWidget()
    calls = []
    widget._gl_ready = True
    widget._renderer._gl = type("GL", (), {"glClearColor": lambda *_args: None, "glClear": lambda *_args: None})()
    monkeypatch.setattr(widget._renderer, "upload", lambda *args: calls.append(args) or True)
    monkeypatch.setattr(widget._renderer, "render", lambda *_args: None)
    scene = SceneMolecule3D(atoms=(SceneAtom3D(1, "C", 0.0, 0.0, 0.0, 1.7, "#404040"),), bonds=())

    widget.set_scene(scene)
    widget._set_selection(PickResult("atom", 1, 1.0))

    assert calls == []
    assert widget._mesh_dirty
    widget.paintGL()
    assert calls == [(scene, 1, None)]
    assert not widget._mesh_dirty


def test_context_destruction_releases_renderer_and_allows_reinitialization(monkeypatch) -> None:
    widget = MolecularViewportWidget()
    context = _ContextSignals()
    initialized = []
    cleaned = []

    class RendererStub:
        _gl = None

        def initialize(self):
            initialized.append(True)
            return True, "ready"

        def upload(self, *_args):
            pass

        def cleanup(self):
            cleaned.append(True)

    widget._renderer = RendererStub()
    monkeypatch.setattr(widget, "context", lambda: context)
    monkeypatch.setattr(widget, "makeCurrent", lambda: None)
    monkeypatch.setattr(widget, "doneCurrent", lambda: None)

    widget.initializeGL()
    context.aboutToBeDestroyed.emit()

    assert initialized == [True]
    assert cleaned == [True]
    assert not widget.gl_ready
    assert widget._mesh_dirty
    widget.initializeGL()
    assert initialized == [True, True]


def test_failed_renderer_initialization_keeps_widget_unready(monkeypatch) -> None:
    widget = MolecularViewportWidget()
    messages = []

    class RendererStub:
        def initialize(self):
            return False, "shader unavailable"

    widget._renderer = RendererStub()
    widget.statusChanged.connect(messages.append)
    widget.initializeGL()

    assert not widget.gl_ready
    assert messages[-1] == "shader unavailable"


def test_auto_fit_uses_resize_until_manual_navigation(monkeypatch) -> None:
    widget = MolecularViewportWidget()
    widget.resize(120, 900)
    scene = SceneMolecule3D(
        atoms=(
            SceneAtom3D(1, "C", -2.0, 0.0, 0.0, 1.7, "#404040"),
            SceneAtom3D(2, "C", 2.0, 0.0, 0.0, 1.7, "#404040"),
        ),
        bonds=(),
    )
    widget.set_scene(scene)
    narrow_distance = widget.camera.distance
    widget.resize(900, 300)
    widget.resizeEvent(QResizeEvent(widget.size(), widget.size()))
    wide_distance = widget.camera.distance
    assert wide_distance < narrow_distance

    widget._last_mouse_position = QPoint(10, 10)
    event = QMouseEvent(
        QEvent.Type.MouseMove,
        QPointF(30, 10),
        Qt.MouseButton.NoButton,
        Qt.MouseButton.LeftButton,
        Qt.KeyboardModifier.NoModifier,
    )
    widget.mouseMoveEvent(event)
    manual_state = (widget.camera.target, widget.camera.distance, widget.camera.yaw, widget.camera.pitch)
    widget.resize(120, 900)
    widget.resizeEvent(QResizeEvent(widget.size(), widget.size()))
    assert (widget.camera.target, widget.camera.distance, widget.camera.yaw, widget.camera.pitch) == manual_state

    replacement_scene = SceneMolecule3D(
        atoms=(SceneAtom3D(3, "O", 100.0, 20.0, 0.0, 1.52, "#FF0D0D"),),
        bonds=(),
    )
    widget.set_scene(replacement_scene)
    widget.resize(900, 300)
    widget.resizeEvent(QResizeEvent(widget.size(), widget.size()))

    assert (widget.camera.target, widget.camera.distance, widget.camera.yaw, widget.camera.pitch) == manual_state
