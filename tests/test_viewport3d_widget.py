from __future__ import annotations

import math

import pytest
from PyQt6.QtCore import QObject, QPoint, QPointF, QEvent, Qt, pyqtSignal
from PyQt6.QtGui import QHideEvent, QKeyEvent, QMouseEvent, QResizeEvent, QShowEvent
from PyQt6.QtWidgets import QApplication

from chemuson.geometry3d import SceneAtom3D, SceneBond3D, SceneMolecule3D
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


def test_hide_and_show_cleanup_renderer_and_reinitialize_same_context(monkeypatch) -> None:
    widget = MolecularViewportWidget()
    context = _ContextSignals()
    initialized = []
    cleaned = []

    class RendererStub:
        _gl = None

        def initialize(self):
            initialized.append(True)
            return True, "ready"

        def cleanup(self):
            cleaned.append(True)

    widget._renderer = RendererStub()
    monkeypatch.setattr(widget, "context", lambda: context)
    monkeypatch.setattr(widget, "makeCurrent", lambda: None)
    monkeypatch.setattr(widget, "doneCurrent", lambda: None)

    widget.initializeGL()
    widget.hideEvent(QHideEvent())

    assert initialized == [True]
    assert cleaned == [True]
    assert not widget.gl_ready

    widget.showEvent(QShowEvent())

    assert initialized == [True, True]
    assert widget.gl_ready


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


def test_edit_mode_drag_moves_atom_and_emits_one_undoable_geometry_edit(monkeypatch) -> None:
    widget = MolecularViewportWidget()
    widget.resize(1000, 1000)
    widget.set_scene(
        SceneMolecule3D(
            atoms=(
                SceneAtom3D(1, "C", 0.0, 0.0, 0.0, 1.7, "#404040"),
                SceneAtom3D(2, "O", 2.0, 0.0, 0.0, 1.52, "#FF0D0D"),
            ),
            bonds=(),
        )
    )
    widget.set_edit_mode(True)
    monkeypatch.setattr(
        "chemuson.gui.viewport3d.widget.pick_scene",
        lambda *_args: PickResult("atom", 1, 1.0),
    )
    previews = []
    commits = []
    widget.geometryPreviewChanged.connect(lambda *args: previews.append(args))
    widget.geometryEditCommitted.connect(lambda *args: commits.append(args))

    press = QMouseEvent(
        QEvent.Type.MouseButtonPress,
        QPointF(500.0, 500.0),
        Qt.MouseButton.LeftButton,
        Qt.MouseButton.LeftButton,
        Qt.KeyboardModifier.NoModifier,
    )
    move = QMouseEvent(
        QEvent.Type.MouseMove,
        QPointF(600.0, 500.0),
        Qt.MouseButton.NoButton,
        Qt.MouseButton.LeftButton,
        Qt.KeyboardModifier.NoModifier,
    )
    release = QMouseEvent(
        QEvent.Type.MouseButtonRelease,
        QPointF(600.0, 500.0),
        Qt.MouseButton.LeftButton,
        Qt.MouseButton.NoButton,
        Qt.KeyboardModifier.NoModifier,
    )
    widget.mousePressEvent(press)
    widget.mouseMoveEvent(move)
    widget.mouseReleaseEvent(release)

    assert len(previews) == 1
    assert len(commits) == 1
    before, after, selected = commits[0]
    assert selected == (1,)
    assert before[1] == (0.0, 0.0, 0.0)
    assert after[1] != before[1]
    assert after[2] == before[2]


def test_control_drag_moves_only_the_connected_fragment(monkeypatch) -> None:
    widget = MolecularViewportWidget()
    widget.resize(1000, 1000)
    widget.set_scene(
        SceneMolecule3D(
            atoms=(
                SceneAtom3D(1, "C", 0.0, 0.0, 0.0, 1.7, "#404040"),
                SceneAtom3D(2, "O", 2.0, 0.0, 0.0, 1.52, "#FF0D0D"),
                SceneAtom3D(3, "N", 6.0, 0.0, 0.0, 1.55, "#3050F8"),
            ),
            bonds=(SceneBond3D(1, 1, 2),),
        )
    )
    widget.set_edit_mode(True)
    monkeypatch.setattr(
        "chemuson.gui.viewport3d.widget.pick_scene",
        lambda *_args: PickResult("atom", 1, 1.0),
    )
    commits = []
    widget.geometryEditCommitted.connect(lambda *args: commits.append(args))

    widget.mousePressEvent(
        QMouseEvent(
            QEvent.Type.MouseButtonPress,
            QPointF(500.0, 500.0),
            Qt.MouseButton.LeftButton,
            Qt.MouseButton.LeftButton,
            Qt.KeyboardModifier.ControlModifier,
        )
    )
    widget.mouseMoveEvent(
        QMouseEvent(
            QEvent.Type.MouseMove,
            QPointF(600.0, 500.0),
            Qt.MouseButton.NoButton,
            Qt.MouseButton.LeftButton,
            Qt.KeyboardModifier.ControlModifier,
        )
    )
    widget.mouseReleaseEvent(
        QMouseEvent(
            QEvent.Type.MouseButtonRelease,
            QPointF(600.0, 500.0),
            Qt.MouseButton.LeftButton,
            Qt.MouseButton.NoButton,
            Qt.KeyboardModifier.ControlModifier,
        )
    )

    assert len(commits) == 1
    before, after, selected = commits[0]
    assert selected == (1, 2)
    assert after[1] != before[1]
    assert after[2] != before[2]
    assert after[3] == before[3]


def test_viewport_dock_exposes_edit_build_live_mm_and_constraint_controls(_qapp):
    from chemuson.gui.viewport3d.dock import MolecularViewportDock

    dock = MolecularViewportDock()
    try:
        assert dock.advanced_controls_button.isCheckable()
        assert dock.advanced_controls_panel.isHidden()
        dock.advanced_controls_button.setChecked(True)
        assert not dock.advanced_controls_panel.isHidden()
        live_values = []
        constraints = []
        dock.liveMMToggled.connect(live_values.append)
        dock.geometryConstraintRequested.connect(lambda *values: constraints.append(values))
        dock.edit_mode_checkbox.setChecked(True)
        dock.tool_combo.setCurrentIndex(1)
        dock.element_combo.setCurrentText("O")
        dock.live_mm_checkbox.setChecked(True)
        dock.constraint_combo.setCurrentIndex(1)
        dock.atom_ids_edit.setText("1, 2, 3")
        dock.target_value.setValue(120.0)
        dock.apply_constraint_button.click()

        assert dock.viewport.edit_mode
        assert dock.viewport.edit_tool == "build"
        assert dock.viewport.build_element == "O"
        assert live_values == [True]
        assert constraints == [("bond_angle", (1, 2, 3), 120.0)]
    finally:
        dock.close()


def test_build_tool_emits_a_world_position_for_empty_space_click(monkeypatch, _qapp):
    widget = MolecularViewportWidget()
    widget.resize(800, 600)
    widget.set_edit_mode(True)
    widget.set_edit_tool("build")
    widget.set_build_element("N")
    monkeypatch.setattr("chemuson.gui.viewport3d.widget.pick_scene", lambda *_args: None)
    requests = []
    widget.buildPointRequested.connect(requests.append)

    widget.mousePressEvent(
        QMouseEvent(
            QEvent.Type.MouseButtonPress,
            QPointF(400.0, 300.0),
            Qt.MouseButton.LeftButton,
            Qt.MouseButton.LeftButton,
            Qt.KeyboardModifier.NoModifier,
        )
    )
    widget.mouseReleaseEvent(
        QMouseEvent(
            QEvent.Type.MouseButtonRelease,
            QPointF(400.0, 300.0),
            Qt.MouseButton.LeftButton,
            Qt.MouseButton.NoButton,
            Qt.KeyboardModifier.NoModifier,
        )
    )

    assert len(requests) == 1
    element, position = requests[0]
    assert element == "N"
    assert len(position) == 3
    assert all(math.isfinite(value) for value in position)


def test_escape_cancels_an_uncommitted_geometry_drag(monkeypatch, _qapp):
    widget = MolecularViewportWidget()
    widget.resize(800, 600)
    before = {1: (0.0, 0.0, 0.0)}
    widget.set_scene(
        SceneMolecule3D(
            atoms=(SceneAtom3D(1, "C", 0.0, 0.0, 0.0, 1.7, "#404040"),),
            bonds=(),
        )
    )
    widget.set_edit_mode(True)
    monkeypatch.setattr(
        "chemuson.gui.viewport3d.widget.pick_scene",
        lambda *_args: PickResult("atom", 1, 1.0),
    )
    commits = []
    cancelled = []
    widget.geometryEditCommitted.connect(lambda *args: commits.append(args))
    widget.geometryEditCancelled.connect(lambda: cancelled.append(True))
    widget.mousePressEvent(
        QMouseEvent(
            QEvent.Type.MouseButtonPress,
            QPointF(400.0, 300.0),
            Qt.MouseButton.LeftButton,
            Qt.MouseButton.LeftButton,
            Qt.KeyboardModifier.NoModifier,
        )
    )
    widget.mouseMoveEvent(
        QMouseEvent(
            QEvent.Type.MouseMove,
            QPointF(480.0, 300.0),
            Qt.MouseButton.NoButton,
            Qt.MouseButton.LeftButton,
            Qt.KeyboardModifier.NoModifier,
        )
    )

    widget.keyPressEvent(
        QKeyEvent(
            QEvent.Type.KeyPress,
            Qt.Key.Key_Escape,
            Qt.KeyboardModifier.NoModifier,
        )
    )

    assert widget.scene.atoms[0].x == pytest.approx(before[1][0])
    assert widget.scene.atoms[0].y == pytest.approx(before[1][1])
    assert widget._geometry_drag_state is None
    assert cancelled == [True]
    assert commits == []
