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


def test_delete_key_rejects_derived_hydrogen_without_mutating_scene(_qapp) -> None:
    widget = MolecularViewportWidget()
    hydrogen = SceneAtom3D(
        "genH1",
        "H",
        1.0,
        0.0,
        0.0,
        1.2,
        "#FFFFFF",
        generated=True,
        parent_atom_id=1,
    )
    scene = SceneMolecule3D(
        atoms=(SceneAtom3D(1, "C", 0.0, 0.0, 0.0, 1.7, "#404040"), hydrogen),
        bonds=(SceneBond3D("genB1", 1, "genH1"),),
    )
    widget.set_scene(scene)
    widget.set_edit_mode(True)
    widget._set_selection(PickResult("atom", "genH1", 1.0))
    messages = []
    widget.statusChanged.connect(messages.append)

    widget.keyPressEvent(
        QKeyEvent(
            QEvent.Type.KeyPress,
            Qt.Key.Key_Delete,
            Qt.KeyboardModifier.NoModifier,
        )
    )

    assert widget.scene is scene
    assert widget.selected_atom_id == "genH1"
    assert messages and "derivado" in messages[-1].lower()


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


def test_main_window_theme_changes_reach_the_viewport_dock(_qapp):
    window = ChemusonWindow()
    try:
        window.current_theme = "dark"
        window._apply_theme()
        assert window.viewport3d_dock._theme_name == "dark"
        assert "#243249" in window.viewport3d_dock.styleSheet()

        window.current_theme = "light"
        window._apply_theme()
        assert window.viewport3d_dock._theme_name == "light"
        assert "#F1F5F9" in window.viewport3d_dock.styleSheet()
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
    assert calls == [(scene, 1, None, "ball_and_stick")]
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
    scene = SceneMolecule3D(
        atoms=(SceneAtom3D(1, "C", 0.0, 0.0, 0.0, 1.7, "#404040"),),
        bonds=(),
    )
    widget.set_scene(scene)
    monkeypatch.setattr(widget, "context", lambda: context)
    monkeypatch.setattr(widget, "makeCurrent", lambda: None)
    monkeypatch.setattr(widget, "doneCurrent", lambda: None)

    widget.initializeGL()
    context.aboutToBeDestroyed.emit()

    assert initialized == [True]
    assert cleaned == [True]
    assert not widget.gl_ready
    assert widget._mesh_dirty
    assert widget.scene is scene
    widget.initializeGL()
    assert initialized == [True, True]
    assert widget.scene is scene


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


def test_measure_tool_reports_distance_and_angle_from_direct_picks(monkeypatch, _qapp):
    widget = MolecularViewportWidget()
    widget.resize(800, 600)
    widget.set_scene(
        SceneMolecule3D(
            atoms=(
                SceneAtom3D(1, "C", 0.0, 0.0, 0.0, 1.7, "#404040"),
                SceneAtom3D(2, "C", 1.0, 0.0, 0.0, 1.7, "#404040"),
                SceneAtom3D(3, "C", 1.0, 1.0, 0.0, 1.7, "#404040"),
            ),
            bonds=(SceneBond3D(1, 1, 2), SceneBond3D(2, 2, 3)),
        )
    )
    widget.set_edit_mode(True)
    widget.set_edit_tool("measure")
    monkeypatch.setattr(
        "chemuson.gui.viewport3d.widget.pick_scene",
        lambda _scene, _camera, _width, _height, x, _y, _representation: PickResult(
            "atom", 1 if x < 200 else 2 if x < 400 else 3, 1.0
        ),
    )
    statuses = []
    widget.statusChanged.connect(statuses.append)

    for x in (100, 300, 500):
        position = QPointF(float(x), 300.0)
        widget.mousePressEvent(
            QMouseEvent(
                QEvent.Type.MouseButtonPress,
                position,
                Qt.MouseButton.LeftButton,
                Qt.MouseButton.LeftButton,
                Qt.KeyboardModifier.NoModifier,
            )
        )
        widget.mouseReleaseEvent(
            QMouseEvent(
                QEvent.Type.MouseButtonRelease,
                position,
                Qt.MouseButton.LeftButton,
                Qt.MouseButton.NoButton,
                Qt.KeyboardModifier.NoModifier,
            )
        )

    assert any("1.00 Å" in message for message in statuses)
    assert statuses[-1].endswith("90.0°")


def test_bond_tool_does_not_create_an_atom_on_empty_space(monkeypatch, _qapp):
    widget = MolecularViewportWidget()
    widget.resize(800, 600)
    widget.set_scene(
        SceneMolecule3D(
            atoms=(SceneAtom3D(1, "C", 0.0, 0.0, 0.0, 1.7, "#404040"),),
            bonds=(),
        )
    )
    widget.set_edit_mode(True)
    widget.set_edit_tool("bond")
    monkeypatch.setattr(
        "chemuson.gui.viewport3d.widget.pick_scene",
        lambda _scene, _camera, _width, _height, x, _y, _representation: (
            PickResult("atom", 1, 1.0) if x < 450 else None
        ),
    )
    point_requests = []
    bond_requests = []
    widget.buildPointRequested.connect(point_requests.append)
    widget.buildBondRequested.connect(bond_requests.append)

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
            QPointF(600.0, 300.0),
            Qt.MouseButton.NoButton,
            Qt.MouseButton.LeftButton,
            Qt.KeyboardModifier.NoModifier,
        )
    )
    widget.mouseReleaseEvent(
        QMouseEvent(
            QEvent.Type.MouseButtonRelease,
            QPointF(600.0, 300.0),
            Qt.MouseButton.LeftButton,
            Qt.MouseButton.NoButton,
            Qt.KeyboardModifier.NoModifier,
        )
    )

    assert point_requests == []
    assert bond_requests == []


def test_tool_keyboard_shortcuts_select_measure_mode(_qapp):
    widget = MolecularViewportWidget()
    widget.keyPressEvent(
        QKeyEvent(
            QEvent.Type.KeyPress,
            Qt.Key.Key_R,
            Qt.KeyboardModifier.NoModifier,
            "r",
        )
    )
    assert widget.edit_mode
    assert widget.edit_tool == "measure"


def test_viewport_dock_exposes_compact_tools_with_shortcut_help(_qapp):
    dock = MolecularViewportDock()
    try:
        assert set(dock.tool_buttons) == {"select", "build", "move", "bond", "measure"}
        assert all(button.toolTip() for button in dock.tool_buttons.values())
        assert "R" in dock.tool_buttons["measure"].toolTip()
        dock.tool_buttons["measure"].click()
        assert dock.viewport.edit_tool == "measure"
        assert dock.edit_mode_checkbox.isChecked()
    finally:
        dock.close()


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


def test_build_drag_from_atom_previews_and_requests_one_bonded_atom(monkeypatch, _qapp):
    widget = MolecularViewportWidget()
    widget.resize(800, 600)
    widget.set_scene(
        SceneMolecule3D(
            atoms=(SceneAtom3D(1, "C", 0.0, 0.0, 0.0, 1.7, "#404040"),),
            bonds=(),
        )
    )
    widget.set_edit_mode(True)
    widget.set_edit_tool("build")
    widget.set_build_element("N")
    monkeypatch.setattr(
        "chemuson.gui.viewport3d.widget.pick_scene",
        lambda _scene, _camera, _width, _height, x, _y, _representation: (
            PickResult("atom", 1, 1.0) if x < 450 else None
        ),
    )
    anchors = []
    requests = []
    previews = []
    widget.buildAnchorSelected.connect(anchors.append)
    widget.buildPointRequested.connect(requests.append)
    widget.buildPreviewChanged.connect(previews.append)
    original_yaw = widget.camera.yaw

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
            QPointF(520.0, 320.0),
            Qt.MouseButton.NoButton,
            Qt.MouseButton.LeftButton,
            Qt.KeyboardModifier.NoModifier,
        )
    )
    assert previews and previews[-1] is not None
    widget.mouseReleaseEvent(
        QMouseEvent(
            QEvent.Type.MouseButtonRelease,
            QPointF(520.0, 320.0),
            Qt.MouseButton.LeftButton,
            Qt.MouseButton.NoButton,
            Qt.KeyboardModifier.NoModifier,
        )
    )

    assert anchors == [1]
    assert len(requests) == 1
    assert requests[0][0] == "N"
    assert len(requests[0][1]) == 3
    assert widget.camera.yaw == original_yaw
    assert previews[-1] is None


def test_dragging_derived_hydrogen_previews_and_requests_parented_edit(monkeypatch, _qapp):
    widget = MolecularViewportWidget()
    widget.resize(800, 600)
    widget.set_scene(
        SceneMolecule3D(
            atoms=(
                SceneAtom3D(1, "C", 0.0, 0.0, 0.0, 1.7, "#404040"),
                SceneAtom3D(
                    "genH1", "H", 1.0, 0.0, 0.0, 1.2, "#FFFFFF",
                    generated=True, parent_atom_id=1,
                ),
            ),
            bonds=(SceneBond3D("genB1", 1, "genH1"),),
        )
    )
    widget.set_edit_mode(True)
    monkeypatch.setattr(
        "chemuson.gui.viewport3d.widget.pick_scene",
        lambda *_args: PickResult("atom", "genH1", 1.0),
    )
    previews = []
    commits = []
    widget.derivedHydrogenPreviewChanged.connect(lambda *args: previews.append(args))
    widget.derivedHydrogenEditRequested.connect(commits.append)

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
            QPointF(500.0, 320.0),
            Qt.MouseButton.NoButton,
            Qt.MouseButton.LeftButton,
            Qt.KeyboardModifier.NoModifier,
        )
    )
    assert previews
    widget.mouseReleaseEvent(
        QMouseEvent(
            QEvent.Type.MouseButtonRelease,
            QPointF(500.0, 320.0),
            Qt.MouseButton.LeftButton,
            Qt.MouseButton.NoButton,
            Qt.KeyboardModifier.NoModifier,
        )
    )

    assert len(commits) == 1
    assert commits[0][0:2] == ("genH1", 1)
    assert commits[0][2] != (1.0, 0.0, 0.0)


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


def test_dock_has_no_active_edit_tool_while_in_navigation_mode(_qapp):
    dock = MolecularViewportDock()
    try:
        assert not dock.viewport.edit_mode
        assert not dock.edit_mode_checkbox.isChecked()
        assert not any(button.isChecked() for button in dock.tool_buttons.values())

        dock.tool_buttons["measure"].click()

        assert dock.viewport.edit_mode
        assert dock.edit_mode_checkbox.isChecked()
        assert dock.tool_buttons["measure"].isChecked()
        dock.edit_mode_checkbox.setChecked(False)
        assert not dock.viewport.edit_mode
        assert not any(button.isChecked() for button in dock.tool_buttons.values())
    finally:
        dock.close()


def test_viewport_dock_theme_is_scoped_and_reapplies_when_floating(_qapp, monkeypatch):
    from PyQt6.QtGui import QPalette
    from PyQt6.QtWidgets import QMainWindow
    from chemuson.geometry3d import SceneAtom3D, SceneMolecule3D
    from chemuson.gui.viewport3d.picking import PickResult
    from chemuson.gui.theme.tokens import get_tokens

    window = QMainWindow()
    dock = MolecularViewportDock(window)
    window.addDockWidget(Qt.DockWidgetArea.RightDockWidgetArea, dock)
    scene = SceneMolecule3D(
        atoms=(SceneAtom3D(1, "C", 0.0, 0.0, 0.0, 1.7, "#909090"),),
        bonds=(),
    )
    dock.set_scene(scene)
    dock.viewport.set_edit_mode(True)
    dock.viewport.set_edit_tool("bond")
    dock.viewport._set_selection(PickResult("atom", 1, 1.0))
    dock.viewport.camera.yaw = 0.75
    applied_themes = []
    set_theme = dock.set_theme

    def record_theme(theme_name):
        applied_themes.append(theme_name)
        set_theme(theme_name)

    monkeypatch.setattr(dock, "set_theme", record_theme)

    try:
        for theme in ("light", "dark"):
            dock.set_theme(theme)
            applied_themes.clear()
            tokens = get_tokens(theme)
            stylesheet = dock.styleSheet()
            assert "QDockWidget#molecularViewportDock" in stylesheet
            assert "QPushButton:disabled" in stylesheet
            assert "QComboBox:disabled" in stylesheet
            assert "QCheckBox::indicator:checked" in stylesheet
            assert dock.palette().color(QPalette.ColorRole.Window).name().upper() == str(
                tokens["bg"]
            ).upper()
            assert dock.palette().color(QPalette.ColorRole.WindowText).name().upper() == str(
                tokens["text1"]
            ).upper()
            dock.setFloating(True)
            assert dock.isFloating()
            dock.setFloating(False)
            assert not dock.isFloating()
            assert len(applied_themes) >= 2
            assert set(applied_themes) == {theme}
            assert dock.styleSheet() == stylesheet
            assert dock.viewport.scene is scene
            assert dock.viewport.edit_mode
            assert dock.viewport.edit_tool == "bond"
            assert dock.viewport.selected_atom_id == 1
            assert dock.viewport.camera.yaw == pytest.approx(0.75)
    finally:
        window.close()
