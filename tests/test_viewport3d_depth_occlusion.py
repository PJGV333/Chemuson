from __future__ import annotations

import math
import os
from pathlib import Path

import pytest
from PyQt6.QtCore import QEventLoop, QTimer, Qt
from PyQt6.QtGui import QOpenGLContext
from PyQt6.QtWidgets import QApplication, QDockWidget, QMainWindow

from chemuson.geometry3d import SceneAtom3D, SceneBond3D, SceneMolecule3D
from chemuson.gui.viewport3d import renderer as viewport_renderer
from chemuson.gui.viewport3d.geometry import ATOM_RADIUS_SCALE
from chemuson.gui.viewport3d.mesh import bond_lane_layout, build_scene_mesh
from chemuson.gui.viewport3d.widget import MolecularViewportWidget


class _RecordingGL:
    def __init__(self) -> None:
        self.depth_test = False
        self.depth_write_mask = False
        self.depth_func = 0x0204
        self.clear_depth = 0.25
        self.clear_states = []

    def glEnable(self, capability: int) -> None:
        if capability == viewport_renderer.GL_DEPTH_TEST:
            self.depth_test = True

    def glDepthMask(self, enabled: bool) -> None:
        self.depth_write_mask = bool(enabled)

    def glDepthFunc(self, function: int) -> None:
        self.depth_func = int(function)

    def glClearDepth(self, depth: float) -> None:
        self.clear_depth = float(depth)

    def glClearColor(self, *_args) -> None:
        pass

    def glClear(self, mask: int) -> None:
        self.clear_states.append(
            (mask, self.depth_test, self.depth_write_mask, self.depth_func, self.clear_depth)
        )


class _RecordingRenderer:
    def __init__(self, gl: _RecordingGL) -> None:
        self._gl = gl
        self.render_states = []

    def upload(self, *_args) -> bool:
        return True

    def render(self, *_args) -> None:
        self.render_states.append(
            (
                self._gl.depth_test,
                self._gl.depth_write_mask,
                self._gl.depth_func,
                self._gl.clear_depth,
            )
        )


def test_paintgl_restores_depth_state_after_each_overlay_frame(monkeypatch) -> None:
    gl = _RecordingGL()
    renderer = _RecordingRenderer(gl)
    widget = MolecularViewportWidget()
    monkeypatch.setattr(widget, "_renderer", renderer)
    widget._gl_ready = True
    widget.set_scene(
        SceneMolecule3D(
            atoms=(SceneAtom3D(1, "C", 0.0, 0.0, 0.0, 1.7, "#FF0000"),),
            bonds=(),
        )
    )

    def painter_leaves_hostile_state() -> None:
        gl.depth_test = False
        gl.depth_write_mask = False
        gl.depth_func = 0x0204
        gl.clear_depth = 0.25

    for overlay in (
        "_paint_geometry_handle",
        "_paint_live_mm_hud",
        "_paint_orientation_axes",
        "_paint_orbit_guide",
    ):
        monkeypatch.setattr(widget, overlay, painter_leaves_hostile_state)

    widget.paintGL()
    widget.paintGL()

    expected = (True, True, viewport_renderer.GL_LESS, 1.0)
    assert renderer.render_states == [expected, expected]
    assert len(gl.clear_states) == 2
    for mask, depth_test, write_mask, depth_func, clear_depth in gl.clear_states:
        assert mask & viewport_renderer.GL_DEPTH_BUFFER_BIT
        assert (depth_test, write_mask, depth_func, clear_depth) == expected


@pytest.mark.parametrize(
    ("order", "aromatic"),
    ((1, False), (2, False), (3, False), (1, True)),
)
def test_ball_and_stick_bond_lanes_start_and_end_at_atom_surfaces(order, aromatic) -> None:
    bond_length = 2.2
    scene = SceneMolecule3D(
        atoms=(
            SceneAtom3D(1, "C", 0.0, 0.0, 0.0, 1.7, "#404040"),
            SceneAtom3D(2, "N", bond_length, 0.0, 0.0, 1.55, "#404040"),
        ),
        bonds=(SceneBond3D(1, 1, 2, order, aromatic),),
    )
    mesh = build_scene_mesh(scene)
    lane_offsets, lane_radius = bond_lane_layout(order, aromatic)
    sphere_vertex_count = 2 * 20 * 10 * 6
    floats_per_vertex = 9
    floats_per_lane = 12 * 6 * floats_per_vertex
    start_radius = 1.7 * ATOM_RADIUS_SCALE
    end_radius = 1.55 * ATOM_RADIUS_SCALE

    for lane_index, offset in enumerate(lane_offsets):
        start = sphere_vertex_count * floats_per_vertex + lane_index * floats_per_lane
        lane = mesh.vertices[start : start + floats_per_lane]
        axial_positions = lane[0::floats_per_vertex]
        radial_extent = abs(offset) + lane_radius
        expected_start_trim = math.sqrt(max(0.0, start_radius**2 - radial_extent**2))
        expected_end_trim = math.sqrt(max(0.0, end_radius**2 - radial_extent**2))
        assert min(axial_positions) == pytest.approx(expected_start_trim)
        assert max(axial_positions) == pytest.approx(bond_length - expected_end_trim)

    assert len(mesh.vertices) == (
        sphere_vertex_count + len(lane_offsets) * 12 * 6
    ) * floats_per_vertex


def test_real_viewport_occludes_bonds_with_depth_after_qpainter_overlays(monkeypatch) -> None:
    app = QApplication.instance()
    if app is None or app.platformName().lower() in {"offscreen", "minimal"}:
        pytest.skip("requires an authorized Qt/OpenGL desktop session")

    main = QMainWindow()
    main.resize(900, 720)
    dock = QDockWidget("Depth-occlusion regression", main)
    widget = MolecularViewportWidget(dock)
    widget.setMinimumSize(640, 520)
    dock.setWidget(widget)
    main.addDockWidget(Qt.DockWidgetArea.RightDockWidgetArea, dock)
    widget.set_scene(
        SceneMolecule3D(
            atoms=(
                SceneAtom3D(1, "C", 0.0, 0.0, 0.0, 1.7, "#FF2020"),
                SceneAtom3D(2, "C", 0.0, 0.0, -3.0, 1.7, "#2040FF"),
                SceneAtom3D(3, "C", 4.0, 0.0, 0.0, 1.7, "#40D080"),
            ),
            bonds=(SceneBond3D(1, 1, 2),),
            source="known-depth-occlusion",
        )
    )
    widget.camera.target = (0.0, 0.0, -1.5)
    widget.camera.yaw = 0.0
    widget.camera.pitch = 0.0
    widget._live_mm_telemetry = {
        "forcefield": "MMFF94",
        "energy": -4.2,
        "delta_energy": 0.01,
        "step": 4,
        "max_steps": 10,
        "fixed_atom_count": 0,
        "constraints": {},
    }
    widget._orbit_guide_active = True
    widget.edit_tool = "geometry"
    widget.geometry_operation = "bond_length"
    widget.selected_atom_ids = [2, 3]
    widget.selected_atom_id = 3

    overlay_counts = {"geometry": 0, "live_mm": 0, "axes": 0, "orbit": 0}
    for method_name, count_key in (
        ("_paint_geometry_handle", "geometry"),
        ("_paint_live_mm_hud", "live_mm"),
        ("_paint_orientation_axes", "axes"),
        ("_paint_orbit_guide", "orbit"),
    ):
        original_overlay = getattr(widget, method_name)

        def counted(original=original_overlay, key=count_key):
            overlay_counts[key] += 1
            return original()

        monkeypatch.setattr(widget, method_name, counted)

    depth_states = []
    original_render = widget._renderer.render

    def observe_render(*args) -> None:
        gl = widget._renderer._gl
        context = QOpenGLContext.currentContext()
        depth_states.append(
            (
                bool(gl.glIsEnabled(viewport_renderer.GL_DEPTH_TEST)),
                bool(gl.glGetBooleanv(0x0B72)),
                int(gl.glGetIntegerv(0x0B74)),
                context.format().depthBufferSize() if context is not None else 0,
            )
        )
        original_render(*args)

    monkeypatch.setattr(widget._renderer, "render", observe_render)
    main.show()
    dock.show()
    widget.show()

    def settle(milliseconds: int = 100) -> None:
        loop = QEventLoop()
        QTimer.singleShot(milliseconds, loop.quit)
        loop.exec()
        app.processEvents()

    evidence_value = os.environ.get("CHEMUSON_DEPTH_EVIDENCE_DIR")
    evidence_dir = Path(evidence_value) if evidence_value else None
    if evidence_dir is not None:
        evidence_dir.mkdir(parents=True, exist_ok=True)

    def assert_front_atom_is_opaque(label: str) -> tuple[int, int, int]:
        image = widget.grabFramebuffer()
        position = widget.camera.world_to_screen(
            (0.0, 0.0, 0.0), widget.width(), widget.height()
        )
        assert position is not None
        ratio = image.devicePixelRatio()
        color = image.pixelColor(round(position[0] * ratio), round(position[1] * ratio))
        pixel = (color.red(), color.green(), color.blue())
        assert pixel[0] > 1.5 * max(pixel[1:])
        if evidence_dir is not None:
            assert image.save(str(evidence_dir / f"{label}.png"))
        return pixel

    try:
        settle(500)
        if not widget.isValid() or widget.context() is None:
            pytest.skip("Qt did not create an OpenGL viewport context")
        assert widget.context().format().depthBufferSize() >= 16

        for _ in range(4):
            widget.update()
            settle()
        assert len(depth_states) >= 4
        assert all(
            enabled and writable and depth_func == viewport_renderer.GL_LESS and depth_bits >= 16
            for enabled, writable, depth_func, depth_bits in depth_states[-4:]
        )
        assert all(count >= 4 for count in overlay_counts.values())
        ball_and_stick_pixel = assert_front_atom_is_opaque("ball_and_stick")

        widget.representation = "space_filling"
        widget._mesh_dirty = True
        widget.update()
        settle(180)
        space_filling_pixel = assert_front_atom_is_opaque("space_filling")

        widget.camera.orbit(18.0, -5.0, widget.height())
        widget.update()
        settle(180)
        rotated_pixel = assert_front_atom_is_opaque("space_filling_rotated")

        dock.setFloating(True)
        dock.show()
        settle(220)
        dock.setFloating(False)
        dock.show()
        settle(220)
        redocked_pixel = assert_front_atom_is_opaque("space_filling_redocked")
        assert all(
            enabled and writable and depth_func == viewport_renderer.GL_LESS and depth_bits >= 16
            for enabled, writable, depth_func, depth_bits in depth_states[-4:]
        )
        assert all(count > 4 for count in overlay_counts.values())
        assert all(
            pixel[0] > 1.5 * max(pixel[1:])
            for pixel in (
                ball_and_stick_pixel,
                space_filling_pixel,
                rotated_pixel,
                redocked_pixel,
            )
        )
    finally:
        main.close()
        app.processEvents()
