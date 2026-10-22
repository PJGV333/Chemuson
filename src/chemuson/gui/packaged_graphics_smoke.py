"""Real-window/OpenGL smoke used by the frozen Linux AppImage validator."""

from __future__ import annotations

import time
from pathlib import Path

from PyQt6.QtCore import QCoreApplication, QEvent, QPointF, QRect, QSettings, Qt, QTimer
from PyQt6.QtGui import QSurfaceFormat
from PyQt6.QtTest import QTest
from PyQt6.QtWidgets import QApplication, QMenu
from rdkit import Chem
from rdkit.Chem import AllChem

from chemuson.core.model import MolGraph
from chemuson.geometry3d import CoordinateSet3D, SceneMolecule3D, scene_molecule_from_graph
from chemuson.gui.main_window import ChemusonWindow
from chemuson.gui.viewport3d.mesh import build_scene_mesh


_SIMPLE_SMILES = "CCO"
_COMPLEX_SMILES = "Cn1c(=O)c2c(ncn2C)n(C)c1=O"


def _wait_until(
    app: QApplication,
    predicate,
    *,
    description: str,
    timeout_s: float = 12.0,
) -> None:
    deadline = time.monotonic() + timeout_s
    while time.monotonic() < deadline:
        app.processEvents()
        if predicate():
            return
        time.sleep(0.025)
    raise RuntimeError(f"Timed out waiting for {description}.")


def _settle(app: QApplication, duration_s: float = 0.12) -> None:
    deadline = time.monotonic() + duration_s
    while time.monotonic() < deadline:
        app.processEvents()
        time.sleep(0.01)


def _changed_pixel_samples(before, after, *, region: QRect | None = None) -> int:
    if before.isNull() or after.isNull():
        return 0
    width = min(before.width(), after.width())
    height = min(before.height(), after.height())
    if region is None:
        region = QRect(0, 0, width, height)
    region = region.intersected(QRect(0, 0, width, height))
    if region.isEmpty():
        return 0
    step_x = max(1, region.width() // 320)
    step_y = max(1, region.height() // 200)
    changed = 0
    for y in range(region.top(), region.bottom() + 1, step_y):
        for x in range(region.left(), region.right() + 1, step_x):
            first = before.pixelColor(x, y)
            second = after.pixelColor(x, y)
            if max(
                abs(first.red() - second.red()),
                abs(first.green() - second.green()),
                abs(first.blue() - second.blue()),
            ) >= 24:
                changed += 1
    return changed


def _rdkit_scene(smiles: str, *, source: str) -> SceneMolecule3D:
    molecule = Chem.MolFromSmiles(smiles)
    if molecule is None:
        raise RuntimeError(f"RDKit could not parse the graphics fixture {source!r}.")
    molecule = Chem.AddHs(molecule)
    if AllChem.EmbedMolecule(molecule, randomSeed=0xC0FFEE, maxAttempts=100) != 0:
        raise RuntimeError(f"RDKit could not embed the graphics fixture {source!r}.")

    graph = MolGraph()
    rdkit_to_graph_id: dict[int, int] = {}
    for atom in molecule.GetAtoms():
        graph_atom = graph.add_atom(
            atom.GetSymbol(),
            0.0,
            0.0,
            is_explicit=True,
            no_implicit=True,
        )
        rdkit_to_graph_id[atom.GetIdx()] = graph_atom.id

    for bond in molecule.GetBonds():
        graph.add_bond(
            rdkit_to_graph_id[bond.GetBeginAtomIdx()],
            rdkit_to_graph_id[bond.GetEndAtomIdx()],
            order=max(1, round(bond.GetBondTypeAsDouble())),
            is_aromatic=bond.GetIsAromatic(),
        )

    conformer = molecule.GetConformer()
    positions = {
        graph_id: (
            float(conformer.GetAtomPosition(rdkit_index).x),
            float(conformer.GetAtomPosition(rdkit_index).y),
            float(conformer.GetAtomPosition(rdkit_index).z),
        )
        for rdkit_index, graph_id in rdkit_to_graph_id.items()
    }
    return scene_molecule_from_graph(
        graph,
        CoordinateSet3D(
            positions=positions,
            source=source,
            method="RDKit ETKDG graphics smoke",
        ),
    )


def _populate_2d_canvas(window: ChemusonWindow, app: QApplication) -> dict[str, int]:
    """Draw a small C-C-O fixture through the editor's real bond tool."""
    canvas = window.canvas
    center = canvas.mapToScene(canvas.viewport().rect().center())
    canvas.set_current_tool("tool_bond")
    canvas.state.default_element = "C"
    start = canvas.mapFromScene(QPointF(center.x() - 35.0, center.y()))
    end = canvas.mapFromScene(QPointF(center.x() + 5.0, center.y()))
    QTest.mousePress(canvas.viewport(), Qt.MouseButton.LeftButton, pos=start)
    QTest.mouseMove(canvas.viewport(), end, 10)
    QTest.mouseRelease(canvas.viewport(), Qt.MouseButton.LeftButton, pos=end)
    app.processEvents()

    first_id = min(canvas.model.atoms)
    first = canvas.model.get_atom(first_id)
    canvas.state.default_element = "O"
    start = canvas.mapFromScene(QPointF(first.x, first.y))
    end = canvas.mapFromScene(QPointF(first.x - 40.0, first.y + 20.0))
    QTest.mousePress(canvas.viewport(), Qt.MouseButton.LeftButton, pos=start)
    QTest.mouseMove(canvas.viewport(), end, 10)
    QTest.mouseRelease(canvas.viewport(), Qt.MouseButton.LeftButton, pos=end)
    app.processEvents()
    canvas.undo_stack.setClean()
    canvas.viewport().update()
    return {"atoms": len(canvas.model.atoms), "bonds": len(canvas.model.bonds)}


def _canvas_fixture_region(window: ChemusonWindow) -> QRect:
    canvas = window.canvas
    points = [
        canvas.mapFromScene(QPointF(atom.x, atom.y))
        for atom in canvas.model.atoms.values()
    ]
    if not points:
        return QRect()
    left = min(point.x() for point in points) - 35
    top = min(point.y() for point in points) - 30
    right = max(point.x() for point in points) + 35
    bottom = max(point.y() for point in points) + 30
    ratio = float(canvas.viewport().devicePixelRatioF())
    pixel_region = QRect(
        round(left * ratio),
        round(top * ratio),
        round((right - left + 1) * ratio),
        round((bottom - top + 1) * ratio),
    )
    image_rect = QRect(
        0,
        0,
        round(canvas.viewport().width() * ratio),
        round(canvas.viewport().height() * ratio),
    )
    return pixel_region.intersected(image_rect)


def _gl_context_report(viewport) -> dict[str, object]:
    if not viewport.isValid() or not viewport.gl_ready:
        raise RuntimeError("The 3D viewport did not obtain a valid OpenGL context.")
    context = viewport.context()
    if context is None or not context.isValid():
        raise RuntimeError("The 3D viewport's QOpenGLContext is invalid.")
    version = tuple(int(part) for part in context.format().version())
    if version < (4, 1):
        raise RuntimeError(f"OpenGL 4.1 Core is required; got {version!r}.")
    if context.format().profile() != QSurfaceFormat.OpenGLContextProfile.CoreProfile:
        raise RuntimeError("The 3D viewport did not receive a CoreProfile context.")

    viewport.makeCurrent()
    try:
        gl = viewport._renderer._gl
        if gl is None:
            raise RuntimeError("The molecular renderer did not initialize OpenGL functions.")
        renderer = str(gl.glGetString(0x1F01) or "")
        gl_version = str(gl.glGetString(0x1F02) or "")
        if not renderer or not gl_version:
            raise RuntimeError("OpenGL returned empty renderer/version strings.")
    finally:
        viewport.doneCurrent()
    return {
        "valid": True,
        "version": list(version),
        "profile": "core",
        "renderer": renderer,
        "gl_version": gl_version,
    }


def _close_window_menus(window: ChemusonWindow, app: QApplication) -> None:
    menus = set(window.findChildren(QMenu))
    for menu_name in ("view_menu", "_menu_aggregator"):
        menu = getattr(window, menu_name, None)
        if isinstance(menu, QMenu):
            menus.add(menu)
    for menu in menus:
        try:
            menu.close()
        except RuntimeError:
            pass
    window.menuBar().close()
    app.processEvents()
    for menu in menus:
        try:
            QCoreApplication.sendPostedEvents(menu, QEvent.Type.DeferredDelete)
        except RuntimeError:
            pass
    app.processEvents()


def _disable_background_geometry_generation(window: ChemusonWindow) -> None:
    """Keep the 2D fixture's delayed conformer job from replacing smoke scenes."""
    timer = getattr(window, "_geometry3d_sync_timer", None)
    if timer is not None:
        timer.stop()
    # This smoke injects its own RDKit scenes after drawing in 2D. The real
    # window's debounce timer would otherwise start a document conformer job and
    # replace those fixtures asynchronously while the OpenGL assertions run.
    setattr(window, "_schedule_geometry3d_sync", _no_op_geometry3d_sync)


def _no_op_geometry3d_sync() -> None:
    return None


def run_packaged_graphics_smoke(*, screenshot_path: Path) -> dict[str, object]:
    """Show the real editor, draw 2D content, and render/close the 3D dock."""
    app = QApplication.instance() or QApplication([])
    app.setApplicationName("Chemuson")
    app.setApplicationDisplayName("ChemUSON")
    app.setQuitOnLastWindowClosed(False)
    settings = QSettings("Chemuson", "Chemuson")
    settings.setValue("ui/onboarding/completed", True)
    settings.sync()

    window = ChemusonWindow()
    _disable_background_geometry_generation(window)
    dock = window.viewport3d_dock
    viewport = dock.viewport
    window.resize(1280, 860)
    window.show()
    outcome: dict[str, object] = {}
    errors: list[Exception] = []

    def perform_smoke() -> None:
        try:
            _wait_until(
                app,
                lambda: window.isVisible() and window.canvas.isVisible(),
                description="the main window and 2D canvas",
            )
            dock.show()
            _wait_until(
                app,
                lambda: viewport.isVisible() and viewport.gl_ready,
                description="the visible 3D dock and initialized OpenGL renderer",
            )
            context_report = _gl_context_report(viewport)

            canvas_before = window.canvas.viewport().grab().toImage()
            canvas_counts = _populate_2d_canvas(window, app)
            _settle(app)
            canvas_after = window.canvas.viewport().grab().toImage()
            canvas_region = _canvas_fixture_region(window)
            canvas_changed = _changed_pixel_samples(
                canvas_before,
                canvas_after,
                region=canvas_region,
            )
            if canvas_counts != {"atoms": 3, "bonds": 2} or canvas_changed < 8:
                debug_dir = screenshot_path.resolve().parent
                debug_dir.mkdir(parents=True, exist_ok=True)
                canvas_before.save(str(debug_dir / "graphics-smoke-canvas-before.png"))
                canvas_after.save(str(debug_dir / "graphics-smoke-canvas-after.png"))
                raise RuntimeError(
                    "The 2D canvas did not visibly render its molecular fixture "
                    f"(counts={canvas_counts}, changed_samples={canvas_changed}, "
                    f"region={canvas_region.getRect()}, "
                    f"sizes={canvas_before.width()}x{canvas_before.height()}/"
                    f"{canvas_after.width()}x{canvas_after.height()})."
                )

            empty_frame = viewport.grabFramebuffer()
            if empty_frame.isNull():
                raise RuntimeError("The OpenGL viewport returned an empty framebuffer.")

            rendered_scenes: list[dict[str, object]] = []
            for label, smiles in (("simple", _SIMPLE_SMILES), ("complex", _COMPLEX_SMILES)):
                scene = _rdkit_scene(smiles, source=f"appimage-graphics-{label}")
                dock.set_scene(scene)
                _wait_until(
                    app,
                    lambda: (
                        viewport.scene_diagnostics.get("upload_status") == "uploaded"
                        and viewport.scene_diagnostics.get("render_status") == "requested"
                    ),
                    description=f"OpenGL upload/render of the {label} molecule",
                )
                vertices_before_capture = viewport._renderer._vertex_count
                framebuffer = viewport.grabFramebuffer()
                vertices_after_capture = viewport._renderer._vertex_count
                changed = _changed_pixel_samples(empty_frame, framebuffer)
                if framebuffer.isNull() or vertices_after_capture <= 0 or changed < 8:
                    scene_with_preview = viewport._scene_with_build_preview()
                    try:
                        cpu_mesh_vertices = build_scene_mesh(
                            scene_with_preview,
                            viewport.selected_atom_id,
                            viewport.selected_bond_id,
                            viewport.representation,
                            viewport.selected_atom_ids,
                        ).vertex_count
                        cpu_mesh_error = None
                    except Exception as exc:
                        cpu_mesh_vertices = None
                        cpu_mesh_error = f"{type(exc).__name__}: {exc}"
                    context = viewport.context()
                    context_valid = context.isValid() if context is not None else False
                    atom_sample = [
                        (atom.id, atom.element, atom.radius, atom.x, atom.y, atom.z)
                        for atom in scene.atoms[:3]
                    ]
                    raise RuntimeError(
                        f"The {label} molecule did not produce visible OpenGL geometry "
                        f"(changed_samples={changed}, vertices={vertices_after_capture}, "
                        f"vertices_before_capture={vertices_before_capture}, "
                        f"cpu_mesh_vertices={cpu_mesh_vertices}, cpu_mesh_error={cpu_mesh_error!r}, "
                        f"scene_atoms={len(scene.atoms)}, scene_bonds={len(scene.bonds)}, "
                        f"atom_sample={atom_sample!r}, diagnostics={viewport.scene_diagnostics!r}, "
                        f"gl_ready={viewport.gl_ready}, renderer_initialized="
                        f"{viewport._renderer.initialized}, mesh_dirty={viewport._mesh_dirty}, "
                        f"context_valid={context_valid}, context_report={context_report!r})."
                    )
                rendered_scenes.append(
                    {
                        "name": label,
                        "atoms": len(scene.atoms),
                        "bonds": len(scene.bonds),
                        "vertices": viewport._renderer._vertex_count,
                        "changed_samples": changed,
                        "diagnostics": viewport.scene_diagnostics,
                    }
                )

            window_image = window.grab().toImage()
            resolved_screenshot = screenshot_path.resolve()
            resolved_screenshot.parent.mkdir(parents=True, exist_ok=True)
            if window_image.isNull() or not window_image.save(str(resolved_screenshot)):
                raise RuntimeError("Could not capture the visible AppImage window.")
            outcome["report"] = {
                "ok": True,
                "qpa": "xcb",
                "window": {
                    "visible": True,
                    "size": [window_image.width(), window_image.height()],
                    "screenshot": str(resolved_screenshot),
                },
                "canvas_2d": {
                    "visible": True,
                    **canvas_counts,
                    "changed_samples": canvas_changed,
                },
                "context": context_report,
                "rendered_scenes": rendered_scenes,
                "screenshot": str(resolved_screenshot),
            }
        except Exception as exc:
            errors.append(exc)
        finally:
            try:
                dock.hide()
                _wait_until(
                    app,
                    lambda: not viewport.isVisible() and not viewport.gl_ready,
                    description="clean OpenGL resource release when closing the 3D dock",
                    timeout_s=8.0,
                )
                _close_window_menus(window, app)
                # The fixture is transient; auto-discard it while still exercising
                # Chemuson's normal worker/Qt teardown path.
                window._confirm_discard_changes = lambda _canvas: True
                window.close()
                _wait_until(
                    app,
                    lambda: not window.isVisible(),
                    description="clean main-window shutdown",
                    timeout_s=8.0,
                )
                if not window._async_shutdown_complete:
                    raise RuntimeError("Background GUI workers did not finish before teardown.")
                window.deleteLater()
                QCoreApplication.sendPostedEvents(window, QEvent.Type.DeferredDelete)
                app.processEvents()
            except Exception as exc:
                errors.append(exc)
            QTimer.singleShot(0, app.quit)

    QTimer.singleShot(0, perform_smoke)
    app.exec()
    if errors:
        raise RuntimeError(f"Packaged graphics smoke failed: {errors[0]}") from errors[0]
    report = outcome.get("report")
    if not isinstance(report, dict):
        raise RuntimeError("Packaged graphics smoke ended without a report.")
    report["dock_closed_cleanly"] = True
    report["window_closed_cleanly"] = True
    return report
