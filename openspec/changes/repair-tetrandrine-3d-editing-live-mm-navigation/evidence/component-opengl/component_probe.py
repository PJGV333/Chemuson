from __future__ import annotations

import json
import math
import os
from pathlib import Path

from PyQt6.QtCore import QEventLoop, QTimer, Qt
from PyQt6.QtGui import QImage, QOpenGLContext
from PyQt6.QtWidgets import QApplication, QDockWidget, QMainWindow

from chemuson.chemio.rdkit_io import molfile_to_molgraph
from chemuson.geometry3d.model import OptimizationSettings, scene_molecule_from_graph
from chemuson.geometry3d.service import conformer_3d_for_graph
from chemuson.gui.viewport3d.mesh import build_scene_mesh
from chemuson.gui.viewport3d.widget import MolecularViewportWidget

ROOT = Path(os.environ["TMPDIR"]) / "chemuson-render-test"
OUT = ROOT / "baseline"
OUT.mkdir(parents=True, exist_ok=True)

graph = molfile_to_molgraph(
    (Path.cwd() / "tests/data/tetrandrina_local_template_repaired.mol").read_text(
        encoding="utf-8"
    )
)
result = conformer_3d_for_graph(
    graph,
    force=True,
    settings=OptimizationSettings(
        max_iters=80,
        steps_per_update=10,
        timeout_s=30.0,
        seed=173,
    ),
)
if not result.ok or result.coordinate_set is None:
    raise RuntimeError(f"conformer generation failed: {result.message}")
scene = scene_molecule_from_graph(graph, result.coordinate_set)
mesh = build_scene_mesh(scene)
app = QApplication([])
main = QMainWindow()
main.resize(1200, 900)
dock = QDockWidget("3D GL probe", main)
widget = MolecularViewportWidget(dock)
dock.setWidget(widget)
main.addDockWidget(Qt.DockWidgetArea.RightDockWidgetArea, dock)
widget.set_scene(scene)
main.show()


def capture(label: str) -> dict[str, object]:
    loop = QEventLoop()
    QTimer.singleShot(250, loop.quit)
    loop.exec()
    app.processEvents()
    image = widget.grabFramebuffer().convertToFormat(QImage.Format.Format_RGB32)
    path = OUT / f"tetrandrina_scene_{label}.png"
    image.save(str(path))
    width, height = image.width(), image.height()
    count = 0
    bounds = [width, height, -1, -1]
    for y in range(height):
        for x in range(width):
            if x < 110 and y >= height - 110:
                continue
            color = image.pixelColor(x, y)
            if sum(
                abs(value - expected)
                for value, expected in zip(
                    (color.red(), color.green(), color.blue()), (14, 17, 20)
                )
            ) > 24:
                count += 1
                bounds = [
                    min(bounds[0], x),
                    min(bounds[1], y),
                    max(bounds[2], x),
                    max(bounds[3], y),
                ]
    widget.makeCurrent()
    try:
        gl = widget._renderer._gl
        if gl is None:
            raise RuntimeError("OpenGL function table is unavailable")
        vbo = widget._renderer._vbo
        diagnostics = {
            "context_present": QOpenGLContext.currentContext() is not None,
            "context_matches_renderer": (
                QOpenGLContext.currentContext() == widget._renderer._context
            ),
            "gl_viewport": gl.glGetIntegerv(0x0BA2),
            "gl_error": int(gl.glGetError()),
            "gl_version": str(gl.glGetString(0x1F02)),
            "vbo_size_bytes": int(vbo.size()) if vbo is not None else None,
            "renderer_vertex_count": widget._renderer._vertex_count,
        }
    finally:
        widget.doneCurrent()
    result_metadata = getattr(result, "metadata", {})
    diagnostics_metadata = (
        result_metadata.get("diagnostics", {})
        if isinstance(result_metadata, dict)
        else {}
    )
    minimization = (
        diagnostics_metadata.get("minimization")
        if isinstance(diagnostics_metadata, dict)
        else None
    )
    record = {
        "label": label,
        "qpa": app.platformName(),
        "logical_size": [widget.width(), widget.height()],
        "dpr": widget.devicePixelRatioF(),
        "scene": {
            "source": scene.source,
            "atoms": len(scene.atoms),
            "bonds": len(scene.bonds),
        },
        "mesh": {
            "vertices": mesh.vertex_count,
            "float_values": len(mesh.vertices),
            "all_finite": all(math.isfinite(value) for value in mesh.vertices),
        },
        "conformer": {
            "method": result.method,
            "energy": result.energy,
            "minimization": minimization,
        },
        "gpu": diagnostics,
        "framebuffer": {
            "path": path.name,
            "size": [width, height],
            "non_background_pixels_excluding_axes_roi": count,
            "bounds": bounds if count else None,
        },
    }
    print(json.dumps(record, ensure_ascii=False, default=str))
    return record


widget.show()
main.show()
app.processEvents()
records = [capture("docked")]
dock.setFloating(True)
dock.show()
records.append(capture("floating"))
dock.setFloating(False)
dock.show()
records.append(capture("redocked"))
(OUT / "tetrandrina_scene_gl.json").write_text(
    json.dumps(records, ensure_ascii=False, indent=2, default=str), encoding="utf-8"
)
main.hide()
app.processEvents()
app.quit()
