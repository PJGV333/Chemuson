from __future__ import annotations

import argparse
from collections import Counter
import hashlib
import json
from pathlib import Path
import subprocess
import sys
import types

from PyQt6.QtCore import QEventLoop, QTimer, Qt
from PyQt6.QtGui import QOpenGLContext, QImage
from PyQt6.QtWidgets import QApplication, QDockWidget, QMainWindow

from chemuson.chemio.rdkit_io import molfile_to_molgraph, smiles_to_molgraph
from chemuson.geometry3d.model import OptimizationSettings, scene_molecule_from_graph
from chemuson.geometry3d.service import conformer_3d_for_graph
from chemuson.gui.viewport3d import renderer as viewport_renderer
from chemuson.gui.viewport3d.widget import MolecularViewportWidget

ROOT = Path(__file__).resolve().parents[5]
EVIDENCE = Path(__file__).resolve().parent
BASELINE_REVISION = "74ebc9a"
GL_DEPTH_WRITEMASK = 0x0B72
GL_DEPTH_FUNC = 0x0B74


def _git_module(revision: str, relative_path: str, module_name: str):
    source = subprocess.check_output(
        ["git", "show", f"{revision}:{relative_path}"],
        cwd=ROOT,
        text=True,
    )
    module = types.ModuleType(module_name)
    module.__package__ = "chemuson.gui.viewport3d"
    sys.modules[module_name] = module
    exec(compile(source, f"{revision}:{relative_path}", "exec"), module.__dict__)
    return module


def _load_baseline_widget():
    old_mesh = _git_module(
        BASELINE_REVISION,
        "src/chemuson/gui/viewport3d/mesh.py",
        "chemuson.gui.viewport3d._depth_baseline_mesh",
    )
    viewport_renderer.build_scene_mesh = old_mesh.build_scene_mesh
    old_widget = _git_module(
        BASELINE_REVISION,
        "src/chemuson/gui/viewport3d/widget.py",
        "chemuson.gui.viewport3d._depth_baseline_widget",
    )
    return old_widget.MolecularViewportWidget


def _molecular_scenes():
    specifications = (
        ("ethane", lambda: smiles_to_molgraph("CC"), 29, 80),
        ("benzene", lambda: smiles_to_molgraph("c1ccccc1"), 29, 80),
        (
            "tetrandrina",
            lambda: molfile_to_molgraph(
                (ROOT / "tests/data/tetrandrina_local_template_repaired.mol").read_text(
                    encoding="utf-8"
                )
            ),
            173,
            80,
        ),
    )
    records = []
    for name, graph_factory, seed, max_iters in specifications:
        graph = graph_factory()
        result = conformer_3d_for_graph(
            graph,
            force=True,
            settings=OptimizationSettings(
                max_iters=max_iters,
                steps_per_update=10,
                timeout_s=30.0,
                seed=seed,
            ),
        )
        if not result.ok or result.coordinate_set is None:
            raise RuntimeError(f"{name} conformer failed: {result.message}")
        scene = scene_molecule_from_graph(graph, result.coordinate_set)
        coordinate_payload = [
            (atom.id, atom.x, atom.y, atom.z)
            for atom in sorted(scene.atoms, key=lambda atom: str(atom.id))
        ]
        coordinate_sha256 = hashlib.sha256(
            json.dumps(
                coordinate_payload,
                separators=(",", ":"),
                allow_nan=False,
            ).encode("utf-8")
        ).hexdigest()
        records.append(
            {
                "name": name,
                "seed": seed,
                "coordinate_sha256": coordinate_sha256,
                "graph_atoms": len(graph.atoms),
                "graph_bonds": len(graph.bonds),
                "scene_atoms": len(scene.atoms),
                "scene_bonds": len(scene.bonds),
                "elements": dict(Counter(atom.element for atom in scene.atoms)),
                "method": result.method,
                "energy": result.energy,
                "scene": scene,
            }
        )
    return records


def _wait(app: QApplication, milliseconds: int = 220) -> None:
    loop = QEventLoop()
    QTimer.singleShot(milliseconds, loop.quit)
    loop.exec()
    app.processEvents()


def run(mode: str) -> list[dict[str, object]]:
    output = EVIDENCE / "captures" / mode
    output.mkdir(parents=True, exist_ok=True)
    widget_type = _load_baseline_widget() if mode == "before" else MolecularViewportWidget
    app = QApplication.instance() or QApplication([])
    if app.platformName().lower() in {"offscreen", "minimal"}:
        raise RuntimeError("The probe requires an authorized Qt/OpenGL desktop session")

    molecular_scenes = _molecular_scenes()
    all_records: list[dict[str, object]] = []
    for item in molecular_scenes:
        main = QMainWindow()
        main.resize(1100, 820)
        dock = QDockWidget(f"{item['name']} · {mode}", main)
        widget = widget_type(dock)
        widget.setMinimumSize(760, 600)
        dock.setWidget(widget)
        main.addDockWidget(Qt.DockWidgetArea.RightDockWidgetArea, dock)
        widget.set_scene(item["scene"])
        samples: list[dict[str, object]] = []
        original_render = widget._renderer.render

        def observe_render(*args) -> None:
            gl = widget._renderer._gl
            context = QOpenGLContext.currentContext()
            state = {
                "depth_test": bool(gl.glIsEnabled(viewport_renderer.GL_DEPTH_TEST)),
                "depth_write_mask": bool(gl.glGetBooleanv(GL_DEPTH_WRITEMASK)),
                "depth_func": int(gl.glGetIntegerv(GL_DEPTH_FUNC)),
                "depth_buffer_bits": (
                    context.format().depthBufferSize() if context is not None else 0
                ),
            }
            original_render(*args)
            state["gl_error_after_draw"] = int(gl.glGetError())
            samples.append(state)

        widget._renderer.render = observe_render
        main.show()
        dock.show()
        widget.show()
        try:
            _wait(app, 500)
            if not widget.isValid() or widget.context() is None:
                raise RuntimeError("Qt did not create a valid MolecularViewportWidget context")
            captures = []
            states = (
                ("docked", None),
                ("rotated", "orbit"),
                ("floating", "float"),
                ("redocked", "redock"),
            )
            for label, action in states:
                if action == "orbit":
                    widget.camera.orbit(32.0, -12.0, max(1, widget.height()))
                    widget._orbit_guide_active = True
                    widget.update()
                elif action == "float":
                    dock.setFloating(True)
                    dock.show()
                elif action == "redock":
                    dock.setFloating(False)
                    dock.show()
                _wait(app, 360 if action else 220)
                image = widget.grabFramebuffer().convertToFormat(QImage.Format.Format_RGB32)
                file_name = f"{item['name']}_{label}.png"
                if not image.save(str(output / file_name)):
                    raise RuntimeError(f"Could not save framebuffer {file_name}")
                context = widget.context()
                captures.append(
                    {
                        "state": label,
                        "image": file_name,
                        "image_size": [image.width(), image.height()],
                        "camera": {
                            "yaw": widget.camera.yaw,
                            "pitch": widget.camera.pitch,
                            "target": widget.camera.target,
                            "distance": widget.camera.distance,
                        },
                        "context_depth_bits": (
                            context.format().depthBufferSize() if context is not None else 0
                        ),
                        "render_samples": samples[-3:],
                    }
                )
            all_records.append(
                {
                    "mode": mode,
                    "baseline_revision": BASELINE_REVISION if mode == "before" else None,
                    "molecule": {
                        key: item[key]
                        for key in (
                            "name",
                            "seed",
                            "coordinate_sha256",
                            "graph_atoms",
                            "graph_bonds",
                            "scene_atoms",
                            "scene_bonds",
                            "elements",
                            "method",
                            "energy",
                        )
                    },
                    "captures": captures,
                }
            )
        finally:
            main.close()
            app.processEvents()
    (output / "results.json").write_text(
        json.dumps(all_records, ensure_ascii=False, indent=2, default=str),
        encoding="utf-8",
    )
    print(json.dumps(all_records, ensure_ascii=False, indent=2, default=str))
    return all_records


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("mode", choices=("before", "after"))
    args = parser.parse_args()
    run(args.mode)
