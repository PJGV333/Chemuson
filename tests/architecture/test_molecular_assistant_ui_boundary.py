from __future__ import annotations

from pathlib import Path

import yaml


ROOT = Path(__file__).resolve().parents[2]
CATALOG = ROOT / "architecture" / "modules.yml"


def _modules() -> dict[str, dict]:
    with CATALOG.open(encoding="utf-8") as stream:
        return {item["id"]: item for item in yaml.safe_load(stream)["modules"]}


def test_m10_controller_consumes_m23_without_reversing_the_boundary():
    modules = _modules()
    controllers = modules["M10"]
    assistant = modules["M23"]

    assert "M23" in controllers["current_dependencies"]
    assert "M23" in controllers["target_dependencies"]
    assert "M10" not in assistant["current_dependencies"]
    assert "M10" in assistant["forbidden_dependencies"]
    assert assistant["current_dependencies"] == ["M00", "M01"]
    assert "M23" not in modules["M02"]["current_dependencies"]
    assert "M23" not in modules["M02"]["target_dependencies"]


def test_only_the_gui_controller_imports_m23_for_the_phase2_flow():
    controller = (
        ROOT
        / "src"
        / "chemuson"
        / "gui"
        / "controllers"
        / "molecular_assistant_controller.py"
    ).read_text(encoding="utf-8")
    main_window = (ROOT / "src/chemuson/gui/main_window.py").read_text(encoding="utf-8")
    shell = (ROOT / "src/chemuson/gui/shell/assembly.py").read_text(encoding="utf-8")
    dialog = (
        ROOT
        / "src"
        / "chemuson"
        / "gui"
        / "dialogs"
        / "molecular_assistant_dialog.py"
    ).read_text(encoding="utf-8")

    assert "chemuson.molecular_assistant" in controller
    assert "MolecularAssistantController" in shell
    assert "chemuson.molecular_assistant" not in main_window
    assert "chemuson.molecular_assistant" not in dialog
