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


def test_reference_orchestration_reuses_name2structure_without_clean2d_or_model_tools():
    root = Path(__file__).resolve().parents[2]
    modules = _modules()
    assert "M16" in modules["M10"]["current_dependencies"]
    controller_path = root / "src/chemuson/gui/controllers/molecular_assistant_controller.py"
    controller = controller_path.read_text(encoding="utf-8")
    provider_path = root / "src/chemuson/molecular_assistant/provider.py"
    provider = provider_path.read_text(encoding="utf-8")

    assert "extract_requested_molecule_name" in controller
    assert "resolve_name_to_structure" in controller
    assert "allow_network=self._allow_external_reference" in controller
    assert "smiles_to_molgraph_isolated" in controller
    assert "chemuson.clean2d" not in controller
    assert '"tools"' not in provider
    assert "browser" not in provider.casefold()
    assert "browsing" not in provider.casefold()


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
