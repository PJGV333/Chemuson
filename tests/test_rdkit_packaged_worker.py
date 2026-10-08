"""RDKit worker protocol and fail-closed frozen-package smoke contracts."""

from __future__ import annotations

import json
import os
import subprocess
import sys
from pathlib import Path
from types import ModuleType
from typing import Any

import pytest

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "packaging" / "release"))

from chemuson.chemio import rdkit_packaged_smoke, rdkit_safe  # noqa: E402
import validate_packaged_rdkit_worker as packaged_validator  # noqa: E402


def _valid_frozen_report(executable: Path) -> dict[str, Any]:
    extensions = {
        name: {
            "path": str(executable.parent / "_MEI123" / f"{name.rsplit('.', 1)[-1]}.so"),
            "exists": True,
            "native": True,
            "inside_bundle": True,
        }
        for name in (
            "rdkit.Chem.rdchem",
            "rdkit.Chem.rdMolDescriptors",
            "rdkit.Chem.rdDistGeom",
        )
    }
    return {
        "ok": True,
        "executable": str(executable.resolve()),
        "frozen": True,
        "meipass": str(executable.parent / "_MEI123"),
        "expected_executable": str(executable.resolve()),
        "parent_rdkit_before": [],
        "parent_rdkit_after": [],
        "parent_gui_before": [],
        "parent_gui_after": [],
        "worker_diagnostics": {
            "ok": True,
            "python_executable": str(executable.resolve()),
            "worker_frozen": True,
            "worker_meipass": str(executable.parent / "_MEI123"),
            "native_extensions": extensions,
        },
        "descriptors": {
            "logp": -0.0014,
            "tpsa": 20.23,
            "hbd": 1,
            "hba": 1,
            "molecular_weight": 46.069,
        },
        "canonical_smiles": "CCO",
        "smiles_input": {"ok": True, "molblock_chars": 120},
        "conformer3d": {
            "atom_count": 3,
            "positions": {"1": [0.0, 0.0, 0.0], "2": [1.0, 0.0, 0.0], "3": [2.0, 0.0, 0.0]},
        },
    }


def test_smoke_allows_qt_runtime_hooks_without_a_gui_application(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    class FakeApplication:
        @staticmethod
        def instance() -> None:
            return None

    qt_core = ModuleType("PyQt6.QtCore")
    qt_gui = ModuleType("PyQt6.QtGui")
    qt_widgets = ModuleType("PyQt6.QtWidgets")
    qt_gui.QGuiApplication = FakeApplication
    qt_widgets.QApplication = FakeApplication
    monkeypatch.setitem(sys.modules, "PyQt6.QtCore", qt_core)
    monkeypatch.setitem(sys.modules, "PyQt6.QtGui", qt_gui)
    monkeypatch.setitem(sys.modules, "PyQt6.QtWidgets", qt_widgets)

    assert rdkit_packaged_smoke._loaded_gui_modules() == []


def test_smoke_detects_active_qt_application_and_chemuson_gui(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    class FakeApplication:
        @staticmethod
        def instance() -> object:
            return object()

    qt_gui = ModuleType("PyQt6.QtGui")
    qt_widgets = ModuleType("PyQt6.QtWidgets")
    qt_gui.QGuiApplication = type(
        "InactiveApplication", (), {"instance": staticmethod(lambda: None)}
    )
    qt_widgets.QApplication = FakeApplication
    monkeypatch.setitem(sys.modules, "PyQt6.QtGui", qt_gui)
    monkeypatch.setitem(sys.modules, "PyQt6.QtWidgets", qt_widgets)
    monkeypatch.setitem(sys.modules, "chemuson.gui", ModuleType("chemuson.gui"))

    assert rdkit_packaged_smoke._loaded_gui_modules() == [
        "PyQt6.QtWidgets.QApplication",
        "chemuson.gui",
    ]


def test_source_worker_imports_rdkit_and_native_extensions() -> None:
    report = rdkit_safe.rdkit_worker_diagnostics(timeout_s=8.0)

    assert report["ok"] is True, report
    assert report["worker_frozen"] is False
    assert report["rdkit_version"]
    for name in (
        "rdkit.Chem.rdchem",
        "rdkit.Chem.rdMolDescriptors",
        "rdkit.Chem.rdDistGeom",
    ):
        extension = report["native_extensions"][name]
        assert extension["exists"] is True
        assert extension["native"] is True


def test_isolated_smiles_and_3d_worker_smoke_is_bounded_and_keeps_parent_rdkit_free() -> None:
    code = """
import sys
from chemuson.core.model import MolGraph
from chemuson.chemio.rdkit_safe import (
    molgraph_to_smiles_isolated, conformer_3d_isolated, text_to_molblock
)
graph = MolGraph()
c1 = graph.add_atom("C", 0.0, 0.0)
c2 = graph.add_atom("C", 40.0, 0.0)
o = graph.add_atom("O", 80.0, 0.0)
graph.add_bond(c1.id, c2.id, order=1)
graph.add_bond(c2.id, o.id, order=1)
smiles, smiles_error = molgraph_to_smiles_isolated(graph, timeout_s=8.0)
coordinates, _metadata, conformer_error = conformer_3d_isolated(graph, timeout_s=15.0)
assert smiles_error is None, smiles_error
assert smiles == "CCO"
smiles_input = text_to_molblock("CCO", fmt="smiles", timeout_s=8.0)
assert smiles_input.get("ok") is True
assert smiles_input.get("molblock")
assert conformer_error is None, conformer_error
assert coordinates is not None and len(coordinates) == 3
assert all(len(position) == 3 for position in coordinates.values())
assert not any(name == "rdkit" or name.startswith("rdkit.") for name in sys.modules)
"""
    environment = os.environ.copy()
    environment["PYTHONPATH"] = os.pathsep.join(
        filter(None, (str(ROOT / "src"), environment.get("PYTHONPATH", "")))
    )
    result = subprocess.run(
        [sys.executable, "-c", code],
        cwd=ROOT,
        env=environment,
        capture_output=True,
        text=True,
        timeout=35,
        check=False,
    )
    assert result.returncode == 0, result.stderr or result.stdout


def test_frozen_worker_uses_private_json_files_and_cleans_temporary_directory(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    executable = tmp_path / "Chemuson.exe"
    executable.touch()
    observed: dict[str, Any] = {}

    def fake_run(command: list[str], **kwargs: Any) -> subprocess.CompletedProcess[str]:
        observed["command"] = command
        observed["temporary"] = Path(command[-1]).parent
        observed["request"] = json.loads(Path(command[-2]).read_text(encoding="utf-8"))
        observed["stdin"] = kwargs["stdin"]
        observed["capture_output"] = kwargs["capture_output"]
        observed["worker_flag"] = kwargs["env"]["CHEMUSON_INTERNAL_RDKIT_WORKER"]
        Path(command[-1]).write_text('{"ok": true, "descriptors": {"hbd": 1}}', encoding="utf-8")
        return subprocess.CompletedProcess(command, 0, "", "")

    monkeypatch.setattr(rdkit_safe.sys, "frozen", True, raising=False)
    monkeypatch.setattr(rdkit_safe.sys, "executable", str(executable))
    monkeypatch.setattr(rdkit_safe.subprocess, "run", fake_run)

    result = rdkit_safe._run_worker({"mode": "graph_descriptors"}, timeout_s=4.0)

    assert result["ok"] is True
    assert result["descriptors"] == {"hbd": 1}
    assert observed["request"] == {"mode": "graph_descriptors"}
    assert observed["command"][1] == "--chemuson-internal-rdkit-worker"
    assert observed["stdin"] == subprocess.DEVNULL
    assert observed["capture_output"] is True
    assert observed["worker_flag"] == "1"
    assert not observed["temporary"].exists()


def test_frozen_worker_timeout_removes_request_files_and_reports_timeout(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    executable = tmp_path / "Chemuson"
    observed: dict[str, Path] = {}

    def fake_run(command: list[str], **kwargs: Any) -> subprocess.CompletedProcess[str]:
        observed["temporary"] = Path(command[-1]).parent
        raise subprocess.TimeoutExpired(command, kwargs["timeout"])

    monkeypatch.setattr(rdkit_safe.sys, "frozen", True, raising=False)
    monkeypatch.setattr(rdkit_safe.sys, "executable", str(executable))
    monkeypatch.setattr(rdkit_safe.subprocess, "run", fake_run)

    result = rdkit_safe._run_worker({"mode": "diagnostics"}, timeout_s=1.0)

    assert result["error"] == "timeout"
    assert not observed["temporary"].exists()


def test_worker_response_parser_distinguishes_invalid_json_and_payload() -> None:
    base = {"python_executable": "Chemuson", "worker_path": "frozen:worker"}

    invalid_json = rdkit_safe._parse_worker_response("not-json", base)
    invalid_payload = rdkit_safe._parse_worker_response("[]", base)

    assert invalid_json["error"] == "invalid_worker_json"
    assert invalid_payload["error"] == "invalid_worker_payload"


def test_frozen_package_validator_accepts_only_complete_smoke_report(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    executable = tmp_path / "Chemuson"
    executable.touch()
    report = _valid_frozen_report(executable)

    def fake_run(command: list[str], **kwargs: Any) -> subprocess.CompletedProcess[str]:
        assert command[0] == str(executable.resolve())
        Path(command[3]).write_text(json.dumps(report), encoding="utf-8")
        return subprocess.CompletedProcess(command, 0, "", "")

    monkeypatch.setattr(packaged_validator.subprocess, "run", fake_run)
    assert packaged_validator.validate(executable)["ok"] is True

    report["worker_diagnostics"]["native_extensions"]["rdkit.Chem.rdchem"]["inside_bundle"] = False
    with pytest.raises(ValueError, match="outside the executable bundle"):
        packaged_validator._validate_report(report, executable)
