"""Source-side ChemName package and GUI-path regression coverage."""

from __future__ import annotations

import re
from contextlib import contextmanager
from pathlib import Path

from chemuson.chemname.packaged_smoke import (
    EXPECTED_MOLECULE_NAMES,
    EXPECTED_TEMPLATE_NAMES,
    run_smoke,
)
from tools.chemname_acceptance import DEFAULT_CASES_PATH, load_cases


def test_source_smoke_names_all_templates_and_representative_molecules() -> None:
    report = run_smoke(require_frozen=False)

    assert report["ok"] is True, report["failures"]
    assert len(report["template_resources"]) == 9
    assert {
        item["relative_path"]: item["name"] for item in report["template_resources"]
    } == EXPECTED_TEMPLATE_NAMES
    assert {item["id"]: item["name"] for item in report["molecule_results"]} == EXPECTED_MOLECULE_NAMES
    assert all(item["exists"] and item["size_bytes"] > 0 for item in report["template_resources"])


def test_acceptance_harness_includes_the_frozen_name_regressions() -> None:
    cases = {case["id"]: case for case in load_cases(DEFAULT_CASES_PATH)}
    smoke = run_smoke(require_frozen=False)
    expected_harness_names = {
        "ethanol": "smiles_ethanol",
        "acetamide": "smiles_acetamide",
        "ethane": "smiles_ethane",
        "cyclohexane": "smiles_cyclohexane",
    }

    assert all(case_id in cases for case_id in expected_harness_names.values())
    for smoke_id, harness_id in expected_harness_names.items():
        regex = re.compile(str(cases[harness_id]["expect"]), re.IGNORECASE)
        result = next(item for item in smoke["molecule_results"] if item["id"] == smoke_id)
        assert regex.search(result["name"]), (harness_id, result["name"])


def test_missing_template_is_reported_with_the_original_exception_cause(
    monkeypatch,
) -> None:
    import chemuson.chemname.special as special

    template_cache = dict(special._TEMPLATE_CACHE)
    molblock_cache = dict(special._TEMPLATE_MOLBLOCK_CACHE)

    @contextmanager
    def missing_template(*_parts: str, **_kwargs: str):
        raise FileNotFoundError("frozen template missing")
        yield Path("unused")

    special._TEMPLATE_CACHE.clear()
    special._TEMPLATE_MOLBLOCK_CACHE.clear()
    monkeypatch.setattr(special, "open_resource_path", missing_template)
    try:
        report = run_smoke(require_frozen=False)
    finally:
        special._TEMPLATE_CACHE.update(template_cache)
        special._TEMPLATE_MOLBLOCK_CACHE.update(molblock_cache)

    assert report["ok"] is False
    first_error = next(
        item for item in report["template_resources"] if item["status"] == "error"
    )
    assert first_error["exception_type"] == "ChemNameInternalError"
    assert first_error["cause_type"] == "FileNotFoundError"
    assert "frozen template missing" in first_error["cause"]


def test_pyinstaller_spec_requires_and_collects_all_nine_templates() -> None:
    root = Path(__file__).resolve().parents[1]
    spec = (root / "chemuson.spec").read_text(encoding="utf-8")

    assert 'EXPECTED_CHEMNAME_TEMPLATES = {' in spec
    assert 'CHEMNAME_TEMPLATE_DIR.glob("*/*.mol")' in spec
    assert 'datas_chemname_templates = [' in spec
    assert '"chemuson/chemname/templates/{path.parent.name}"' in spec
    assert 'datas = datas_c + datas_chemname_templates + datas_icons + datas_qt' in spec
    assert '"chemuson.chemname.packaged_smoke"' in spec


def test_flatpak_manifest_smokes_the_installed_package_data_from_app() -> None:
    root = Path(__file__).resolve().parents[1]
    manifest = (root / "packaging/flatpak/io.github.PJGV333.Chemuson.yml").read_text(
        encoding="utf-8"
    )
    pyproject = (root / "pyproject.toml").read_text(encoding="utf-8")

    assert '"**/*.mol"' in pyproject
    assert "pip3 install --prefix=/app --no-build-isolation --no-deps ." in manifest
    assert "packaged_smoke import run_smoke" in manifest
    assert "expected_package_root='/app'" in manifest
