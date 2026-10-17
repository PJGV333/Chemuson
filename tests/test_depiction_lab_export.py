from __future__ import annotations

import hashlib
import json
import os
from pathlib import Path
import subprocess
import sys

from PyQt6.QtGui import QImage

from tools.visual.depiction_lab_app.exporting import export_case_bundle
from tools.visual.depiction_lab_app.parameters import DepictionStyleParameters
from tools.visual.depiction_lab_app.profiles import encode_profile, load_profile
from tools.visual.depiction_lab_app.scenes import CASES

ROOT = Path(__file__).resolve().parents[1]


def _sha(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def _isolated_env(config: Path) -> dict[str, str]:
    env = os.environ.copy()
    env["QT_QPA_PLATFORM"] = "offscreen"
    env["XDG_CONFIG_HOME"] = str(config)
    python_paths = [str(ROOT), str(ROOT / "src")]
    if env.get("PYTHONPATH"):
        python_paths.append(env["PYTHONPATH"])
    env["PYTHONPATH"] = os.pathsep.join(python_paths)
    return env


def test_case_export_produces_matched_png_svg_profile_and_provenance(tmp_path) -> None:
    current = DepictionStyleParameters.current()
    values = current.to_mapping()
    values["solid_wedge_width_ratio"] = 0.34
    working = DepictionStyleParameters.from_mapping(values)
    result = export_case_bundle(
        tmp_path / "case-evidence",
        "A1",
        current,
        working,
        theme="dark",
        profile_name="Visual candidate",
    )
    output = tmp_path / "case-evidence"
    current_png = output / "current-A1-dark.png"
    working_png = output / "working-A1-dark.png"
    comparison_png = output / "comparison-A1-dark.png"
    assert {Path(path).name for path in result["artifacts"]} >= {
        current_png.name,
        working_png.name,
        comparison_png.name,
        "working-profile.json",
        "report.md",
    }
    assert QImage(str(current_png)).size().width() == 1200
    assert QImage(str(working_png)).size().height() == 820
    assert QImage(str(comparison_png)).size().width() == 2400
    assert _sha(current_png) != _sha(working_png)
    name, restored = load_profile(output / "working-profile.json")
    assert name == "Visual candidate"
    assert restored == working
    report = (output / "report.md").read_text(encoding="utf-8")
    assert "Renderer source SHA at export" in report
    assert "Git working tree at export" in report
    assert "}\n```\n\n## Changed values" in report
    assert "A1" in report and "solid_wedge_width_ratio" in report
    if "exported through" in result["svg"]:
        for svg in (output / "current-A1-dark.svg", output / "working-A1-dark.svg"):
            assert svg.is_file()
            assert all(
                line == line.rstrip(b" \\t")
                for line in svg.read_bytes().splitlines()
            )


def test_full_contact_sheet_export_is_bounded_and_reproducible(tmp_path) -> None:
    profile_path = tmp_path / "candidate.json"
    parameters = DepictionStyleParameters.current().to_mapping()
    parameters["solid_wedge_width_ratio"] = 0.32
    profile_path.write_text(
        encode_profile(DepictionStyleParameters.from_mapping(parameters), name="Candidate"),
        encoding="utf-8",
    )
    output = tmp_path / "gallery"
    command = [
        sys.executable,
        str(ROOT / "tools/visual/depiction_lab.py"),
        "--export-gallery",
        str(output),
        "--profile",
        str(profile_path),
        "--theme",
        "light",
    ]
    result = subprocess.run(
        command,
        cwd=ROOT,
        env=_isolated_env(tmp_path / "config"),
        capture_output=True,
        text=True,
        timeout=90.0,
        check=False,
    )
    assert result.returncode == 0, result.stdout + result.stderr
    payload = json.loads(result.stdout)
    assert payload["case_count"] == len(CASES) >= 12
    contact = Path(payload["contact_sheet"])
    assert contact.is_file()
    image = QImage(str(contact))
    assert image.width() > 2000 and image.height() > 1000
    report = Path(payload["report"]).read_text(encoding="utf-8")
    assert "Renderer source SHA at export" in report
    assert "Git working tree at export" in report
    for case in CASES:
        assert f"`{case.case_id}`" in report
    assert "solid_wedge_width_ratio" in report
    assert "B4: `geometric fixture`; `non-editorial`; `excluded from preset selection / style calibration`" in report
    assert "outside the approved A1–B3 editorial calibration cohort" in report


def test_cli_loads_named_editorial_preset_for_reproducible_export(tmp_path) -> None:
    from tools.visual.depiction_lab_app.named_presets import (
        EDITORIAL_AB3,
        EDITORIAL_AB3_BASED_ON,
        EDITORIAL_AB3_NAME,
    )

    output = tmp_path / "editorial-ab3-case"
    command = [
        sys.executable,
        str(ROOT / "tools/visual/depiction_lab.py"),
        "--export-case",
        "C3",
        "--output",
        str(output),
        "--preset",
        EDITORIAL_AB3_NAME,
        "--theme",
        "light",
    ]
    result = subprocess.run(
        command,
        cwd=ROOT,
        env=_isolated_env(tmp_path / "editorial-config"),
        capture_output=True,
        text=True,
        timeout=90.0,
        check=False,
    )
    assert result.returncode == 0, result.stdout + result.stderr
    payload = json.loads(result.stdout)
    assert payload["case_id"] == "C3"
    name, parameters = load_profile(output / "working-profile.json")
    raw_profile = json.loads((output / "working-profile.json").read_text(encoding="utf-8"))
    assert name == EDITORIAL_AB3_NAME
    assert parameters == EDITORIAL_AB3.parameters
    assert raw_profile["based_on"] == EDITORIAL_AB3_BASED_ON
    report = (output / "report.md").read_text(encoding="utf-8")
    assert "B4 remains a geometric fixture and is excluded" in report
    assert "outside the approved A1–B3 editorial calibration cohort" in report
