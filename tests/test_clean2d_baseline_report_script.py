from __future__ import annotations

import json
from pathlib import Path

from scripts.clean2d_baseline_report import EXIT_CHANGED, EXIT_ERROR, EXIT_OK, main


def test_write_creates_json_stable_report(tmp_path: Path) -> None:
    output = tmp_path / "report.json"

    assert main(["write", "--output", str(output)]) == EXIT_OK
    data = json.loads(output.read_text(encoding="utf-8"))

    assert data["schema"] == "chemuson.clean2d.baseline-report"
    assert data["version"] == 1
    assert data["records"]
    assert data["summary"]["case_count"] == len(data["records"])


def test_compare_report_with_itself_returns_zero(tmp_path: Path, capsys) -> None:
    report = _write_report(tmp_path / "report.json")

    assert main(["compare", "--left", str(report), "--right", str(report)]) == EXIT_OK
    diff = json.loads(capsys.readouterr().out)

    assert diff["equivalent"] is True


def test_compare_changed_report_returns_one(tmp_path: Path, capsys) -> None:
    left = _write_report(tmp_path / "left.json")
    right = _changed_report(left, tmp_path / "right.json")

    assert main(["compare", "--left", str(left), "--right", str(right)]) == EXIT_CHANGED
    diff = json.loads(capsys.readouterr().out)

    assert diff["equivalent"] is False
    assert diff["changed_cases"]


def test_review_equivalent_reports_returns_zero(tmp_path: Path, capsys) -> None:
    report = _write_report(tmp_path / "report.json")

    assert main(["review", "--left", str(report), "--right", str(report)]) == EXIT_OK
    summary = json.loads(capsys.readouterr().out)

    assert summary["item_count"] == 0
    assert summary["observational_only"] is True


def test_review_changed_report_returns_one_and_summary(tmp_path: Path, capsys) -> None:
    left = _write_report(tmp_path / "left.json")
    right = _changed_report(left, tmp_path / "right.json")

    assert main(["review", "--left", str(left), "--right", str(right)]) == EXIT_CHANGED
    summary = json.loads(capsys.readouterr().out)

    assert summary["item_count"] == 1
    assert summary["observational_only"] is True
    assert summary["items"][0]["field"] == "result_state"


def test_invalid_json_returns_error(tmp_path: Path, capsys) -> None:
    invalid = tmp_path / "invalid.json"
    invalid.write_text("not-json", encoding="utf-8")
    report = _write_report(tmp_path / "report.json")

    assert main(["compare", "--left", str(invalid), "--right", str(report)]) == EXIT_ERROR
    assert "error:" in capsys.readouterr().err


def test_missing_file_returns_error(tmp_path: Path, capsys) -> None:
    report = _write_report(tmp_path / "report.json")

    assert main(["review", "--left", str(tmp_path / "missing.json"), "--right", str(report)]) == EXIT_ERROR
    assert "error:" in capsys.readouterr().err


def test_script_is_developer_test_only() -> None:
    import scripts.clean2d_baseline_report as script

    assert "/scripts/" in script.__file__


def _write_report(path: Path) -> Path:
    assert main(["write", "--output", str(path)]) == EXIT_OK
    return path


def _changed_report(source: Path, target: Path) -> Path:
    data = json.loads(source.read_text(encoding="utf-8"))
    record = data["records"][0]
    record["result_state"] = "failed-controlled" if record["result_state"] != "failed-controlled" else "applied"
    target.write_text(json.dumps(data, allow_nan=False, sort_keys=True), encoding="utf-8")
    return target
