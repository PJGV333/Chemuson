"""Pruebas de telemetría local de actualización."""

import json


from chemuson.update.telemetry import UpdateTelemetryLogger


def test_update_telemetry_writes_jsonl(tmp_path) -> None:
    logger = UpdateTelemetryLogger(log_dir=str(tmp_path))
    logger.log_event(
        "check_result",
        channel="stable",
        mode="notify",
        available=False,
        error_type="TimeoutError",
    )
    log_path = tmp_path / "events.jsonl"
    assert log_path.exists()
    lines = [line for line in log_path.read_text(encoding="utf-8").splitlines() if line.strip()]
    assert len(lines) == 1
    payload = json.loads(lines[0])
    assert payload["event"] == "check_result"
    assert payload["channel"] == "stable"
    assert payload["mode"] == "notify"
    assert payload["available"] is False
    assert payload["error_type"] == "TimeoutError"


def test_update_telemetry_skips_sensitive_like_paths_and_urls(tmp_path) -> None:
    logger = UpdateTelemetryLogger(log_dir=str(tmp_path))
    logger.log_event(
        "download_failed",
        reason_code="network_error",
        artifact_path="/home/user/secret/file.exe",
        source_url="https://github.com/private/token",
        context="windows_installer",
    )
    payload = json.loads((tmp_path / "events.jsonl").read_text(encoding="utf-8").strip())
    assert payload["event"] == "download_failed"
    assert payload["reason_code"] == "network_error"
    assert payload["context"] == "windows_installer"
    assert "artifact_path" not in payload
    assert "source_url" not in payload
