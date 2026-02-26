"""Telemetría local mínima del subsistema de actualización.

No envía datos a red. Solo guarda eventos técnicos básicos en disco.
"""

from __future__ import annotations

import json
import os
from datetime import datetime, timezone
from pathlib import Path

_ALLOWED_FIELDS = {
    "event",
    "channel",
    "mode",
    "available",
    "current_version",
    "latest_version",
    "source",
    "reason_code",
    "error_type",
    "force",
    "interactive",
    "action",
    "result",
    "context",
}


def _now_iso() -> str:
    return datetime.now(timezone.utc).isoformat()


class UpdateTelemetryLogger:
    """Logger de eventos de update en JSON Lines (local-only)."""

    def __init__(self, log_dir: str | None = None) -> None:
        base_dir = (
            os.path.abspath(log_dir)
            if log_dir
            else os.path.join(os.path.expanduser("~"), ".chemuson", "update_logs")
        )
        self.log_dir = base_dir
        self.log_path = os.path.join(base_dir, "events.jsonl")

    def _sanitize(self, payload: dict) -> dict:
        clean: dict = {}
        for key, value in payload.items():
            if key not in _ALLOWED_FIELDS:
                continue
            if isinstance(value, bool):
                clean[key] = value
                continue
            if isinstance(value, (int, float)):
                clean[key] = value
                continue
            text = str(value or "").strip()
            # Evita rutas/URLs y recorta tamaño.
            if "://" in text or text.startswith("/") or ":\\" in text:
                continue
            clean[key] = text[:120]
        return clean

    def log_event(self, event: str, **fields) -> None:
        """Registra un evento no sensible.

        Si ocurre error de escritura, se ignora para no afectar UX.
        """
        try:
            os.makedirs(self.log_dir, exist_ok=True)
            payload = {"timestamp": _now_iso(), "event": str(event or "").strip()[:80]}
            payload.update(self._sanitize(fields))
            with Path(self.log_path).open("a", encoding="utf-8") as fh:
                fh.write(json.dumps(payload, ensure_ascii=False) + "\n")
        except Exception:
            return
