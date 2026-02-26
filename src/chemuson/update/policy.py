"""Políticas de canal, frecuencia y elegibilidad de updates."""

from __future__ import annotations

from datetime import datetime, timedelta, timezone

from chemuson.update.semver import channel_accepts_version, compare_versions, is_newer_version
from chemuson.update.types import UpdateSettings


def now_utc() -> datetime:
    """Entrega fecha actual en UTC."""
    return datetime.now(timezone.utc)


def parse_iso_utc(value: str) -> datetime | None:
    """Parsea ISO8601 y devuelve datetime con tz, o `None`."""
    raw = str(value or "").strip()
    if not raw:
        return None
    try:
        parsed = datetime.fromisoformat(raw.replace("Z", "+00:00"))
    except Exception:
        return None
    if parsed.tzinfo is None:
        return parsed.replace(tzinfo=timezone.utc)
    return parsed.astimezone(timezone.utc)


def should_check_now(settings: UpdateSettings, now: datetime | None = None) -> bool:
    """Valida frecuencia mínima de chequeo según configuración."""
    if not settings.enabled:
        return False
    current = now or now_utc()
    last = parse_iso_utc(settings.last_check_iso)
    if last is None:
        return True
    interval = max(1, int(settings.check_interval_hours))
    return current >= (last + timedelta(hours=interval))


def mark_checked(settings: UpdateSettings, now: datetime | None = None) -> UpdateSettings:
    """Actualiza timestamp de último chequeo en formato ISO."""
    settings.last_check_iso = (now or now_utc()).isoformat()
    return settings


def can_offer_update(
    *,
    channel: str,
    current_version: str,
    candidate_version: str,
) -> bool:
    """Determina si una versión candidata es elegible para el canal."""
    return channel_accepts_version(channel, candidate_version) and is_newer_version(
        candidate_version,
        current_version,
    )


def is_downgrade_suspected(
    observed_version: str,
    highest_seen_version: str,
) -> bool:
    """Detecta posible replay/downgrade respecto a la versión remota máxima conocida."""
    observed = str(observed_version or "").strip()
    highest = str(highest_seen_version or "").strip()
    if not observed or not highest:
        return False
    try:
        return compare_versions(observed, highest) < 0
    except Exception:
        return False
