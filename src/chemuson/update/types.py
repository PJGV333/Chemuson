"""Tipos compartidos para el subsistema de actualización."""

from __future__ import annotations

from dataclasses import dataclass, field
from enum import Enum
from typing import Optional


class UpdateChannel(str, Enum):
    """Canales de publicación soportados."""

    STABLE = "stable"
    BETA = "beta"


class UpdateMode(str, Enum):
    """Modo de notificación/ejecución de updates."""

    NOTIFY = "notify"
    SILENT = "silent"


def coerce_update_channel(value: UpdateChannel | str) -> UpdateChannel:
    """Convierte string/enum a `UpdateChannel`."""
    if isinstance(value, UpdateChannel):
        return value
    return UpdateChannel(str(value or "").strip().lower())


def coerce_update_mode(value: UpdateMode | str) -> UpdateMode:
    """Convierte string/enum a `UpdateMode`."""
    if isinstance(value, UpdateMode):
        return value
    return UpdateMode(str(value or "").strip().lower())


@dataclass(slots=True)
class UpdateSettings:
    """Configuración persistente de actualización."""

    enabled: bool = False
    channel: UpdateChannel = UpdateChannel.STABLE
    mode: UpdateMode = UpdateMode.NOTIFY
    check_interval_hours: int = 24
    last_check_iso: str = ""
    highest_seen_version: str = ""
    require_sha256: bool = True
    require_signature: bool = False


@dataclass(slots=True)
class ReleaseAsset:
    """Asset individual de GitHub Release."""

    name: str
    url: str
    size: int = 0
    content_type: str = ""
    sha256_url: str = ""
    signature_url: str = ""


@dataclass(slots=True)
class ReleaseInfo:
    """Metadatos de release relevantes para update."""

    tag_name: str
    version: str
    prerelease: bool
    published_at: str
    html_url: str
    assets: list[ReleaseAsset] = field(default_factory=list)


@dataclass(slots=True)
class UpdateCandidate:
    """Release + asset objetivo para actualizar."""

    release: ReleaseInfo
    artifact: ReleaseAsset


@dataclass(slots=True)
class UpdateCheckResult:
    """Resultado de chequeo de update."""

    available: bool
    current_version: str
    latest_version: str = ""
    channel: str = ""
    reason: str = ""
    candidate: Optional[UpdateCandidate] = None


@dataclass(slots=True)
class DownloadedUpdate:
    """Rutas locales de un update descargado."""

    candidate: UpdateCandidate
    artifact_path: str
    checksum_path: str = ""
    signature_path: str = ""
