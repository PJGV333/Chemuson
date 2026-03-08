"""Proveedor de releases desde GitHub para auto-update."""

from __future__ import annotations

import json
import os
import platform
import re
import time
from dataclasses import dataclass
from datetime import datetime, timedelta, timezone
from pathlib import Path
from typing import Any, Optional
from urllib.request import Request, urlopen

from chemuson.update.semver import parse_semver
from chemuson.update.types import (
    ReleaseAsset,
    ReleaseInfo,
    UpdateChannel,
    coerce_update_channel,
)

_ASSET_VERSION_RE = re.compile(r"v?(\d+\.\d+\.\d+(?:-[0-9A-Za-z.-]+)?)")
_OWNER_REPO_RE = re.compile(r"^[A-Za-z0-9_.-]+$")


@dataclass(slots=True)
class AssetMatchRules:
    """Reglas de matching de nombres de artefacto por plataforma."""

    preferred_tokens: tuple[str, ...]
    required_extensions: tuple[str, ...]


def detect_platform_tag() -> str:
    """Devuelve clave de plataforma compatible con parser de assets."""
    sys_name = platform.system().lower()
    machine = platform.machine().lower()
    if machine in {"amd64", "x86_64", "x64"}:
        arch = "x86_64"
    elif machine in {"arm64", "aarch64"}:
        arch = "arm64"
    else:
        arch = machine.replace(" ", "_")
    if "windows" in sys_name:
        return f"windows-{arch}"
    if "linux" in sys_name:
        return f"linux-{arch}"
    if "darwin" in sys_name or "mac" in sys_name:
        return f"macos-{arch}"
    return f"{sys_name}-{arch}"


def _version_from_tag_or_name(tag_name: str, name: str = "") -> str:
    """Extrae versión semántica desde tag o nombre de release."""
    for raw in (tag_name, name):
        if not raw:
            continue
        try:
            return parse_semver(raw).normalized()
        except ValueError:
            pass
        match = _ASSET_VERSION_RE.search(raw)
        if match:
            try:
                return parse_semver(match.group(1)).normalized()
            except ValueError:
                continue
    return ""


def _asset_rules_for_platform(platform_tag: str) -> AssetMatchRules:
    """Reglas de selección según plataforma destino."""
    normalized = platform_tag.lower()
    if normalized.startswith("windows-"):
        return AssetMatchRules(
            preferred_tokens=("setup", "installer", "portable", "win"),
            required_extensions=(".exe", ".msi"),
        )
    if normalized.startswith("linux-"):
        return AssetMatchRules(
            preferred_tokens=("appimage", "linux"),
            required_extensions=(".appimage",),
        )
    if normalized.startswith("macos-"):
        return AssetMatchRules(
            preferred_tokens=("macos", "darwin", "universal", "arm64", "x86_64"),
            required_extensions=(".dmg", ".pkg", ".zip"),
        )
    return AssetMatchRules(preferred_tokens=(), required_extensions=())


def _is_installer_asset_name(name: str) -> bool:
    """Indica si un nombre de asset parece ser instalador de escritorio."""
    lowered = str(name or "").strip().lower()
    return lowered.endswith(".msi") or "setup" in lowered or "installer" in lowered


def _is_portable_asset_name(name: str) -> bool:
    """Indica si un asset parece ser un binario portable."""
    lowered = str(name or "").strip().lower()
    if lowered.endswith(".appimage"):
        return True
    if lowered.endswith(".exe") and not _is_installer_asset_name(lowered):
        return True
    return "portable" in lowered


class GitHubReleasesProvider:
    """Proveedor que consulta releases públicos vía API GitHub."""

    def __init__(
        self,
        owner: str,
        repo: str,
        timeout: float = 8.0,
        api_base: str = "https://api.github.com",
        cache_dir: str | None = None,
        cache_max_age_hours: int = 168,
        allow_cached_fallback: bool = True,
        retries: int = 2,
        retry_backoff_sec: float = 0.35,
        max_payload_bytes: int = 4 * 1024 * 1024,
    ) -> None:
        self.owner = str(owner).strip()
        self.repo = str(repo).strip()
        if not _OWNER_REPO_RE.match(self.owner or ""):
            raise ValueError("Owner de GitHub inválido para proveedor de updates.")
        if not _OWNER_REPO_RE.match(self.repo or ""):
            raise ValueError("Repo de GitHub inválido para proveedor de updates.")
        self.timeout = float(timeout)
        self.api_base = api_base.rstrip("/")
        if not self.api_base.lower().startswith("https://"):
            raise ValueError("api_base debe usar HTTPS.")
        self.cache_dir = (
            os.path.abspath(cache_dir)
            if cache_dir
            else os.path.join(os.path.expanduser("~"), ".chemuson", "update_cache")
        )
        self.cache_max_age_hours = max(1, int(cache_max_age_hours))
        self.allow_cached_fallback = bool(allow_cached_fallback)
        self.retries = max(0, int(retries))
        self.retry_backoff_sec = max(0.0, float(retry_backoff_sec))
        self.max_payload_bytes = max(1024, int(max_payload_bytes))
        self.last_fetch_source = "none"
        self.last_fetch_error = ""

    def _api_url(self, suffix: str) -> str:
        url = f"{self.api_base}/repos/{self.owner}/{self.repo}{suffix}"
        if not url.lower().startswith("https://"):
            raise ValueError("URL de API insegura.")
        return url

    def _cache_path(self) -> str:
        """Devuelve ruta local de caché del payload crudo de releases."""
        safe_owner = re.sub(r"[^A-Za-z0-9_.-]+", "_", self.owner or "owner")
        safe_repo = re.sub(r"[^A-Za-z0-9_.-]+", "_", self.repo or "repo")
        name = f"github_releases_{safe_owner}_{safe_repo}.json"
        return os.path.join(self.cache_dir, name)

    def _write_cache(self, payload: list[dict[str, Any]]) -> None:
        """Persiste payload remoto para fallback offline seguro."""
        os.makedirs(self.cache_dir, exist_ok=True)
        wrapper = {
            "fetched_at": datetime.now(timezone.utc).isoformat(),
            "payload": payload,
        }
        path = Path(self._cache_path())
        tmp_path = path.with_suffix(".tmp")
        tmp_path.write_text(json.dumps(wrapper, ensure_ascii=False), encoding="utf-8")
        os.replace(str(tmp_path), str(path))

    def _read_cache(self) -> list[dict[str, Any]] | None:
        """Carga caché si está vigente; devuelve `None` cuando no aplica."""
        path = Path(self._cache_path())
        if not path.exists():
            return None
        try:
            raw = json.loads(path.read_text(encoding="utf-8"))
        except Exception:
            return None
        if not isinstance(raw, dict):
            return None
        fetched_at = str(raw.get("fetched_at") or "").strip()
        payload = raw.get("payload")
        if not isinstance(payload, list):
            return None
        try:
            stamp = datetime.fromisoformat(fetched_at.replace("Z", "+00:00"))
            if stamp.tzinfo is None:
                stamp = stamp.replace(tzinfo=timezone.utc)
            now = datetime.now(timezone.utc)
            if now > (stamp.astimezone(timezone.utc) + timedelta(hours=self.cache_max_age_hours)):
                return None
        except Exception:
            return None
        clean_payload: list[dict[str, Any]] = []
        for item in payload:
            if isinstance(item, dict):
                clean_payload.append(item)
        return clean_payload

    def _read_json(self, url: str) -> Any:
        if not str(url or "").lower().startswith("https://"):
            raise ValueError("Solo se permiten endpoints HTTPS.")
        attempts = self.retries + 1
        last_error: Exception | None = None
        for attempt in range(attempts):
            try:
                request = Request(
                    url=url,
                    headers={
                        "Accept": "application/vnd.github+json",
                        "User-Agent": "chemuson-updater",
                    },
                )
                with urlopen(request, timeout=self.timeout) as response:
                    status = int(getattr(response, "status", 200))
                    if status < 200 or status >= 300:
                        raise RuntimeError(f"HTTP status inesperado: {status}")
                    payload = response.read(self.max_payload_bytes + 1)
                if len(payload) > self.max_payload_bytes:
                    raise ValueError("Payload de API excede tamaño máximo permitido.")
                decoded = payload.decode("utf-8")
                parsed = json.loads(decoded)
                if not isinstance(parsed, (list, dict)):
                    raise ValueError("Payload de API no es JSON esperado.")
                return parsed
            except Exception as exc:
                last_error = exc
                if attempt >= attempts - 1:
                    raise
                time.sleep(self.retry_backoff_sec * (2 ** attempt))
        if last_error is not None:
            raise last_error
        raise RuntimeError("Error inesperado leyendo API JSON.")

    def fetch_releases_payload(self) -> list[dict[str, Any]]:
        """Devuelve payload crudo de releases, excluyendo drafts."""
        self.last_fetch_source = "none"
        self.last_fetch_error = ""
        try:
            payload = self._read_json(self._api_url("/releases?per_page=30"))
            if not isinstance(payload, list):
                payload = []
            releases = []
            for item in payload:
                if not isinstance(item, dict):
                    continue
                if bool(item.get("draft", False)):
                    continue
                releases.append(item)
            self._write_cache(releases)
            self.last_fetch_source = "remote"
            return releases
        except Exception as exc:
            self.last_fetch_error = exc.__class__.__name__
            if self.allow_cached_fallback:
                cached = self._read_cache()
                if cached is not None:
                    self.last_fetch_source = "cache"
                    return cached
            self.last_fetch_source = "error"
            raise

    def _parse_assets(self, assets_payload: list[dict[str, Any]]) -> list[ReleaseAsset]:
        assets: list[ReleaseAsset] = []
        asset_by_name: dict[str, ReleaseAsset] = {}
        for raw in assets_payload:
            name = str(raw.get("name") or "").strip()
            url = str(raw.get("browser_download_url") or "").strip()
            if not name or not url:
                continue
            if not url.lower().startswith("https://"):
                continue
            asset = ReleaseAsset(
                name=name,
                url=url,
                size=int(raw.get("size") or 0),
                content_type=str(raw.get("content_type") or ""),
            )
            assets.append(asset)
            asset_by_name[name] = asset

        # Enlaza sidecars por convención: <asset>.sha256 y <asset>.sig
        for asset in assets:
            hash_name = f"{asset.name}.sha256"
            sig_name = f"{asset.name}.sig"
            hash_asset = asset_by_name.get(hash_name)
            sig_asset = asset_by_name.get(sig_name)
            if hash_asset is not None:
                asset.sha256_url = hash_asset.url
            if sig_asset is not None:
                asset.signature_url = sig_asset.url
        return assets

    def parse_release(self, payload: dict[str, Any]) -> Optional[ReleaseInfo]:
        """Convierte payload de API en `ReleaseInfo`."""
        tag_name = str(payload.get("tag_name") or "").strip()
        version = _version_from_tag_or_name(tag_name, str(payload.get("name") or ""))
        if not version:
            return None
        assets = self._parse_assets(list(payload.get("assets") or []))
        html_url = str(payload.get("html_url") or "").strip()
        if html_url and not html_url.lower().startswith("https://"):
            html_url = ""
        return ReleaseInfo(
            tag_name=tag_name,
            version=version,
            prerelease=bool(payload.get("prerelease", False)),
            published_at=str(payload.get("published_at") or ""),
            html_url=html_url,
            assets=assets,
        )

    def list_releases(self) -> list[ReleaseInfo]:
        """Lista releases parseadas y ordenadas por fecha descendente."""
        releases: list[ReleaseInfo] = []
        for payload in self.fetch_releases_payload():
            parsed = self.parse_release(payload)
            if parsed is not None:
                releases.append(parsed)

        def _key(release: ReleaseInfo) -> datetime:
            raw = (release.published_at or "").replace("Z", "+00:00")
            try:
                return datetime.fromisoformat(raw)
            except Exception:
                return datetime.min

        releases.sort(key=_key, reverse=True)
        return releases

    def latest_for_channel(self, channel: UpdateChannel | str) -> Optional[ReleaseInfo]:
        """Obtiene release candidata según canal `stable`/`beta`."""
        ch = coerce_update_channel(channel)
        for release in self.list_releases():
            if ch == UpdateChannel.STABLE and release.prerelease:
                continue
            return release
        return None

    def find_platform_asset(
        self,
        release: ReleaseInfo,
        platform_tag: Optional[str] = None,
        preferred_flavor: Optional[str] = None,
    ) -> Optional[ReleaseAsset]:
        """Selecciona el asset más apropiado para la plataforma."""
        if not release.assets:
            return None
        tag = (platform_tag or detect_platform_tag()).lower()
        rules = _asset_rules_for_platform(tag)
        flavor = str(preferred_flavor or "").strip().lower()
        if flavor not in {"installer", "portable"}:
            flavor = ""
        candidates: list[tuple[int, ReleaseAsset]] = []

        for asset in release.assets:
            name = asset.name.lower()
            if name.endswith(".sha256") or name.endswith(".sig"):
                continue
            if rules.required_extensions and not any(
                name.endswith(ext) for ext in rules.required_extensions
            ):
                continue
            score = 0
            if any(token in name for token in rules.preferred_tokens):
                score += 20
            if "portable" in name:
                score += 5
            if "setup" in name or "installer" in name:
                score += 5
            if tag.startswith("windows-") and ("win" in name or "windows" in name):
                score += 10
            if tag.startswith("linux-") and "linux" in name:
                score += 10
            arch = tag.split("-", 1)[1] if "-" in tag else ""
            if arch and arch in name:
                score += 10
            if flavor == "installer":
                if _is_installer_asset_name(name):
                    score += 30
                if _is_portable_asset_name(name):
                    score -= 10
            elif flavor == "portable":
                if _is_portable_asset_name(name):
                    score += 30
                if _is_installer_asset_name(name):
                    score -= 10
            candidates.append((score, asset))

        if not candidates:
            return None
        candidates.sort(key=lambda item: item[0], reverse=True)
        return candidates[0][1]
