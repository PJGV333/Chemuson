"""Genera manifiesto de canal (`stable.json` / `beta.json`) para updater."""

from __future__ import annotations

import argparse
import hashlib
import hmac
import json
import os
from datetime import datetime, timezone
from pathlib import Path


def sha256_of(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as fh:
        while True:
            chunk = fh.read(8192)
            if not chunk:
                break
            digest.update(chunk)
    return digest.hexdigest()


def artifact_key(name: str) -> str:
    lowered = name.lower()
    if lowered.endswith(".flatpak"):
        return "linux-x86_64-flatpak-bundle"
    if lowered.endswith(".appimage"):
        return "linux-x86_64-appimage"
    if lowered.endswith("-setup.exe") or "setup" in lowered or "installer" in lowered:
        return "windows-x86_64-installer"
    if lowered.endswith(".exe"):
        return "windows-x86_64-portable"
    return name


def build_manifest(channel: str, version: str, base_url: str, artifacts_dir: Path) -> dict:
    artifacts = {}
    for path in sorted(artifacts_dir.iterdir()):
        if not path.is_file():
            continue
        if path.name.endswith(".sha256") or path.name.endswith(".sig"):
            continue
        if path.name.endswith(".zsync") or path.name.endswith(".updateinfo"):
            continue
        if path.name.endswith(".update.json") or path.name.endswith(".flatpakref"):
            continue
        if path.name in {"checksums.txt", "stable.json", "beta.json"}:
            continue
        artifacts[artifact_key(path.name)] = {
            "version": version,
            "url": f"{base_url.rstrip('/')}/{path.name}",
            "sha256": sha256_of(path),
        }
    return {
        "channel": channel,
        "latest": version,
        "published_at": datetime.now(timezone.utc).isoformat(),
        "artifacts": artifacts,
    }


def sign_manifest(manifest: dict, key: str) -> dict:
    payload = json.dumps(manifest, sort_keys=True, separators=(",", ":")).encode("utf-8")
    signature = hmac.new(key.encode("utf-8"), payload, hashlib.sha256).hexdigest()
    signed = dict(manifest)
    signed["signature"] = {"algorithm": "hmac-sha256", "value": signature}
    return signed


def main() -> None:
    parser = argparse.ArgumentParser(description="Genera manifest de updates por canal.")
    parser.add_argument("--channel", required=True, choices=["stable", "beta"])
    parser.add_argument("--version", required=True, help="Versión semántica.")
    parser.add_argument("--artifacts-dir", required=True, help="Directorio de artifacts.")
    parser.add_argument("--base-url", required=True, help="Base URL de descarga de assets.")
    parser.add_argument("--output", required=True, help="Ruta destino del JSON.")
    parser.add_argument("--key", default="", help="Clave HMAC opcional.")
    args = parser.parse_args()

    manifest = build_manifest(
        channel=args.channel,
        version=args.version,
        base_url=args.base_url,
        artifacts_dir=Path(args.artifacts_dir).resolve(),
    )
    key = args.key or os.getenv("CHEMUSON_SIGN_KEY", "")
    if key:
        manifest = sign_manifest(manifest, key)
    output = Path(args.output).resolve()
    output.write_text(json.dumps(manifest, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    print(f"Wrote manifest: {output}")


if __name__ == "__main__":
    main()
