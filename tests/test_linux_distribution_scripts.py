"""Pruebas smoke de scripts de distribución Linux."""

from __future__ import annotations

import importlib.util
import json
import os
import subprocess
from pathlib import Path


def _load_manifest_module():
    script_path = (
        Path(__file__).resolve().parent.parent
        / "packaging"
        / "release"
        / "generate_channel_manifest.py"
    )
    spec = importlib.util.spec_from_file_location("generate_channel_manifest", script_path)
    assert spec is not None and spec.loader is not None
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def test_generate_channel_manifest_ignores_sidecars_and_includes_flatpak(tmp_path) -> None:
    module = _load_manifest_module()
    artifacts_dir = tmp_path / "artifacts"
    artifacts_dir.mkdir(parents=True, exist_ok=True)

    (artifacts_dir / "Chemuson-v1.2.3-linux-x86_64.flatpak").write_bytes(b"flatpak")
    (artifacts_dir / "Chemuson-v1.2.3-linux-x86_64.AppImage").write_bytes(b"appimage")
    (artifacts_dir / "Chemuson-v1.2.3-linux-x86_64.AppImage.updateinfo").write_text(
        "gh-releases-zsync|PJGV333|Chemuson|latest|Chemuson.zsync",
        encoding="utf-8",
    )
    (artifacts_dir / "Chemuson-v1.2.3-linux-x86_64.AppImage.update.json").write_text(
        "{}",
        encoding="utf-8",
    )
    (artifacts_dir / "Chemuson-v1.2.3-linux-x86_64.AppImage.zsync").write_text(
        "zsync",
        encoding="utf-8",
    )

    manifest = module.build_manifest(
        channel="stable",
        version="1.2.3",
        base_url="https://example.invalid/download",
        artifacts_dir=artifacts_dir,
    )

    artifacts = manifest.get("artifacts", {})
    assert "linux-x86_64-flatpak-bundle" in artifacts
    assert "linux-x86_64-appimage" in artifacts
    assert all(not key.endswith(".updateinfo") for key in artifacts.keys())
    assert all(not key.endswith(".update.json") for key in artifacts.keys())
    assert all(not key.endswith(".zsync") for key in artifacts.keys())


def test_build_appimage_script_writes_update_metadata(tmp_path) -> None:
    repo_root = Path(__file__).resolve().parent.parent
    script_path = repo_root / "packaging" / "linux" / "build_appimage.sh"
    dist_dir = tmp_path / "dist"
    out_dir = tmp_path / "dist-appimage"
    dist_dir.mkdir(parents=True, exist_ok=True)

    # Binario ejecutable mínimo para activar ruta fallback del script.
    app_bin = dist_dir / "Chemuson"
    app_bin.write_text("#!/usr/bin/env bash\necho chemuson\n", encoding="utf-8")
    app_bin.chmod(0o755)

    cmd = [
        "bash",
        str(script_path),
        "1.2.3",
        str(dist_dir),
        str(out_dir),
        "PJGV333",
        "Chemuson",
        "beta",
        "v1.2.3-beta.1",
    ]
    subprocess.run(cmd, check=True, cwd=str(repo_root))

    appimage = out_dir / "Chemuson-v1.2.3-linux-x86_64.AppImage"
    updateinfo = out_dir / "Chemuson-v1.2.3-linux-x86_64.AppImage.updateinfo"
    updatejson = out_dir / "Chemuson-v1.2.3-linux-x86_64.AppImage.update.json"

    assert appimage.exists()
    assert os.access(appimage, os.X_OK)
    assert updateinfo.exists()
    assert updatejson.exists()

    update_info_text = updateinfo.read_text(encoding="utf-8").strip()
    assert "gh-releases-zsync|PJGV333|Chemuson|prerelease|" in update_info_text

    payload = json.loads(updatejson.read_text(encoding="utf-8"))
    assert payload["channel"] == "beta"
    assert payload["tag"] == "v1.2.3-beta.1"
