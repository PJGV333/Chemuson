"""Pruebas smoke de scripts de distribución Linux."""

from __future__ import annotations

import base64
import importlib.util
import json
import os
import subprocess
from pathlib import Path

import pytest


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


def _load_flatpak_remote_module():
    script_path = (
        Path(__file__).resolve().parent.parent
        / "packaging"
        / "release"
        / "generate_flatpak_remote_files.py"
    )
    spec = importlib.util.spec_from_file_location(
        "generate_flatpak_remote_files", script_path
    )
    assert spec is not None and spec.loader is not None
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def _load_flatpak_pages_index_module():
    script_path = (
        Path(__file__).resolve().parent.parent
        / "packaging"
        / "release"
        / "generate_flatpak_pages_index.py"
    )
    spec = importlib.util.spec_from_file_location(
        "generate_flatpak_pages_index", script_path
    )
    assert spec is not None and spec.loader is not None
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def _load_flatpak_validate_module():
    script_path = (
        Path(__file__).resolve().parent.parent
        / "packaging"
        / "release"
        / "validate_flatpak_remote_artifacts.py"
    )
    spec = importlib.util.spec_from_file_location(
        "validate_flatpak_remote_artifacts", script_path
    )
    assert spec is not None and spec.loader is not None
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def _write_flatpak_remote_configs(
    root: Path,
    *,
    channel: str,
    gpg_key: str = "",
) -> None:
    repo_lines = [
        "[Flatpak Repo]",
        "Title=Chemuson",
        "Url=https://example.invalid/repo/",
        f"DefaultBranch={channel}",
    ]
    ref_lines = [
        "[Flatpak Ref]",
        "Title=Chemuson",
        "Name=io.github.PJGV333.Chemuson",
        f"Branch={channel}",
        "IsRuntime=false",
        "Url=https://example.invalid/repo/",
        "RuntimeRepo=https://dl.flathub.org/repo/flathub.flatpakrepo",
    ]
    if gpg_key:
        repo_lines.append(f"GPGKey={gpg_key}")
        ref_lines.append(f"GPGKey={gpg_key}")

    (root / f"Chemuson-{channel}.flatpakrepo").write_text(
        "\n".join(repo_lines) + "\n",
        encoding="utf-8",
    )
    (root / f"Chemuson-{channel}.flatpakref").write_text(
        "\n".join(ref_lines) + "\n",
        encoding="utf-8",
    )


def test_generate_channel_manifest_ignores_sidecars_and_includes_flatpak(tmp_path) -> None:
    module = _load_manifest_module()
    artifacts_dir = tmp_path / "artifacts"
    artifacts_dir.mkdir(parents=True, exist_ok=True)

    (artifacts_dir / "Chemuson-v1.2.3-linux-x86_64.flatpak").write_bytes(b"flatpak")
    (artifacts_dir / "Chemuson-v1.2.3-linux-x86_64.AppImage").write_bytes(b"appimage")
    (artifacts_dir / "Chemuson-stable.flatpakrepo").write_text(
        "[Flatpak Repo]\nTitle=Chemuson\n",
        encoding="utf-8",
    )
    (artifacts_dir / "Chemuson-stable.flatpakref").write_text(
        "[Flatpak Ref]\nTitle=Chemuson\n",
        encoding="utf-8",
    )
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


def test_generate_flatpak_remote_files_builds_repo_and_ref_payloads() -> None:
    module = _load_flatpak_remote_module()

    repo_text = module.build_flatpak_repo_config(
        title="Chemuson (beta)",
        repo_url="https://pjgv333.github.io/Chemuson/flatpak/beta/repo",
        homepage="https://github.com/PJGV333/Chemuson",
        comment="Canal oficial beta",
        description="Repositorio oficial beta",
        icon_url="https://pjgv333.github.io/Chemuson/flatpak/icon.svg",
        default_branch="beta",
    )
    ref_text = module.build_flatpak_ref_config(
        title="Chemuson (beta)",
        app_id="io.github.PJGV333.Chemuson",
        branch="beta",
        repo_url="https://pjgv333.github.io/Chemuson/flatpak/beta/repo",
        runtime_repo="https://dl.flathub.org/repo/flathub.flatpakrepo",
        suggest_remote_name="chemuson-beta",
        gpg_key="ZmFrZS1rZXk=",
    )

    assert "[Flatpak Repo]" in repo_text
    assert "Title=Chemuson (beta)" in repo_text
    assert "Url=https://pjgv333.github.io/Chemuson/flatpak/beta/repo/" in repo_text
    assert "DefaultBranch=beta" in repo_text
    assert "[Flatpak Ref]" in ref_text
    assert "Branch=beta" in ref_text
    assert "Name=io.github.PJGV333.Chemuson" in ref_text
    assert "RuntimeRepo=https://dl.flathub.org/repo/flathub.flatpakrepo" in ref_text
    assert "SuggestRemoteName=chemuson-beta" in ref_text
    assert "GPGKey=ZmFrZS1rZXk=" in ref_text


def test_validate_flatpak_remote_artifacts_accepts_build_output_and_pages_payload(
    tmp_path,
) -> None:
    module = _load_flatpak_validate_module()

    build_root = tmp_path / "dist-flatpak"
    repo_root = build_root / "repo"
    (repo_root / "objects" / "00").mkdir(parents=True, exist_ok=True)
    (repo_root / "refs" / "heads").mkdir(parents=True, exist_ok=True)
    _write_flatpak_remote_configs(build_root, channel="stable")
    (repo_root / "config").write_text("config", encoding="utf-8")
    (repo_root / "summary").write_text("summary", encoding="utf-8")
    (repo_root / "objects" / "00" / "payload.filez").write_text(
        "payload",
        encoding="utf-8",
    )
    (repo_root / "refs" / "heads" / "app").write_text("ref", encoding="utf-8")

    module.validate_build_output(root=build_root, basename="Chemuson", channel="stable")

    payload_root = tmp_path / "flatpak-remote"
    channel_root = payload_root / "stable"
    repo_payload_root = channel_root / "repo"
    (repo_payload_root / "objects" / "00").mkdir(parents=True, exist_ok=True)
    (repo_payload_root / "refs" / "heads").mkdir(parents=True, exist_ok=True)
    (payload_root / "icon.svg").write_text("<svg />", encoding="utf-8")
    _write_flatpak_remote_configs(channel_root, channel="stable")
    (repo_payload_root / "config").write_text("config", encoding="utf-8")
    (repo_payload_root / "summary").write_text("summary", encoding="utf-8")
    (repo_payload_root / "objects" / "00" / "payload.filez").write_text(
        "payload",
        encoding="utf-8",
    )
    (repo_payload_root / "refs" / "heads" / "app").write_text(
        "ref",
        encoding="utf-8",
    )

    module.validate_publish_payload(
        root=payload_root,
        basename="Chemuson",
        channel="stable",
    )


def test_validate_flatpak_remote_artifacts_requires_commitmeta_for_signed_repo(
    tmp_path,
    monkeypatch,
) -> None:
    module = _load_flatpak_validate_module()
    gpg_key = base64.b64encode(b"fake-gpg-key").decode("ascii")
    commit = "a" * 64

    build_root = tmp_path / "dist-flatpak"
    repo_root = build_root / "repo"
    ref_path = repo_root / "refs" / "heads" / "app" / "io.github.PJGV333.Chemuson" / "x86_64"
    (repo_root / "objects" / "aa").mkdir(parents=True, exist_ok=True)
    ref_path.mkdir(parents=True, exist_ok=True)
    _write_flatpak_remote_configs(build_root, channel="stable", gpg_key=gpg_key)
    (repo_root / "config").write_text("config", encoding="utf-8")
    (repo_root / "summary").write_text("summary", encoding="utf-8")
    (repo_root / "summary.sig").write_text("sig", encoding="utf-8")
    (repo_root / "objects" / "aa" / "payload.filez").write_text("payload", encoding="utf-8")
    (ref_path / "stable").write_text(f"{commit}\n", encoding="utf-8")

    def fake_run(args, check, capture_output, text):
        assert check is True
        if args[:2] == ["ostree", "rev-parse"]:
            return subprocess.CompletedProcess(args, 0, stdout=f"{commit}\n", stderr="")
        raise AssertionError(f"Unexpected command: {args}")

    monkeypatch.setattr(module.subprocess, "run", fake_run)

    with pytest.raises(FileNotFoundError):
        module.validate_build_output(root=build_root, basename="Chemuson", channel="stable")


def test_validate_flatpak_remote_artifacts_verifies_signed_app_ref(
    tmp_path,
    monkeypatch,
) -> None:
    module = _load_flatpak_validate_module()
    gpg_key = base64.b64encode(b"fake-gpg-key").decode("ascii")
    commit = "b" * 64
    commands: list[list[str]] = []

    build_root = tmp_path / "dist-flatpak"
    repo_root = build_root / "repo"
    ref_path = repo_root / "refs" / "heads" / "app" / "io.github.PJGV333.Chemuson" / "x86_64"
    commitmeta_path = repo_root / "objects" / "bb" / f"{'b' * 62}.commitmeta"
    payload_path = repo_root / "objects" / "bb" / "payload.filez"
    ref_path.mkdir(parents=True, exist_ok=True)
    commitmeta_path.parent.mkdir(parents=True, exist_ok=True)
    _write_flatpak_remote_configs(build_root, channel="stable", gpg_key=gpg_key)
    (repo_root / "config").write_text("config", encoding="utf-8")
    (repo_root / "summary").write_text("summary", encoding="utf-8")
    (repo_root / "summary.sig").write_text("sig", encoding="utf-8")
    payload_path.write_text("payload", encoding="utf-8")
    commitmeta_path.write_text("signed", encoding="utf-8")
    (ref_path / "stable").write_text(f"{commit}\n", encoding="utf-8")

    def fake_run(args, check, capture_output, text):
        assert check is True
        commands.append(args)
        if args[:2] == ["ostree", "rev-parse"]:
            return subprocess.CompletedProcess(args, 0, stdout=f"{commit}\n", stderr="")
        return subprocess.CompletedProcess(args, 0, stdout="", stderr="")

    monkeypatch.setattr(module.subprocess, "run", fake_run)

    module.validate_build_output(root=build_root, basename="Chemuson", channel="stable")

    assert any(args[:2] == ["ostree", "pull"] for args in commands)


def test_generate_flatpak_pages_index_only_links_existing_channels(tmp_path) -> None:
    module = _load_flatpak_pages_index_module()

    stable_dir = tmp_path / "flatpak" / "stable" / "repo"
    stable_dir.mkdir(parents=True, exist_ok=True)
    (tmp_path / "flatpak" / "stable" / "Chemuson-stable.flatpakref").write_text(
        "ref",
        encoding="utf-8",
    )
    (tmp_path / "flatpak" / "stable" / "Chemuson-stable.flatpakrepo").write_text(
        "repo",
        encoding="utf-8",
    )
    (stable_dir / "summary").write_text("summary", encoding="utf-8")

    channels = module.collect_channels(root=tmp_path, basename="Chemuson")
    html = module.build_index_html(channels)

    assert len(channels) == 1
    assert "./flatpak/stable/Chemuson-stable.flatpakref" in html
    assert "./flatpak/stable/Chemuson-stable.flatpakrepo" in html
    assert "./flatpak/beta/Chemuson-beta.flatpakref" not in html


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
