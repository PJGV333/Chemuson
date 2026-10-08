"""Pruebas smoke de scripts de distribución Linux."""

from __future__ import annotations

import base64
import importlib.util
import json
import os
import subprocess
from pathlib import Path

import pytest

from chemuson import __version__ as APP_VERSION


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


def _write_fake_appimagetool(tmp_path: Path) -> Path:
    tool = tmp_path / "fake-appimagetool"
    tool.write_text(
        "#!/usr/bin/env python3\n"
        "import json, os, pathlib, sys\n"
        "args = sys.argv[1:]\n"
        "if args == ['--version']:\n"
        "    print('appimagetool, continuous build (commit 5735cc5)')\n"
        "    raise SystemExit(0)\n"
        "appdir, output = pathlib.Path(args[-2]), pathlib.Path(args[-1])\n"
        "for required in ('AppRun', 'io.github.PJGV333.Chemuson.desktop', '.DirIcon', 'usr/bin/Chemuson'):\n"
        "    if not (appdir / required).exists(): raise SystemExit(f'missing {required}')\n"
        "header = bytearray(64)\n"
        "header[:8] = b'\\x7fELF\\x02\\x01\\x01\\x00'\n"
        "header[8:11] = b'AI\\x02'\n"
        "header[16:18] = (2).to_bytes(2, 'little')\n"
        "header[18:20] = (62).to_bytes(2, 'little')\n"
        "output.parent.mkdir(parents=True, exist_ok=True)\n"
        "output.write_bytes(header + b'fake-squashfs-payload')\n"
        "output.chmod(0o755)\n"
        "if os.environ.get('APPIMAGE_TOOL_CAPTURE'):\n"
        "    pathlib.Path(os.environ['APPIMAGE_TOOL_CAPTURE']).write_text(json.dumps(args))\n",
        encoding="utf-8",
    )
    tool.chmod(0o755)
    return tool


def _write_fake_zsyncmake(tmp_path: Path) -> Path:
    bindir = tmp_path / "bin"
    bindir.mkdir(exist_ok=True)
    tool = bindir / "zsyncmake"
    tool.write_text(
        "#!/usr/bin/env python3\n"
        "import pathlib, sys\n"
        "args = sys.argv[1:]\n"
        "out = pathlib.Path(args[args.index('-o') + 1])\n"
        "url = args[args.index('-u') + 1]\n"
        "out.write_text('zsync: 0.6.2\\nURL: ' + url + '\\n')\n",
        encoding="utf-8",
    )
    tool.chmod(0o755)
    return bindir


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
        source_sha="a" * 40,
    )

    artifacts = manifest.get("artifacts", {})
    assert manifest["version"] == "1.2.3"
    assert manifest["tag"] == "v1.2.3"
    assert manifest["source_sha"] == "a" * 40
    assert "linux-x86_64-flatpak-bundle" in artifacts
    assert "linux-x86_64-appimage" in artifacts
    assert all(not key.endswith(".updateinfo") for key in artifacts.keys())
    assert all(not key.endswith(".update.json") for key in artifacts.keys())
    assert all(not key.endswith(".zsync") for key in artifacts.keys())

    with pytest.raises(ValueError, match="belongs to channel"):
        module.build_manifest(
            channel="beta",
            version="1.2.3",
            base_url="https://example.invalid/download",
            artifacts_dir=artifacts_dir,
        )


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


def test_build_appimage_script_embeds_existing_update_information(tmp_path) -> None:
    repo_root = Path(__file__).resolve().parent.parent
    script_path = repo_root / "packaging" / "linux" / "build_appimage.sh"
    dist_dir = tmp_path / "dist"
    out_dir = tmp_path / "dist-appimage"
    dist_dir.mkdir(parents=True, exist_ok=True)
    app_bin = dist_dir / "Chemuson"
    app_bin.write_text("#!/usr/bin/env bash\necho chemuson\n", encoding="utf-8")
    app_bin.chmod(0o755)
    fake_tool = _write_fake_appimagetool(tmp_path)
    fake_zsync_bin = _write_fake_zsyncmake(tmp_path)
    capture = tmp_path / "appimagetool-args.json"
    env = os.environ.copy()
    env["APPIMAGETOOL_BIN"] = str(fake_tool)
    env["APPIMAGE_TOOL_CAPTURE"] = str(capture)
    env["PATH"] = f"{fake_zsync_bin}{os.pathsep}{env['PATH']}"
    source_sha = subprocess.run(
        ["git", "rev-parse", "HEAD"], cwd=repo_root, check=True, capture_output=True, text=True
    ).stdout.strip()

    subprocess.run(
        [
            "bash", str(script_path), "1.2.3-beta.1", str(dist_dir), str(out_dir),
            "PJGV333", "Chemuson", "beta", "v1.2.3-beta.1", source_sha,
        ],
        check=True,
        cwd=str(repo_root),
        env=env,
    )

    appimage = out_dir / "Chemuson-v1.2.3-beta.1-linux-x86_64.AppImage"
    updateinfo = Path(f"{appimage}.updateinfo")
    updatejson = Path(f"{appimage}.update.json")
    zsync = Path(f"{appimage}.zsync")
    assert appimage.read_bytes()[8:11] == b"AI\x02"
    assert os.access(appimage, os.X_OK)
    update_information = (
        "gh-releases-zsync|PJGV333|Chemuson|prerelease|"
        "Chemuson-v1.2.3-beta.1-linux-x86_64.AppImage.zsync"
    )
    assert updateinfo.read_text(encoding="utf-8").strip() == update_information
    assert update_information in json.loads(capture.read_text(encoding="utf-8"))
    assert json.loads(updatejson.read_text(encoding="utf-8"))["appimage_update_information"] == update_information
    assert "URL: https://github.com/PJGV333/Chemuson/releases/download/v1.2.3-beta.1/" in zsync.read_text(encoding="utf-8")


def test_build_appimage_preview_creates_type2_without_public_updater_metadata(tmp_path) -> None:
    repo_root = Path(__file__).resolve().parent.parent
    script_path = repo_root / "packaging" / "linux" / "build_appimage.sh"
    dist_dir = tmp_path / "dist"
    out_dir = tmp_path / "dist-preview"
    dist_dir.mkdir(parents=True, exist_ok=True)
    app_bin = dist_dir / "Chemuson"
    app_bin.write_text("#!/usr/bin/env bash\\necho chemuson\\n", encoding="utf-8")
    app_bin.chmod(0o755)
    fake_tool = _write_fake_appimagetool(tmp_path)
    sha = subprocess.run(
        ["git", "rev-parse", "HEAD"],
        cwd=repo_root,
        check=True,
        capture_output=True,
        text=True,
    ).stdout.strip()
    env = os.environ.copy()
    env["APPIMAGETOOL_BIN"] = str(fake_tool)

    subprocess.run(
        [
            "bash", str(script_path), APP_VERSION, str(dist_dir), str(out_dir), "", "",
            "beta", "", sha, "preview", "release/v0.3.0-beta.1-prep",
        ],
        check=True,
        cwd=str(repo_root),
        env=env,
    )

    artifact = out_dir / f"Chemuson-v{APP_VERSION}-preview-{sha[:8]}-linux-x86_64.AppImage"
    assert artifact.read_bytes()[8:11] == b"AI\x02"
    assert os.access(artifact, os.X_OK)
    assert not Path(f"{artifact}.updateinfo").exists()
    assert not Path(f"{artifact}.update.json").exists()
    assert not Path(f"{artifact}.zsync").exists()


def test_preview_flatpak_refuses_public_remote_and_signing_credentials(tmp_path) -> None:
    repo_root = Path(__file__).resolve().parent.parent
    script_path = repo_root / "packaging" / "linux" / "build_flatpak.sh"
    sha = subprocess.run(
        ["git", "rev-parse", "HEAD"],
        cwd=repo_root,
        check=True,
        capture_output=True,
        text=True,
    ).stdout.strip()
    env = os.environ.copy()
    env["CHEMUSON_FLATPAK_REPO_URL"] = "https://example.invalid/flatpak/beta/repo/"
    output_dir = tmp_path / "dist-flatpak-preview"

    result = subprocess.run(
        [
            "bash",
            str(script_path),
            APP_VERSION,
            "beta",
            str(output_dir),
            "packaging/flatpak/io.github.PJGV333.Chemuson.yml",
            "preview",
            sha,
            "release/v0.3.0-beta.1-prep",
        ],
        cwd=repo_root,
        env=env,
        capture_output=True,
        text=True,
    )

    assert result.returncode == 2
    assert "cannot use public remote URLs or signing credentials" in result.stderr
    assert not output_dir.exists()


def test_build_appimage_rejects_mismatched_channel(tmp_path) -> None:
    repo_root = Path(__file__).resolve().parent.parent
    script_path = repo_root / "packaging" / "linux" / "build_appimage.sh"
    dist_dir = tmp_path / "dist"
    dist_dir.mkdir(parents=True, exist_ok=True)
    app_bin = dist_dir / "Chemuson"
    app_bin.write_text("#!/usr/bin/env bash\\necho chemuson\\n", encoding="utf-8")
    app_bin.chmod(0o755)

    result = subprocess.run(
        [
            "bash",
            str(script_path),
            "1.2.3-beta.1",
            str(dist_dir),
            str(tmp_path / "dist-appimage"),
            "PJGV333",
            "Chemuson",
            "stable",
            "v1.2.3-beta.1",
            subprocess.run(
                ["git", "rev-parse", "HEAD"], cwd=repo_root, check=True, capture_output=True, text=True
            ).stdout.strip(),
        ],
        cwd=str(repo_root),
        capture_output=True,
        text=True,
    )
    assert result.returncode != 0
    assert "belongs to channel" in result.stderr
