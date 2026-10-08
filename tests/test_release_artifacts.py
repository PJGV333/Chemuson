"""Version/channel/provenance checks for assembled release artifacts."""

from __future__ import annotations

import hashlib
import json
import subprocess
import sys
from pathlib import Path

import pytest

RELEASE_DIR = Path(__file__).resolve().parent.parent / "packaging" / "release"
sys.path.insert(0, str(RELEASE_DIR))

from validate_release_artifacts import validate_release_artifacts  # noqa: E402
from write_build_provenance import build_provenance  # noqa: E402


def _fake_type2_header() -> bytes:
    header = bytearray(64)
    header[:8] = b"\x7fELF\x02\x01\x01\x00"
    header[8:11] = b"AI\x02"
    header[16:18] = (2).to_bytes(2, "little")
    header[18:20] = (62).to_bytes(2, "little")
    return bytes(header)


def _make_release_dir(root: Path, *, version: str = "0.3.0-beta.1", channel: str = "beta") -> str:
    tag = f"v{version}"
    sha = "a" * 40
    appimage_name = f"Chemuson-v{version}-linux-x86_64.AppImage"
    update_track = "prerelease" if channel == "beta" else "latest"
    update_information = f"gh-releases-zsync|PJGV333|Chemuson|{update_track}|{appimage_name}.zsync"
    names = [
        f"Chemuson-v{version}-windows-x86_64-portable.exe",
        f"Chemuson-v{version}-windows-x86_64-setup.exe",
        f"Chemuson-v{version}-linux-x86_64.AppImage",
        f"Chemuson-v{version}-linux-x86_64.AppImage.updateinfo",
        f"Chemuson-v{version}-linux-x86_64.AppImage.update.json",
        f"Chemuson-v{version}-linux-x86_64.AppImage.zsync",
        f"Chemuson-v{version}-linux-x86_64.flatpak",
    ]
    for name in names:
        (root / name).write_bytes(b"artifact")
    (root / appimage_name).write_bytes(_fake_type2_header())
    (root / f"{appimage_name}.updateinfo").write_text(update_information + "\n", encoding="utf-8")
    (root / f"{appimage_name}.zsync").write_text(
        f"zsync: 0.6.2\nURL: https://github.com/PJGV333/Chemuson/releases/download/{tag}/{appimage_name}\n",
        encoding="utf-8",
    )
    (root / f"{appimage_name}.update.json").write_text(
        json.dumps({
            "version": version,
            "channel": channel,
            "tag": tag,
            "source_sha": sha,
            "appimage_update_information": update_information,
        }),
        encoding="utf-8",
    )
    (root / f"Chemuson-{channel}.flatpakref").write_text(
        "[Flatpak Ref]\n"
        f"Branch={channel}\n"
        f"Url=https://example.invalid/flatpak/{channel}/repo/\n",
        encoding="utf-8",
    )
    (root / f"Chemuson-{channel}.flatpakrepo").write_text(
        "[Flatpak Repo]\n"
        f"DefaultBranch={channel}\n"
        f"Url=https://example.invalid/flatpak/{channel}/repo/\n",
        encoding="utf-8",
    )
    (root / "build-provenance.json").write_text(
        json.dumps(build_provenance(version=version, channel=channel, tag=tag, source_sha=sha)),
        encoding="utf-8",
    )
    return sha


def test_build_provenance_and_release_artifact_set_match_exact_release(tmp_path) -> None:
    version = "0.3.0-beta.1"
    tag = f"v{version}"
    sha = _make_release_dir(tmp_path)

    validated = validate_release_artifacts(
        root=tmp_path,
        version=version,
        channel="beta",
        tag=tag,
        source_sha=sha,
    )

    assert len(validated) == 10
    record = json.loads((tmp_path / "build-provenance.json").read_text(encoding="utf-8"))
    assert record["tag"] == tag
    assert record["source_sha"] == sha


def test_artifact_validator_rejects_channel_metadata_mismatch(tmp_path) -> None:
    version = "0.3.0-beta.1"
    sha = _make_release_dir(tmp_path)
    update_path = tmp_path / f"Chemuson-v{version}-linux-x86_64.AppImage.update.json"
    update_path.write_text(
        json.dumps({"version": version, "channel": "stable", "tag": f"v{version}"}),
        encoding="utf-8",
    )

    with pytest.raises(ValueError, match="channel"):
        validate_release_artifacts(
            root=tmp_path,
            version=version,
            channel="beta",
            tag=f"v{version}",
            source_sha=sha,
        )


def test_artifact_validator_rejects_wrong_flatpak_channel(tmp_path) -> None:
    version = "0.3.0-beta.1"
    sha = _make_release_dir(tmp_path)
    ref_path = tmp_path / "Chemuson-beta.flatpakref"
    ref_path.write_text(
        "[Flatpak Ref]\nBranch=stable\nUrl=https://example.invalid/flatpak/stable/repo/\n",
        encoding="utf-8",
    )

    with pytest.raises(ValueError, match="Flatpak ref branch"):
        validate_release_artifacts(
            root=tmp_path,
            version=version,
            channel="beta",
            tag=f"v{version}",
            source_sha=sha,
        )


def test_artifact_validator_rejects_foreign_versioned_asset(tmp_path) -> None:
    version = "0.3.0-beta.1"
    sha = _make_release_dir(tmp_path)
    (tmp_path / "Chemuson-v0.2.5-linux-x86_64.AppImage").write_bytes(b"stale")

    with pytest.raises(ValueError, match="foreign version"):
        validate_release_artifacts(
            root=tmp_path,
            version=version,
            channel="beta",
            tag=f"v{version}",
            source_sha=sha,
        )


def test_release_checksum_generator_covers_packages_and_provenance(tmp_path) -> None:
    package = tmp_path / "chemuson-package.exe"
    provenance = tmp_path / "build-provenance.json"
    package.write_bytes(b"package")
    provenance.write_text('{"source_sha":"' + "a" * 40 + '"}\n', encoding="utf-8")

    subprocess.run(
        [
            sys.executable,
            str(RELEASE_DIR / "generate_checksums.py"),
            "--dir",
            str(tmp_path),
        ],
        check=True,
    )

    lines = (tmp_path / "checksums.txt").read_text(encoding="utf-8")
    for path in (package, provenance):
        expected = hashlib.sha256(path.read_bytes()).hexdigest()
        assert f"{expected}  {path.name}" in lines
        assert (tmp_path / f"{path.name}.sha256").read_text(encoding="utf-8") == (
            f"{expected}  {path.name}\n"
        )


def test_provenance_rejects_version_or_sha_mismatch() -> None:
    with pytest.raises(ValueError, match="does not match"):
        build_provenance(
            version="0.3.0-beta.1",
            channel="beta",
            tag="v0.3.0-beta.2",
            source_sha="a" * 40,
        )
    with pytest.raises(ValueError, match="full Git SHAs"):
        build_provenance(
            version="0.3.0-beta.1",
            channel="beta",
            tag="v0.3.0-beta.1",
            source_sha="short",
        )
