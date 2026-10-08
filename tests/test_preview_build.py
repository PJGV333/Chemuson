"""Exact-source identity and provenance tests for preview build helpers."""

from __future__ import annotations

import json
import subprocess
import sys
from pathlib import Path

import pytest

RELEASE_DIR = Path(__file__).resolve().parent.parent / "packaging" / "release"
sys.path.insert(0, str(RELEASE_DIR))

import preview_build  # noqa: E402


def _git_repo(root: Path) -> tuple[str, str]:
    (root / "src/chemuson").mkdir(parents=True)
    (root / "src/chemuson/_version.py").write_text(
        '__version__ = "0.3.0-beta.1"\n', encoding="utf-8"
    )
    subprocess.run(["git", "init", "-q", str(root)], check=True)
    subprocess.run(["git", "-C", str(root), "add", "."], check=True)
    subprocess.run(
        [
            "git",
            "-C",
            str(root),
            "-c",
            "user.name=Preview Test",
            "-c",
            "user.email=preview@example.invalid",
            "commit",
            "-qm",
            "baseline",
        ],
        check=True,
    )
    sha = subprocess.run(
        ["git", "-C", str(root), "rev-parse", "HEAD"],
        check=True,
        capture_output=True,
        text=True,
    ).stdout.strip()
    return sha, "0.3.0-beta.1"


def test_preview_identity_requires_exact_sha_version_and_preparation_branch(tmp_path) -> None:
    sha, version = _git_repo(tmp_path)
    assert preview_build.validate_preview_identity(
        version=version,
        git_sha=sha,
        source_branch="release/v0.3.0-beta.1-prep",
        repo_root=tmp_path,
    ) == sha

    with pytest.raises(ValueError, match="differs from canonical"):
        preview_build.validate_preview_identity(
            version="0.3.0-beta.2",
            git_sha=sha,
            source_branch="release/v0.3.0-beta.1-prep",
            repo_root=tmp_path,
        )
    with pytest.raises(ValueError, match="does not match the expected Git SHA"):
        preview_build.validate_preview_identity(
            version=version,
            git_sha="b" * 40,
            source_branch="release/v0.3.0-beta.1-prep",
            repo_root=tmp_path,
        )
    with pytest.raises(ValueError, match="restricted to"):
        preview_build.validate_preview_identity(
            version=version,
            git_sha=sha,
            source_branch="main",
            repo_root=tmp_path,
        )


def test_preview_provenance_checksums_and_package_names_are_verified(tmp_path) -> None:
    sha, version = _git_repo(tmp_path)
    stage = tmp_path / "preview-assets"
    stage.mkdir()
    package = "Chemuson-preview-abcdef12-windows-portable.exe"
    (stage / package).write_bytes(b"preview package")

    record = preview_build.write_group(
        output_dir=stage,
        version=version,
        git_sha=sha,
        source_branch="release/v0.3.0-beta.1-prep",
        operating_system="Windows runner",
        files=[package],
        repo_root=tmp_path,
    )

    assert record["version"] == version
    assert record["build_type"] == "preview"
    assert record["git_sha"] == sha
    assert record["source_branch"] == "release/v0.3.0-beta.1-prep"
    assert record["build_status"] == "success"
    assert record["publication"] is False
    saved = json.loads((stage / "preview-provenance.json").read_text(encoding="utf-8"))
    assert saved == record
    checksums = (stage / "checksums.sha256").read_text(encoding="utf-8")
    assert f"{record['artifacts'][package]['sha256']}  {package}" in checksums
    assert "  preview-provenance.json" in checksums


def test_preview_manifest_rejects_release_like_or_missing_package(tmp_path) -> None:
    sha, version = _git_repo(tmp_path)
    stage = tmp_path / "preview-assets"
    stage.mkdir()
    official_name = "Chemuson-v0.3.0-beta.1-linux-x86_64.AppImage"
    (stage / official_name).write_bytes(b"not distinguishable")

    with pytest.raises(ValueError, match="distinct from releases"):
        preview_build.write_group(
            output_dir=stage,
            version=version,
            git_sha=sha,
            source_branch="release/v0.3.0-beta.1-prep",
            operating_system="Linux runner",
            files=[official_name],
            repo_root=tmp_path,
        )
    with pytest.raises(ValueError, match="missing or empty"):
        preview_build.write_group(
            output_dir=stage,
            version=version,
            git_sha=sha,
            source_branch="release/v0.3.0-beta.1-prep",
            operating_system="Linux runner",
            files=["Chemuson-preview-abcdef12-missing.flatpak"],
            repo_root=tmp_path,
        )


def test_preview_prepare_rejects_tag_events(tmp_path, monkeypatch) -> None:
    sha, _ = _git_repo(tmp_path)
    monkeypatch.setenv("GITHUB_REF_TYPE", "tag")
    monkeypatch.setenv("GITHUB_REF_NAME", "v0.3.0-beta.1")
    monkeypatch.setenv("GITHUB_SHA", sha)

    with pytest.raises(ValueError, match="never a tag"):
        preview_build.prepare_from_environment(tmp_path)
