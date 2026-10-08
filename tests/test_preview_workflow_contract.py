"""Static security contract for the Actions-only preview workflow."""

from __future__ import annotations

from pathlib import Path

import yaml

ROOT = Path(__file__).resolve().parent.parent
PREVIEW_PATH = ROOT / ".github/workflows/build-preview.yml"
RELEASE_PATH = ROOT / ".github/workflows/release.yml"
EXPECTED_ARTIFACTS = {
    "chemuson-preview-windows-portable",
    "chemuson-preview-windows-installer",
    "chemuson-preview-linux-appimage",
    "chemuson-preview-linux-flatpak",
}


def _workflow() -> tuple[str, dict]:
    text = PREVIEW_PATH.read_text(encoding="utf-8")
    parsed = yaml.load(text, Loader=yaml.BaseLoader)
    assert isinstance(parsed, dict)
    return text, parsed


def test_preview_triggers_only_manually_or_on_prep_branches() -> None:
    _text, workflow = _workflow()
    triggers = workflow["on"]
    assert "workflow_dispatch" in triggers
    assert triggers["push"]["branches"] == ["release/**-prep"]


def test_preview_permissions_are_read_only_and_no_secrets_are_injected() -> None:
    text, workflow = _workflow()
    assert workflow["permissions"] == {"contents": "read"}
    for job in workflow["jobs"].values():
        if "permissions" in job:
            assert job["permissions"] == {"contents": "read"}
    assert "contents: write" not in text
    assert "secrets." not in text
    assert "GITHUB_TOKEN" not in text
    assert "CHEMUSON_SIGN_KEY" not in text
    assert "FLATPAK_GPG" not in text
    assert "WINDOWS_CODESIGN" not in text


def test_preview_has_no_release_tag_gh_pages_or_public_channel_operations() -> None:
    text, _workflow_data = _workflow()
    forbidden = (
        "softprops/action-gh-release",
        "gh release",
        "gh-pages",
        "git push",
        "git tag",
        "refs/tags/",
        "generate_channel_manifest.py",
        "generate_flatpak_pages_index.py",
        "generate_flatpak_remote_files.py",
        "stable.json",
        "beta.json",
        "flatpak/stable/",
        "flatpak/beta/",
        "publish_flatpak_remote",
        "CHEMUSON_FLATPAK_REPO_URL",
        "releases/download/",
    )
    lowered = text.lower()
    assert not [token for token in forbidden if token.lower() in lowered]

    appimage_builder = (ROOT / "packaging/linux/build_appimage.sh").read_text(encoding="utf-8")
    preview_appimage_path = appimage_builder.split('if [[ "${BUILD_TYPE}" == "preview" ]]', 1)[1]
    assert "no public updater metadata" in appimage_builder.lower()
    assert "exit 0" in preview_appimage_path
    assert appimage_builder.index('if [[ "${BUILD_TYPE}" == "preview" ]]') < appimage_builder.index("# Metadata AppImageUpdate")
    flatpak_builder = (ROOT / "packaging/linux/build_flatpak.sh").read_text(encoding="utf-8")
    assert "Preview Flatpak builds cannot use public remote URLs or signing credentials." in flatpak_builder
    assert 'BRANCH="preview-${SOURCE_SHA:0:8}"' in flatpak_builder
    assert '"preview" "$PREVIEW_SHA" "$PREVIEW_BRANCH"' in text


def test_all_build_jobs_use_the_prevalidated_event_sha_and_version() -> None:
    text, workflow = _workflow()
    jobs = workflow["jobs"]
    assert jobs["prepare"]["steps"][0]["with"]["ref"] == "${{ github.sha }}"
    for name in ("build_windows", "build_linux_appimage", "build_linux_flatpak"):
        job = jobs[name]
        assert job["needs"] == "prepare"
        checkout = job["steps"][0]
        assert checkout["with"]["ref"] == "${{ needs.prepare.outputs.sha }}"
        assert "preview_build.py verify" in text
    assert "preview_build.py prepare" in text
    assert "needs.prepare.outputs.version" in text
    assert "needs.prepare.outputs.branch" in text


def test_four_preview_artifact_groups_require_existing_nonempty_packages() -> None:
    text, workflow = _workflow()
    uploads = {
        step["with"]["name"]: step["with"]
        for job in workflow["jobs"].values()
        for step in job["steps"]
        if step.get("uses", "").startswith("actions/upload-artifact@")
    }
    assert EXPECTED_ARTIFACTS <= set(uploads)
    for name in EXPECTED_ARTIFACTS:
        assert uploads[name]["if-no-files-found"] == "error"
        assert uploads[name]["path"].endswith("/*")
    assert "test -s" in text
    assert "Test-Path" in text
    assert "preview_build.py manifest" in text
    assert "--files" in text


def test_preview_metadata_uses_canonical_version_and_official_workflow_stays_separate() -> None:
    text, _workflow_data = _workflow()
    helper = (ROOT / "packaging/release/preview_build.py").read_text(encoding="utf-8")
    assert "read_canonical_version" in helper
    assert '"build_type": "preview"' in helper
    assert '"publication": False' in helper
    assert "build_status" in helper
    assert "Chemuson-v$env:PREVIEW_VERSION-preview-$env:PREVIEW_SHORT_SHA" in text
    assert "Chemuson-v${PREVIEW_VERSION}-preview-${PREVIEW_SHORT_SHA}" in text

    release_text = RELEASE_PATH.read_text(encoding="utf-8")
    assert "build-preview.yml" not in release_text
    assert "build_windows:" in release_text
    assert "build_linux:" in release_text
    assert "build_flatpak:" in release_text
