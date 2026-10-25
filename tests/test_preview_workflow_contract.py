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
        "python packaging/release/generate_channel_manifest.py",
        "python packaging/release/generate_flatpak_pages_index.py",
        "python packaging/release/generate_flatpak_remote_files.py",
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
    appimage_validator = (ROOT / "packaging/release/validate_appimage.py").read_text(encoding="utf-8")
    assert '"${APPIMAGETOOL}" "${APPDIR}" "${APPIMAGE_PATH}"' in appimage_builder
    assert "--updateinformation \"${UPDATE_INFO}\"" in appimage_builder
    assert 'if build_type == "preview":' in appimage_validator
    preview_validation = appimage_validator.split('if build_type == "preview":', 1)[1].split("else:", 1)[0]
    assert "embedded_update_info" in preview_validation
    assert "updateinfo_path" in preview_validation and "update_json_path" in preview_validation
    flatpak_builder = (ROOT / "packaging/linux/build_flatpak.sh").read_text(encoding="utf-8")
    assert "Preview Flatpak builds cannot use public remote URLs or signing credentials." in flatpak_builder
    assert 'BRANCH="preview-${SOURCE_SHA:0:8}"' in flatpak_builder
    assert '"preview" "$PREVIEW_SHA" "$PREVIEW_BRANCH"' in text


def test_all_build_jobs_use_the_prevalidated_event_sha_and_version() -> None:
    text, workflow = _workflow()
    jobs = workflow["jobs"]
    assert jobs["prepare"]["steps"][0]["with"]["ref"] == "${{ github.sha }}"
    for name in (
        "preflight_rehearsal",
        "release_gate_rehearsal",
        "build_windows",
        "build_linux_appimage",
        "build_linux_flatpak",
    ):
        job = jobs[name]
        assert job["needs"] == "prepare"
        checkout = job["steps"][0]
        assert checkout["with"]["ref"] == "${{ needs.prepare.outputs.sha }}"
        assert "preview_build.py verify" in text
    assert "preview_build.py prepare" in text
    assert "needs.prepare.outputs.version" in text
    assert "needs.prepare.outputs.branch" in text


def _normalize(run: str) -> str:
    return " ".join(str(run).split())


def test_preflight_rehearsal_replays_the_official_validator_fail_closed() -> None:
    text, workflow = _workflow()
    job = workflow["jobs"]["preflight_rehearsal"]
    assert job["runs-on"] == "ubuntu-latest"
    assert job["permissions"] == {"contents": "read"}
    steps = job["steps"]
    checkout = steps[0]
    assert checkout["with"]["fetch-depth"] == "0"
    python_setup = next(step for step in steps if step.get("name") == "Setup Python")
    assert python_setup["with"]["python-version"] == "3.11"

    run_bodies = [_normalize(step.get("run", "")) for step in steps]
    joined = "\n".join(run_bodies)

    # The rehearsal must run the exact official validator, never a reimplementation.
    assert "python packaging/release/validate_release_tag.py" in joined
    # It must first prove the remote tag is absent, then simulate a tag event fully.
    assert 'git ls-remote --tags origin "refs/tags/v${PREVIEW_VERSION}"' in joined
    assert '"ref_type":"tag","deleted":false' in joined
    assert 'export GITHUB_REF_NAME="${TAG}"' in joined
    assert 'export GITHUB_REF_TYPE="tag"' in joined
    assert 'export GITHUB_SHA="${PREVIEW_SHA}"' in joined
    assert 'export GITHUB_EVENT_PATH="${EVENT_PATH}"' in joined
    assert 'export GITHUB_OUTPUT="${OUTPUT_FILE}"' in joined
    assert "validate_release_tag.REPO_ROOT" in joined
    assert '"git", "rev-parse", "--show-toplevel"' in joined
    # It must use an authorized read-only token and enforce exit code and outputs.
    preflight_step = next(
        step for step in steps if "validate_release_tag.py" in step.get("run", "")
    )
    assert preflight_step["env"]["RELEASE_PREFLIGHT_TOKEN"] == "${{ github.token }}"
    assert "set -euo pipefail" in preflight_step["run"]
    verify = next(step for step in steps if "match the candidate identity" in step.get("name", ""))
    assert f'grep -qxF "tag=v' in verify["run"]
    assert 'grep -qxF "version=${PREVIEW_VERSION}"' in verify["run"]
    assert 'grep -qxF "channel=' in verify["run"]
    assert 'grep -qxF "sha=${PREVIEW_SHA}"' in verify["run"]
    # Local-only annotated tag: created for the exact SHA and always deleted afterwards.
    assert 'git tag -a "v${PREVIEW_VERSION}"' in joined
    cleanup = steps[-1]
    assert cleanup["if"] == "always()"
    assert 'git tag -d "v${PREVIEW_VERSION}"' in cleanup["run"]
    assert 'test -z "$(git tag -l "v${PREVIEW_VERSION}")"' in cleanup["run"]
    assert "git push" not in text


def test_release_gate_rehearsal_mirrors_official_gate_commands_exactly() -> None:
    _text, preview = _workflow()
    release_text = RELEASE_PATH.read_text(encoding="utf-8")
    release = yaml.load(release_text, Loader=yaml.BaseLoader)
    gate_steps = {
        step.get("name"): step for step in release["jobs"]["release_gate"]["steps"]
    }
    rehearsal_steps = {
        step.get("name"): step
        for step in preview["jobs"]["release_gate_rehearsal"]["steps"]
    }
    mirrored = (
        "Setup Python",
        "Install Qt runtime and test dependencies",
        "Validate exact baseline exception catalog",
        "Run Python ChemName acceptance",
        "Compile source, tests and packaging tools",
        "Run scoped release Ruff checks",
        "Run architecture tests on validated tag SHA",
        "Run bounded release, version, packaging and workflow contracts",
    )
    for name in mirrored:
        assert name in rehearsal_steps, f"gate rehearsal is missing step: {name}"
        assert _normalize(rehearsal_steps[name].get("run", "")) == _normalize(
            gate_steps[name].get("run", "")
        ), f"gate rehearsal step differs from official gate: {name}"
    assert (
        rehearsal_steps["Setup Python"]["with"]["python-version"]
        == gate_steps["Setup Python"]["with"]["python-version"]
    )
    assert preview["jobs"]["release_gate_rehearsal"]["env"]["QT_QPA_PLATFORM"] == "offscreen"
    assert preview["jobs"]["release_gate_rehearsal"]["permissions"] == {"contents": "read"}
    assert "environment" not in preview["jobs"]["release_gate_rehearsal"]


def test_preview_report_gates_every_rehearsal_and_build() -> None:
    _text, workflow = _workflow()
    report = workflow["jobs"]["preview_report"]
    assert report["needs"] == [
        "prepare",
        "preflight_rehearsal",
        "release_gate_rehearsal",
        "build_windows",
        "build_linux_appimage",
        "build_linux_flatpak",
    ]
    step_run = report["steps"][0]["run"]
    for key in ("release_preflight_rehearsal", "release_gate_bounded_rehearsal"):
        assert f'"{key}"' in step_run
    assert 'all(value == "success" for value in statuses.values())' in step_run


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


def test_preview_and_release_appimage_builds_preflight_qt_xcb_loader() -> None:
    preview_text, preview = _workflow()
    release_text = RELEASE_PATH.read_text(encoding="utf-8")
    release = yaml.load(release_text, Loader=yaml.BaseLoader)
    runtime_packages = (
        "libxcb-keysyms1",
        "libxcb-icccm4",
        "libxcb-shape0",
        "xauth",
        "xvfb",
    )
    for workflow_text, job_name, install_name in (
        (preview_text, "build_linux_appimage", "Install Qt runtime libraries"),
        (release_text, "build_linux", "Install AppImage and Qt runtime tools"),
    ):
        workflow = preview if job_name == "build_linux_appimage" else release
        steps = workflow["jobs"][job_name]["steps"]
        install_step = next(step for step in steps if step.get("name") == install_name)
        assert all(package in install_step["run"] for package in runtime_packages)
        diagnostic_index = next(
            index for index, step in enumerate(steps)
            if step.get("name") == "Diagnose Qt xcb plugin loader before AppImage graphics smoke"
        )
        validation_index = next(
            index for index, step in enumerate(steps)
            if "validate_appimage.py" in step.get("run", "")
        )
        assert diagnostic_index < validation_index
        assert "diagnose_qt_xcb.py" in steps[diagnostic_index]["run"]
        assert "validate_appimage.py" in workflow_text

    diagnostic = (ROOT / "packaging/release/diagnose_qt_xcb.py").read_text(encoding="utf-8")
    for required in ("readelf", "ldd", "QT_DEBUG_PLUGINS", "LD_DEBUG", "xvfb-run", "platform=xcb"):
        assert required in diagnostic


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
    assert "tools/chemname_acceptance.py" in text
    assert "validate_packaged_chemname.py --executable \"dist/Chemuson.exe\"" in text
    assert "validate_packaged_chemname.py --executable dist/Chemuson" in text
    assert "validate_packaged_chemname.py --executable \"dist/Chemuson.exe\"" in release_text
    assert "validate_packaged_chemname.py --executable dist/Chemuson" in release_text
    assert "validate_packaged_rdkit_worker.py --executable \"dist/Chemuson.exe\"" in text
    assert "validate_packaged_rdkit_worker.py --executable dist/Chemuson" in text
    assert "validate_packaged_rdkit_worker.py --executable \"dist/Chemuson.exe\"" in release_text
    assert "validate_packaged_rdkit_worker.py --executable dist/Chemuson" in release_text
    workflow = _workflow()[1]
    windows_build = workflow["jobs"]["build_windows"]
    validation_step = next(
        step for step in windows_build["steps"]
        if "validate_packaged_rdkit_worker.py" in step.get("run", "")
    )
    assert validation_step["run"].count("if ($LASTEXITCODE -ne 0)") == 3
    assert "Packaged RDKit validation failed" in validation_step["run"]
    assert "Packaged ChemName validation failed" in validation_step["run"]
