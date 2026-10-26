"""Focused CI contract for the Flatpak host-toolchain preparation.

The beta.4 blocker was not a ChemUSON regression: the bare apt step that
installs the Flatpak toolchain downloads ~99 MB and, when the assigned
`azure.archive.ubuntu.com` node is degraded, it can take 37-44 minutes and
consume the whole job budget. These tests pin the repaired contract: a shared,
bounded, cache-backed preparation that keeps apt as the only package source,
never weakens signature verification, and leaves the scientific packaging
inputs untouched.
"""

from __future__ import annotations

import re
import shutil
import subprocess
from pathlib import Path

import yaml

ROOT = Path(__file__).resolve().parent.parent
PREVIEW_PATH = ROOT / ".github/workflows/build-preview.yml"
RELEASE_PATH = ROOT / ".github/workflows/release.yml"
INSTALLER = ROOT / "packaging/linux/install_flatpak_tools.sh"
LOCK_FILE = ROOT / "packaging/flatpak/toolchain.lock"
BUILDER = ROOT / "packaging/linux/build_flatpak.sh"
MANIFEST = ROOT / "packaging/flatpak/io.github.PJGV333.Chemuson.yml"

OFFICIAL_APT_HOSTS = {
    "azure.archive.ubuntu.com",
    "archive.ubuntu.com",
    "us.archive.ubuntu.com",
    "security.ubuntu.com",
    "ports.ubuntu.com",
}
CACHE_ACTION = "actions/cache@55cc8345863c7cc4c66a329aec7e433d2d1c52a9"
UNSAFE_APT_TOKENS = (
    "--allow-unauthenticated",
    "AllowInsecureRepositories",
    "AllowDowngradeToInsecureRepositories",
    "allow-change-held-packages",
    "add-apt-repository",
    "trusted=yes",
    "apt-transport-",
    "keyserver",
    "--allow-insecure",
)


def _yaml(path: Path) -> dict:
    parsed = yaml.load(path.read_text(encoding="utf-8"), Loader=yaml.BaseLoader)
    assert isinstance(parsed, dict)
    return parsed


def _lock_values() -> dict[str, str]:
    values: dict[str, str] = {}
    for line in LOCK_FILE.read_text(encoding="utf-8").splitlines():
        if "=" in line and not line.startswith("#"):
            key, value = line.split("=", 1)
            values[key.strip()] = value.strip()
    return values


def _steps(job: dict) -> list[dict]:
    return job["steps"]


def _step_index(steps: list[dict], needle: str) -> int:
    for index, step in enumerate(steps):
        body = f"{step.get('name', '')} {step.get('uses', '')} {step.get('run', '')}"
        if needle in body:
            return index
    raise AssertionError(f"step containing {needle!r} not found")


def test_installer_exists_and_is_valid_bash() -> None:
    assert INSTALLER.is_file()
    text = INSTALLER.read_text(encoding="utf-8")
    assert text.startswith("#!/usr/bin/env bash")
    assert "set -euo pipefail" in text
    bash = shutil.which("bash")
    if bash is None:  # pragma: no cover - runners and dev shells always ship bash
        return
    result = subprocess.run(
        [bash, "-n", str(INSTALLER)], capture_output=True, text=True, timeout=60
    )
    assert result.returncode == 0, result.stderr


def test_installer_keeps_apt_as_the_only_signed_package_source() -> None:
    text = INSTALLER.read_text(encoding="utf-8")
    for token in UNSAFE_APT_TOKENS:
        assert token not in text, f"weakened apt acquisition: {token}"
    assert "sudo apt-get" not in text  # privilege comes from the SUDO variable
    assert "${SUDO} apt-get" in text
    for required in (
        "apt-get update",
        "--download-only",
        "--no-download",
        "dpkg-query",
        "toolchain.lock",
    ):
        assert required in text


def test_installer_bounds_every_network_phase() -> None:
    text = INSTALLER.read_text(encoding="utf-8")
    for required in (
        "Acquire::Retries",
        "Acquire::http::Timeout",
        "Acquire::https::Timeout",
        'timeout "${UPDATE_TIMEOUT}"',
        'timeout "${DOWNLOAD_TIMEOUT}"',
        'timeout "${INSTALL_TIMEOUT}"',
        "--max-time",
    ):
        assert required in text
    lock = _lock_values()
    for key in (
        "apt_mirror_probe_seconds",
        "apt_update_timeout_seconds",
        "apt_download_timeout_seconds",
        "apt_install_timeout_seconds",
        "apt_acquire_retries",
        "apt_acquire_timeout_seconds",
    ):
        assert lock[key].isdigit(), f"missing numeric bound: {key}"
        assert int(lock[key]) > 0


def test_installer_only_reaches_official_ubuntu_archives() -> None:
    lock = _lock_values()
    allowed = set(lock["allowed_apt_hosts"].split())
    assert allowed and allowed <= OFFICIAL_APT_HOSTS
    text = INSTALLER.read_text(encoding="utf-8")
    hosts = set(re.findall(r"[a-z0-9.-]+\.(?:ubuntu|flathub)\.com", text))
    assert hosts <= OFFICIAL_APT_HOSTS | {"dl.flathub.org"}, hosts


def test_installer_verifies_the_origin_files_it_changes() -> None:
    """An origin rewrite that silently fails would leave the throttled mirror.

    The runner points apt at its own `/etc/apt/apt-mirrors.txt`. If that file is
    not writable the script must say so and pin the main stanza of the deb822
    sources file to a probed official origin, while leaving the security stanza
    and every Signed-By entry untouched.
    """
    text = INSTALLER.read_text(encoding="utf-8")
    assert 'head -n 1 "${MIRROR_LIST}"' in text
    assert "sources-direct" in text
    assert "if (direct && $2 ~ /^mirror\\+file:/)" in text
    assert 'is not writable; refusing to install from an unverified origin.' in text
    assert 'if ($2 ~ /security\\.ubuntu\\.com/ && host != "security.ubuntu.com") { print; next }' in text
    assert "Signed-By" in text


def test_installer_fails_closed_and_reports_timings() -> None:
    text = INSTALLER.read_text(encoding="utf-8")
    for required in (
        "no official Ubuntu archive answered",
        "could not read the version of",
        "is older than the locked minimum",
        "not registered as installed",
        "Flatpak toolchain preparation timings",
        "GITHUB_STEP_SUMMARY",
    ):
        assert required in text
    assert re.search(r"exit 3", text) and re.search(r"exit 4", text)


def test_lock_file_pins_tools_runtime_and_job_image() -> None:
    lock = _lock_values()
    assert lock["runner_image"] == "ubuntu-24.04"
    for pkg in ("flatpak", "flatpak-builder", "ostree", "gnupg"):
        assert lock[f"min_{pkg.replace('-', '_')}" if f"min_{pkg.replace('-', '_')}" in lock else f"min_{pkg}"]
    manifest_text = MANIFEST.read_text(encoding="utf-8")
    assert lock["manifest"] == "packaging/flatpak/io.github.PJGV333.Chemuson.yml"
    assert f"runtime: {lock['runtime']}" in manifest_text
    assert f'runtime-version: "{lock["runtime_version"]}"' in manifest_text
    assert f"sdk: {lock['sdk']}" in manifest_text
    assert lock["runtime_repo"] in BUILDER.read_text(encoding="utf-8")
    assert "flatpak flatpak-builder ostree" == lock["apt_packages"]
    assert "gnupg" in lock["apt_release_packages"]


def _install_step(job: dict) -> dict:
    steps = _steps(job)
    return steps[_step_index(steps, "install_flatpak_tools.sh")]


def test_preview_flatpak_job_uses_bounded_cached_preparation() -> None:
    workflow = _yaml(PREVIEW_PATH)
    job = workflow["jobs"]["build_linux_flatpak"]
    assert job["runs-on"] == "ubuntu-24.04"
    assert job["timeout-minutes"] == "45"
    assert job["permissions"] == {"contents": "read"}
    steps = _steps(job)
    assert steps[0]["uses"].startswith("actions/checkout@")
    install = _install_step(job)
    install_index = _step_index(steps, "install_flatpak_tools.sh")
    cache = _step_index(steps, CACHE_ACTION)
    assert 0 < cache < install_index, "the package cache must be restored before installing"
    assert install["env"]["CHEMUSON_FLATPAK_APT_PACKAGES"] == "flatpak flatpak-builder ostree"
    assert install["env"]["CHEMUSON_FLATPAK_APT_CACHE_DIR"] == "${{ runner.temp }}/chemuson-flatpak-apt"
    cache_with = steps[cache]["with"]
    assert cache_with["path"] == install["env"]["CHEMUSON_FLATPAK_APT_CACHE_DIR"]
    assert "hashFiles('packaging/flatpak/toolchain.lock')" in cache_with["key"]
    assert "restore-keys" in cache_with
    assert install["run"].strip() == "bash packaging/linux/install_flatpak_tools.sh"
    assert "sudo apt-get install -y flatpak" not in PREVIEW_PATH.read_text(encoding="utf-8")


def test_release_flatpak_jobs_use_the_same_preparation() -> None:
    workflow = _yaml(RELEASE_PATH)
    build = workflow["jobs"]["build_flatpak"]
    publish = workflow["jobs"]["publish_flatpak_remote"]
    assert build["runs-on"] == "ubuntu-24.04"
    assert build["timeout-minutes"] == "45"
    assert build["needs"] == "release_gate"
    assert publish["runs-on"] == "ubuntu-24.04"
    assert publish["timeout-minutes"] == "25"
    for job, packages in (
        (build, "flatpak flatpak-builder gnupg ostree"),
        (publish, "ostree gnupg"),
    ):
        install = _install_step(job)
        assert install["env"]["CHEMUSON_FLATPAK_APT_PACKAGES"] == packages
        steps = _steps(job)
        cache = _step_index(steps, CACHE_ACTION)
        assert 0 < cache < _step_index(steps, "install_flatpak_tools.sh")
        assert steps[cache]["with"]["path"] == install["env"]["CHEMUSON_FLATPAK_APT_CACHE_DIR"]
    text = RELEASE_PATH.read_text(encoding="utf-8")
    assert "sudo apt-get install -y flatpak" not in text
    assert "sudo apt-get install -y ostree" not in text


def test_job_level_env_never_uses_the_runner_context() -> None:
    """GitHub rejects `${{ runner.* }}` in job-level env; it is a step context.

    A rejected workflow file produces a run with zero jobs and the file path as
    its name, which silently prevents any package from being built. The same
    parser rejects any expression inside a step `uses:` reference.
    """
    step_only_contexts = ("runner.", "steps.", "job.", "env.", "inputs.")
    for path in (PREVIEW_PATH, RELEASE_PATH, ROOT / ".github/workflows/test.yml"):
        workflow = _yaml(path)
        for job_name, job in workflow["jobs"].items():
            env = job.get("env")
            if isinstance(env, dict):
                for key, value in env.items():
                    for context in step_only_contexts:
                        assert context not in str(value), f"{path.name}:{job_name}.env.{key}"
            for step in _steps(job):
                reference = step.get("uses")
                if reference is not None:
                    assert "${{" not in str(reference), f"{path.name}:{job_name} uses: {reference}"


def test_release_safety_controls_stay_reachable() -> None:
    text = RELEASE_PATH.read_text(encoding="utf-8")
    for required in (
        "validate_release_tag.py",
        "generate_checksums.py",
        "validate_flatpak_remote_artifacts.py",
        "Validate local Flatpak commit metadata",
        "Verify GitHub Pages publication",
        "overwrite_files: false",
        "fail_on_unmatched_files: true",
    ):
        assert required in text
    preview_text = PREVIEW_PATH.read_text(encoding="utf-8")
    assert "secrets." not in preview_text
    assert "GITHUB_TOKEN" not in preview_text
    assert "git push" not in preview_text


def test_builder_bounds_flatpak_phases_and_reports_timings() -> None:
    text = BUILDER.read_text(encoding="utf-8")
    for required in (
        'BUILDER_TIMEOUT_SECONDS="${CHEMUSON_FLATPAK_BUILDER_TIMEOUT_SECONDS:-1500}"',
        'REPO_UPDATE_TIMEOUT_SECONDS="${CHEMUSON_FLATPAK_REPO_UPDATE_TIMEOUT_SECONDS:-300}"',
        'BUNDLE_TIMEOUT_SECONDS="${CHEMUSON_FLATPAK_BUNDLE_TIMEOUT_SECONDS:-300}"',
        'timeout "${BUILDER_TIMEOUT_SECONDS}" flatpak-builder',
        'timeout "${REPO_UPDATE_TIMEOUT_SECONDS}" flatpak build-update-repo',
        'timeout "${BUNDLE_TIMEOUT_SECONDS}" flatpak build-bundle',
        '-eq 124',
        "wall-clock limit",
        "Flatpak packaging timings",
    ):
        assert required in text
    assert 'BRANCH="preview-${SOURCE_SHA:0:8}"' in text
    assert "Preview Flatpak builds cannot use public remote URLs or signing credentials." in text


def test_scientific_packaging_inputs_remain_pinned() -> None:
    manifest = MANIFEST.read_text(encoding="utf-8")
    for module in (
        "python3-pyqt6-qt6",
        "python3-pyqt6-sip",
        "python3-pyqt6",
        "python3-numpy",
        "python3-rdkit",
        "chemuson",
    ):
        assert f"name: {module}" in manifest
    sha_lines = re.findall(r"sha256:\s[0-9a-f]{64}", manifest)
    assert len(sha_lines) == 5, sha_lines
    assert re.findall(r"url: https://files\.pythonhosted\.org/", manifest)
    for required in (
        "--share=network",
        "--filesystem=home/.chemuson:create",
        "--filesystem=home/.config/Chemuson:create",
        "io.github.PJGV333.Chemuson.metainfo.xml",
        "io.github.PJGV333.Chemuson.svg",
        "chemuson.chemname.packaged_smoke",
    ):
        assert required in manifest
