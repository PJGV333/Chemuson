"""Release tag policy and fail-closed GitHub preflight tests."""

from __future__ import annotations

import importlib.util
import json
import os
import subprocess
import sys
import urllib.error
from pathlib import Path

import pytest

RELEASE_DIR = Path(__file__).resolve().parent.parent / "packaging" / "release"
sys.path.insert(0, str(RELEASE_DIR))

import release_policy  # noqa: E402
import validate_release_tag  # noqa: E402


def test_release_tag_policy_maps_only_supported_tags() -> None:
    stable = release_policy.parse_release_tag("v0.3.0")
    beta = release_policy.parse_release_tag("v0.3.0-beta.1")
    rc = release_policy.parse_release_tag("v0.3.0-rc.1")

    assert (stable.version, stable.channel, stable.prerelease) == ("0.3.0", "stable", False)
    assert (beta.version, beta.channel, beta.prerelease) == ("0.3.0-beta.1", "beta", True)
    assert (rc.version, rc.channel, rc.prerelease) == ("0.3.0-rc.1", "beta", True)


@pytest.mark.parametrize(
    "tag",
    [
        "0.3.0",
        "v0.3.0-dev",
        "v0.3.0-alpha.1",
        "v0.3.0-beta.0",
        "v0.3.0-beta.01",
        "v01.3.0",
        "v0.3.0+build.1",
        "v0.3.0-rc.1.2",
        "v0.3.0-unknown",
    ],
)
def test_release_tag_policy_rejects_unsupported_tags(tag: str) -> None:
    with pytest.raises(ValueError):
        release_policy.parse_release_tag(tag)


def test_channel_mismatch_is_rejected() -> None:
    with pytest.raises(ValueError, match="belongs to channel"):
        release_policy.validate_version_channel("0.3.0-beta.1", "stable")
    with pytest.raises(ValueError, match="belongs to channel"):
        release_policy.validate_version_channel("0.3.0", "beta")


def test_prepared_canonical_version_matches_appstream_and_dynamic_package_metadata() -> None:
    root = Path(__file__).resolve().parent.parent
    source_version = validate_release_tag.read_source_version(root / "src/chemuson/_version.py")
    appstream_version = validate_release_tag.read_appstream_version(
        root / "packaging/flatpak/io.github.PJGV333.Chemuson.metainfo.xml"
    )
    assert source_version == appstream_version
    pyproject = (root / "pyproject.toml").read_text(encoding="utf-8")
    assert 'dynamic = ["version"]' in pyproject
    assert 'version = {attr = "chemuson._version.__version__"}' in pyproject


def test_tag_canonical_version_and_appstream_must_match() -> None:
    result = release_policy.validate_release_metadata(
        tag="v0.3.0-beta.1",
        source_version="0.3.0-beta.1",
        appstream_version="0.3.0-beta.1",
    )
    assert result.channel == "beta"

    with pytest.raises(ValueError, match="canonical application version"):
        release_policy.validate_release_metadata(
            tag="v0.3.0-beta.1",
            source_version="0.3.0-dev",
            appstream_version="0.3.0-beta.1",
        )
    with pytest.raises(ValueError, match="AppStream version"):
        release_policy.validate_release_metadata(
            tag="v0.3.0-beta.1",
            source_version="0.3.0-beta.1",
            appstream_version="0.2.5",
        )


def test_tag_event_checkout_and_git_sha_must_match() -> None:
    sha = "a" * 40
    assert release_policy.validate_commit_identity(sha, sha.upper(), sha) == sha
    with pytest.raises(ValueError, match="do not match"):
        release_policy.validate_commit_identity(sha, "b" * 40)
    with pytest.raises(ValueError, match="full Git SHAs"):
        release_policy.validate_commit_identity("short", "short")


def test_release_absence_accepts_only_github_404(monkeypatch) -> None:
    class _Response:
        status = 404

        def __enter__(self):
            return self

        def __exit__(self, *_args):
            return False

        def getcode(self):
            return self.status

    monkeypatch.setattr(validate_release_tag.urllib.request, "urlopen", lambda *_a, **_k: _Response())
    validate_release_tag.assert_release_absent("PJGV333/Chemuson", "v0.3.0", "token")

    def _existing(*_args, **_kwargs):
        response = _Response()
        response.status = 200
        return response

    monkeypatch.setattr(validate_release_tag.urllib.request, "urlopen", _existing)
    with pytest.raises(ValueError, match="already exists"):
        validate_release_tag.assert_release_absent("PJGV333/Chemuson", "v0.3.0", "token")


def test_release_absence_fails_closed_on_api_errors(monkeypatch) -> None:
    def _forbidden(*_args, **_kwargs):
        raise urllib.error.HTTPError("https://api.github.com", 403, "forbidden", {}, None)

    monkeypatch.setattr(validate_release_tag.urllib.request, "urlopen", _forbidden)
    with pytest.raises(RuntimeError, match="failing closed"):
        validate_release_tag.assert_release_absent("PJGV333/Chemuson", "v0.3.0", "token")

    with pytest.raises(ValueError, match="GITHUB_TOKEN"):
        validate_release_tag.assert_release_absent("PJGV333/Chemuson", "v0.3.0", "")


def test_release_absence_rejects_invalid_repo() -> None:
    with pytest.raises(ValueError, match="GITHUB_REPOSITORY"):
        validate_release_tag.assert_release_absent("https://example.com", "v0.3.0", "token")


def _git(cwd: Path, *args: str) -> str:
    result = subprocess.run(
        ["git", *args],
        cwd=cwd,
        check=True,
        capture_output=True,
        text=True,
    )
    return result.stdout.strip()


def test_module_repo_root_matches_git_toplevel() -> None:
    """Regression for the beta.2 failure: REPO_ROOT must be the checkout root.

    The preflight resolved REPO_ROOT one level above the repository, so the
    official tag run executed `git rev-parse HEAD` outside any Git work tree
    and failed closed with exit status 128.
    """
    toplevel = _git(RELEASE_DIR, "rev-parse", "--show-toplevel")
    assert validate_release_tag.REPO_ROOT == Path(toplevel).resolve()
    assert validate_release_tag.REPO_ROOT / "packaging/release/validate_release_tag.py" == (
        Path(__file__).resolve().parent.parent / "packaging/release/validate_release_tag.py"
    )


def _make_rehearsal_repository(tmp_path: Path, version: str) -> tuple[Path, str]:
    """Copy the real preflight scripts into a throwaway annotated-tag checkout."""
    scripts = tmp_path / "packaging/release"
    scripts.mkdir(parents=True)
    for name in ("validate_release_tag.py", "release_policy.py"):
        (scripts / name).write_text((RELEASE_DIR / name).read_text(encoding="utf-8"), encoding="utf-8")
    version_file = tmp_path / "src/chemuson/_version.py"
    version_file.parent.mkdir(parents=True)
    version_file.write_text(f'__version__ = "{version}"\n', encoding="utf-8")
    metainfo = tmp_path / "packaging/flatpak/io.github.PJGV333.Chemuson.metainfo.xml"
    metainfo.parent.mkdir(parents=True)
    metainfo.write_text(
        f'<component><releases><release version="{version}" date="2026-01-01" /></releases></component>\n',
        encoding="utf-8",
    )
    identity = ["-c", "user.name=Release Regression", "-c", "user.email=release@example.invalid"]
    _git(tmp_path, "init", "-q", "--initial-branch=candidate")
    _git(tmp_path, *identity, "add", "-A")
    _git(tmp_path, *identity, "commit", "-q", "-m", "candidate commit")
    sha = _git(tmp_path, "rev-parse", "HEAD")
    _git(tmp_path, *identity, "tag", "-a", f"v{version}", "-m", f"annotated rehearsal tag v{version}")
    return tmp_path, sha


def _rehearsal_environment(root: Path, version: str, sha: str, *, event_sha: str | None = None) -> dict[str, str]:
    event_path = root / "tag-event.json"
    event_path.write_text(
        json.dumps({"ref": f"refs/tags/v{version}", "ref_type": "tag", "deleted": False}),
        encoding="utf-8",
    )
    env = {
        "GITHUB_REF_NAME": f"v{version}",
        "GITHUB_REF_TYPE": "tag",
        "GITHUB_SHA": event_sha or sha,
        "GITHUB_EVENT_PATH": str(event_path),
        "GITHUB_OUTPUT": str(root / "preflight-output.txt"),
        "GITHUB_REPOSITORY": "PJGV333/Chemuson",
        "PATH": os.environ.get("PATH", ""),
        "GIT_CONFIG_GLOBAL": "/dev/null",
        "GIT_CONFIG_SYSTEM": "/dev/null",
    }
    return env


def _load_rehearsal_module(root: Path):
    script = root / "packaging/release/validate_release_tag.py"
    spec = importlib.util.spec_from_file_location("rehearsal_validate_release_tag", script)
    assert spec is not None and spec.loader is not None
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def test_rehearsal_module_repo_root_resolves_to_rehearsal_checkout(tmp_path: Path) -> None:
    root, _sha = _make_rehearsal_repository(tmp_path, "0.3.0-beta.3")
    module = _load_rehearsal_module(root)
    assert module.REPO_ROOT == root.resolve()
    assert module.REPO_ROOT == Path(_git(root, "rev-parse", "--show-toplevel")).resolve()


def test_verify_tag_checkout_runs_real_git_logic_on_rehearsal_checkout(tmp_path: Path) -> None:
    """Exercise verify_tag_checkout itself against a real annotated tag checkout."""
    root, sha = _make_rehearsal_repository(tmp_path, "0.3.0-beta.3")
    module = _load_rehearsal_module(root)
    assert module.verify_tag_checkout(
        repo_root=module.REPO_ROOT,
        tag="v0.3.0-beta.3",
        event_sha=sha,
        ref_type="tag",
    ) == sha
    with pytest.raises(ValueError, match="do not match"):
        module.verify_tag_checkout(
            repo_root=module.REPO_ROOT,
            tag="v0.3.0-beta.3",
            event_sha="f" * 40,
            ref_type="tag",
        )
    with pytest.raises(ValueError, match="tag ref"):
        module.verify_tag_checkout(
            repo_root=module.REPO_ROOT,
            tag="v0.3.0-beta.3",
            event_sha=sha,
            ref_type="branch",
        )
    with pytest.raises(subprocess.CalledProcessError):
        module.verify_tag_checkout(
            repo_root=module.REPO_ROOT,
            tag="v0.3.0-beta.4",
            event_sha=sha,
            ref_type="tag",
        )


def test_full_preflight_main_succeeds_on_simulated_tag_run(tmp_path: Path, monkeypatch) -> None:
    """Run the real main() wiring: path resolution, git identity, metadata, outputs."""
    root, sha = _make_rehearsal_repository(tmp_path, "0.3.0-beta.3")
    module = _load_rehearsal_module(root)
    looked_up: list[tuple[str, str, str]] = []
    monkeypatch.setattr(
        module,
        "assert_release_absent",
        lambda repository, tag, token, **_kwargs: looked_up.append((repository, tag, token)),
    )
    env = _rehearsal_environment(root, "0.3.0-beta.3", sha)
    env["RELEASE_PREFLIGHT_TOKEN"] = "read-only-token"
    for key, value in env.items():
        monkeypatch.setenv(key, value)
    monkeypatch.chdir(root)
    module.main()
    outputs = (root / "preflight-output.txt").read_text(encoding="utf-8")
    assert outputs.splitlines() == [
        "tag=v0.3.0-beta.3",
        "version=0.3.0-beta.3",
        "channel=beta",
        f"sha={sha}",
    ]
    assert looked_up == [("PJGV333/Chemuson", "v0.3.0-beta.3", "read-only-token")]


def test_preflight_subprocess_fails_closed_offline_after_all_local_gates(tmp_path: Path) -> None:
    """Offline end-to-end: every local gate must pass before the network lookup stops it."""
    root, sha = _make_rehearsal_repository(tmp_path, "0.3.0-beta.3")
    env = os.environ.copy()
    env.update(_rehearsal_environment(root, "0.3.0-beta.3", sha))
    result = subprocess.run(
        [sys.executable, "packaging/release/validate_release_tag.py"],
        cwd=root,
        env=env,
        capture_output=True,
        text=True,
    )
    assert result.returncode != 0
    combined = result.stdout + result.stderr
    assert "GITHUB_TOKEN is required" in combined
    for leak in ("rev-parse", "128", "canonical application version", "AppStream version"):
        assert leak not in combined


def test_preflight_subprocess_fails_closed_on_sha_mismatch_without_token(tmp_path: Path) -> None:
    """The SHA immutability gate runs before the network and cannot be bypassed."""
    root, sha = _make_rehearsal_repository(tmp_path, "0.3.0-beta.3")
    env = os.environ.copy()
    env.update(_rehearsal_environment(root, "0.3.0-beta.3", sha, event_sha="f" * 40))
    result = subprocess.run(
        [sys.executable, "packaging/release/validate_release_tag.py"],
        cwd=root,
        env=env,
        capture_output=True,
        text=True,
    )
    assert result.returncode != 0
    assert "do not match" in result.stdout + result.stderr
    assert "GITHUB_TOKEN" not in result.stdout + result.stderr
