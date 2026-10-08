"""Release tag policy and fail-closed GitHub preflight tests."""

from __future__ import annotations

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
