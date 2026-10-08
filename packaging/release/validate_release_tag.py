"""Fail-closed GitHub tag preflight for the release workflow."""

from __future__ import annotations

import json
import os
import re
import subprocess
import sys
import urllib.error
import urllib.parse
import urllib.request
import xml.etree.ElementTree as ET
from pathlib import Path

RELEASE_DIR = Path(__file__).resolve().parent
REPO_ROOT = RELEASE_DIR.parents[2]
if str(RELEASE_DIR) not in sys.path:
    sys.path.insert(0, str(RELEASE_DIR))

from release_policy import (  # noqa: E402
    parse_release_tag,
    validate_commit_identity,
    validate_release_metadata,
)

_VERSION_RE = re.compile(r'^\s*__version__\s*=\s*"([^"]+)"\s*$')


def read_source_version(path: Path) -> str:
    for line in path.read_text(encoding="utf-8").splitlines():
        match = _VERSION_RE.fullmatch(line)
        if match:
            return match.group(1)
    raise ValueError(f"No canonical __version__ found in {path}.")


def read_appstream_version(path: Path) -> str:
    root = ET.parse(path).getroot()
    releases = root.find("releases")
    if releases is None:
        raise ValueError("AppStream metainfo has no releases element.")
    current = releases.find("release")
    if current is None:
        raise ValueError("AppStream metainfo has no current release entry.")
    value = str(current.attrib.get("version", "")).strip()
    if not value:
        raise ValueError("Current AppStream release entry has no version.")
    return value


def _git_output(repo_root: Path, *args: str) -> str:
    result = subprocess.run(
        ["git", *args],
        cwd=repo_root,
        check=True,
        capture_output=True,
        text=True,
        timeout=10,
    )
    return result.stdout.strip()


def verify_tag_checkout(
    *, repo_root: Path, tag: str, event_sha: str, ref_type: str
) -> str:
    if ref_type != "tag":
        raise ValueError("Release workflow must run from a Git tag ref.")
    parse_release_tag(tag)
    checkout_sha = _git_output(repo_root, "rev-parse", "HEAD")
    tag_sha = _git_output(repo_root, "rev-parse", f"refs/tags/{tag}^{{commit}}")
    return validate_commit_identity(event_sha, tag_sha, checkout_sha)


def assert_release_absent(
    repository: str,
    tag: str,
    token: str,
    *,
    opener=None,
) -> None:
    if not re.fullmatch(r"[A-Za-z0-9_.-]+/[A-Za-z0-9_.-]+", repository or ""):
        raise ValueError("GITHUB_REPOSITORY is missing or invalid.")
    if not token:
        raise ValueError("GITHUB_TOKEN is required to check for an existing release.")
    quote_tag = urllib.parse.quote(tag, safe="")
    url = f"https://api.github.com/repos/{repository}/releases/tags/{quote_tag}"
    request = urllib.request.Request(
        url,
        headers={
            "Accept": "application/vnd.github+json",
            "Authorization": f"Bearer {token}",
            "User-Agent": "Chemuson-release-preflight",
            "X-GitHub-Api-Version": "2022-11-28",
        },
    )
    call = opener or urllib.request.urlopen
    try:
        response = call(request, timeout=15)
    except urllib.error.HTTPError as exc:
        if exc.code == 404:
            return
        raise RuntimeError(
            f"GitHub release lookup failed with HTTP {exc.code}; failing closed."
        ) from exc
    except Exception as exc:
        raise RuntimeError(
            f"GitHub release lookup failed ({exc.__class__.__name__}); failing closed."
        ) from exc

    with response:
        status_value = getattr(response, "status", None)
        status = int(status_value if status_value is not None else response.getcode())
    if status == 404:
        return
    if 200 <= status < 300:
        raise ValueError(f"A GitHub Release already exists for {tag}; refusing overwrite.")
    raise RuntimeError(f"GitHub release lookup returned HTTP {status}; failing closed.")


def _reject_deleted_event(event_path: str) -> None:
    if not event_path:
        raise ValueError("GITHUB_EVENT_PATH is required.")
    payload = json.loads(Path(event_path).read_text(encoding="utf-8"))
    if payload.get("deleted") is True:
        raise ValueError("Deleted tag events cannot publish a release.")


def main() -> None:
    tag = os.environ.get("GITHUB_REF_NAME", "").strip()
    ref_type = os.environ.get("GITHUB_REF_TYPE", "").strip()
    event_sha = os.environ.get("GITHUB_SHA", "").strip()
    _reject_deleted_event(os.environ.get("GITHUB_EVENT_PATH", ""))
    verify_tag_checkout(
        repo_root=REPO_ROOT,
        tag=tag,
        event_sha=event_sha,
        ref_type=ref_type,
    )

    version_file = Path(
        os.environ.get("CHEMUSON_VERSION_FILE", REPO_ROOT / "src/chemuson/_version.py")
    )
    metainfo_file = Path(
        os.environ.get(
            "CHEMUSON_METAINFO_FILE",
            REPO_ROOT / "packaging/flatpak/io.github.PJGV333.Chemuson.metainfo.xml",
        )
    )
    ref = validate_release_metadata(
        tag=tag,
        source_version=read_source_version(version_file),
        appstream_version=read_appstream_version(metainfo_file),
    )
    sha = validate_commit_identity(event_sha, _git_output(REPO_ROOT, "rev-parse", "HEAD"))
    assert_release_absent(
        os.environ.get("GITHUB_REPOSITORY", ""),
        tag,
        os.environ.get("RELEASE_PREFLIGHT_TOKEN", ""),
    )

    output_path = os.environ.get("GITHUB_OUTPUT", "").strip()
    if not output_path:
        raise ValueError("GITHUB_OUTPUT is required.")
    with Path(output_path).open("a", encoding="utf-8") as stream:
        stream.write(f"tag={ref.tag}\n")
        stream.write(f"version={ref.version}\n")
        stream.write(f"channel={ref.channel}\n")
        stream.write(f"sha={sha}\n")
    print(f"Validated {ref.tag} ({ref.channel}) at {sha}; no existing release.")


if __name__ == "__main__":
    try:
        main()
    except Exception as exc:
        raise SystemExit(f"Release preflight failed closed: {exc}") from exc
