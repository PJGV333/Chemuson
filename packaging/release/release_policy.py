"""Strict release tag/version/channel policy (standard library only)."""

from __future__ import annotations

import argparse
import re
from dataclasses import dataclass

_VERSION_RE = re.compile(
    r"^(?P<major>0|[1-9][0-9]*)\."
    r"(?P<minor>0|[1-9][0-9]*)\."
    r"(?P<patch>0|[1-9][0-9]*)"
    r"(?:-(?P<kind>beta|rc)\.(?P<sequence>[1-9][0-9]*))?$"
)
_SHA_RE = re.compile(r"^(?:[0-9a-f]{40}|[0-9a-f]{64})$", re.IGNORECASE)


@dataclass(frozen=True, slots=True)
class ReleaseRef:
    tag: str
    version: str
    channel: str
    prerelease: bool


def parse_release_tag(tag: str) -> ReleaseRef:
    """Parse only stable, beta.N, or rc.N tags accepted for publication."""
    normalized_tag = str(tag or "").strip()
    if not normalized_tag.startswith("v"):
        raise ValueError("Release tag must start with 'v'.")
    version = normalized_tag[1:]
    match = _VERSION_RE.fullmatch(version)
    if match is None:
        raise ValueError(f"Unsupported or invalid release tag: {tag!r}")
    prerelease = match.group("kind") is not None
    return ReleaseRef(
        tag=normalized_tag,
        version=version,
        channel="beta" if prerelease else "stable",
        prerelease=prerelease,
    )


def validate_version_channel(version: str, channel: str) -> ReleaseRef:
    """Require a release version and its only valid distribution channel."""
    ref = parse_release_tag(f"v{str(version or '').strip()}")
    if str(channel or "").strip() != ref.channel:
        raise ValueError(
            f"Version {ref.version} belongs to channel {ref.channel!r}, "
            f"not {channel!r}."
        )
    return ref


def validate_release_metadata(
    *, tag: str, source_version: str, appstream_version: str
) -> ReleaseRef:
    """Require tag, canonical application version, and current AppStream entry to match."""
    ref = parse_release_tag(tag)
    if str(source_version or "").strip() != ref.version:
        raise ValueError(
            f"Tag {ref.tag} does not match canonical application version "
            f"{source_version!r}."
        )
    if str(appstream_version or "").strip() != ref.version:
        raise ValueError(
            f"Tag {ref.tag} does not match current AppStream version "
            f"{appstream_version!r}."
        )
    return ref


def validate_commit_identity(*shas: str) -> str:
    """Require valid, identical full Git SHAs for event, tag and checkout."""
    normalized = [str(sha or "").strip().lower() for sha in shas]
    if len(normalized) < 2 or any(not _SHA_RE.fullmatch(sha) for sha in normalized):
        raise ValueError("All release source identifiers must be full Git SHAs.")
    if len(set(normalized)) != 1:
        raise ValueError("Tag, event and checked-out source SHAs do not match.")
    return normalized[0]


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--version", required=True)
    parser.add_argument("--channel", required=True, choices=("beta", "stable"))
    parser.add_argument("--tag", help="Optional exact v-prefixed tag to verify.")
    args = parser.parse_args()
    ref = validate_version_channel(args.version, args.channel)
    if args.tag and args.tag != ref.tag:
        raise SystemExit(f"Tag {args.tag!r} does not match {ref.tag!r}.")
    print(f"Validated {ref.tag} for {ref.channel}")


if __name__ == "__main__":
    main()
