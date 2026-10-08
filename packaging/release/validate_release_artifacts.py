"""Validate the exact cross-platform artifact set before release assembly."""

from __future__ import annotations

import argparse
import configparser
import json
from pathlib import Path

from release_policy import validate_commit_identity, validate_version_channel


def _require_nonempty(root: Path, name: str) -> Path:
    path = root / name
    if not path.is_file() or path.stat().st_size == 0:
        raise ValueError(f"Required release artifact is missing or empty: {name}")
    return path


def _read_ini(path: Path) -> configparser.ConfigParser:
    parser = configparser.ConfigParser(interpolation=None)
    try:
        parser.read_string(path.read_text(encoding="utf-8"))
    except (configparser.Error, UnicodeDecodeError) as exc:
        raise ValueError(f"Invalid channel config: {path.name}") from exc
    return parser


def validate_release_artifacts(
    *, root: Path, version: str, channel: str, tag: str, source_sha: str
) -> list[Path]:
    ref = validate_version_channel(version, channel)
    if tag != ref.tag:
        raise ValueError(f"Expected release tag {ref.tag}, got {tag!r}.")
    validate_commit_identity(source_sha, source_sha)

    names = [
        f"Chemuson-v{version}-windows-x86_64-portable.exe",
        f"Chemuson-v{version}-windows-x86_64-setup.exe",
        f"Chemuson-v{version}-linux-x86_64.AppImage",
        f"Chemuson-v{version}-linux-x86_64.AppImage.updateinfo",
        f"Chemuson-v{version}-linux-x86_64.AppImage.update.json",
        f"Chemuson-v{version}-linux-x86_64.flatpak",
        f"Chemuson-{channel}.flatpakref",
        f"Chemuson-{channel}.flatpakrepo",
        "build-provenance.json",
    ]
    paths = [_require_nonempty(root, name) for name in names]

    update_path = root / f"Chemuson-v{version}-linux-x86_64.AppImage.update.json"
    update = json.loads(update_path.read_text(encoding="utf-8"))
    expected_update = {
        "version": version,
        "channel": channel,
        "tag": tag,
        "source_sha": source_sha.lower(),
    }
    for key, expected in expected_update.items():
        if update.get(key) != expected:
            raise ValueError(
                f"AppImage updater metadata {key!r} must equal {expected!r}."
            )

    provenance = json.loads((root / "build-provenance.json").read_text(encoding="utf-8"))
    expected_provenance = {
        "schema_version": 1,
        "application": "Chemuson",
        "version": version,
        "tag": tag,
        "channel": channel,
        "source_sha": source_sha.lower(),
    }
    for key, expected in expected_provenance.items():
        if provenance.get(key) != expected:
            raise ValueError(f"Build provenance {key!r} must equal {expected!r}.")

    flatpak_ref = _read_ini(root / f"Chemuson-{channel}.flatpakref")
    flatpak_repo = _read_ini(root / f"Chemuson-{channel}.flatpakrepo")
    if flatpak_ref.get("Flatpak Ref", "Branch", fallback="") != channel:
        raise ValueError("Flatpak ref branch does not match release channel.")
    if flatpak_repo.get("Flatpak Repo", "DefaultBranch", fallback="") != channel:
        raise ValueError("Flatpak repo default branch does not match release channel.")
    for parser, section in ((flatpak_ref, "Flatpak Ref"), (flatpak_repo, "Flatpak Repo")):
        url = parser.get(section, "Url", fallback="")
        if f"/flatpak/{channel}/repo/" not in url:
            raise ValueError("Flatpak remote URL does not target the release channel.")

    # Catch versioned Chemuson assets accidentally copied from a different build.
    expected_prefixes = (
        f"Chemuson-v{version}-windows-",
        f"Chemuson-v{version}-linux-",
    )
    for path in root.iterdir():
        if path.is_file() and path.name.startswith("Chemuson-v"):
            if not path.name.startswith(expected_prefixes):
                raise ValueError(
                    f"Artifact {path.name!r} carries a foreign version; "
                    f"expected {version!r}."
                )
    return paths


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--dir", required=True, type=Path)
    parser.add_argument("--version", required=True)
    parser.add_argument("--channel", required=True, choices=("beta", "stable"))
    parser.add_argument("--tag", required=True)
    parser.add_argument("--source-sha", required=True)
    args = parser.parse_args()
    validated = validate_release_artifacts(
        root=args.dir.resolve(),
        version=args.version,
        channel=args.channel,
        tag=args.tag,
        source_sha=args.source_sha,
    )
    print(f"Validated {len(validated)} release artifacts for {args.tag}.")


if __name__ == "__main__":
    main()
