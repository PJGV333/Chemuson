"""Write a minimal release provenance record for the validated source commit."""

from __future__ import annotations

import argparse
import json
from pathlib import Path

from release_policy import validate_commit_identity, validate_version_channel


def build_provenance(*, version: str, channel: str, tag: str, source_sha: str) -> dict[str, object]:
    ref = validate_version_channel(version, channel)
    if ref.tag != tag:
        raise ValueError(f"Tag {tag!r} does not match version {version!r}.")
    sha = validate_commit_identity(source_sha, source_sha)
    return {
        "schema_version": 1,
        "application": "Chemuson",
        "version": ref.version,
        "tag": ref.tag,
        "channel": ref.channel,
        "source_sha": sha,
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--version", required=True)
    parser.add_argument("--channel", required=True, choices=("beta", "stable"))
    parser.add_argument("--tag", required=True)
    parser.add_argument("--source-sha", required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    record = build_provenance(
        version=args.version,
        channel=args.channel,
        tag=args.tag,
        source_sha=args.source_sha,
    )
    output = args.output.resolve()
    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text(json.dumps(record, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(f"Wrote release provenance: {output}")


if __name__ == "__main__":
    main()
