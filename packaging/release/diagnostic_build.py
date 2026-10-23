"""Validation and provenance helpers for diagnostic package builds.

Diagnostic builds exercise the full packaging pipeline on repair/feature
branches without touching the preview/release channel contract: they keep
every integrity control (canonical version, exact source SHA, clean worktree,
SHA-256 provenance manifest, packaged validators) but they are never
published, never embed updater information, and their artifact names always
carry a ``diagnostic`` marker so they can never impersonate a release.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import os
import platform
import re
import subprocess
from datetime import datetime, timezone
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parents[2]
_VERSION_RE = re.compile(
    r"^(?:0|[1-9][0-9]*)\.(?:0|[1-9][0-9]*)\.(?:0|[1-9][0-9]*)"
    r"(?:-[0-9A-Za-z][0-9A-Za-z.-]*)?(?:\+[0-9A-Za-z][0-9A-Za-z.-]*)?$"
)
_SHA_RE = re.compile(r"^(?:[0-9a-f]{40}|[0-9a-f]{64})$", re.IGNORECASE)
_DIAGNOSTIC_BRANCH_RE = re.compile(
    r"^(?:ai|chore|docs|feat|fix|hotfix|integration|investigation|spike|stabilize)/[A-Za-z0-9._/-]+$"
)
_RELEASE_BRANCH_RE = re.compile(r"^release/")


def read_canonical_version(repo_root: Path = REPO_ROOT) -> str:
    version_path = repo_root / "src/chemuson/_version.py"
    pattern = re.compile(r'^\s*__version__\s*=\s*"([^"]+)"\s*$')
    for line in version_path.read_text(encoding="utf-8").splitlines():
        match = pattern.fullmatch(line)
        if match:
            version = match.group(1)
            if not _VERSION_RE.fullmatch(version):
                raise ValueError(f"Unsupported diagnostic version in canonical source: {version!r}")
            return version
    raise ValueError(f"No canonical __version__ found in {version_path}.")


def _git_head(repo_root: Path) -> str:
    result = subprocess.run(
        ["git", "rev-parse", "HEAD"],
        cwd=repo_root,
        check=True,
        capture_output=True,
        text=True,
        timeout=10,
    )
    return result.stdout.strip().lower()


def _git_worktree_is_clean(repo_root: Path) -> bool:
    result = subprocess.run(
        ["git", "status", "--porcelain"],
        cwd=repo_root,
        check=True,
        capture_output=True,
        text=True,
        timeout=30,
    )
    return not result.stdout.strip()


def validate_diagnostic_identity(
    *,
    version: str,
    git_sha: str,
    source_branch: str,
    repo_root: Path = REPO_ROOT,
) -> str:
    canonical_version = read_canonical_version(repo_root)
    if version != canonical_version:
        raise ValueError(
            f"Requested diagnostic version {version!r} differs from canonical "
            f"{canonical_version!r}."
        )
    if not _VERSION_RE.fullmatch(version):
        raise ValueError(f"Invalid diagnostic version: {version!r}")
    normalized_sha = str(git_sha or "").strip().lower()
    if not _SHA_RE.fullmatch(normalized_sha):
        raise ValueError("Diagnostic source must be identified by a full Git SHA.")
    if _git_head(repo_root) != normalized_sha:
        raise ValueError("Checked-out diagnostic source does not match the expected Git SHA.")
    if not _git_worktree_is_clean(repo_root):
        raise ValueError(
            "Diagnostic builds require a clean worktree so the recorded SHA "
            "describes every byte of the build."
        )
    branch = str(source_branch or "")
    if _RELEASE_BRANCH_RE.match(branch):
        raise ValueError(
            "Diagnostic builds must not run from release/** branches; use the "
            "preview channel for release preparation."
        )
    if not _DIAGNOSTIC_BRANCH_RE.fullmatch(branch):
        raise ValueError(
            "Diagnostic builds are restricted to ai/, chore/, docs/, feat/, fix/, "
            "hotfix/, integration/, investigation/, spike/ or stabilize/ branches."
        )
    return normalized_sha


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for chunk in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def build_provenance(
    *,
    version: str,
    git_sha: str,
    source_branch: str,
    operating_system: str,
    files: list[str],
    repo_root: Path = REPO_ROOT,
    artifact_dir: Path | None = None,
    checks: dict[str, str] | None = None,
) -> dict[str, object]:
    sha = validate_diagnostic_identity(
        version=version,
        git_sha=git_sha,
        source_branch=source_branch,
        repo_root=repo_root,
    )
    if not files:
        raise ValueError("A diagnostic artifact group must contain at least one package file.")
    artifacts: dict[str, dict[str, object]] = {}
    for name in files:
        if Path(name).name != name or "diagnostic" not in name.lower():
            raise ValueError(
                f"Diagnostic package name must carry the 'diagnostic' marker: {name!r}"
            )
        path = (artifact_dir / name) if artifact_dir is not None else (repo_root / name)
        if not path.is_file() or path.stat().st_size <= 0:
            raise ValueError(f"Diagnostic artifact is missing or empty: {path}")
        artifacts[name] = {"size_bytes": path.stat().st_size, "sha256": _sha256(path)}
    return {
        "version": version,
        "build_type": "diagnostic",
        "git_sha": sha,
        "source_branch": source_branch,
        "build_time_utc": datetime.now(timezone.utc).isoformat(),
        "operating_system": {
            "platform": operating_system,
            "system": platform.system(),
            "release": platform.release(),
            "machine": platform.machine(),
            "runner_os": os.environ.get("RUNNER_OS", ""),
            "runner_arch": os.environ.get("RUNNER_ARCH", ""),
            "runner_image": os.environ.get("ImageOS", ""),
            "runner_image_version": os.environ.get("ImageVersion", ""),
        },
        "build_status": "success",
        "publication": False,
        "distribution": "none",
        "not_for_publication": True,
        "checks": dict(checks or {}),
        "artifacts": artifacts,
    }


def write_group(
    *,
    output_dir: Path,
    version: str,
    git_sha: str,
    source_branch: str,
    operating_system: str,
    files: list[str],
    repo_root: Path = REPO_ROOT,
    checks: dict[str, str] | None = None,
) -> dict[str, object]:
    output_dir.mkdir(parents=True, exist_ok=True)
    root = repo_root.resolve()
    resolved_dir = output_dir.resolve()
    if not resolved_dir.is_relative_to(root):
        raise ValueError("Diagnostic artifacts must be staged inside the checked-out repository.")
    record = build_provenance(
        version=version,
        git_sha=git_sha,
        source_branch=source_branch,
        operating_system=operating_system,
        files=files,
        repo_root=root,
        artifact_dir=resolved_dir,
        checks=checks,
    )
    provenance_path = resolved_dir / "diagnostic-provenance.json"
    provenance_path.write_text(
        json.dumps(record, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    checksum_lines = [
        f"{record['artifacts'][name]['sha256']}  {name}\n" for name in files
    ]
    checksum_lines.append(f"{_sha256(provenance_path)}  {provenance_path.name}\n")
    (resolved_dir / "checksums.sha256").write_text("".join(checksum_lines), encoding="utf-8")
    return record


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    subparsers = parser.add_subparsers(dest="command", required=True)

    verify_parser = subparsers.add_parser("verify", help="Verify an exact checked-out source SHA.")
    verify_parser.add_argument("--repo-root", type=Path, default=REPO_ROOT)
    verify_parser.add_argument("--version", required=True)
    verify_parser.add_argument("--git-sha", required=True)
    verify_parser.add_argument("--source-branch", required=True)

    manifest_parser = subparsers.add_parser(
        "manifest", help="Validate and checksum one diagnostic artifact group."
    )
    manifest_parser.add_argument("--output-dir", required=True, type=Path)
    manifest_parser.add_argument("--repo-root", type=Path, default=REPO_ROOT)
    manifest_parser.add_argument("--version", required=True)
    manifest_parser.add_argument("--git-sha", required=True)
    manifest_parser.add_argument("--source-branch", required=True)
    manifest_parser.add_argument("--operating-system", required=True)
    manifest_parser.add_argument("--files", nargs="+", required=True)

    args = parser.parse_args()
    if args.command == "verify":
        validate_diagnostic_identity(
            version=args.version,
            git_sha=args.git_sha,
            source_branch=args.source_branch,
            repo_root=args.repo_root.resolve(),
        )
        print(f"Verified diagnostic source {args.git_sha} on branch {args.source_branch}")
    else:
        record = write_group(
            output_dir=args.output_dir,
            version=args.version,
            git_sha=args.git_sha,
            source_branch=args.source_branch,
            operating_system=args.operating_system,
            files=args.files,
            repo_root=args.repo_root.resolve(),
        )
        print(f"Validated and checksummed {len(record['artifacts'])} diagnostic package(s).")


if __name__ == "__main__":
    try:
        main()
    except Exception as exc:
        raise SystemExit(f"Diagnostic build validation failed: {exc}") from exc
