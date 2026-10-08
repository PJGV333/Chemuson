"""Validate the auditable, exact-identity release baseline exception catalog."""

from __future__ import annotations

import argparse
import json
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parents[2]
_REQUIRED_ENTRY_FIELDS = {
    "id",
    "kind",
    "case",
    "status",
    "owner",
    "tracking",
    "evidence",
    "retire_when",
}
_ALLOWED_KINDS = {"pytest_node", "ordered_reproducer", "static_finding"}


def _repo_file(repo_root: Path, value: str, field: str) -> Path:
    candidate = (repo_root / value).resolve()
    if not candidate.is_relative_to(repo_root.resolve()) or not candidate.is_file():
        raise ValueError(f"{field} must reference an existing repository file: {value!r}")
    return candidate


def _validate_case(kind: str, case: str) -> None:
    if not case or any(char in case for char in "*?[]"):
        raise ValueError("Exception cases must use exact identities, not wildcards.")
    if kind == "pytest_node" and (not case.startswith("tests/") or "::" not in case):
        raise ValueError("pytest_node cases must be exact pytest node IDs.")
    if kind == "ordered_reproducer":
        parts = case.split(" -> ")
        if len(parts) < 2 or any(not p.startswith("tests/") or "::" not in p for p in parts):
            raise ValueError("ordered_reproducer must list exact pytest nodes in order.")
    if kind == "static_finding" and not case.startswith("ruff "):
        raise ValueError("static_finding cases must identify the exact Ruff finding.")


def validate_catalog(data: object, repo_root: Path = REPO_ROOT) -> int:
    if not isinstance(data, dict) or data.get("schema_version") != 1:
        raise ValueError("Unsupported baseline exception catalog schema.")
    policy = data.get("policy")
    expected_policy = {
        "automatic_skips": False,
        "blanket_ignores": False,
        "unknown_failures_block": True,
        "suite_crash_is_pass": False,
    }
    if policy != expected_policy:
        raise ValueError("Baseline exception policy must remain fail-closed.")
    entries = data.get("exceptions")
    if not isinstance(entries, list):
        raise ValueError("Catalog exceptions must be a list.")
    seen_ids: set[str] = set()
    for entry in entries:
        if not isinstance(entry, dict) or set(entry) != _REQUIRED_ENTRY_FIELDS:
            raise ValueError("Every exception must contain exactly the required audit fields.")
        if any(not isinstance(entry[field], str) or not entry[field].strip() for field in _REQUIRED_ENTRY_FIELDS):
            raise ValueError("Exception audit fields must be non-empty strings.")
        exception_id = entry["id"]
        if exception_id in seen_ids:
            raise ValueError(f"Duplicate baseline exception id: {exception_id}")
        seen_ids.add(exception_id)
        kind = entry["kind"]
        if kind not in _ALLOWED_KINDS:
            raise ValueError(f"Unknown baseline exception kind: {kind}")
        if entry["status"] != "open":
            raise ValueError("This catalog records open baseline exceptions only.")
        _validate_case(kind, entry["case"])
        _repo_file(repo_root, entry["tracking"], "tracking")
        _repo_file(repo_root, entry["evidence"], "evidence")
    return len(entries)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--catalog",
        type=Path,
        default=REPO_ROOT / "docs/release/KNOWN_BASELINE_EXCEPTIONS.json",
    )
    args = parser.parse_args()
    data = json.loads(args.catalog.read_text(encoding="utf-8"))
    count = validate_catalog(data, REPO_ROOT)
    print(f"Validated {count} exact, open baseline exception records.")


if __name__ == "__main__":
    main()
