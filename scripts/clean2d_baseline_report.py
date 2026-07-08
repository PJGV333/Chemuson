#!/usr/bin/env python3
from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path
from typing import Any, Sequence


ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from tests.clean2d_regression.reports import (  # noqa: E402
    baseline_report_to_json,
    build_baseline_report,
    compare_baseline_reports,
)
from tests.clean2d_regression.review_policy import (  # noqa: E402
    classify_baseline_diff,
    summarize_baseline_diff_review,
)


EXIT_OK = 0
EXIT_CHANGED = 1
EXIT_ERROR = 2


def main(argv: Sequence[str] | None = None) -> int:
    parser = _parser()
    try:
        args = parser.parse_args(argv)
        if args.command == "write":
            return _write(args.output)
        if args.command == "compare":
            return _compare(args.left, args.right)
        if args.command == "review":
            return _review(args.left, args.right)
        parser.print_help(sys.stderr)
        return EXIT_ERROR
    except (OSError, ValueError, json.JSONDecodeError) as exc:
        print(f"error: {exc}", file=sys.stderr)
        return EXIT_ERROR


def _parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description="Clean 2D baseline report developer utility")
    subparsers = parser.add_subparsers(dest="command", required=True)

    write = subparsers.add_parser("write", help="write a Clean 2D baseline report")
    write.add_argument("--output", required=True, type=Path)

    compare = subparsers.add_parser("compare", help="compare two Clean 2D baseline reports")
    compare.add_argument("--left", required=True, type=Path)
    compare.add_argument("--right", required=True, type=Path)

    review = subparsers.add_parser("review", help="review differences between two Clean 2D baseline reports")
    review.add_argument("--left", required=True, type=Path)
    review.add_argument("--right", required=True, type=Path)
    return parser


def _write(output: Path) -> int:
    report = build_baseline_report(metadata={"source": "clean2d_baseline_report.py"})
    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text(baseline_report_to_json(report), encoding="utf-8")
    return EXIT_OK


def _compare(left: Path, right: Path) -> int:
    diff = compare_baseline_reports(_read_json(left), _read_json(right))
    print(json.dumps(diff, allow_nan=False, sort_keys=True))
    return EXIT_OK if diff["equivalent"] else EXIT_CHANGED


def _review(left: Path, right: Path) -> int:
    left_report = _read_json(left)
    right_report = _read_json(right)
    diff = compare_baseline_reports(left_report, right_report)
    items = classify_baseline_diff(diff, report=left_report)
    summary = summarize_baseline_diff_review(items)
    print(json.dumps(summary, allow_nan=False, sort_keys=True))
    return EXIT_OK if summary["item_count"] == 0 else EXIT_CHANGED


def _read_json(path: Path) -> dict[str, Any]:
    data = json.loads(path.read_text(encoding="utf-8"))
    if not isinstance(data, dict):
        raise ValueError(f"expected JSON object in {path}")
    return data


if __name__ == "__main__":
    raise SystemExit(main())
