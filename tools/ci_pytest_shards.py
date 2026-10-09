#!/usr/bin/env python3
"""Plan, supervise, and verify exact pytest shards for GitHub Actions."""

from __future__ import annotations

import argparse
import hashlib
import json
import os
import signal
import subprocess
import sys
import threading
import time
import xml.etree.ElementTree as ET
from collections import Counter, defaultdict
from pathlib import Path, PurePosixPath
from typing import Any, Sequence

ROOT = Path(__file__).resolve().parents[1]
DEFAULT_MANIFEST = ROOT / "tests" / "ci" / "expected_pytest_nodeids.txt"
PLAN_SCHEMA = 1
DEFAULT_SHARD_COUNT = 8
MAX_TIMEOUT_SECONDS = 600


def _digest_nodeids(nodeids: Sequence[str]) -> str:
    payload = "\n".join(sorted(nodeids)) + "\n"
    return hashlib.sha256(payload.encode("utf-8")).hexdigest()


def _read_expected_nodeids(path: Path) -> list[str]:
    try:
        nodeids = path.read_text(encoding="utf-8").splitlines()
    except OSError as exc:
        raise ValueError(f"Cannot read expected node ID manifest {path}: {exc}") from exc
    if not nodeids:
        raise ValueError(f"Expected node ID manifest is empty: {path}")
    if nodeids != sorted(nodeids):
        raise ValueError(f"Expected node ID manifest must be sorted: {path}")
    duplicates = [nodeid for nodeid, count in Counter(nodeids).items() if count > 1]
    if duplicates:
        raise ValueError(f"Expected node ID manifest has duplicates: {duplicates[:10]}")
    if any(not nodeid.startswith("tests/") or "::" not in nodeid for nodeid in nodeids):
        raise ValueError(f"Expected node ID manifest contains an invalid pytest ID: {path}")
    return nodeids


class _NodeIdCollector:
    def __init__(self) -> None:
        self.nodeids: list[str] = []

    def pytest_collection_finish(self, session: Any) -> None:
        self.nodeids = [item.nodeid for item in session.items]


def _collect_nodeids() -> list[str]:
    if str(ROOT) not in sys.path:
        sys.path.insert(0, str(ROOT))
    import pytest

    collector = _NodeIdCollector()
    from contextlib import redirect_stderr, redirect_stdout
    from io import StringIO

    stdout = StringIO()
    stderr = StringIO()
    with redirect_stdout(stdout), redirect_stderr(stderr):
        exit_code = pytest.main(
            ["--collect-only", "-q", "--disable-warnings"],
            plugins=[collector],
        )
    if int(exit_code) != 0:
        diagnostic = "\n".join(
            part for part in (stdout.getvalue(), stderr.getvalue()) if part
        )
        raise ValueError(
            f"pytest collection failed with exit code {int(exit_code)}:\n{diagnostic[-12000:]}"
        )
    if not collector.nodeids:
        raise ValueError("pytest collection produced no test node IDs")
    duplicates = [nodeid for nodeid, count in Counter(collector.nodeids).items() if count > 1]
    if duplicates:
        raise ValueError(f"pytest collection produced duplicate IDs: {duplicates[:10]}")
    return collector.nodeids


def _check_collected_against_manifest(
    nodeids: Sequence[str], manifest_path: Path
) -> list[str]:
    expected = _read_expected_nodeids(manifest_path)
    collected = set(nodeids)
    expected_set = set(expected)
    missing = sorted(expected_set - collected)
    unexpected = sorted(collected - expected_set)
    if missing or unexpected:
        sections = [
            f"pytest collection does not match {manifest_path}: "
            f"expected={len(expected)}, collected={len(nodeids)}"
        ]
        if missing:
            sections.append(f"missing ({len(missing)}): {missing[:20]}")
        if unexpected:
            sections.append(f"unexpected/new ({len(unexpected)}): {unexpected[:20]}")
        sections.append(
            "Review the change and update the manifest explicitly with "
            "`python tools/ci_pytest_shards.py update-expected --confirm`."
        )
        raise ValueError("\n".join(sections))
    return expected


def _partition_by_file(nodeids: Sequence[str], shard_count: int) -> list[list[str]]:
    if shard_count < 1:
        raise ValueError("shard count must be positive")
    by_file: dict[str, list[tuple[int, str]]] = {}
    for ordinal, nodeid in enumerate(nodeids):
        test_file = nodeid.split("::", 1)[0]
        by_file.setdefault(test_file, []).append((ordinal, nodeid))

    loads = [0] * shard_count
    shards_with_order: list[list[tuple[int, str]]] = [[] for _ in range(shard_count)]
    files = sorted(
        by_file.items(),
        key=lambda pair: (-len(pair[1]), pair[1][0][0]),
    )
    for _test_file, members in files:
        shard_index = min(range(shard_count), key=lambda index: (loads[index], index))
        shards_with_order[shard_index].extend(members)
        loads[shard_index] += len(members)
    return [
        [nodeid for _ordinal, nodeid in sorted(members)]
        for members in shards_with_order
    ]


def _validate_plan(plan: dict[str, Any], expected: Sequence[str]) -> None:
    if plan.get("schema_version") != PLAN_SCHEMA:
        raise ValueError(f"Unsupported shard plan schema: {plan.get('schema_version')!r}")
    shards = plan.get("shards")
    shard_count = plan.get("shard_count")
    total = plan.get("total_tests")
    if not isinstance(shard_count, int) or shard_count < 1:
        raise ValueError("Shard plan has an invalid shard_count")
    if not isinstance(shards, list) or len(shards) != shard_count:
        raise ValueError("Shard plan does not contain every shard")
    if total != len(expected) or plan.get("expected_nodeids_sha256") != _digest_nodeids(expected):
        raise ValueError("Shard plan collection count/hash does not match the committed manifest")

    flattened: list[str] = []
    distribution: list[int] = []
    for index, shard in enumerate(shards):
        if shard.get("index") != index:
            raise ValueError(f"Shard plan index is missing or duplicated at {index}")
        nodeids = shard.get("nodeids")
        if not isinstance(nodeids, list) or any(not isinstance(item, str) for item in nodeids):
            raise ValueError(f"Shard {index} has an invalid node ID list")
        if shard.get("count") != len(nodeids):
            raise ValueError(f"Shard {index} count does not match its assignment")
        distribution.append(len(nodeids))
        flattened.extend(nodeids)
    actual_counts = Counter(flattened)
    expected_counts = Counter(expected)
    if actual_counts != expected_counts:
        duplicates = sorted(nodeid for nodeid, count in actual_counts.items() if count > 1)
        missing = sorted(expected_counts.keys() - actual_counts.keys())
        unexpected = sorted(actual_counts.keys() - expected_counts.keys())
        raise ValueError(
            "Shard assignment does not cover the manifest exactly once: "
            f"duplicates={duplicates[:10]}, missing={missing[:10]}, "
            f"unexpected={unexpected[:10]}"
        )
    if plan.get("distribution") != distribution:
        raise ValueError("Shard plan distribution does not match assigned node IDs")


def _write_json(path: Path, data: dict[str, Any]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(
        json.dumps(data, ensure_ascii=False, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )


def create_plan(args: argparse.Namespace) -> int:
    manifest_path = args.expected.resolve()
    nodeids = _collect_nodeids()
    expected = _check_collected_against_manifest(nodeids, manifest_path)
    assignments = _partition_by_file(nodeids, args.shards)
    plan = {
        "schema_version": PLAN_SCHEMA,
        "total_tests": len(nodeids),
        "expected_nodeids_sha256": _digest_nodeids(expected),
        "shard_count": args.shards,
        "distribution": [len(shard) for shard in assignments],
        "shards": [
            {"index": index, "count": len(shard), "nodeids": shard}
            for index, shard in enumerate(assignments)
        ],
    }
    _validate_plan(plan, expected)
    output_path = args.output.resolve()
    _write_json(output_path, plan)
    matrix = {"shard": list(range(args.shards))}
    print(
        f"pytest collection: {len(nodeids)} node IDs; "
        f"{args.shards} shards; distribution={plan['distribution']}; "
        f"sha256={plan['expected_nodeids_sha256']}"
    )
    print(f"plan: {output_path}")
    if args.github_output:
        with args.github_output.open("a", encoding="utf-8") as output:
            output.write(f"shard-matrix={json.dumps(matrix, separators=(',', ':'))}\n")
            output.write(f"total-tests={len(nodeids)}\n")
            output.write(
                "distribution="
                + json.dumps(plan["distribution"], separators=(",", ":"))
                + "\n"
            )
    return 0


def update_expected(args: argparse.Namespace) -> int:
    if not args.confirm:
        raise ValueError("Refusing to update the manifest without --confirm")
    nodeids = _collect_nodeids()
    path = args.expected.resolve()
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text("\n".join(sorted(nodeids)) + "\n", encoding="utf-8")
    print(f"wrote {len(nodeids)} sorted unique node IDs to {path}")
    return 0


def _stop_process_group(process: subprocess.Popen[str], grace_seconds: float = 2.0) -> None:
    if os.name != "posix":
        if process.poll() is None:
            process.terminate()
        try:
            process.wait(timeout=grace_seconds)
        except subprocess.TimeoutExpired:
            process.kill()
            process.wait()
        return

    process_group = process.pid
    try:
        os.killpg(process_group, signal.SIGTERM)
    except ProcessLookupError:
        return
    deadline = time.monotonic() + grace_seconds
    while time.monotonic() < deadline:
        try:
            os.killpg(process_group, 0)
        except ProcessLookupError:
            return
        time.sleep(0.05)
    try:
        os.killpg(process_group, signal.SIGKILL)
    except ProcessLookupError:
        pass


def _split_pytest_nodeid(nodeid: str) -> list[str]:
    """Split `::` separators without treating parameter values as path parts."""
    parts: list[str] = []
    start = 0
    bracket_depth = 0
    index = 0
    while index < len(nodeid):
        character = nodeid[index]
        if character == "[":
            bracket_depth += 1
        elif character == "]" and bracket_depth:
            bracket_depth -= 1
        elif nodeid.startswith("::", index) and bracket_depth == 0:
            parts.append(nodeid[start:index])
            start = index + 2
            index += 1
        index += 1
    parts.append(nodeid[start:])
    return parts


def _junit_nodeids(junit_path: Path, assigned: Sequence[str]) -> tuple[list[str], list[str]]:
    root = ET.parse(junit_path).getroot()
    by_case: dict[tuple[str, str], list[str]] = defaultdict(list)
    for nodeid in assigned:
        parts = _split_pytest_nodeid(nodeid)
        test_path = PurePosixPath(parts[0]).with_suffix("").as_posix()
        module_name = test_path.replace("/", ".")
        class_name = module_name
        if len(parts) > 2:
            class_name += "." + ".".join(parts[1:-1])
        by_case[(class_name, parts[-1])].append(nodeid)

    actual: list[str] = []
    unmatched: list[str] = []
    for case in root.iter("testcase"):
        key = (case.attrib.get("classname", ""), case.attrib.get("name", ""))
        candidates = by_case.get(key)
        if not candidates:
            unmatched.append(f"{key[0]}::{key[1]}")
            continue
        actual.append(candidates.pop(0))
    remaining = [nodeid for candidates in by_case.values() for nodeid in candidates]
    if remaining:
        unmatched.extend(f"not executed: {nodeid}" for nodeid in remaining)
    return actual, unmatched


def _junit_test_count(path: Path) -> int:
    root = ET.parse(path).getroot()
    return sum(1 for _case in root.iter("testcase"))


def run_shard(args: argparse.Namespace) -> int:
    manifest_path = args.expected.resolve()
    expected = _read_expected_nodeids(manifest_path)
    try:
        plan = json.loads(args.plan.read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError) as exc:
        raise ValueError(f"Cannot read shard plan {args.plan}: {exc}") from exc
    _validate_plan(plan, expected)
    if not 0 <= args.shard_index < plan["shard_count"]:
        raise ValueError(f"Invalid shard index {args.shard_index}")
    if not 1 <= args.timeout_seconds <= MAX_TIMEOUT_SECONDS:
        raise ValueError(f"timeout must be between 1 and {MAX_TIMEOUT_SECONDS} seconds")

    shard = plan["shards"][args.shard_index]
    nodeids: list[str] = shard["nodeids"]
    output_dir = args.output_dir.resolve()
    output_dir.mkdir(parents=True, exist_ok=True)
    prefix = f"pytest-shard-{args.shard_index}"
    log_path = output_dir / f"{prefix}.log"
    junit_path = output_dir / f"{prefix}.xml"
    report_path = output_dir / f"{prefix}.json"
    command = [
        sys.executable,
        "-m",
        "pytest",
        "-vv",
        "-o",
        "faulthandler_timeout=120",
        f"--junitxml={junit_path}",
        *nodeids,
    ]
    environment = os.environ.copy()
    environment.pop("PYTEST_ADDOPTS", None)
    environment["PYTHONFAULTHANDLER"] = "1"
    environment.setdefault("QT_QPA_PLATFORM", "offscreen")

    started = time.monotonic()
    last_test: dict[str, str | None] = {"nodeid": None}
    timed_out = False
    launch_error: str | None = None
    process: subprocess.Popen[str] | None = None
    reader: threading.Thread | None = None
    returncode: int | None = None
    with log_path.open("w", encoding="utf-8", errors="replace") as logfile:
        try:
            process = subprocess.Popen(
                command,
                cwd=ROOT,
                env=environment,
                stdout=subprocess.PIPE,
                stderr=subprocess.STDOUT,
                text=True,
                bufsize=1,
                start_new_session=(os.name == "posix"),
            )
        except OSError as exc:
            launch_error = f"{type(exc).__name__}: {exc}"
            print(f"could not start pytest shard {args.shard_index}: {launch_error}")
        else:
            assert process.stdout is not None

            def stream_output() -> None:
                for line in process.stdout:
                    logfile.write(line)
                    logfile.flush()
                    sys.stdout.write(line)
                    sys.stdout.flush()
                    stripped = line.strip()
                    if stripped.startswith("tests/") and "::" in stripped:
                        last_test["nodeid"] = stripped.split()[0]

            reader = threading.Thread(target=stream_output, daemon=True)
            reader.start()
            deadline = started + args.timeout_seconds
            while process.poll() is None:
                if time.monotonic() >= deadline:
                    timed_out = True
                    print(
                        f"pytest shard {args.shard_index} exceeded "
                        f"{args.timeout_seconds}s; terminating its process group",
                        flush=True,
                    )
                    _stop_process_group(process)
                    break
                time.sleep(0.1)
            try:
                returncode = process.wait(timeout=3)
            except subprocess.TimeoutExpired:
                _stop_process_group(process, grace_seconds=1.0)
                returncode = process.wait(timeout=3)
            # Also reap subprocesses a test may have left behind after pytest exits.
            _stop_process_group(process, grace_seconds=0.5)
            reader.join(timeout=5)
            if reader.is_alive() and process.stdout is not None:
                process.stdout.close()
                reader.join(timeout=1)

    elapsed = round(time.monotonic() - started, 3)
    signal_number = -returncode if returncode is not None and returncode < 0 else None
    try:
        signal_name = signal.Signals(signal_number).name if signal_number is not None else None
    except ValueError:
        signal_name = f"SIG{signal_number}" if signal_number is not None else None

    junit_count: int | None = None
    junit_matches_assignment = False
    junit_error: str | None = None
    if junit_path.is_file():
        try:
            junit_count = _junit_test_count(junit_path)
            actual_nodeids, unmatched = _junit_nodeids(junit_path, nodeids)
            junit_matches_assignment = not unmatched and Counter(actual_nodeids) == Counter(nodeids)
            if unmatched:
                junit_error = "; ".join(unmatched[:20])
        except (ET.ParseError, OSError, ValueError) as exc:
            junit_error = f"{type(exc).__name__}: {exc}"
    else:
        junit_error = "JUnit report was not created"

    report = {
        "schema_version": PLAN_SCHEMA,
        "shard_index": args.shard_index,
        "shard_count": plan["shard_count"],
        "total_tests": plan["total_tests"],
        "expected_nodeids_sha256": plan["expected_nodeids_sha256"],
        "assigned_test_count": len(nodeids),
        "assigned_nodeids_sha256": _digest_nodeids(nodeids),
        "timeout_seconds": args.timeout_seconds,
        "process_exit_code": returncode,
        "timed_out": timed_out,
        "signal_number": signal_number,
        "signal_name": signal_name,
        "launch_error": launch_error,
        "elapsed_seconds": elapsed,
        "last_test_started": last_test["nodeid"],
        "junit_file": junit_path.name,
        "junit_test_count": junit_count,
        "junit_matches_assignment": junit_matches_assignment,
        "junit_error": junit_error,
        "log_file": log_path.name,
    }
    _write_json(report_path, report)
    print(
        f"shard {args.shard_index}: exit={returncode}, timeout={timed_out}, "
        f"signal={signal_name}, tests={junit_count}/{len(nodeids)}, "
        f"elapsed={elapsed:.3f}s, last={last_test['nodeid']}",
        flush=True,
    )
    return 0 if (
        returncode == 0
        and not timed_out
        and junit_count == len(nodeids)
        and junit_matches_assignment
        and launch_error is None
    ) else 1


def _load_json(path: Path) -> dict[str, Any]:
    try:
        value = json.loads(path.read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError) as exc:
        raise ValueError(f"Cannot read JSON {path}: {exc}") from exc
    if not isinstance(value, dict):
        raise ValueError(f"Expected JSON object in {path}")
    return value


def verify_results(args: argparse.Namespace) -> int:
    errors: list[str] = []
    if args.plan_status != "success":
        errors.append(f"pytest plan job status is {args.plan_status!r}")
    if args.shards_status != "success":
        errors.append(f"pytest matrix job status is {args.shards_status!r}")

    manifest_path = args.expected.resolve()
    expected = _read_expected_nodeids(manifest_path)
    plan_files = sorted(args.artifacts.rglob("pytest-plan.json"))
    if len(plan_files) != 1:
        errors.append(f"expected exactly one pytest-plan.json artifact, found {len(plan_files)}")
        plan = None
    else:
        try:
            plan = _load_json(plan_files[0])
            _validate_plan(plan, expected)
        except ValueError as exc:
            errors.append(str(exc))
            plan = None

    reports = sorted(args.artifacts.rglob("pytest-shard-*.json"))
    if plan is None:
        expected_indices = set(range(args.shard_count))
        planned_assignments: dict[int, list[str]] = {}
    else:
        if plan.get("shard_count") != args.shard_count:
            errors.append(
                f"plan has {plan.get('shard_count')} shards, expected {args.shard_count}"
            )
        expected_indices = set(range(plan.get("shard_count", 0)))
        planned_assignments = {
            shard["index"]: shard["nodeids"] for shard in plan.get("shards", [])
        }

    reports_by_index: dict[int, tuple[Path, dict[str, Any]]] = {}
    for path in reports:
        try:
            report = _load_json(path)
            index = report.get("shard_index")
        except ValueError as exc:
            errors.append(str(exc))
            continue
        if not isinstance(index, int) or index in reports_by_index:
            errors.append(f"invalid or duplicate shard report index in {path}")
            continue
        reports_by_index[index] = (path, report)
    missing_reports = sorted(expected_indices - reports_by_index.keys())
    extra_reports = sorted(reports_by_index.keys() - expected_indices)
    if missing_reports:
        errors.append(f"missing shard result reports: {missing_reports}")
    if extra_reports:
        errors.append(f"unexpected shard result reports: {extra_reports}")

    for index in sorted(reports_by_index):
        path, report = reports_by_index[index]
        assigned = planned_assignments.get(index)
        print(
            f"shard {index}: exit={report.get('process_exit_code')}, "
            f"timeout={report.get('timed_out')}, signal={report.get('signal_name')}, "
            f"tests={report.get('junit_test_count')}/{report.get('assigned_test_count')}, "
            f"elapsed={report.get('elapsed_seconds')}s, "
            f"last={report.get('last_test_started')}"
        )
        if report.get("schema_version") != PLAN_SCHEMA:
            errors.append(f"shard {index} report has an unsupported schema")
        if assigned is None:
            errors.append(f"shard {index} has no assignment in the plan")
            continue
        if report.get("assigned_test_count") != len(assigned):
            errors.append(f"shard {index} assigned test count does not match the plan")
        if report.get("assigned_nodeids_sha256") != _digest_nodeids(assigned):
            errors.append(f"shard {index} assigned node IDs do not match the plan")
        if report.get("expected_nodeids_sha256") != _digest_nodeids(expected):
            errors.append(f"shard {index} used a different collection manifest")
        if report.get("process_exit_code") != 0:
            errors.append(
                f"shard {index} pytest exit code is {report.get('process_exit_code')!r}; "
                f"log={report.get('log_file')}, last={report.get('last_test_started')}"
            )
        if report.get("timed_out") is not False:
            errors.append(f"shard {index} timed out")
        if not isinstance(report.get("timeout_seconds"), int) or report["timeout_seconds"] > 300:
            errors.append(f"shard {index} external timeout exceeds 300 seconds")
        if report.get("signal_name") is not None:
            errors.append(f"shard {index} exited by {report.get('signal_name')}")
        if report.get("launch_error") is not None:
            errors.append(f"shard {index} launch error: {report.get('launch_error')}")

        log_name = report.get("log_file")
        log_path = path.parent / log_name if isinstance(log_name, str) else None
        if log_path is None or not log_path.is_file():
            errors.append(f"shard {index} log artifact is missing")

        junit_name = report.get("junit_file")
        junit_path = path.parent / junit_name if isinstance(junit_name, str) else None
        if junit_path is None or not junit_path.is_file():
            errors.append(f"shard {index} JUnit artifact is missing")
            continue
        try:
            actual_nodeids, unmatched = _junit_nodeids(junit_path, assigned)
        except (ET.ParseError, OSError, ValueError) as exc:
            errors.append(f"shard {index} JUnit parse error: {exc}")
            continue
        junit_count = len(actual_nodeids) + len(unmatched)
        if junit_count != len(assigned) or unmatched or Counter(actual_nodeids) != Counter(assigned):
            errors.append(
                f"shard {index} JUnit coverage mismatch: expected={len(assigned)}, "
                f"matched={len(actual_nodeids)}, unmatched={unmatched[:10]}"
            )
        if report.get("junit_test_count") != junit_count:
            errors.append(
                f"shard {index} report/JUnit count mismatch: "
                f"{report.get('junit_test_count')} != {junit_count}"
            )
        if report.get("junit_matches_assignment") is not True:
            errors.append(f"shard {index} runner marked its JUnit assignment invalid")

    if errors:
        print("pytest shard coverage verification FAILED:")
        for error in errors:
            print(f"- {error}")
        return 1
    print(
        f"pytest shard coverage verified: {len(expected)} unique node IDs, "
        f"{args.shard_count} complete successful shards"
    )
    return 0


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    subparsers = parser.add_subparsers(dest="command", required=True)

    plan_parser = subparsers.add_parser("plan", help="collect and partition the full suite")
    plan_parser.add_argument("--shards", type=int, default=DEFAULT_SHARD_COUNT)
    plan_parser.add_argument("--expected", type=Path, default=DEFAULT_MANIFEST)
    plan_parser.add_argument("--output", type=Path, required=True)
    plan_parser.add_argument("--github-output", type=Path)
    plan_parser.set_defaults(handler=create_plan)

    run_parser = subparsers.add_parser("run", help="run one bounded pytest shard")
    run_parser.add_argument("--plan", type=Path, required=True)
    run_parser.add_argument("--shard-index", type=int, required=True)
    run_parser.add_argument("--timeout-seconds", type=int, default=300)
    run_parser.add_argument("--expected", type=Path, default=DEFAULT_MANIFEST)
    run_parser.add_argument("--output-dir", type=Path, required=True)
    run_parser.set_defaults(handler=run_shard)

    verify_parser = subparsers.add_parser("verify", help="verify all shard artifacts and outcomes")
    verify_parser.add_argument("--artifacts", type=Path, required=True)
    verify_parser.add_argument("--expected", type=Path, default=DEFAULT_MANIFEST)
    verify_parser.add_argument("--shard-count", type=int, default=DEFAULT_SHARD_COUNT)
    verify_parser.add_argument("--plan-status", required=True)
    verify_parser.add_argument("--shards-status", required=True)
    verify_parser.set_defaults(handler=verify_results)

    update_parser = subparsers.add_parser(
        "update-expected", help="explicitly regenerate the reviewed node ID manifest"
    )
    update_parser.add_argument("--expected", type=Path, default=DEFAULT_MANIFEST)
    update_parser.add_argument("--confirm", action="store_true", required=True)
    update_parser.set_defaults(handler=update_expected)
    return parser


def main(argv: Sequence[str] | None = None) -> int:
    args = build_parser().parse_args(argv)
    try:
        return args.handler(args)
    except (OSError, ValueError) as exc:
        print(f"pytest shard tooling error: {exc}", file=sys.stderr)
        return 2


if __name__ == "__main__":
    raise SystemExit(main())
