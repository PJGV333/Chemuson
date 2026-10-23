from __future__ import annotations

import json
from pathlib import Path
from types import SimpleNamespace

from tools.ci_pytest_timing_plugin import ShardTimingPlugin
from tools.ci_pytest_shards import (
    DEFAULT_DURATION_PROFILE,
    _load_duration_profile,
    _partition_by_file,
    _read_timing_records,
    _validate_plan,
)


def test_duration_profile_is_pinned_to_the_timeout_baseline() -> None:
    profile, digest = _load_duration_profile(DEFAULT_DURATION_PROFILE)

    assert profile["source"]["run_id"] == "39697199469"
    assert profile["source"]["commit"] == "8c28d47c0518f4ca40619839ae3504df8cf87096"
    assert profile["source"]["supplemental_run_ids"] == [
        "39679607676",
        "38996089270",
    ]
    assert profile["fallback_seconds_per_test"] == 0.2
    assert len(digest) == 64

    manifest_files = {
        nodeid.split("::", 1)[0]
        for nodeid in (Path(__file__).parent / "ci/expected_pytest_nodeids.txt")
        .read_text(encoding="utf-8")
        .splitlines()
    }
    assert set(profile["files"]) == manifest_files
    assert len(profile["files"]) == 254
    assert profile["files"]["tests/test_clean2d_baseline_reports.py"]["seconds_per_test"] > 5.0
    assert profile["files"]["tests/test_beta_save_shortcut_menu_access.py"]["seconds_per_test"] > 4.64
    assert profile["files"]["tests/test_ui_polish.py"]["seconds_per_test"] > 2.98
    assert profile["files"]["tests/test_ui_svg_icons.py"]["seconds_per_test"] > 1.21
    assert profile["files"]["tests/test_beta_save_shortcut_menu_access.py"]["observations"] == 40
    assert profile["files"]["tests/test_ui_polish.py"]["observations"] == 76


def test_weighted_lpt_is_deterministic_and_keeps_each_file_indivisible() -> None:
    nodeids = [
        *(f"tests/test_slow_a.py::test_{index}" for index in range(4)),
        *(f"tests/test_slow_b.py::test_{index}" for index in range(3)),
        *(f"tests/test_fast_c.py::test_{index}" for index in range(2)),
        *(f"tests/test_fast_d.py::test_{index}" for index in range(2)),
        "tests/test_fast_e.py::test_one",
    ]
    profile = {
        "fallback_seconds_per_test": 0.5,
        "files": {
            "tests/test_slow_a.py": {"seconds_per_test": 10.0},
            "tests/test_slow_b.py": {"seconds_per_test": 9.0},
            "tests/test_fast_c.py": {"seconds_per_test": 1.0},
            "tests/test_fast_d.py": {"seconds_per_test": 1.0},
        },
    }

    shards, estimates = _partition_by_file(nodeids, 2, profile)
    repeated, repeated_estimates = _partition_by_file(nodeids, 2, profile)

    assert shards == repeated
    assert estimates == repeated_estimates == [40.0, 31.5]
    assert set(shards[0]).isdisjoint(shards[1])
    assert sorted(node for shard in shards for node in shard) == sorted(nodeids)
    assert any(node.startswith("tests/test_slow_a.py::") for node in shards[0])
    assert any(node.startswith("tests/test_slow_b.py::") for node in shards[1])
    for test_file in ("test_slow_a.py", "test_slow_b.py", "test_fast_c.py", "test_fast_d.py", "test_fast_e.py"):
        owners = [index for index, shard in enumerate(shards) if any(test_file in node for node in shard)]
        assert len(owners) == 1


def test_plan_requires_profile_provenance_and_consistent_estimates() -> None:
    expected = ["tests/test_one.py::test_a", "tests/test_two.py::test_b"]
    plan = {
        "schema_version": 2,
        "total_tests": 2,
        "expected_nodeids_sha256": "unused-by-this-assertion",
        "duration_profile": {"path": "tests/ci/profile.json", "sha256": "a" * 64},
        "shard_count": 2,
        "distribution": [1, 1],
        "estimated_duration_seconds": [1.0, 2.0],
        "shards": [
            {"index": 0, "count": 1, "nodeids": [expected[0]], "estimated_duration_seconds": 1.0},
            {"index": 1, "count": 1, "nodeids": [expected[1]], "estimated_duration_seconds": 2.0},
        ],
    }
    from tools.ci_pytest_shards import _digest_nodeids

    plan["expected_nodeids_sha256"] = _digest_nodeids(expected)
    _validate_plan(plan, expected)
    plan["shards"][1]["estimated_duration_seconds"] = 3.0
    try:
        _validate_plan(plan, expected)
    except ValueError as exc:
        assert "duration estimate" in str(exc)
    else:
        raise AssertionError("inconsistent duration estimate was accepted")


def test_timing_plugin_flushes_complete_test_phases_incrementally(tmp_path: Path) -> None:
    output = tmp_path / "pytest-shard-0-timings.jsonl"
    plugin = ShardTimingPlugin(output)
    nodeid = "tests/test_ui.py::test_render"
    for phase, duration in (("setup", 0.1), ("call", 0.4), ("teardown", 0.05)):
        plugin.pytest_runtest_logreport(
            SimpleNamespace(
                nodeid=nodeid,
                when=phase,
                duration=duration,
                failed=False,
                skipped=False,
            )
        )
        if phase != "teardown":
            assert output.read_text(encoding="utf-8") == ""

    records, per_file = _read_timing_records(output, [nodeid])
    assert records[0]["duration_seconds"] == 0.55
    assert per_file == {"tests/test_ui.py": 0.55}
    assert json.loads(output.read_text(encoding="utf-8")) == records[0]


def test_timing_reader_rejects_duplicate_or_unassigned_node_ids(tmp_path: Path) -> None:
    path = tmp_path / "timings.jsonl"
    record = {
        "nodeid": "tests/test_ui.py::test_render",
        "test_file": "tests/test_ui.py",
        "duration_seconds": 0.5,
        "outcome": "passed",
        "phases": {"call": 0.5},
    }
    path.write_text(json.dumps(record) + "\n" + json.dumps(record) + "\n", encoding="utf-8")

    try:
        _read_timing_records(path, [record["nodeid"]])
    except ValueError as exc:
        assert "duplicate timing record" in str(exc)
    else:
        raise AssertionError("duplicate timing record was accepted")


def test_timing_reader_preserves_flushed_records_before_a_truncated_final_line(tmp_path: Path) -> None:
    path = tmp_path / "timings.jsonl"
    record = {
        "nodeid": "tests/test_ui.py::test_render",
        "test_file": "tests/test_ui.py",
        "duration_seconds": 0.5,
        "outcome": "passed",
        "phases": {"call": 0.5},
    }
    path.write_text(json.dumps(record) + "\n{\"nodeid\":", encoding="utf-8")

    records, per_file = _read_timing_records(path, [record["nodeid"]])
    assert [item["nodeid"] for item in records] == [record["nodeid"]]
    assert per_file == {"tests/test_ui.py": 0.5}
