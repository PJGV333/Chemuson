from __future__ import annotations

import json

from tests.clean2d_regression.baselines import (
    build_baseline_record,
    canonicalize_baseline_record,
)
from tests.clean2d_regression.cases import get_regression_cases
from tests.clean2d_regression.metadata import (
    SIZE_CLASSES,
    derive_topology_metadata,
    validate_case_taxonomy,
)

EXPECTED_VECTOR_KEYS = {
    "bond_length_error",
    "bond_length_variance",
    "bond_angle_penalty",
    "atom_collision_count",
    "label_collision_count",
    "bond_crossing_count",
    "ring_distortion",
    "ring_degeneracy",
    "rigid_block_distortion",
    "branch_separation",
    "connector_congestion",
    "compactness",
    "whitespace_balance",
    "global_extent",
    "symmetry_preservation",
    "runtime_ms",
    "candidate_count",
    "candidate_sources",
    "result_state",
    "stable_reason",
}


def test_campaign1_cases_have_explicit_stable_taxonomy() -> None:
    cases = get_regression_cases()
    assert len({case.name for case in cases}) == len(cases)
    for case in cases:
        validate_case_taxonomy(case)
        assert case.size_class in SIZE_CLASSES
        assert case.family
        assert case.tags


def test_campaign1_topology_metadata_is_reproducible_and_json_safe() -> None:
    for case in get_regression_cases():
        first = derive_topology_metadata(case.builder())
        second = derive_topology_metadata(case.builder())
        assert first == second
        json.dumps(first, allow_nan=False, sort_keys=True)
        assert first["atom_count"] > 0
        assert first["heavy_atom_count"] <= first["atom_count"]
        assert first["bond_count"] > 0
        assert first["ring_count"] >= 0
        assert first["connected_components"] > 0


def test_campaign1_baseline_record_contains_complete_observability_vector() -> None:
    record = canonicalize_baseline_record(build_baseline_record(get_regression_cases()[0]))

    assert set(record["metric_vector"]) == EXPECTED_VECTOR_KEYS
    assert record["topology"]["atom_count"] > 0
    assert record["metric_vector"]["candidate_count"] == len(record["candidate_sources"])
    assert record["metric_vector"]["candidate_sources"] == record["candidate_sources"]
    assert record["metric_vector"]["result_state"] == record["result_state"]
    assert record["metric_vector"]["stable_reason"] == record["stable_reason"]
    assert record["runtime_ms"] >= 0
    json.dumps(record, allow_nan=False, sort_keys=True)


def test_campaign1_runtime_is_evidence_but_not_report_identity() -> None:
    first = canonicalize_baseline_record(build_baseline_record(get_regression_cases()[0]))
    second = canonicalize_baseline_record(first)
    second["runtime_ms"] = first["runtime_ms"] + 1000.0
    second["metric_vector"]["runtime_ms"] = second["runtime_ms"]

    from tests.clean2d_regression.baselines import baseline_records_equivalent

    assert baseline_records_equivalent(first, second)
