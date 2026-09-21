from __future__ import annotations

import subprocess
from pathlib import Path

ROOT = Path(__file__).parents[2]
CHANGE = ROOT / "openspec" / "changes" / "establish-clean2d-campaign-policy"
SPEC = CHANGE / "specs" / "clean-2d-campaign-policy" / "spec.md"
ROADMAP = ROOT / "docs" / "clean2d" / "CAMPAIGN.md"


NINE_CAMPAIGNS = (
    "Benchmark & observability",
    "Topology/decomposition",
    "Medium molecule assembly",
    "Rigid/fused/multiring layout",
    "Flexible connectors and branch routing",
    "Macrocycles and large structures",
    "Global candidate search/ranking",
    "Local polish",
    "Production acceptance",
)


def _text(path: Path) -> str:
    return path.read_text(encoding="utf-8")


def test_master_openspec_and_roadmap_exist() -> None:
    assert all((CHANGE / name).is_file() for name in ("proposal.md", "design.md", "tasks.md"))
    assert SPEC.is_file()
    assert ROADMAP.is_file()


def test_roadmap_contains_nine_campaigns() -> None:
    roadmap = _text(ROADMAP).lower()
    for campaign in NINE_CAMPAIGNS:
        assert campaign.lower() in roadmap


def test_policy_distinguishes_hard_gates_and_soft_metrics() -> None:
    policy = _text(SPEC).lower()
    assert "hard gates" in policy
    assert "soft metrics" in policy
    assert "shall not compensate" in policy


def test_policy_contains_complexity_posture_and_preserve_only() -> None:
    policy = _text(SPEC).lower()
    for term in ("simple", "medium", "large", "complex-scale", "preserve-only", "no-op"):
        assert term in policy


def test_policy_prohibits_molecule_specific_routing() -> None:
    policy = _text(SPEC).lower()
    assert "molecule-specific production routing is prohibited" in policy
    assert "known case ids" in policy
    assert "reusable signals" in policy


def test_policy_requires_baseline_before_algorithm_changes_and_independent_openspecs() -> None:
    policy = _text(SPEC).lower()
    assert "capture a baseline" in policy
    assert "separate openspec" in policy
    assert "before/after diff" in policy


def test_policy_declares_metric_vector_and_corpus_taxonomy() -> None:
    policy = _text(SPEC).lower()
    for metric in (
        "bond_length_error",
        "bond_crossing_count",
        "ring_degeneracy",
        "connector_congestion",
        "whitespace_balance",
    ):
        assert metric in policy
    for tag in ("macrocycle", "multiblock", "stereo-sensitive", "selection-boundary"):
        assert f"`{tag}`" in policy


def test_change_does_not_modify_production_clean2d_or_architecture_catalog() -> None:
    result = subprocess.run(
        ["git", "status", "--porcelain", "--untracked-files=all"],
        cwd=ROOT,
        check=True,
        capture_output=True,
        text=True,
    )
    changed = [line[3:] for line in result.stdout.splitlines() if len(line) >= 4]
    forbidden_prefixes = (
        "src/chemuson/clean2d/",
        "src/chemuson/gui/",
        "src/chemuson/core/",
        "src/chemuson/chemio/",
        "architecture/",
    )
    assert not [path for path in changed if path.startswith(forbidden_prefixes)]
