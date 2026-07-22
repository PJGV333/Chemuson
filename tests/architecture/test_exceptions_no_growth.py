"""Freeze the reviewed temporary-exception identities without importing production."""

from __future__ import annotations

from collections import Counter
from dataclasses import dataclass
from pathlib import Path, PurePosixPath
import re
from typing import Mapping

import pytest
import yaml


REPO_ROOT = Path(__file__).resolve().parent.parent.parent
CATALOG_PATH = REPO_ROOT / "architecture" / "modules.yml"


@dataclass(frozen=True, order=True)
class ExceptionIdentity:
    source_id: str
    target_id: str
    file: str
    import_path: str
    type_checking_only: bool


@dataclass(frozen=True)
class ExceptionGrowthAudit:
    current_identities: tuple[ExceptionIdentity, ...]
    unexpected_identities: tuple[ExceptionIdentity, ...]
    removed_baseline_identities: tuple[ExceptionIdentity, ...]
    duplicate_identities: tuple[tuple[ExceptionIdentity, int], ...]
    baseline_duplicate_identities: tuple[tuple[ExceptionIdentity, int], ...]
    current_counts_by_module: tuple[tuple[str, int], ...]
    baseline_counts_by_module: tuple[tuple[str, int], ...]
    modules_with_growth: tuple[tuple[str, int, int], ...]


def normalize_exception_path(value: str) -> str:
    if not isinstance(value, str):
        raise TypeError("exception file must be a string")
    return PurePosixPath(value.strip().replace("\\", "/")).as_posix()


def normalize_exception_import(value: str) -> str:
    if not isinstance(value, str):
        raise TypeError("exception import_path must be a string")
    return re.sub(r"\s+", " ", value).strip()


def exception_identity_from_mapping(exc: Mapping[str, object]) -> ExceptionIdentity:
    source_id = exc["source_id"]
    target_id = exc["target_id"]
    type_checking_only = exc["type_checking_only"]
    if not isinstance(source_id, str) or not isinstance(target_id, str):
        raise TypeError("exception module IDs must be strings")
    if not isinstance(type_checking_only, bool):
        raise TypeError("exception type_checking_only must be a bool")
    return ExceptionIdentity(
        source_id=source_id,
        target_id=target_id,
        file=normalize_exception_path(exc["file"]),
        import_path=normalize_exception_import(exc["import_path"]),
        type_checking_only=type_checking_only,
    )


def _identity(
    source_id: str,
    target_id: str,
    file: str,
    import_path: str,
    type_checking_only: bool,
) -> ExceptionIdentity:
    return exception_identity_from_mapping(
        {
            "source_id": source_id,
            "target_id": target_id,
            "file": file,
            "import_path": import_path,
            "type_checking_only": type_checking_only,
        }
    )


# Expanding this baseline requires an explicit architectural decision and reviewed change.
FROZEN_EXCEPTION_BASELINE_ROWS: tuple[ExceptionIdentity, ...] = (
    _identity("M01", "M09", "src/chemuson/chemio/persistence.py", "from chemuson.gui.canvas import ChemusonCanvas", True),
)
ELIMINATED_M01_M02_EXCEPTIONS = (
    _identity("M01", "M02", "src/chemuson/chemio/depiction_candidates.py", "from chemuson.clean2d.geometry import count_crossings, cycle_basis, segments_intersect", False),
    _identity("M01", "M02", "src/chemuson/chemio/depiction_candidates.py", "from chemuson.clean2d.safety import has_cycles, min_nonbonded_distance, ring_degeneracy_score", False),
    _identity("M01", "M02", "src/chemuson/chemio/rdkit_io.py", "from chemuson.clean2d.geometry import apply_coords_in_place", False),
    _identity("M01", "M02", "src/chemuson/chemio/rdkit_io.py", "from chemuson.clean2d.scaffold_depiction import scaffold_depiction_candidates", False),
    _identity("M01", "M02", "src/chemuson/chemio/rdkit_io.py", "from chemuson.clean2d.block_unwrap import block_unwrap_layout", False),
)
FROZEN_EXCEPTION_BASELINE = frozenset(FROZEN_EXCEPTION_BASELINE_ROWS)
ELIMINATED_M15_EXCEPTIONS = (
    _identity("M15", "M01", "src/chemuson/utils/autosave.py", "from chemuson.chemio.persistence import PersistenceManager", False),
    _identity("M15", "M09", "src/chemuson/utils/autosave.py", "from chemuson.gui.canvas import ChemusonCanvas", True),
)
ELIMINATED_M03_EXCEPTIONS = (
    _identity("M03", "M04", "src/chemuson/chemcalc/formula.py", "from chemuson.chemname.molview import MolView", False),
    _identity("M03", "M04", "src/chemuson/chemcalc/valence.py", "from chemuson.chemname.molview import MolView", False),
)


def _counts_by_source(rows: tuple[ExceptionIdentity, ...]) -> tuple[tuple[str, int], ...]:
    return tuple(sorted(Counter(row.source_id for row in rows).items()))


def _duplicates(rows: tuple[ExceptionIdentity, ...]) -> tuple[tuple[ExceptionIdentity, int], ...]:
    return tuple(sorted((identity, count) for identity, count in Counter(rows).items() if count > 1))


def audit_exception_growth(
    modules: list[dict],
    baseline_rows: tuple[ExceptionIdentity, ...] = FROZEN_EXCEPTION_BASELINE_ROWS,
) -> ExceptionGrowthAudit:
    current_rows = tuple(
        exception_identity_from_mapping(exc)
        for module in modules
        for exc in module["temporary_exceptions"]
    )
    current_set = frozenset(current_rows)
    baseline_set = frozenset(baseline_rows)
    current_counts = _counts_by_source(current_rows)
    baseline_counts = _counts_by_source(baseline_rows)
    current_by_source = dict(current_counts)
    baseline_by_source = dict(baseline_counts)
    growth = tuple(
        (source_id, current_count, baseline_by_source.get(source_id, 0))
        for source_id, current_count in sorted(current_by_source.items())
        if current_count > baseline_by_source.get(source_id, 0)
    )
    return ExceptionGrowthAudit(
        current_identities=tuple(sorted(current_rows)),
        unexpected_identities=tuple(sorted(current_set - baseline_set)),
        removed_baseline_identities=tuple(sorted(baseline_set - current_set)),
        duplicate_identities=_duplicates(current_rows),
        baseline_duplicate_identities=_duplicates(baseline_rows),
        current_counts_by_module=current_counts,
        baseline_counts_by_module=baseline_counts,
        modules_with_growth=growth,
    )


def _mapping(identity: ExceptionIdentity, **metadata: object) -> dict[str, object]:
    return {
        "source_id": identity.source_id,
        "target_id": identity.target_id,
        "file": identity.file,
        "import_path": identity.import_path,
        "type_checking_only": identity.type_checking_only,
        "reason": metadata.get("reason", "reviewed debt"),
        "debt_ref": metadata.get("debt_ref", "test-debt"),
        "elimination_condition": metadata.get("elimination_condition", "remove dependency"),
    }


def _modules_for(*identities: ExceptionIdentity, metadata: dict[int, dict[str, object]] | None = None) -> list[dict]:
    by_source: dict[str, list[dict[str, object]]] = {}
    for index, identity in enumerate(identities):
        by_source.setdefault(identity.source_id, []).append(_mapping(identity, **(metadata or {}).get(index, {})))
    return [{"id": source_id, "temporary_exceptions": entries} for source_id, entries in by_source.items()]


def _changed(identity: ExceptionIdentity, field: str, value: object) -> ExceptionIdentity:
    values = {
        "source_id": identity.source_id,
        "target_id": identity.target_id,
        "file": identity.file,
        "import_path": identity.import_path,
        "type_checking_only": identity.type_checking_only,
    }
    values[field] = value
    return _identity(**values)  # type: ignore[arg-type]


def test_frozen_baseline_is_complete_and_exact() -> None:
    assert len(FROZEN_EXCEPTION_BASELINE_ROWS) == 1
    assert len(FROZEN_EXCEPTION_BASELINE) == 1
    assert _duplicates(FROZEN_EXCEPTION_BASELINE_ROWS) == ()
    assert _counts_by_source(FROZEN_EXCEPTION_BASELINE_ROWS) == (("M01", 1),)
    assert sum(not row.type_checking_only for row in FROZEN_EXCEPTION_BASELINE_ROWS) == 0
    assert sum(row.type_checking_only for row in FROZEN_EXCEPTION_BASELINE_ROWS) == 1
    for row in FROZEN_EXCEPTION_BASELINE_ROWS:
        assert row.file == normalize_exception_path(row.file)
        assert row.import_path == normalize_exception_import(row.import_path)
        assert not PurePosixPath(row.file).is_absolute()
        assert "*" not in row.file


def test_real_catalog_matches_frozen_exception_baseline() -> None:
    modules = yaml.safe_load(CATALOG_PATH.read_text(encoding="utf-8"))["modules"]
    audit = audit_exception_growth(modules)
    runtime = sum(not row.type_checking_only for row in audit.current_identities)
    type_checking = sum(row.type_checking_only for row in audit.current_identities)
    problems: list[str] = []
    for identity in audit.unexpected_identities:
        problems.append(f"{identity} | not present in frozen exception baseline")
    for identity, count in audit.duplicate_identities:
        problems.append(f"{identity} | count={count} | duplicate normalized temporary_exception")
    for source_id, current_count, baseline_count in audit.modules_with_growth:
        problems.append(f"{source_id} | current={current_count} | baseline={baseline_count} | temporary_exception count exceeds frozen baseline")
    if audit.baseline_duplicate_identities:
        problems.append(f"frozen baseline duplicates: {audit.baseline_duplicate_identities}")
    problems.extend(
        [
            f"current counts: {audit.current_counts_by_module}",
            f"baseline counts: {audit.baseline_counts_by_module}",
            f"removed baseline identities: {audit.removed_baseline_identities}",
        ]
    )
    assert not audit.unexpected_identities and not audit.duplicate_identities and not audit.modules_with_growth and len(audit.current_identities) == 1 and audit.current_counts_by_module == (("M01", 1),) and runtime == 0 and type_checking == 1, "\n".join(problems)


def test_eliminated_m01_m02_exceptions_cannot_reappear() -> None:
    audit = audit_exception_growth(
        _modules_for(*FROZEN_EXCEPTION_BASELINE_ROWS, *ELIMINATED_M01_M02_EXCEPTIONS)
    )

    assert audit.unexpected_identities == tuple(sorted(ELIMINATED_M01_M02_EXCEPTIONS))
    assert audit.modules_with_growth == (("M01", 6, 1),)
    assert set(audit.unexpected_identities) == set(ELIMINATED_M01_M02_EXCEPTIONS)


@pytest.mark.parametrize("eliminated", ELIMINATED_M01_M02_EXCEPTIONS)
def test_single_eliminated_m01_m02_exception_is_rejected(eliminated: ExceptionIdentity) -> None:
    audit = audit_exception_growth(_modules_for(*FROZEN_EXCEPTION_BASELINE_ROWS, eliminated))

    assert audit.unexpected_identities == (eliminated,)
    assert audit.modules_with_growth == (("M01", 2, 1),)


def test_replacing_persistence_exception_with_m01_m02_debt_is_rejected() -> None:
    eliminated = ELIMINATED_M01_M02_EXCEPTIONS[0]
    audit = audit_exception_growth(_modules_for(eliminated))

    assert audit.unexpected_identities == (eliminated,)
    assert audit.removed_baseline_identities == FROZEN_EXCEPTION_BASELINE_ROWS
    assert len(audit.current_identities) == len(FROZEN_EXCEPTION_BASELINE_ROWS)
    assert not audit.modules_with_growth


def test_normalized_eliminated_m01_m02_identity_is_still_rejected() -> None:
    eliminated = ELIMINATED_M01_M02_EXCEPTIONS[0]
    normalized_variant = _identity(
        eliminated.source_id,
        eliminated.target_id,
        eliminated.file.replace("/", "\\"),
        "from chemuson.clean2d.geometry import count_crossings,\n  cycle_basis, segments_intersect",
        eliminated.type_checking_only,
    )
    audit = audit_exception_growth(_modules_for(*FROZEN_EXCEPTION_BASELINE_ROWS, normalized_variant))

    assert audit.unexpected_identities == (eliminated,)
    assert audit.modules_with_growth == (("M01", 2, 1),)


def test_eliminated_m15_exceptions_cannot_reappear() -> None:
    audit = audit_exception_growth(_modules_for(*FROZEN_EXCEPTION_BASELINE_ROWS, *ELIMINATED_M15_EXCEPTIONS))

    assert audit.unexpected_identities == ELIMINATED_M15_EXCEPTIONS
    assert audit.modules_with_growth == (("M15", 2, 0),)


def test_eliminated_m03_exceptions_cannot_reappear() -> None:
    audit = audit_exception_growth(_modules_for(*FROZEN_EXCEPTION_BASELINE_ROWS, *ELIMINATED_M03_EXCEPTIONS))

    assert audit.unexpected_identities == ELIMINATED_M03_EXCEPTIONS
    assert audit.modules_with_growth == (("M03", 2, 0),)


def test_equal_empty_and_reduced_catalogs_are_no_growth() -> None:
    equal = audit_exception_growth(_modules_for(*FROZEN_EXCEPTION_BASELINE_ROWS))
    empty = audit_exception_growth([])
    reduced = audit_exception_growth(_modules_for(*FROZEN_EXCEPTION_BASELINE_ROWS[:-1]))
    assert not equal.unexpected_identities and not equal.duplicate_identities and not equal.modules_with_growth and not equal.removed_baseline_identities
    assert not empty.unexpected_identities and len(empty.removed_baseline_identities) == 1
    assert not reduced.unexpected_identities and len(reduced.removed_baseline_identities) == 1 and not reduced.modules_with_growth


def test_new_identities_and_replacement_are_detected() -> None:
    new_m01 = _identity("M01", "M02", "src/new.py", "from chemuson.clean2d import new", False)
    new_m02 = _identity("M02", "M01", "src/new.py", "from chemuson.chemio import new", False)
    audit = audit_exception_growth(_modules_for(*FROZEN_EXCEPTION_BASELINE_ROWS, new_m01, new_m02))
    replacement = audit_exception_growth(_modules_for(*FROZEN_EXCEPTION_BASELINE_ROWS[1:], new_m01))
    assert audit.unexpected_identities == (new_m01, new_m02)
    assert audit.modules_with_growth == (("M01", 2, 1), ("M02", 1, 0))
    assert replacement.unexpected_identities == (new_m01,)
    assert replacement.removed_baseline_identities == (FROZEN_EXCEPTION_BASELINE_ROWS[0],)
    assert not replacement.modules_with_growth


@pytest.mark.parametrize(
    ("field", "value"),
    [
        ("source_id", "M02"),
        ("target_id", "M03"),
        ("file", "src/chemuson/chemio/other.py"),
        ("import_path", "from chemuson.clean2d.geometry import other"),
        ("type_checking_only", False),
    ],
)
def test_each_identity_field_change_is_unexpected(field: str, value: object) -> None:
    original = FROZEN_EXCEPTION_BASELINE_ROWS[0]
    changed = _changed(original, field, value)
    audit = audit_exception_growth(_modules_for(*FROZEN_EXCEPTION_BASELINE_ROWS[1:], changed))
    assert audit.unexpected_identities == (changed,), field
    assert audit.removed_baseline_identities == (original,)


@pytest.mark.parametrize(
    "metadata",
    [
        {0: {"reason": "new reason"}},
        {0: {"debt_ref": "new-ref"}},
        {0: {"elimination_condition": "new condition"}},
        {0: {"reason": "new reason", "debt_ref": "new-ref", "elimination_condition": "new condition"}},
    ],
)
def test_metadata_changes_do_not_change_identity(metadata: dict[int, dict[str, object]]) -> None:
    audit = audit_exception_growth(_modules_for(*FROZEN_EXCEPTION_BASELINE_ROWS, metadata=metadata))
    assert not audit.unexpected_identities and not audit.duplicate_identities


@pytest.mark.parametrize(
    ("field", "value"),
    [
        ("import_path", "\tfrom chemuson.gui.canvas\n import ChemusonCanvas "),
        ("file", r"src\chemuson\chemio\persistence.py"),
        ("file", "./src/chemuson/chemio/persistence.py"),
        ("file", " src/chemuson/chemio/persistence.py "),
    ],
)
def test_normalized_identity_forms_match_baseline(field: str, value: str) -> None:
    changed = _changed(FROZEN_EXCEPTION_BASELINE_ROWS[0], field, value)
    audit = audit_exception_growth(_modules_for(*FROZEN_EXCEPTION_BASELINE_ROWS[1:], changed))
    assert not audit.unexpected_identities and not audit.duplicate_identities


@pytest.mark.parametrize(
    "duplicate",
    [
        FROZEN_EXCEPTION_BASELINE_ROWS[0],
        _changed(FROZEN_EXCEPTION_BASELINE_ROWS[0], "import_path", "from chemuson.gui.canvas\n import ChemusonCanvas"),
        _changed(FROZEN_EXCEPTION_BASELINE_ROWS[0], "file", r"src\chemuson\chemio\persistence.py"),
    ],
)
def test_duplicates_after_normalization_are_detected(duplicate: ExceptionIdentity) -> None:
    audit = audit_exception_growth(_modules_for(*FROZEN_EXCEPTION_BASELINE_ROWS, duplicate))
    assert audit.duplicate_identities == ((FROZEN_EXCEPTION_BASELINE_ROWS[0], 2),)
    assert audit.modules_with_growth == (("M01", 2, 1),)


def test_multiple_growth_order_baseline_duplicates_and_metadata_duplicates() -> None:
    new_rows = (
        _identity("M02", "M01", "src/two.py", "from chemuson.chemio import two", False),
        _identity("M02", "M01", "src/one.py", "from chemuson.chemio import one", False),
        _identity("M04", "M03", "src/three.py", "from chemuson.chemcalc import three", False),
    )
    audit = audit_exception_growth(_modules_for(*FROZEN_EXCEPTION_BASELINE_ROWS, *reversed(new_rows)))
    duplicate_metadata = audit_exception_growth(_modules_for(FROZEN_EXCEPTION_BASELINE_ROWS[0], FROZEN_EXCEPTION_BASELINE_ROWS[0], metadata={1: {"reason": "different"}}))
    duplicate_baseline = audit_exception_growth([], (FROZEN_EXCEPTION_BASELINE_ROWS[0], FROZEN_EXCEPTION_BASELINE_ROWS[0]))
    assert audit.unexpected_identities == tuple(sorted(new_rows))
    assert audit.modules_with_growth == (("M02", 2, 0), ("M04", 1, 0))
    assert duplicate_metadata.duplicate_identities == ((FROZEN_EXCEPTION_BASELINE_ROWS[0], 2),)
    assert duplicate_baseline.baseline_duplicate_identities == ((FROZEN_EXCEPTION_BASELINE_ROWS[0], 2),)


def test_order_and_simultaneous_removal_addition_are_deterministic() -> None:
    added = _identity("M02", "M01", "src/new.py", "from chemuson.chemio import new", False)
    rows = FROZEN_EXCEPTION_BASELINE_ROWS[1:] + (added,)
    forward = audit_exception_growth(_modules_for(*rows))
    reverse = audit_exception_growth(_modules_for(*reversed(rows)))
    assert forward == reverse
    assert forward.unexpected_identities == (added,)
    assert forward.removed_baseline_identities == (FROZEN_EXCEPTION_BASELINE_ROWS[0],)
