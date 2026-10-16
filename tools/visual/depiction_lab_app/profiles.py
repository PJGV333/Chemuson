"""Versioned JSON profiles isolated from normal ChemUSON preferences."""
from __future__ import annotations

import json
import os
from pathlib import Path
import re
from typing import Any, Mapping

from PyQt6.QtCore import QStandardPaths

from .parameters import (
    SCHEMA,
    SCHEMA_VERSION,
    DepictionStyleParameters,
    ParameterValidationError,
)

PROFILE_KEYS = {"schema", "version", "name", "based_on", "parameters"}
V0_PARAMETER_KEYS = {
    "stroke_ratio": "bond_stroke_ratio",
    "double_offset_ratio": "double_bond_offset_ratio",
    "inner_trim_ratio": "double_bond_inner_trim_ratio",
    "wedge_width_ratio": "solid_wedge_width_ratio",
    "hash_count": "hashed_wedge_count",
    "hash_min_width_ratio": "hashed_wedge_min_width_ratio",
    "hash_max_width_ratio": "hashed_wedge_max_width_ratio",
    "hash_stroke_ratio": "hashed_wedge_stroke_ratio",
    "cap_style": "cap_style",
    "join_style": "join_style",
}


class ProfileError(ValueError):
    """Profile file is invalid or incompatible."""


def profile_payload(
    parameters: DepictionStyleParameters,
    *,
    name: str = "Working",
    based_on: str = "current-renderer",
) -> dict[str, Any]:
    if not isinstance(name, str) or not name.strip():
        raise ProfileError("profile name must be a non-empty string")
    if not isinstance(based_on, str) or not based_on.strip():
        raise ProfileError("based_on must be a non-empty string")
    return {
        "schema": SCHEMA,
        "version": SCHEMA_VERSION,
        "name": str(name),
        "based_on": str(based_on),
        "parameters": parameters.to_mapping(),
    }


def _reject_constant(value: str) -> None:
    raise ProfileError(f"non-finite JSON constant is not allowed: {value}")


def _reject_duplicate_keys(pairs: list[tuple[str, Any]]) -> dict[str, Any]:
    result: dict[str, Any] = {}
    for key, value in pairs:
        if key in result:
            raise ProfileError(f"duplicate JSON field: {key}")
        result[key] = value
    return result


def _migrate_v0(payload: Mapping[str, Any]) -> dict[str, Any]:
    unknown_fields = sorted(set(payload) - PROFILE_KEYS)
    if unknown_fields:
        raise ProfileError(f"unknown version 0 profile field(s): {', '.join(unknown_fields)}")
    old_parameters = payload.get("parameters")
    if not isinstance(old_parameters, Mapping):
        raise ProfileError("version 0 parameters must be an object")
    unknown = sorted(set(old_parameters) - set(V0_PARAMETER_KEYS))
    if unknown:
        raise ProfileError(f"unknown version 0 parameter(s): {', '.join(unknown)}")
    values = DepictionStyleParameters.current().to_mapping()
    for old_name, value in old_parameters.items():
        values[V0_PARAMETER_KEYS[old_name]] = value
    name = payload.get("name", "Migrated profile")
    based_on = payload.get("based_on", "current-renderer")
    if not isinstance(name, str) or not name.strip():
        raise ProfileError("version 0 profile name must be a non-empty string")
    if not isinstance(based_on, str) or not based_on.strip():
        raise ProfileError("version 0 based_on must be a non-empty string")
    return profile_payload(
        DepictionStyleParameters.from_mapping(values),
        name=name,
        based_on=based_on,
    )


MIGRATIONS = {0: _migrate_v0}


def decode_profile_text(text: str) -> tuple[str, DepictionStyleParameters]:
    try:
        payload = json.loads(
            text,
            parse_constant=_reject_constant,
            object_pairs_hook=_reject_duplicate_keys,
        )
    except ProfileError:
        raise
    except (json.JSONDecodeError, TypeError) as exc:
        raise ProfileError(f"invalid profile JSON: {exc}") from exc
    if not isinstance(payload, dict):
        raise ProfileError("profile root must be a JSON object")
    if payload.get("schema") != SCHEMA:
        raise ProfileError(f"unsupported profile schema: {payload.get('schema')!r}")
    version = payload.get("version")
    if isinstance(version, bool) or not isinstance(version, int):
        raise ProfileError("profile version must be an integer")
    if version in MIGRATIONS:
        payload = MIGRATIONS[version](payload)
        version = payload["version"]
    if version != SCHEMA_VERSION:
        raise ProfileError(f"unsupported profile version: {version!r}")
    unknown = sorted(set(payload) - PROFILE_KEYS)
    missing = sorted(PROFILE_KEYS - set(payload))
    if unknown:
        raise ProfileError(f"unknown profile field(s): {', '.join(unknown)}")
    if missing:
        raise ProfileError(f"missing profile field(s): {', '.join(missing)}")
    if not isinstance(payload["name"], str) or not payload["name"].strip():
        raise ProfileError("profile name must be a non-empty string")
    if not isinstance(payload["based_on"], str) or not payload["based_on"].strip():
        raise ProfileError("based_on must be a non-empty string")
    try:
        parameters = DepictionStyleParameters.from_mapping(payload["parameters"])
    except ParameterValidationError as exc:
        raise ProfileError(str(exc)) from exc
    return payload["name"].strip(), parameters


def encode_profile(
    parameters: DepictionStyleParameters,
    *,
    name: str = "Working",
    based_on: str = "current-renderer",
) -> str:
    # Revalidate dataclass contents too; callers may construct it directly.
    validated = DepictionStyleParameters.from_mapping(parameters.to_mapping())
    return json.dumps(
        profile_payload(validated, name=name, based_on=based_on),
        ensure_ascii=False,
        indent=2,
        sort_keys=True,
        allow_nan=False,
    ) + "\n"


def profile_directory() -> Path:
    """Return the lab-only profile directory, never the normal settings store."""
    configured = QStandardPaths.writableLocation(
        QStandardPaths.StandardLocation.AppConfigLocation
    )
    base = Path(configured) if configured else Path.home() / ".config" / "chemuson"
    return base / "depiction-lab" / "profiles"


def save_profile(
    path: Path,
    parameters: DepictionStyleParameters,
    *,
    name: str = "Working",
) -> Path:
    target = Path(path)
    target.parent.mkdir(parents=True, exist_ok=True)
    text = encode_profile(parameters, name=name)
    temporary = target.with_name(target.name + ".tmp")
    temporary.write_text(text, encoding="utf-8")
    os.replace(temporary, target)
    return target


def load_profile(path: Path) -> tuple[str, DepictionStyleParameters]:
    try:
        text = Path(path).read_text(encoding="utf-8")
    except OSError as exc:
        raise ProfileError(f"cannot read profile {path}: {exc}") from exc
    return decode_profile_text(text)


def safe_profile_filename(name: str) -> str:
    stem = re.sub(r"[^A-Za-z0-9._-]+", "-", name.strip()).strip(".-_")
    if not stem:
        raise ProfileError("profile name must contain letters or numbers")
    return f"{stem}.json"


def compare_profiles(
    left: DepictionStyleParameters,
    right: DepictionStyleParameters,
) -> dict[str, dict[str, Any]]:
    return left.changed_values(right)
