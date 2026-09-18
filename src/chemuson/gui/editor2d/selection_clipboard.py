"""Deterministic selection clipboard policy and payload codec."""

from __future__ import annotations

import json
from collections.abc import Iterable

MIME_SELECTION = "application/x-chemuson-selection"
MIME_TEXT_ITEMS = "application/x-chemuson-text-items"
MIME_MDL_MOLFILE = "chemical/x-mdl-molfile"
MIME_PNG = "image/png"
MIME_SVG = "image/svg+xml"


def mime_has_pasteable_format(mime: object) -> bool:
    """Return whether MIME data has a format accepted by ChemUSON."""
    return bool(
        mime.hasFormat(MIME_SELECTION)
        or mime.hasFormat(MIME_TEXT_ITEMS)
        or mime.hasFormat(MIME_MDL_MOLFILE)
        or mime.hasUrls()
        or mime.hasText()
        or mime.hasFormat(MIME_PNG)
        or mime.hasImage()
        or mime.hasFormat(MIME_SVG)
    )


def encode_selection_payload(payload: dict) -> bytes:
    """Encode the existing selection payload as UTF-8 JSON bytes."""
    return json.dumps(payload).encode("utf-8")


def decode_selection_payload(data: bytes) -> dict | None:
    """Decode a selection payload, returning None for malformed/non-object data."""
    try:
        payload = json.loads(bytes(data).decode("utf-8"))
    except (TypeError, UnicodeDecodeError, json.JSONDecodeError):
        return None
    return payload if isinstance(payload, dict) else None


def is_large_clipboard_structure(
    atom_count: int | None,
    bond_count: int | None,
    *,
    atom_threshold: int,
    bond_threshold: int,
) -> bool:
    """Apply the existing lightweight-export thresholds to structure counts."""
    if atom_count is None or bond_count is None:
        return False
    return atom_count >= atom_threshold or bond_count >= bond_threshold


def bond_copy_priority(bond: object) -> int:
    """Return the existing priority used when duplicate bond pairs are copied."""
    style = getattr(bond, "style", None)
    style_name = getattr(style, "value", style)
    style_bonus = 5 if style_name == "coordination" else 0
    aromatic_bonus = 50 if bond.is_aromatic else 0
    display_bonus = int(bond.display_order or 0)
    return int(bond.order or 1) * 10 + style_bonus + aromatic_bonus + display_bonus


def unique_bonds_for_copy(bonds: Iterable[object]) -> list[object]:
    """Deduplicate copied bonds by atom pair, retaining the highest priority."""
    unique: dict[tuple[int, int], object] = {}
    for bond in sorted(bonds, key=lambda item: item.id):
        pair = (min(int(bond.a1_id), int(bond.a2_id)), max(int(bond.a1_id), int(bond.a2_id)))
        existing = unique.get(pair)
        if existing is None or bond_copy_priority(bond) > bond_copy_priority(existing):
            unique[pair] = bond
    return list(unique.values())
