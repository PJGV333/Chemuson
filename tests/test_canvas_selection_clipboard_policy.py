"""Functional contracts for selection clipboard policy and codec."""

from __future__ import annotations

from dataclasses import dataclass
from types import SimpleNamespace

from chemuson.gui.canvas.selection_clipboard import (
    MIME_MDL_MOLFILE,
    MIME_PNG,
    MIME_SELECTION,
    MIME_SVG,
    MIME_TEXT_ITEMS,
    bond_copy_priority,
    decode_selection_payload,
    encode_selection_payload,
    is_large_clipboard_structure,
    mime_has_pasteable_format,
    unique_bonds_for_copy,
)


@dataclass
class Mime:
    formats: set[str]
    urls: bool = False
    text: bool = False
    image: bool = False

    def hasFormat(self, value: str) -> bool:
        return value in self.formats

    def hasUrls(self) -> bool:
        return self.urls

    def hasText(self) -> bool:
        return self.text

    def hasImage(self) -> bool:
        return self.image


def test_mime_policy_preserves_existing_formats() -> None:
    assert {MIME_SELECTION, MIME_TEXT_ITEMS, MIME_MDL_MOLFILE, MIME_PNG, MIME_SVG} == {
        "application/x-chemuson-selection",
        "application/x-chemuson-text-items",
        "chemical/x-mdl-molfile",
        "image/png",
        "image/svg+xml",
    }
    assert mime_has_pasteable_format(Mime({MIME_SELECTION}))
    assert mime_has_pasteable_format(Mime({MIME_TEXT_ITEMS}))
    assert mime_has_pasteable_format(Mime({MIME_MDL_MOLFILE}))
    assert mime_has_pasteable_format(Mime(set(), urls=True))
    assert mime_has_pasteable_format(Mime(set(), text=True))
    assert mime_has_pasteable_format(Mime({MIME_PNG}))
    assert mime_has_pasteable_format(Mime(set(), image=True))
    assert mime_has_pasteable_format(Mime({MIME_SVG}))
    assert not mime_has_pasteable_format(Mime(set()))


def test_payload_codec_roundtrips_json_and_rejects_invalid_data() -> None:
    payload = {"atoms": [{"id": 1}], "offset": [20.0, 20.0]}

    encoded = encode_selection_payload(payload)

    assert decode_selection_payload(encoded) == payload
    assert decode_selection_payload(b"not-json") is None
    assert decode_selection_payload(b"[]") is None


def test_large_clipboard_policy_preserves_threshold_boundaries() -> None:
    assert not is_large_clipboard_structure(17, 19, atom_threshold=18, bond_threshold=20)
    assert is_large_clipboard_structure(18, 0, atom_threshold=18, bond_threshold=20)
    assert is_large_clipboard_structure(0, 20, atom_threshold=18, bond_threshold=20)
    assert not is_large_clipboard_structure(None, None, atom_threshold=18, bond_threshold=20)


def test_duplicate_bond_policy_keeps_highest_priority_pair() -> None:
    plain = SimpleNamespace(id=1, a1_id=1, a2_id=2, order=1, style="plain", is_aromatic=False, display_order=0)
    aromatic = SimpleNamespace(id=2, a1_id=2, a2_id=1, order=1, style="plain", is_aromatic=True, display_order=0)

    assert bond_copy_priority(aromatic) > bond_copy_priority(plain)
    assert unique_bonds_for_copy([plain, aromatic]) == [aromatic]
