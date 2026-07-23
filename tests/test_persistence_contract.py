"""Characterization tests for the CMSN persistence document surface."""

from __future__ import annotations

import json
from pathlib import Path

import pytest

from chemuson.chemio.persistence import PersistenceManager
from chemuson.core.model import MolGraph


class FakePersistenceDocument:
    """Pure document double with no GUI or Qt dependency."""

    def __init__(self) -> None:
        self.model = MolGraph()
        self.canvas_data = {"settings": {"show_grid": True}}
        self.loaded_canvas_data: dict[str, object] | None = None
        self.rebuild_count = 0
        self.rebuild_snapshot: tuple[int, int, int, int, dict[str, object] | None] | None = None

    def get_persistence_data(self) -> dict[str, object]:
        return self.canvas_data

    def load_persistence_data(self, data: dict[str, object]) -> None:
        self.loaded_canvas_data = data

    def rebuild_persistence_view(self) -> None:
        self.rebuild_count += 1
        self.rebuild_snapshot = (
            len(self.model.atoms),
            len(self.model.bonds),
            self.model._next_atom_id,
            self.model._next_bond_id,
            self.loaded_canvas_data,
        )


def _document_with_bond() -> FakePersistenceDocument:
    document = FakePersistenceDocument()
    atom_a = document.model.add_atom("C", 1.0, 2.0, formal_charge=1)
    atom_b = document.model.add_atom("O", 3.0, 4.0)
    document.model.add_bond(atom_a.id, atom_b.id, order=2)
    return document


def test_save_to_dict_preserves_minimal_payload_without_mutating_document() -> None:
    document = _document_with_bond()

    payload = PersistenceManager.save_to_dict(document)

    assert payload["application"] == "Chemuson"
    assert payload["version"] == PersistenceManager.VERSION
    assert len(payload["model"]["atoms"]) == 2
    assert len(payload["model"]["bonds"]) == 1
    assert payload["canvas"] == document.canvas_data
    assert len(document.model.atoms) == 2
    assert len(document.model.bonds) == 1


def test_load_restores_model_canvas_ids_then_rebuilds_once() -> None:
    source = _document_with_bond()
    payload = PersistenceManager.save_to_dict(source)
    target = FakePersistenceDocument()
    target.model.add_atom("N", 99.0, 99.0)

    PersistenceManager.load_from_dict(payload, target)

    assert target.loaded_canvas_data == source.canvas_data
    assert target.rebuild_count == 1
    assert target.rebuild_snapshot == (2, 1, 3, 2, source.canvas_data)
    assert target.model.get_atom(1).formal_charge == 1


def test_invalid_application_does_not_mutate_or_rebuild() -> None:
    document = _document_with_bond()

    with pytest.raises(ValueError, match="^Not a valid Chemuson file$"):
        PersistenceManager.load_from_dict({"application": "Other"}, document)

    assert len(document.model.atoms) == 2
    assert document.rebuild_count == 0


def test_save_to_file_writes_json_and_removes_temporary_file(tmp_path: Path) -> None:
    document = _document_with_bond()
    filepath = tmp_path / "document.cmsn"

    PersistenceManager.save_to_file(str(filepath), document)

    assert filepath.exists()
    assert not (tmp_path / "document.cmsn.tmp").exists()
    assert json.loads(filepath.read_text(encoding="utf-8")) == PersistenceManager.save_to_dict(document)


def test_save_to_file_does_not_replace_when_dump_fails(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    document = _document_with_bond()
    replaced = False

    def fail_dump(*args: object, **kwargs: object) -> None:
        raise OSError("write failed")

    def record_replace(*args: object) -> None:
        nonlocal replaced
        replaced = True

    monkeypatch.setattr(json, "dump", fail_dump)
    monkeypatch.setattr("chemuson.chemio.persistence.os.replace", record_replace)

    with pytest.raises(OSError, match="write failed"):
        PersistenceManager.save_to_file(str(tmp_path / "document.cmsn"), document)

    assert not replaced


def test_load_from_file_delegates_to_dictionary_load(tmp_path: Path) -> None:
    source = _document_with_bond()
    filepath = tmp_path / "document.cmsn"
    filepath.write_text(json.dumps(PersistenceManager.save_to_dict(source)), encoding="utf-8")
    target = FakePersistenceDocument()

    PersistenceManager.load_from_file(str(filepath), target)

    assert len(target.model.atoms) == 2
    assert target.rebuild_count == 1
