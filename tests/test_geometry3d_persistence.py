import math

from chemuson.chemio.persistence import PersistenceManager
from chemuson.geometry3d.model import CoordinateSet3D
from chemuson.gui.canvas import ChemusonCanvas


def _canvas_with_atom() -> tuple[ChemusonCanvas, int]:
    canvas = ChemusonCanvas()
    atom = canvas.model.add_atom("C", 15.0, 25.0)
    return canvas, atom.id


def test_geometry3d_round_trips_without_changing_2d_model():
    canvas, atom_id = _canvas_with_atom()
    canvas.coordinate_set_3d = CoordinateSet3D(
        {atom_id: (1.25, -2.5, 3.75)},
        source="rdkit",
        method="MMFF94",
        energy=-4.25,
        metadata={"forcefield": "MMFF94", "energy_unit": "kcal/mol", "generated_atoms": [{"id": "H1"}]},
    )

    payload = PersistenceManager.save_to_dict(canvas)
    saved_geometry = payload["canvas"]["geometry3d"]
    restored = ChemusonCanvas()
    PersistenceManager.load_from_dict(payload, restored)

    assert saved_geometry["schema_version"] == 1
    assert saved_geometry["positions"] == {str(atom_id): [1.25, -2.5, 3.75]}
    assert restored.coordinate_set_3d is not None
    assert restored.coordinate_set_3d.positions == {atom_id: (1.25, -2.5, 3.75)}
    assert restored.coordinate_set_3d.energy == -4.25
    assert restored.coordinate_set_3d.metadata.get("generated_atoms") is None
    assert restored.model.atoms[atom_id].x == 15.0
    assert restored.model.atoms[atom_id].y == 25.0


def test_geometry3d_legacy_files_load_without_inventing_coordinates():
    canvas, _atom_id = _canvas_with_atom()
    payload = PersistenceManager.save_to_dict(canvas)
    payload["canvas"].pop("geometry3d", None)
    restored = ChemusonCanvas()

    PersistenceManager.load_from_dict(payload, restored)

    assert restored.coordinate_set_3d is None
    assert restored.geometry3d_load_error is None
    assert len(restored.model.atoms) == 1


def test_invalid_geometry3d_is_dropped_without_losing_molecule():
    canvas, atom_id = _canvas_with_atom()
    payload = PersistenceManager.save_to_dict(canvas)
    payload["canvas"]["geometry3d"] = {
        "schema_version": 1,
        "source": "test",
        "method": "manual",
        "positions": {str(atom_id): [math.nan, 0.0, 0.0]},
    }
    restored = ChemusonCanvas()

    PersistenceManager.load_from_dict(payload, restored)

    assert restored.coordinate_set_3d is None
    assert restored.geometry3d_load_error == "invalid_geometry3d"
    assert len(restored.model.atoms) == 1


def test_geometry3d_save_filters_orphan_ids_retained_for_undo():
    canvas, atom_id = _canvas_with_atom()
    canvas.coordinate_set_3d = CoordinateSet3D(
        {atom_id: (0.0, 0.0, 0.0), 999: (8.0, 9.0, 10.0)},
        source="test",
        method="manual",
    )

    payload = PersistenceManager.save_to_dict(canvas)

    assert payload["canvas"]["geometry3d"]["positions"] == {str(atom_id): [0.0, 0.0, 0.0]}
