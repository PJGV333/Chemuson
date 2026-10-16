from chemuson.geometry3d.model import CoordinateSet3D
from chemuson.gui.canvas import ChemusonCanvas
from chemuson.gui.commands import SetGeometry3DCommand


def test_geometry3d_command_is_one_undoable_action_without_mutating_2d():
    canvas = ChemusonCanvas()
    atom = canvas.model.add_atom("C", 12.0, 24.0)
    before = CoordinateSet3D({atom.id: (0.0, 0.0, 0.0)}, "manual", "test")
    after = CoordinateSet3D({atom.id: (1.5, -0.5, 0.25)}, "manual", "drag")
    canvas.set_geometry3d_coordinates(before)
    revision_before = canvas.document_revision

    canvas.undo_stack.push(SetGeometry3DCommand(canvas, after))

    assert canvas.undo_stack.count() == 1
    assert canvas.coordinate_set_3d == after
    assert canvas.document_revision > revision_before
    assert (canvas.model.atoms[atom.id].x, canvas.model.atoms[atom.id].y) == (12.0, 24.0)

    canvas.undo_stack.undo()
    assert canvas.coordinate_set_3d == before
    assert (canvas.model.atoms[atom.id].x, canvas.model.atoms[atom.id].y) == (12.0, 24.0)

    canvas.undo_stack.redo()
    assert canvas.coordinate_set_3d == after
    assert (canvas.model.atoms[atom.id].x, canvas.model.atoms[atom.id].y) == (12.0, 24.0)


def test_geometry_state_identity_and_revision_are_per_document():
    first = ChemusonCanvas()
    second = ChemusonCanvas()

    assert first.document_token != second.document_token
    first_revision = first.document_revision
    first.set_geometry3d_coordinates(CoordinateSet3D({1: (0.0, 0.0, 0.0)}, "test", "manual"))

    assert first.document_revision > first_revision
    assert second.coordinate_set_3d is None
    assert second.document_revision == 0


def test_command_can_undo_first_geometry_assignment_to_none():
    canvas = ChemusonCanvas()
    atom = canvas.model.add_atom("N", 0.0, 0.0)
    coordinates = CoordinateSet3D({atom.id: (0.0, 1.0, 2.0)}, "manual", "build")

    canvas.undo_stack.push(SetGeometry3DCommand(canvas, coordinates))
    assert canvas.coordinate_set_3d == coordinates
    canvas.undo_stack.undo()

    assert canvas.coordinate_set_3d is None
