from __future__ import annotations

from unittest.mock import MagicMock, patch


from chemuson.clean2d import Clean2DCandidate, Clean2DMode, Clean2DResult
from chemuson.core.model import MolGraph
from chemuson.gui.controllers.clean2d_controller import Clean2DController


def _make_chain_graph() -> tuple[MolGraph, list[int]]:
    graph = MolGraph()
    a1 = graph.add_atom("C", 0.0, 0.0, atom_id=1)
    a2 = graph.add_atom("C", 40.0, 0.0, atom_id=2)
    a3 = graph.add_atom("C", 80.0, 0.0, atom_id=3)
    b12 = graph.add_bond(a1.id, a2.id, order=1)
    b23 = graph.add_bond(a2.id, a3.id, order=1)
    return graph, [b12.id, b23.id]


def _mock_window(graph: MolGraph) -> MagicMock:
    window = MagicMock()
    window.canvas.model = graph
    window.canvas.state.bond_length = 40.0
    window.canvas.undo_stack = MagicMock()
    window.canvas._update_selection_overlay = MagicMock()
    window.canvas._build_selection_graph = MagicMock(return_value=graph)
    window.statusBar.return_value = MagicMock()
    return window


def test_controller_rebuilds_authoritative_bonds_from_model() -> None:
    graph, bond_ids = _make_chain_graph()
    controller = Clean2DController()
    window = _mock_window(graph)
    partial_bond = graph.get_bond(bond_ids[0])
    full_bonds = [graph.get_bond(bond_ids[0]), graph.get_bond(bond_ids[1])]

    window.canvas._selected_structure_ids.return_value = ({1, 2}, [])
    window.canvas._structural_bonds_for_clean2d = MagicMock(return_value=full_bonds)
    result = Clean2DResult(
        mode=Clean2DMode.QUICK,
        atom_ids={1, 2},
        before_coords={1: (0.0, 0.0), 2: (40.0, 0.0)},
        candidates=(),
        selected=None,
        message="sin candidato",
    )

    with patch("chemuson.gui.controllers.clean2d_controller.run_clean2d_engine", return_value=result):
        controller.run_clean_2d(window, 1.0, 200, "(test)")

    window.canvas._build_selection_graph.assert_called_once()
    called_atom_ids, called_bonds = window.canvas._build_selection_graph.call_args.args
    assert called_atom_ids == {1, 2}
    assert called_bonds == full_bonds
    assert partial_bond in called_bonds


def test_validate_canvas_bond_integrity_rejects_boundary_bond_break() -> None:
    graph, bond_ids = _make_chain_graph()
    controller = Clean2DController()
    canvas = _mock_window(graph).canvas
    canvas._structural_bonds_for_clean2d = MagicMock(return_value=[graph.get_bond(bond_ids[0])])
    before = {1: (0.0, 0.0), 2: (40.0, 0.0)}
    after = {1: (0.0, 0.0), 2: (160.0, 0.0)}

    validation = controller.validate_canvas_bond_integrity_before_apply(
        canvas,
        {1, 2},
        before,
        after,
        40.0,
        "quick",
        candidate_source="internal_templates",
    )

    assert not validation.accepted
    assert validation.reason == "reparacion_rechazada_por_integridad_de_enlaces"
    assert validation.diagnostics["stable_rejection_reason"] == "boundary-bond-risk"
    assert validation.bond_integrity_regressions
    assert validation.real_bond_count == 2


def test_run_clean_2d_rejects_candidate_that_breaks_boundary_bond() -> None:
    graph, bond_ids = _make_chain_graph()
    controller = Clean2DController()
    window = _mock_window(graph)
    window.canvas._selected_structure_ids.return_value = ({1, 2}, [])
    window.canvas._structural_bonds_for_clean2d = MagicMock(return_value=[graph.get_bond(bond_ids[0])])

    candidate = Clean2DCandidate(
        source="internal_templates",
        coords={1: (0.0, 0.0), 2: (160.0, 0.0)},
        message="Estructura 2D limpiada (motor interno)",
    )
    result = Clean2DResult(
        mode=Clean2DMode.QUICK,
        atom_ids={1, 2},
        before_coords={1: (0.0, 0.0), 2: (40.0, 0.0)},
        candidates=(candidate,),
        selected=candidate,
        message=candidate.message,
    )

    attempt = controller._apply_engine_result(
        window,
        graph,
        {1, 2},
        [graph.get_bond(bond_ids[0])],
        {1: (0.0, 0.0), 2: (40.0, 0.0)},
        result,
        40.0,
        mode="quick",
    )

    window.canvas.undo_stack.push.assert_not_called()
    assert attempt.result_state == "rejected"
    assert attempt.stable_reason == "boundary-bond-risk"


def test_run_clean_2d_allows_safe_boundary_candidate() -> None:
    graph, bond_ids = _make_chain_graph()
    controller = Clean2DController()
    window = _mock_window(graph)
    window.canvas._selected_structure_ids.return_value = ({1, 2}, [])
    window.canvas._structural_bonds_for_clean2d = MagicMock(return_value=[graph.get_bond(bond_ids[0])])

    candidate = Clean2DCandidate(
        source="internal_templates",
        coords={1: (0.0, 0.0), 2: (42.0, 0.0)},
        message="Estructura 2D limpiada (motor interno)",
    )
    result = Clean2DResult(
        mode=Clean2DMode.QUICK,
        atom_ids={1, 2},
        before_coords={1: (0.0, 0.0), 2: (40.0, 0.0)},
        candidates=(candidate,),
        selected=candidate,
        message=candidate.message,
    )

    with patch("chemuson.gui.controllers.clean2d_controller.run_clean2d_engine", return_value=result):
        controller.run_clean_2d(window, 1.0, 200, "(test)")

    window.canvas.undo_stack.push.assert_called_once()
