from __future__ import annotations

import threading
import time

import pytest
from PyQt6.QtCore import QPointF
from PyQt6.QtWidgets import QApplication

from chemuson.core.model import BondStyle, MolGraph
from chemuson.gui.controllers import molecular_assistant_controller as assistant_controller_module
from chemuson.gui.dialogs import MolecularAssistantDialog
from chemuson.gui.main_window import ChemusonWindow
from chemuson.molecular_assistant import (
    MolecularAssistantResult,
    MolecularAssistantStatus,
)


@pytest.fixture(autouse=True)
def _isolated_config_home(tmp_path, monkeypatch):
    monkeypatch.setenv("XDG_CONFIG_HOME", str(tmp_path))


def _wait_for(predicate, timeout_s: float = 6.0) -> bool:
    deadline = time.monotonic() + timeout_s
    while time.monotonic() < deadline:
        QApplication.processEvents()
        if predicate():
            return True
        time.sleep(0.002)
    QApplication.processEvents()
    return bool(predicate())


def _graph(elements: tuple[str, ...]) -> MolGraph:
    graph = MolGraph()
    atom_ids = [
        graph.add_atom(element, float(index * 40), 0.0, atom_id=index + 1).id
        for index, element in enumerate(elements)
    ]
    for left, right in zip(atom_ids, atom_ids[1:]):
        graph.add_bond(left, right, order=1)
    return graph


def _success(graph: MolGraph, smiles: str) -> MolecularAssistantResult:
    return MolecularAssistantResult(
        status=MolecularAssistantStatus.SUCCESS,
        provider_id="test-provider",
        model_id="offline-test-model",
        proposed_smiles=smiles,
        graph=graph,
        validation_passed=True,
    )


def _make_window_with_two_molecules() -> tuple[ChemusonWindow, set[int], set[int]]:
    window = ChemusonWindow()
    window._properties_update_timer.stop()
    canvas = window.canvas
    canvas._insert_molgraph(_graph(("C", "C", "O")), select_inserted=True)
    source_atom_ids = set(canvas.state.selected_atoms)
    canvas._insert_molgraph_at(_graph(("N", "C")), QPointF(650.0, 450.0))
    other_atom_ids = set(canvas.model.atoms) - source_atom_ids
    canvas._select_inserted_items(source_atom_ids)
    canvas.undo_stack.setClean()
    window._properties_update_timer.stop()
    return window, source_atom_ids, other_atom_ids


def _finish_window(window: ChemusonWindow) -> None:
    for index in range(window.tabs.count()):
        window.tabs.widget(index).undo_stack.setClean()
    window.close()
    window.deleteLater()
    QApplication.processEvents()


def test_transform_action_is_discoverable_and_rejects_partial_or_ambiguous_selection():
    window, source_atom_ids, other_atom_ids = _make_window_with_two_molecules()
    try:
        action = window.action_ai_molecular_transform
        assert action.text() == "Transformar molécula seleccionada con IA..."
        assert window._command_registry.find_by_action(action) is not None
        window.command_palette._apply_filter("transformar")
        assert any(item.action is action for item in window.command_palette._filtered)
        structure_action = next(
            item for item in window.menuBar().actions() if item.text() == "Estructura"
        )
        assert action in structure_action.menu().actions()

        canvas = window.canvas
        canvas._select_inserted_items({min(source_atom_ids)})
        window._on_ai_molecular_transform()
        assert not window.findChildren(MolecularAssistantDialog)
        assert "molécula conectada" in window.statusBar().currentMessage()
        assert window._molecular_assistant_controller.active_jobs() == ()

        canvas._select_inserted_items(source_atom_ids | other_atom_ids)
        window._on_ai_molecular_transform()
        assert not window.findChildren(MolecularAssistantDialog)
        assert window._molecular_assistant_controller.active_jobs() == ()

        canvas._select_inserted_items(source_atom_ids)
        external_bond_id = next(
            bond_id
            for bond_id, bond in canvas.model.bonds.items()
            if bond.a1_id in other_atom_ids and bond.a2_id in other_atom_ids
        )
        canvas.state.selected_bonds.add(external_bond_id)
        window._on_ai_molecular_transform()
        assert not window.findChildren(MolecularAssistantDialog)
        assert window._molecular_assistant_controller.active_jobs() == ()
        canvas.state.selected_bonds.clear()

        canvas.scene.clearSelection()
        QApplication.processEvents()
        window._on_ai_molecular_transform()
        assert not window.findChildren(MolecularAssistantDialog)
        assert window._molecular_assistant_controller.active_jobs() == ()
    finally:
        _finish_window(window)


def test_transform_reviews_source_and_proposal_then_replaces_in_one_undo_step(
    monkeypatch,
):
    window, source_atom_ids, other_atom_ids = _make_window_with_two_molecules()
    proposed_graph = _graph(("N", "C", "Cl"))
    gui_thread_id = threading.get_ident()
    observations = {}

    def fake_source_export(graph, *, timeout_s):
        observations["export_thread"] = threading.get_ident()
        observations["export_timeout"] = timeout_s
        observations["export_atoms"] = len(graph.atoms)
        source_copy_atom = graph.atoms[min(source_atom_ids)]
        observations["source_stereo"] = source_copy_atom.stereo_cip
        observations["source_groups"] = source_copy_atom.r_group_substituents
        return "CCO"

    def fake_generator(request, config):
        observations["generator_thread"] = threading.get_ident()
        observations["transformation_request"] = request
        observations["api_key"] = config.api_key
        return _success(proposed_graph, "NCCCl")

    monkeypatch.setattr(
        assistant_controller_module,
        "molgraph_to_smiles_isolated_or_error",
        fake_source_export,
    )
    window._molecular_assistant_controller._generator = fake_generator
    canvas = window.canvas
    source_atom = canvas.model.get_atom(min(source_atom_ids))
    source_atom.stereo_cip = "R"
    source_atom.stereo_axial = "R_a"
    source_atom.stereo_helical = "P"
    source_atom.stereo_si_re = "si"
    source_atom.group_h_cap = 2
    source_atom.r_group_substituents = ("Me", "Et")
    source_bond = canvas.model.get_bond(min(canvas.model.bonds))
    source_bond.style = BondStyle.FLEX
    source_bond.stereo_ez = "E"
    source_bond.stereo_axial = "R_a"
    source_bond.stereo_endo_exo = "endo"
    source_bond.stereo_helical = "P"
    source_bond.flex_curve_1 = 0.25
    source_bond.flex_curve_2 = -0.25
    source_bond.pi_offset_sign = -1
    unrelated_signature = window._assistant_graph_signature(
        window._build_assistant_component_graph(canvas, other_atom_ids)
    )
    before_signature = window._assistant_graph_signature(canvas.model)
    before_index = canvas.undo_stack.index()
    before_count = canvas.undo_stack.count()
    source_atoms = [canvas.model.get_atom(atom_id) for atom_id in source_atom_ids]
    source_center = QPointF(
        (min(atom.x for atom in source_atoms) + max(atom.x for atom in source_atoms))
        / 2.0,
        (min(atom.y for atom in source_atoms) + max(atom.y for atom in source_atoms))
        / 2.0,
    )

    try:
        window._on_ai_molecular_transform()
        dialog = window.findChildren(MolecularAssistantDialog)[-1]
        assert dialog.windowTitle() == "Transformar molécula con IA"
        assert dialog.insert_button.text() == "Reemplazar molécula seleccionada"
        dialog.provider_combo.setCurrentIndex(
            dialog.provider_combo.findData("llama-cpp")
        )
        dialog.model_edit.setText("offline-test-model")
        dialog.api_key_edit.setText("transient-transform-key")
        dialog.description_edit.setPlainText("Sustituye el oxígeno terminal por cloro")
        dialog.generate_button.click()

        assert _wait_for(lambda: dialog.insert_button.isVisible())
        assert observations["export_thread"] != gui_thread_id
        assert observations["generator_thread"] == observations["export_thread"]
        assert observations["export_timeout"] == 8.0
        assert observations["export_atoms"] == len(source_atom_ids)
        assert observations["source_stereo"] == "R"
        assert observations["source_groups"] == ("Me", "Et")
        transformation_request = observations["transformation_request"]
        assert transformation_request.source_smiles == "CCO"
        assert transformation_request.instruction == "Sustituye el oxígeno terminal por cloro"
        assert "transient-transform-key" not in repr(transformation_request)
        assert dialog.source_smiles_preview.toPlainText() == "CCO"
        assert dialog.proposal_smiles_label.text() == "Estructura propuesta (SMILES)"
        assert dialog.proposal_smiles_label.isVisible()
        assert dialog.smiles_preview.toPlainText() == "NCCCl"
        assert window._assistant_graph_signature(canvas.model) == before_signature
        assert canvas.undo_stack.index() == before_index

        dialog.insert_button.click()
        assert canvas.undo_stack.count() == before_count + 1
        assert canvas.undo_stack.index() == before_index + 1
        assert set(canvas.model.atoms) >= other_atom_ids
        assert not source_atom_ids.intersection(canvas.model.atoms)
        assert len(canvas.model.atoms) == len(other_atom_ids) + len(
            proposed_graph.atoms
        )
        assert (
            window._assistant_graph_signature(
                window._build_assistant_component_graph(canvas, other_atom_ids)
            )
            == unrelated_signature
        )
        assert set(canvas.state.selected_atoms)
        inserted_atoms = [
            canvas.model.get_atom(atom_id) for atom_id in canvas.state.selected_atoms
        ]
        assert (
            (
                min(atom.x for atom in inserted_atoms)
                + max(atom.x for atom in inserted_atoms)
            )
            / 2.0
        ) == pytest.approx(source_center.x())
        assert (
            (
                min(atom.y for atom in inserted_atoms)
                + max(atom.y for atom in inserted_atoms)
            )
            / 2.0
        ) == pytest.approx(source_center.y())
        after_signature = window._assistant_graph_signature(canvas.model)

        canvas.undo_stack.undo()
        assert window._assistant_graph_signature(canvas.model) == before_signature
        assert canvas.undo_stack.index() == before_index
        assert canvas.undo_stack.isClean()

        canvas.undo_stack.redo()
        assert window._assistant_graph_signature(canvas.model) == after_signature
        assert canvas.undo_stack.index() == before_index + 1
    finally:
        _finish_window(window)


def test_transform_insert_variant_preserves_source_and_has_exact_undo_redo(monkeypatch):
    window, source_atom_ids, other_atom_ids = _make_window_with_two_molecules()
    canvas = window.canvas
    proposal = _graph(("N", "C", "Cl"))
    monkeypatch.setattr(
        assistant_controller_module,
        "molgraph_to_smiles_isolated_or_error",
        lambda *_args, **_kwargs: "CCO",
    )

    def fake_generator(request, _config):
        assert request.source_smiles == "CCO"
        assert request.instruction == "Añade un átomo de cloro"
        return _success(proposal, "NCCCl")

    window._molecular_assistant_controller._generator = fake_generator
    before_signature = window._assistant_graph_signature(canvas.model)
    source_signature = window._assistant_graph_signature(
        window._build_assistant_component_graph(canvas, source_atom_ids)
    )
    other_signature = window._assistant_graph_signature(
        window._build_assistant_component_graph(canvas, other_atom_ids)
    )
    before_atom_ids = set(canvas.model.atoms)
    before_index = canvas.undo_stack.index()
    before_count = canvas.undo_stack.count()

    try:
        window._on_ai_molecular_transform()
        dialog = window.findChildren(MolecularAssistantDialog)[-1]
        dialog.provider_combo.setCurrentIndex(dialog.provider_combo.findData("llama-cpp"))
        dialog.model_edit.setText("offline-test-model")
        dialog.description_edit.setPlainText("Añade un átomo de cloro")
        dialog.generate_button.click()
        assert _wait_for(lambda: dialog.insert_variant_button.isVisible())
        assert dialog.insert_variant_button.text() == "Insertar variante"
        assert dialog.insert_button.text() == "Reemplazar molécula seleccionada"

        dialog.insert_variant_button.click()
        assert canvas.undo_stack.count() == before_count + 1
        assert canvas.undo_stack.index() == before_index + 1
        assert source_atom_ids.issubset(canvas.model.atoms)
        assert other_atom_ids.issubset(canvas.model.atoms)
        assert window._assistant_graph_signature(
            window._build_assistant_component_graph(canvas, source_atom_ids)
        ) == source_signature
        assert window._assistant_graph_signature(
            window._build_assistant_component_graph(canvas, other_atom_ids)
        ) == other_signature
        variant_atom_ids = set(canvas.model.atoms) - before_atom_ids
        assert len(variant_atom_ids) == len(proposal.atoms)
        source_right = max(canvas.model.get_atom(atom_id).x for atom_id in source_atom_ids)
        variant_left = min(canvas.model.get_atom(atom_id).x for atom_id in variant_atom_ids)
        assert variant_left > source_right
        after_signature = window._assistant_graph_signature(canvas.model)

        canvas.undo_stack.undo()
        assert window._assistant_graph_signature(canvas.model) == before_signature
        assert canvas.undo_stack.index() == before_index
        assert source_atom_ids.issubset(canvas.model.atoms)

        canvas.undo_stack.redo()
        assert window._assistant_graph_signature(canvas.model) == after_signature
        assert canvas.undo_stack.index() == before_index + 1
        assert window._assistant_graph_signature(
            window._build_assistant_component_graph(canvas, source_atom_ids)
        ) == source_signature
    finally:
        _finish_window(window)


def test_transform_export_failure_decline_and_stale_source_fail_closed(monkeypatch):
    window, source_atom_ids, _other_atom_ids = _make_window_with_two_molecules()
    calls = []
    proposed_graph = _graph(("N", "C"))
    window._molecular_assistant_controller._generator = lambda *_args: (
        calls.append("provider") or _success(proposed_graph, "CN")
    )
    canvas = window.canvas
    original_signature = window._assistant_graph_signature(canvas.model)

    try:
        monkeypatch.setattr(
            assistant_controller_module,
            "molgraph_to_smiles_isolated_or_error",
            lambda *_args, **_kwargs: (_ for _ in ()).throw(
                RuntimeError("private detail")
            ),
        )
        window._on_ai_molecular_transform()
        export_failure_dialog = window.findChildren(MolecularAssistantDialog)[-1]
        export_failure_dialog.provider_combo.setCurrentIndex(
            export_failure_dialog.provider_combo.findData("llama-cpp")
        )
        export_failure_dialog.model_edit.setText("offline-test-model")
        export_failure_dialog.description_edit.setPlainText("Cambiar un grupo")
        export_failure_dialog.generate_button.click()
        assert _wait_for(
            lambda: (
                "exportar la molécula original" in export_failure_dialog.status_label.text()
            )
        )
        assert calls == []
        assert "private detail" not in export_failure_dialog.status_label.text()
        assert window._assistant_graph_signature(canvas.model) == original_signature
        export_failure_dialog.close()
        assert _wait_for(
            lambda: not window._molecular_assistant_controller.active_jobs()
        )

        monkeypatch.setattr(
            assistant_controller_module,
            "molgraph_to_smiles_isolated_or_error",
            lambda *_args, **_kwargs: "CCO",
        )
        window._on_ai_molecular_transform()
        declined_dialog = window.findChildren(MolecularAssistantDialog)[-1]
        declined_dialog.provider_combo.setCurrentIndex(
            declined_dialog.provider_combo.findData("llama-cpp")
        )
        declined_dialog.model_edit.setText("offline-test-model")
        declined_dialog.description_edit.setPlainText("Cambiar un grupo")
        declined_dialog.generate_button.click()
        assert _wait_for(lambda: declined_dialog.insert_button.isVisible())
        before_decline = window._assistant_graph_signature(canvas.model)
        before_decline_index = canvas.undo_stack.index()
        before_decline_selection = set(canvas.state.selected_atoms)
        before_decline_clean = canvas.undo_stack.isClean()
        declined_dialog.reject()
        QApplication.processEvents()
        assert window._assistant_graph_signature(canvas.model) == before_decline
        assert canvas.undo_stack.index() == before_decline_index
        assert set(canvas.state.selected_atoms) == before_decline_selection
        assert canvas.undo_stack.isClean() == before_decline_clean
        assert _wait_for(
            lambda: not window._molecular_assistant_controller.active_jobs()
        )

        window._on_ai_molecular_transform()
        stale_dialog = window.findChildren(MolecularAssistantDialog)[-1]
        stale_dialog.provider_combo.setCurrentIndex(
            stale_dialog.provider_combo.findData("llama-cpp")
        )
        stale_dialog.model_edit.setText("offline-test-model")
        stale_dialog.description_edit.setPlainText("Cambiar un grupo")
        stale_dialog.generate_button.click()
        assert _wait_for(lambda: stale_dialog.insert_button.isVisible())
        calls.clear()
        before_replace_index = canvas.undo_stack.index()
        canvas._select_inserted_items(source_atom_ids - {min(source_atom_ids)})
        stale_dialog.insert_button.click()
        assert "cambió" in stale_dialog.status_label.text()
        assert window._assistant_graph_signature(canvas.model) == original_signature
        assert canvas.undo_stack.index() == before_replace_index

        canvas._select_inserted_items(source_atom_ids)
        stale_atom = canvas.model.get_atom(min(source_atom_ids))
        stale_atom.x += 5.0
        stale_atom.stereo_cip = "S"
        stale_signature = window._assistant_graph_signature(canvas.model)
        stale_dialog.insert_button.click()
        assert "cambió" in stale_dialog.status_label.text()
        assert window._assistant_graph_signature(canvas.model) == stale_signature
        assert canvas.undo_stack.index() == before_replace_index

        replacement_canvas = window._create_document_tab(make_current=True)
        stale_dialog.insert_button.click()
        assert "Activa el documento" in stale_dialog.status_label.text()
        assert window._assistant_graph_signature(canvas.model) == stale_signature
        window._tab_manager.discard_canvas(canvas)
        stale_dialog.insert_button.click()
        assert "ya no está disponible" in stale_dialog.status_label.text()
        assert window._assistant_graph_signature(replacement_canvas.model) == ((), ())
    finally:
        _finish_window(window)


def test_transform_dialog_close_suppresses_late_worker_result(monkeypatch):
    window, _source_atom_ids, _other_atom_ids = _make_window_with_two_molecules()
    started = threading.Event()
    release = threading.Event()
    proposed_graph = _graph(("N", "C"))

    def delayed_generator(*_args):
        started.set()
        release.wait(4.0)
        return _success(proposed_graph, "CN")

    window._molecular_assistant_controller._generator = delayed_generator
    monkeypatch.setattr(
        assistant_controller_module,
        "molgraph_to_smiles_isolated_or_error",
        lambda *_args, **_kwargs: "CCO",
    )
    canvas = window.canvas
    before_signature = window._assistant_graph_signature(canvas.model)
    before_index = canvas.undo_stack.index()

    try:
        window._on_ai_molecular_transform()
        dialog = window.findChildren(MolecularAssistantDialog)[-1]
        dialog.provider_combo.setCurrentIndex(
            dialog.provider_combo.findData("llama-cpp")
        )
        dialog.model_edit.setText("offline-test-model")
        dialog.description_edit.setPlainText("Cambiar un grupo")
        dialog.generate_button.click()
        assert _wait_for(started.is_set)
        assert window._molecular_assistant_controller.active_jobs()

        dialog.reject()
        release.set()
        assert _wait_for(
            lambda: not window._molecular_assistant_controller.active_jobs()
        )
        assert window._molecular_assistant_transform_jobs == {}
        assert window._molecular_assistant_results == {}
        assert window._assistant_graph_signature(canvas.model) == before_signature
        assert canvas.undo_stack.index() == before_index
    finally:
        release.set()
        _finish_window(window)
