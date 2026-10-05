from __future__ import annotations

import threading
import time

import pytest
from PyQt6.QtTest import QSignalSpy
from PyQt6.QtWidgets import QApplication

from chemuson.core.model import MolGraph
from chemuson.gui.dialogs import MolecularAssistantDialog
from chemuson.gui.main_window import ChemusonWindow
from chemuson.molecular_assistant import (
    MolecularAssistantResult,
    MolecularAssistantStatus,
)


@pytest.fixture(scope="module", autouse=True)
def _qapp():
    return QApplication.instance() or QApplication([])


@pytest.fixture(autouse=True)
def _isolated_config_home(tmp_path, monkeypatch):
    monkeypatch.setenv("XDG_CONFIG_HOME", str(tmp_path))


def _success_result(graph: MolGraph | None = None) -> MolecularAssistantResult:
    graph = graph or _single_carbon_graph()
    return MolecularAssistantResult(
        status=MolecularAssistantStatus.SUCCESS,
        provider_id="test-provider",
        model_id="test-model",
        proposed_smiles="C",
        graph=graph,
        validation_passed=True,
    )


def _single_carbon_graph() -> MolGraph:
    graph = MolGraph()
    graph.add_atom("C", 0.0, 0.0)
    return graph


def _graph_signature(canvas) -> tuple:
    atoms = tuple(
        sorted(
            (
                atom.id,
                atom.element,
                atom.x,
                atom.y,
                atom.formal_charge,
                atom.isotope,
            )
            for atom in canvas.model.atoms.values()
        )
    )
    bonds = tuple(
        sorted(
            (bond.a1_id, bond.a2_id, bond.order, bond.is_aromatic)
            for bond in canvas.model.bonds.values()
        )
    )
    return atoms, bonds


def _editor_snapshot(canvas) -> tuple:
    return (
        _graph_signature(canvas),
        tuple(id(item) for item in canvas.scene.selectedItems()),
        canvas.undo_stack.index(),
        canvas.undo_stack.isClean(),
    )


def _open_dialog(window: ChemusonWindow) -> MolecularAssistantDialog:
    window._on_ai_molecular_assistant()
    QApplication.processEvents()
    dialogs = window.findChildren(MolecularAssistantDialog)
    assert dialogs
    return dialogs[-1]


def _seed_canvas(window: ChemusonWindow) -> None:
    window.canvas._insert_molgraph(_single_carbon_graph(), select_inserted=True)
    window.canvas.undo_stack.setClean()


def _wait_for(predicate, timeout_s: float = 5.0) -> bool:
    deadline = time.monotonic() + timeout_s
    while time.monotonic() < deadline:
        QApplication.processEvents()
        if predicate():
            return True
        time.sleep(0.005)
    QApplication.processEvents()
    return bool(predicate())


def test_ai_action_is_discoverable_in_structure_menu_and_command_palette():
    window = ChemusonWindow()
    try:
        action = window.action_ai_molecular_assistant
        assert action.text() == "Dibujar estructura con IA..."
        entry = window._command_registry.find_by_action(action)
        assert entry is not None and entry.section == "Estructura"
        structure_action = next(
            item for item in window.menuBar().actions() if item.text() == "Estructura"
        )
        assert action in structure_action.menu().actions()
        window.command_palette._apply_filter("generar")
        assert any(item.action is action for item in window.command_palette._filtered)
    finally:
        window.close()


def test_dialog_is_modeless_masks_key_and_rejects_missing_request_fields():
    window = ChemusonWindow()
    try:
        dialog = _open_dialog(window)
        assert not dialog.isModal()
        assert dialog.api_key_edit.echoMode() == dialog.api_key_edit.EchoMode.Password
        emitted = []
        dialog.generation_requested.connect(lambda *args: emitted.append(args))

        dialog.generate_button.click()
        assert emitted == []
        assert "descripción" in dialog.status_label.text()

        dialog.description_edit.setPlainText("Dibuja cafeína")
        dialog.generate_button.click()
        assert emitted == []
        assert "endpoint" in dialog.status_label.text()
    finally:
        window.close()


def test_controller_runs_generation_on_worker_thread_and_keeps_key_transient():
    gui_thread_id = threading.get_ident()
    observed = {}
    result = _success_result()

    def fake_generator(description, config):
        observed["thread_id"] = threading.get_ident()
        observed["description"] = description
        observed["api_key"] = config.api_key
        observed["config_repr"] = repr(config)
        return result

    from chemuson.gui.controllers import MolecularAssistantController

    controller = MolecularAssistantController(generator=fake_generator)
    try:
        finished = QSignalSpy(controller.job_finished)
        completed = []
        controller.job_finished.connect(
            lambda received_id, value: completed.append((received_id, value))
        )
        job_id = controller.start_job(
            "Dibuja cafeína",
            base_url="https://provider.example/v1",
            model="test-model",
            api_key="transient-secret",
        )
        assert job_id is not None
        assert finished.wait(5000)
        assert _wait_for(lambda: not controller.active_jobs())

        assert completed == [(job_id, result)]
        assert observed["thread_id"] != gui_thread_id
        assert observed["description"] == "Dibuja cafeína"
        assert observed["api_key"] == "transient-secret"
        assert "transient-secret" not in observed["config_repr"]
    finally:
        for job_id in controller.active_jobs():
            controller.abandon_job(job_id)


def test_controller_rejects_bad_endpoint_without_starting_a_job():
    from chemuson.gui.controllers import MolecularAssistantController

    calls = []
    controller = MolecularAssistantController(
        generator=lambda *_args: calls.append(True) or _success_result()
    )
    assert controller.start_job(
        "Dibuja cafeína",
        base_url="javascript:alert(1)",
        model="test-model",
    ) is None
    assert controller.active_jobs() == ()
    assert calls == []


def test_abandoned_job_suppresses_late_success_result():
    started = threading.Event()
    release = threading.Event()

    def slow_generator(_description, _config):
        started.set()
        release.wait(timeout=3.0)
        return _success_result()

    from chemuson.gui.controllers import MolecularAssistantController

    controller = MolecularAssistantController(generator=slow_generator)
    finished = QSignalSpy(controller.job_finished)
    job_id = controller.start_job(
        "Dibuja cafeína",
        base_url="https://provider.example/v1",
        model="test-model",
    )
    assert job_id is not None
    assert started.wait(timeout=2.0)
    controller.abandon_job(job_id)
    release.set()
    assert _wait_for(lambda: not controller.active_jobs())
    assert len(finished) == 0


@pytest.mark.parametrize(
    ("status", "reason_code", "validation_passed"),
    [
        (MolecularAssistantStatus.INVALID_STRUCTURE, "invalid_smiles", False),
        (MolecularAssistantStatus.INVALID_REQUEST, "empty_prompt", None),
        (MolecularAssistantStatus.VALIDATION_ERROR, "parser_timeout", None),
        (MolecularAssistantStatus.PROVIDER_ERROR, "network_error", None),
        (MolecularAssistantStatus.MALFORMED_RESPONSE, "invalid_json", None),
        (MolecularAssistantStatus.CANCELLED, "cancelled", None),
    ],
)
def test_failure_result_preserves_canvas_selection_undo_and_dirty_state(
    status, reason_code, validation_passed
):
    window = ChemusonWindow()
    try:
        _seed_canvas(window)
        before = _editor_snapshot(window.canvas)
        dialog = _open_dialog(window)
        dialog.set_job_id(71)
        window._molecular_assistant_dialogs[71] = (dialog, window.canvas)
        failure = MolecularAssistantResult(
            status=status,
            provider_id="test-provider",
            model_id="test-model",
            proposed_smiles="invalid" if status is MolecularAssistantStatus.INVALID_STRUCTURE else None,
            graph=None,
            validation_passed=validation_passed,
            reason_code=reason_code,
        )

        window._on_molecular_assistant_job_finished(71, failure)

        assert _editor_snapshot(window.canvas) == before
        assert not dialog.insert_button.isVisible()
        assert dialog.status_label.text() == (
            f"No se generó una estructura válida. Estado: {status.value}; "
            f"motivo: {reason_code}."
        )
    finally:
        window.close()


def test_success_is_previewed_before_normal_canvas_insertion(monkeypatch):
    window = ChemusonWindow()
    try:
        canvas = window.canvas
        dialog = _open_dialog(window)
        dialog.set_job_id(72)
        window._molecular_assistant_dialogs[72] = (dialog, canvas)
        result = _success_result()
        insert_calls = []
        monkeypatch.setattr(
            canvas,
            "_insert_molgraph",
            lambda graph, select_inserted=False: insert_calls.append(
                (graph, select_inserted)
            ),
        )

        window._on_molecular_assistant_job_finished(72, result)

        assert insert_calls == []
        assert dialog.insert_button.isVisible()
        assert "test-provider" in dialog.provenance_label.text()
        assert "no demuestra" in dialog.preview_group.text()

        dialog.insert_button.click()

        assert insert_calls == [(result.graph, True)]
        assert 72 not in window._molecular_assistant_dialogs
        assert 72 not in window._molecular_assistant_results
    finally:
        window.close()


def test_declining_preview_and_late_result_after_close_do_not_mutate_canvas():
    window = ChemusonWindow()
    try:
        _seed_canvas(window)
        canvas = window.canvas
        before = _editor_snapshot(canvas)
        dialog = _open_dialog(window)
        dialog.set_job_id(73)
        window._molecular_assistant_dialogs[73] = (dialog, canvas)
        window._on_molecular_assistant_job_finished(73, _success_result())
        dialog.reject()
        QApplication.processEvents()
        assert _editor_snapshot(canvas) == before
        assert 73 not in window._molecular_assistant_results

        late_dialog = _open_dialog(window)
        late_dialog.set_job_id(74)
        window._molecular_assistant_dialogs[74] = (late_dialog, canvas)
        late_dialog.reject()
        QApplication.processEvents()
        window._on_molecular_assistant_job_finished(74, _success_result())
        assert _editor_snapshot(canvas) == before
    finally:
        window.close()


def test_canvas_graph_insertion_is_one_normal_undo_redo_operation():
    from chemuson.gui.canvas import ChemusonCanvas

    canvas = ChemusonCanvas()
    try:
        canvas._insert_molgraph(_single_carbon_graph(), select_inserted=True)
        inserted = _graph_signature(canvas)

        assert len(canvas.model.atoms) == 1
        assert canvas.undo_stack.index() == 1
        assert canvas.scene.selectedItems()

        canvas.undo_stack.undo()
        assert not canvas.model.atoms
        assert canvas.undo_stack.isClean()

        canvas.undo_stack.redo()
        assert _graph_signature(canvas) == inserted
        assert canvas.undo_stack.index() == 1
        assert all(item.scene() is canvas.scene for item in canvas.scene.selectedItems())
    finally:
        canvas.deleteLater()
        QApplication.processEvents()
