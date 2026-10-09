from __future__ import annotations

import threading
import time

import pytest
from PyQt6.QtCore import QSettings, QStandardPaths
from PyQt6.QtTest import QSignalSpy, QTest
from PyQt6.QtWidgets import QApplication, QMessageBox

from chemuson.core.model import MolGraph
from chemuson.gui.dialogs import MolecularAssistantDialog
from chemuson.gui.main_window import ChemusonWindow
from chemuson.name2structure import (
    MolecularIdentityStatus,
    MolecularIdentityVerification,
    NameToStructureResult,
)
from chemuson.gui.controllers.molecular_assistant_controller import (
    MolecularAssistantResolution,
    MolecularResolutionMethod,
    StructureOrigin,
)
from chemuson.molecular_assistant import (
    OPENAI_COMPATIBLE_PROFILES,
    MolecularAssistant,
    MolecularAssistantResult,
    MolecularAssistantStatus,
    ProviderResponse,
)
from chemuson.molecular_assistant import service as service_module


@pytest.fixture(scope="module", autouse=True)
def _qapp():
    return QApplication.instance() or QApplication([])


@pytest.fixture(autouse=True)
def _isolated_config_home(tmp_path, monkeypatch):
    config_location = QStandardPaths.writableLocation(
        QStandardPaths.StandardLocation.ConfigLocation
    )
    monkeypatch.setenv("XDG_CONFIG_HOME", str(tmp_path))
    QSettings.setPath(
        QSettings.Format.NativeFormat,
        QSettings.Scope.UserScope,
        str(tmp_path),
    )
    try:
        settings = QSettings("Chemuson", "Chemuson")
        settings.clear()
        settings.sync()
        yield
    finally:
        QSettings.setPath(
            QSettings.Format.NativeFormat,
            QSettings.Scope.UserScope,
            config_location,
        )


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
        assert dialog.resolution_method == "ai_reference"
        assert dialog.identity_verification_enabled is True
        assert dialog.allow_external_identity_reference is False
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


def test_provider_profiles_fill_editable_endpoint_and_clear_transient_keys():
    dialog = MolecularAssistantDialog(profiles=OPENAI_COMPATIBLE_PROFILES)
    try:
        emitted = []
        dialog.generation_requested.connect(lambda *args: emitted.append(args))
        dialog.description_edit.setPlainText("Dibuja etanol")
        dialog.model_edit.setText("loaded-model-id")
        dialog.api_key_edit.setText("key-for-another-profile")

        openai_index = dialog.provider_combo.findData("openai")
        dialog.provider_combo.setCurrentIndex(openai_index)
        assert dialog.base_url_edit.text() == "https://api.openai.com/v1"
        assert dialog.model_edit.text() == ""
        assert dialog.api_key_label.text() == "API key (requerida)"
        dialog.model_edit.setText("openai-model")
        assert dialog.api_key_edit.text() == ""
        dialog.generate_button.click()
        assert emitted == []
        assert "requiere una API key" in dialog.status_label.text()

        dialog.api_key_edit.setText("transient-openai-key")
        dialog.generate_button.click()
        assert emitted == [
            (
                "Dibuja etanol",
                "openai",
                "https://api.openai.com/v1",
                "openai-model",
                "transient-openai-key",
                False,
                60,
                4096,
            )
        ]

        local_index = dialog.provider_combo.findData("lm-studio")
        dialog.provider_combo.setCurrentIndex(local_index)
        assert dialog.base_url_edit.text() == "http://127.0.0.1:1234/v1"
        assert dialog.api_key_label.text() == "API key (opcional)"
        assert dialog.api_key_edit.text() == ""
        assert dialog.model_edit.text() == ""
    finally:
        dialog.close()


def test_dialog_elapsed_timer_and_human_facing_provider_failures():
    dialog = MolecularAssistantDialog(profiles=OPENAI_COMPATIBLE_PROFILES)
    try:
        dialog.timeout_spin.setValue(180)
        dialog.set_pending()
        assert dialog.status_label.text().startswith("Generando y validando…")
        QTest.qWait(1050)
        assert " s" in dialog.status_label.text()

        dialog.show_failure("provider_error", "timeout")
        assert "(180 s)" in dialog.status_label.text()
        assert "timeout" not in dialog.status_label.text()

        dialog.show_failure(
            "malformed_response",
            "invalid_json",
            structured_output_requested=True,
            structured_output_native=False,
            structured_output_fallback_used=True,
            format_repair_used=True,
            format_repair_succeeded=False,
        )
        assert "formato estructurado requerido" in dialog.status_label.text()
        assert "invalid_json" not in dialog.status_label.text()
        assert "contrato textual" in dialog.output_diagnostic_label.text()
        assert "Reintento de formato: fallido" in dialog.output_diagnostic_label.text()
        assert "secret" not in dialog.output_diagnostic_label.text()

        dialog.show_failure("provider_error", "network_error")
        assert not dialog.output_diagnostic_label.isVisible()
        assert "contactar el endpoint configurado" in dialog.status_label.text()
    finally:
        dialog.close()


def test_controller_runs_generation_on_worker_thread_and_keeps_key_transient():
    gui_thread_id = threading.get_ident()
    observed = {}
    result = _success_result()

    def fake_generator(request, config):
        observed["thread_id"] = threading.get_ident()
        observed["description"] = request.description
        observed["api_key"] = config.api_key
        observed["config_repr"] = repr(config)
        observed["timeout_s"] = config.timeout_s
        observed["max_tokens"] = config.max_tokens
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
            timeout_s=180,
            max_tokens=2048,
        )
        assert job_id is not None
        assert finished.wait(5000)
        assert _wait_for(lambda: not controller.active_jobs())

        assert len(completed) == 1
        assert completed[0][0] == job_id
        assert completed[0][1].ai_result is result
        assert completed[0][1].default_candidate == "ai"
        assert observed["thread_id"] != gui_thread_id
        assert observed["description"] == "Dibuja cafeína"
        assert observed["api_key"] == "transient-secret"
        assert observed["timeout_s"] == 180
        assert observed["max_tokens"] == 2048
        assert "transient-secret" not in observed["config_repr"]
    finally:
        for job_id in controller.active_jobs():
            controller.abandon_job(job_id)


def test_controller_runs_identity_verification_in_the_worker_and_relays_separately():
    gui_thread_id = threading.get_ident()
    events = []
    identity = MolecularIdentityVerification(
        MolecularIdentityStatus.VERIFIED,
        requested_name="ethanol",
        reference_identifier="fake:ethanol",
    )

    def generator(_request, _config):
        events.append(("generation", threading.get_ident()))
        return _success_result()

    def resolver(name, *, allow_network):
        events.append(("reference", threading.get_ident(), name, allow_network))
        return NameToStructureResult(
            name,
            _single_carbon_graph(),
            "fake-reference",
            1.0,
            smiles="C",
            resolved_name=name,
        )

    def verifier(description, result, reference):
        assert description == "Draw ethanol"
        assert result.status is MolecularAssistantStatus.SUCCESS
        assert reference.source == "fake-reference"
        events.append(("identity", threading.get_ident()))
        return identity

    from chemuson.gui.controllers import MolecularAssistantController

    controller = MolecularAssistantController(
        generator=generator,
        identity_verifier=verifier,
        reference_resolver=resolver,
    )
    try:
        finished = QSignalSpy(controller.job_finished)
        identity_ready = QSignalSpy(controller.identity_ready)
        controller.start_job(
            "Draw ethanol",
            base_url="http://127.0.0.1:8081/v1",
            model="local-model",
            resolution_method=MolecularResolutionMethod.AI_REFERENCE,
        )
        assert finished.wait(5000)
        assert len(identity_ready) == 1
        assert [event[0] for event in events] == ["generation", "reference", "identity"]
        assert all(event[1] != gui_thread_id for event in events)
        assert events[1][2:] == ("ethanol", False)
        assert identity_ready[0][1] is identity
    finally:
        for job_id in controller.active_jobs():
            controller.abandon_job(job_id)


def test_successful_format_repair_still_runs_identity_verification(monkeypatch):
    monkeypatch.setattr(
        service_module,
        "smiles_to_molgraph_isolated",
        lambda _smiles, *, timeout_s: (_single_carbon_graph(), None),
    )

    class RepairProvider:
        provider_id = "repair-provider"
        model_id = "repair-model"

        def __init__(self):
            self.count = 0

        def generate(self, _request):
            self.count += 1
            content = "not JSON" if self.count == 1 else '{"smiles":"CCO"}'
            return ProviderResponse(content, self.model_id)

    def generator(request, _config):
        return MolecularAssistant(RepairProvider()).generate(request)

    checked = []

    def resolver(name, *, allow_network):
        assert name == "ethanol"
        assert allow_network is False
        return NameToStructureResult(
            name,
            _single_carbon_graph(),
            "fake-reference",
            1.0,
            smiles="C",
            resolved_name=name,
        )

    def verifier(_description, result, reference):
        checked.append((result, reference))
        assert result.format_repair_used is True
        assert result.format_repair_succeeded is True
        return MolecularIdentityVerification(
            MolecularIdentityStatus.VERIFIED,
            requested_name="ethanol",
            reference_identifier="fake:ethanol",
        )

    from chemuson.gui.controllers import MolecularAssistantController

    controller = MolecularAssistantController(
        generator=generator,
        identity_verifier=verifier,
        reference_resolver=resolver,
    )
    try:
        finished = QSignalSpy(controller.job_finished)
        identity_ready = QSignalSpy(controller.identity_ready)
        controller.start_job(
            "Draw ethanol",
            base_url="https://provider.example/v1",
            model="repair-model",
            resolution_method=MolecularResolutionMethod.AI_REFERENCE,
        )
        assert finished.wait(5000)
        assert _wait_for(lambda: not controller.active_jobs())
        assert len(identity_ready) == 1
        assert len(checked) == 1
        assert checked[0][1].source == "fake-reference"
        assert checked[0][0].status is MolecularAssistantStatus.SUCCESS
        assert identity_ready[0][1].status is MolecularIdentityStatus.VERIFIED
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


def test_controller_enforces_profile_key_before_worker_and_sets_profile_identity():
    observed = []
    from chemuson.gui.controllers import MolecularAssistantController

    controller = MolecularAssistantController(
        generator=lambda _description, config: observed.append(config) or _success_result()
    )
    assert controller.start_job(
        "Dibuja cafeína",
        base_url="https://api.openai.com/v1",
        model="test-model",
        provider_id="openai",
    ) is None
    assert controller.start_job(
        "Dibuja cafeína",
        base_url="https://provider.example/v1",
        model="test-model",
        provider_id="unknown-profile",
        api_key="secret",
    ) is None
    assert controller.active_jobs() == ()
    assert observed == []

    finished = QSignalSpy(controller.job_finished)
    job_id = controller.start_job(
        "Dibuja cafeína",
        base_url="https://api.openai.com/v1",
        model="test-model",
        provider_id="openai",
        api_key="openai-secret",
    )
    assert job_id is not None
    assert finished.wait(5000)
    assert _wait_for(lambda: not controller.active_jobs())
    assert observed[0].provider_id == "openai"
    assert "openai-secret" not in repr(observed[0])


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
        window._register_molecular_assistant_dialog_job(dialog, 71, window.canvas)
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
        assert dialog.status_label.text()
        assert reason_code not in dialog.status_label.text()
        if reason_code == "invalid_json":
            assert "formato estructurado" in dialog.status_label.text()
    finally:
        window.close()


def test_success_is_previewed_before_normal_canvas_insertion(monkeypatch):
    window = ChemusonWindow()
    try:
        canvas = window.canvas
        dialog = _open_dialog(window)
        window._register_molecular_assistant_dialog_job(dialog, 72, canvas)
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


def test_identity_mismatch_requires_explicit_override_confirmation(monkeypatch):
    window = ChemusonWindow()
    try:
        canvas = window.canvas
        before = _editor_snapshot(canvas)
        dialog = _open_dialog(window)
        window._register_molecular_assistant_dialog_job(dialog, 75, canvas)
        identity = MolecularIdentityVerification(
            MolecularIdentityStatus.MISMATCH,
            requested_name="colesterol",
            reference_identifier="pubchem:cholesterol",
        )
        ai_result = _success_result()
        reference_graph, reference_error = service_module.smiles_to_molgraph_isolated(
            "CCO", timeout_s=5.0
        )
        assert reference_error is None and reference_graph is not None
        reference_result = NameToStructureResult(
            "colesterol",
            reference_graph,
            "pubchem",
            0.9,
            smiles="CCO",
            resolved_name="cholesterol",
        )
        result = MolecularAssistantResolution(
            MolecularResolutionMethod.AI_REFERENCE,
            ai_result,
            reference_result,
            identity,
            "reference",
            StructureOrigin.AI_MISMATCH_REFERENCE,
        )
        window._molecular_assistant_identity_results[75] = identity
        window._molecular_assistant_results[75] = result

        window._on_molecular_assistant_job_finished(75, result)
        assert "no coincide" in dialog.identity_label.text()
        assert "y referencia no coinciden" in dialog.provenance_label.text()
        assert "test-provider/test-model" in dialog.provenance_label.text()
        assert dialog.insert_button.text() == "Usar referencia PubChem"
        assert dialog.use_ai_proposal_button.isVisible()
        assert dialog.reference_smiles_preview.toPlainText() == "CCO"
        assert _editor_snapshot(canvas) == before

        answers = [QMessageBox.StandardButton.No, QMessageBox.StandardButton.Yes]
        prompts = []

        def confirm(*args):
            prompts.append(args[2])
            return answers.pop(0)

        monkeypatch.setattr(QMessageBox, "warning", confirm)
        dialog.use_ai_proposal_button.click()
        assert _editor_snapshot(canvas) == before
        assert "no confirmaste" in dialog.status_label.text()

        dialog.use_ai_proposal_button.click()
        assert len(canvas.model.atoms) == 1
        assert len(prompts) == 2
        assert all("colesterol" in prompt for prompt in prompts)
    finally:
        window.canvas.undo_stack.setClean()
        window.close()


def test_mismatch_reference_candidate_uses_normal_undo_redo_without_override(monkeypatch):
    window = ChemusonWindow()
    try:
        canvas = window.canvas
        dialog = _open_dialog(window)
        window._register_molecular_assistant_dialog_job(dialog, 76, canvas)
        ai_result = _success_result()
        reference_graph = service_module.smiles_to_molgraph_isolated
        reference, error = reference_graph("CCO", timeout_s=5.0)
        assert error is None and reference is not None
        identity = MolecularIdentityVerification(
            MolecularIdentityStatus.MISMATCH,
            requested_name="ethanol",
            reference_identifier="pubchem:ethanol",
        )
        outcome = MolecularAssistantResolution(
            MolecularResolutionMethod.AI_REFERENCE,
            ai_result,
            NameToStructureResult("ethanol", reference, "pubchem", 0.9, "CCO", "ethanol"),
            identity,
            "reference",
            StructureOrigin.AI_MISMATCH_REFERENCE,
        )
        window._molecular_assistant_identity_results[76] = identity
        window._molecular_assistant_results[76] = outcome
        window._on_molecular_assistant_job_finished(76, outcome)
        monkeypatch.setattr(
            QMessageBox,
            "warning",
            lambda *_args: pytest.fail("reference selection must not ask for AI override"),
        )

        dialog.insert_button.click()
        assert len(canvas.model.atoms) == len(reference.atoms)
        assert canvas.undo_stack.index() == 1
        canvas.undo_stack.undo()
        assert not canvas.model.atoms
        canvas.undo_stack.redo()
        assert len(canvas.model.atoms) == len(reference.atoms)
    finally:
        window.canvas.undo_stack.setClean()
        window.close()


def test_reference_fallback_preview_and_insertion_use_normal_undo_path():
    window = ChemusonWindow()
    try:
        canvas = window.canvas
        dialog = _open_dialog(window)
        window._register_molecular_assistant_dialog_job(dialog, 77, canvas)
        reference_graph, error = service_module.smiles_to_molgraph_isolated(
            "CCO", timeout_s=5.0
        )
        assert error is None and reference_graph is not None
        ai_failure = MolecularAssistantResult(
            status=MolecularAssistantStatus.PROVIDER_ERROR,
            provider_id="llama-cpp",
            model_id="qwen-local",
            reason_code="generation_exhausted",
            finish_reason="length",
            completion_tokens=4096,
            reasoning_tokens=4096,
        )
        reference_result = NameToStructureResult(
            "tetrandrine",
            reference_graph,
            "pubchem",
            0.95,
            smiles="CCO",
            resolved_name="Tetrandrine",
        )
        identity = MolecularIdentityVerification(
            MolecularIdentityStatus.NOT_APPLICABLE,
            requested_name="tetrandrine",
            reason_code="proposal_unavailable",
        )
        outcome = MolecularAssistantResolution(
            MolecularResolutionMethod.AI_REFERENCE,
            ai_failure,
            reference_result,
            identity,
            "reference",
            StructureOrigin.REFERENCE,
        )
        window._molecular_assistant_identity_results[77] = identity
        window._molecular_assistant_results[77] = outcome

        window._on_molecular_assistant_job_finished(77, outcome)
        assert "agotó el límite de generación" in dialog.status_label.text()
        assert "PubChem" in dialog.provenance_label.text()
        assert dialog.smiles_preview.toPlainText() == "CCO"
        assert dialog.insert_button.text() == "Insertar referencia química"

        dialog.insert_button.click()
        assert len(canvas.model.atoms) == len(reference_graph.atoms)
        assert canvas.undo_stack.index() == 1
        canvas.undo_stack.undo()
        assert not canvas.model.atoms
        canvas.undo_stack.redo()
        assert len(canvas.model.atoms) == len(reference_graph.atoms)
    finally:
        window.canvas.undo_stack.setClean()
        window.close()


def test_declining_preview_and_late_result_after_close_do_not_mutate_canvas():
    window = ChemusonWindow()
    try:
        _seed_canvas(window)
        canvas = window.canvas
        before = _editor_snapshot(canvas)
        dialog = _open_dialog(window)
        window._register_molecular_assistant_dialog_job(dialog, 73, canvas)
        window._on_molecular_assistant_job_finished(73, _success_result())
        dialog.reject()
        QApplication.processEvents()
        assert _editor_snapshot(canvas) == before
        assert 73 not in window._molecular_assistant_results

        late_dialog = _open_dialog(window)
        window._register_molecular_assistant_dialog_job(late_dialog, 74, canvas)
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
