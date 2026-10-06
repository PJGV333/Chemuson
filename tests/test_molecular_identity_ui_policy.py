from __future__ import annotations

import pytest
from PyQt6.QtTest import QSignalSpy
from PyQt6.QtWidgets import QApplication

from chemuson.gui.dialogs import MolecularAssistantDialog
from chemuson.molecular_assistant import OPENAI_COMPATIBLE_PROFILES


@pytest.fixture(scope="module", autouse=True)
def _qapp():
    return QApplication.instance() or QApplication([])


def test_resolution_methods_and_external_reference_permission_are_explicit():
    dialog = MolecularAssistantDialog(profiles=OPENAI_COMPATIBLE_PROFILES)
    try:
        assert dialog.resolution_method == "ai_reference"
        assert dialog.allow_external_reference is False
        assert "pubchem" in dialog.external_identity_check.text().lower()
        assert "sólo se envía el nombre químico" in dialog.external_identity_check.toolTip().lower()
        assert dialog.external_identity_check.isEnabled()

        dialog.external_identity_check.setChecked(True)
        assert dialog.allow_external_reference is True
        dialog.resolution_combo.setCurrentIndex(dialog.resolution_combo.findData("ai"))
        assert dialog.resolution_method == "ai"
        assert not dialog.external_identity_check.isEnabled()
        assert dialog.allow_external_reference is True
        dialog.resolution_combo.setCurrentIndex(dialog.resolution_combo.findData("reference"))
        assert dialog.resolution_method == "reference"
        assert dialog.allow_external_reference is True
    finally:
        dialog.close()


def test_reference_only_dialog_does_not_require_provider_endpoint_or_model():
    dialog = MolecularAssistantDialog(profiles=OPENAI_COMPATIBLE_PROFILES)
    spy = QSignalSpy(dialog.generation_requested)
    try:
        dialog.description_edit.setPlainText("Dibuja tetrandrina")
        dialog.api_key_edit.setText("must-not-reach-reference-connector")
        dialog.base_url_edit.clear()
        dialog.model_edit.clear()
        dialog.resolution_combo.setCurrentIndex(dialog.resolution_combo.findData("reference"))
        dialog._submit()

        assert len(spy) == 1
        assert spy[0][2] == ""
        assert spy[0][3] == ""
        assert spy[0][4] == ""

        dialog.resolution_combo.setCurrentIndex(dialog.resolution_combo.findData("ai"))
        dialog._submit()
        assert len(spy) == 1
        assert "endpoint base y el modelo" in dialog.status_label.text()
    finally:
        dialog.close()


def test_double_failure_keeps_ai_and_reference_reasons_separate():
    dialog = MolecularAssistantDialog(profiles=OPENAI_COMPATIBLE_PROFILES)
    try:
        dialog.show_failure(
            "provider_error",
            "timeout",
            reference_failure="reference_not_found",
        )
        assert "no respondió antes del límite" in dialog.status_label.text()
        assert "No se encontró una referencia química utilizable" in dialog.status_label.text()
    finally:
        dialog.close()


def test_reference_fallback_preview_keeps_ai_failure_and_provenance_distinct():
    dialog = MolecularAssistantDialog(profiles=OPENAI_COMPATIBLE_PROFILES)
    try:
        dialog.show_preview(
            provider_id="reference",
            model_id=None,
            smiles="CCO",
            identity_status="not_applicable",
            identity_reason_code="proposal_unavailable",
            structure_origin="reference",
            reference_source="pubchem",
            reference_resolved_name="Tetrandrine",
            reference_from_cache=True,
            ai_failure_reason="generation_exhausted",
            completion_tokens=4096,
            reasoning_tokens=4096,
            ai_provider_id="llama-cpp",
            ai_model_id="qwen-local",
        )

        assert "Fuente principal: referencia química" in dialog.provenance_label.text()
        assert "PubChem" in dialog.provenance_label.text()
        assert "Tetrandrine" in dialog.provenance_label.text()
        assert "IA: llama-cpp/qwen-local" in dialog.provenance_label.text()
        assert "La IA agotó el límite de generación" in dialog.status_label.text()
        assert dialog.smiles_preview.toPlainText() == "CCO"
        assert dialog.insert_button.text() == "Insertar referencia química"
        assert "completion tokens: 4096" in dialog.output_diagnostic_label.text()
        assert "reasoning tokens: 4096" in dialog.output_diagnostic_label.text()
    finally:
        dialog.close()


def test_mismatch_preview_shows_both_smiles_and_recommends_reference():
    dialog = MolecularAssistantDialog(profiles=OPENAI_COMPATIBLE_PROFILES)
    try:
        dialog.show_preview(
            provider_id="llama-cpp",
            model_id="model",
            smiles="CCN",
            ai_smiles="CCN",
            reference_smiles="CCO",
            reference_source="pubchem",
            reference_resolved_name="caffeine",
            structure_origin="ai_mismatch_reference",
            identity_status="mismatch",
            requested_name="caffeine",
        )

        assert dialog.smiles_preview.toPlainText() == "CCN"
        assert dialog.reference_smiles_preview.toPlainText() == "CCO"
        assert "IA (llama-cpp/model)" in dialog.provenance_label.text()
        assert not dialog.reference_smiles_preview.isHidden()
        assert dialog.insert_button.text() == "Usar referencia PubChem"
        assert dialog.close_button.text() == "Cancelar"
        assert not dialog.use_ai_proposal_button.isHidden()
        assert dialog._insert_candidate == "reference"
    finally:
        dialog.close()


def test_dialog_explains_offline_unverified_and_disabled_identity_states():
    dialog = MolecularAssistantDialog(profiles=OPENAI_COMPATIBLE_PROFILES)
    try:
        dialog.show_preview(
            provider_id="test-provider",
            model_id="test-model",
            smiles="CCO",
            identity_status="unverified",
            identity_reason_code="reference_not_found_offline",
        )
        assert "no se consultó la red" in dialog.identity_label.text()

        dialog.show_preview(
            provider_id="reference",
            model_id=None,
            smiles="CCO",
            identity_status="not_applicable",
            identity_reason_code="reference_selected",
            structure_origin="reference",
            reference_source="pubchem",
            reference_resolved_name="ethanol",
        )
        assert "no es una propuesta generada por IA" in dialog.identity_label.text()
        assert "PubChem" in dialog.provenance_label.text()
    finally:
        dialog.close()
