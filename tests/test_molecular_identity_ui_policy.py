from __future__ import annotations

import pytest
from PyQt6.QtWidgets import QApplication

from chemuson.gui.dialogs import MolecularAssistantDialog
from chemuson.molecular_assistant import OPENAI_COMPATIBLE_PROFILES


@pytest.fixture(scope="module", autouse=True)
def _qapp():
    return QApplication.instance() or QApplication([])


def test_identity_policy_defaults_offline_and_external_lookup_is_explicit():
    dialog = MolecularAssistantDialog(profiles=OPENAI_COMPATIBLE_PROFILES)
    try:
        assert dialog.identity_verification_enabled is True
        assert dialog.allow_external_identity_reference is False
        assert "extern" in dialog.external_identity_check.text().lower()
        assert "extern" in dialog.external_identity_check.toolTip().lower()

        dialog.external_identity_check.setChecked(True)
        assert dialog.allow_external_identity_reference is True
        dialog.identity_verification_check.setChecked(False)
        assert not dialog.external_identity_check.isEnabled()
        assert dialog.allow_external_identity_reference is False
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
        assert "no se consultaron fuentes externas" in dialog.identity_label.text()

        dialog.show_preview(
            provider_id="test-provider",
            model_id="test-model",
            smiles="CCO",
            identity_status="not_applicable",
            identity_reason_code="verification_disabled",
        )
        assert "verificación está desactivada" in dialog.identity_label.text()
    finally:
        dialog.close()
