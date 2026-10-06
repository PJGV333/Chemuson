"""Small modeless UI for requesting and reviewing an AI structure proposal."""

from __future__ import annotations

from collections.abc import Mapping, Sequence
from typing import Protocol

from PyQt6.QtCore import QElapsedTimer, Qt, QTimer, pyqtSignal
from PyQt6.QtWidgets import (
    QCheckBox,
    QComboBox,
    QDialog,
    QFormLayout,
    QGroupBox,
    QHBoxLayout,
    QLabel,
    QLineEdit,
    QPlainTextEdit,
    QPushButton,
    QSpinBox,
    QVBoxLayout,
)

class ProviderProfileView(Protocol):
    """Read-only profile fields injected by the M23-owning controller."""

    profile_id: str
    display_name: str
    default_base_url: str
    api_key_required: bool


class MolecularAssistantDialog(QDialog):
    """Collect one explicit request/configuration and review before insertion."""

    generation_requested = pyqtSignal(str, str, str, str, str, bool, int, int)
    insert_requested = pyqtSignal()
    insert_variant_requested = pyqtSignal()

    def __init__(
        self,
        parent=None,
        *,
        profiles: Sequence[ProviderProfileView],
        transform_mode: bool = False,
        profile_preferences: Mapping[str, Mapping[str, object]] | None = None,
    ) -> None:
        super().__init__(parent)
        self._profiles_by_id = {profile.profile_id: profile for profile in profiles}
        self._profile_preferences = {
            str(profile_id): dict(values)
            for profile_id, values in (profile_preferences or {}).items()
        }
        self._active_profile_id: str | None = None
        self._transform_mode = bool(transform_mode)
        self.setWindowTitle(
            "Transformar molécula con IA"
            if self._transform_mode
            else "Generar estructura con IA"
        )
        self.setModal(False)
        self.setMinimumWidth(520)
        self.setAttribute(Qt.WidgetAttribute.WA_DeleteOnClose, True)
        self._job_id: int | None = None

        self.description_edit = QPlainTextEdit(self)
        self.description_edit.setPlaceholderText(
            "Describe cómo transformar la molécula…"
            if self._transform_mode
            else "Dibuja cafeína…"
        )
        self.description_edit.document().setMaximumBlockCount(100)
        self.description_edit.setMinimumHeight(76)

        self.provider_combo = QComboBox(self)
        for profile in profiles:
            self.provider_combo.addItem(profile.display_name, profile.profile_id)
        self.base_url_edit = QLineEdit(self)
        self.base_url_edit.setPlaceholderText("https://servidor/v1 o http://localhost:1234/v1")
        self.model_edit = QLineEdit(self)
        self.model_edit.setPlaceholderText("ID del modelo expuesto por el endpoint")
        self.api_key_label = QLabel("API key (opcional)", self)
        self.api_key_edit = QLineEdit(self)
        self.api_key_edit.setEchoMode(QLineEdit.EchoMode.Password)
        self.timeout_spin = QSpinBox(self)
        self.timeout_spin.setRange(10, 600)
        self.timeout_spin.setValue(60)
        self.timeout_spin.setSuffix(" s")
        self.max_tokens_spin = QSpinBox(self)
        self.max_tokens_spin.setRange(64, 8192)
        self.max_tokens_spin.setValue(4096)
        self.json_output_check = QCheckBox("Solicitar salida JSON estructurada si el endpoint la admite")

        advanced_group = QGroupBox("Opciones avanzadas", self)
        advanced_form = QFormLayout(advanced_group)
        advanced_form.addRow("Timeout (s)", self.timeout_spin)
        advanced_form.addRow("Máximo de tokens de salida", self.max_tokens_spin)
        advanced_form.addRow("", self.json_output_check)

        form = QFormLayout()
        form.addRow("Descripción", self.description_edit)
        form.addRow("Proveedor", self.provider_combo)
        form.addRow("Endpoint base", self.base_url_edit)
        form.addRow("Modelo", self.model_edit)
        form.addRow(self.api_key_label, self.api_key_edit)

        self.status_label = QLabel(self)
        self.status_label.setWordWrap(True)
        self.status_label.setText("La estructura se validará antes de poder insertarla.")

        self.preview_group = QLabel(
            "La aceptación del parser confirma que la estructura se pudo interpretar; "
            "no demuestra que corresponda científicamente a la descripción.",
            self,
        )
        self.preview_group.setWordWrap(True)
        self.preview_group.setVisible(False)

        self.provenance_label = QLabel(self)
        self.provenance_label.setTextInteractionFlags(Qt.TextInteractionFlag.TextSelectableByMouse)
        self.identity_label = QLabel(self)
        self.identity_label.setWordWrap(True)
        self.identity_label.setVisible(False)
        self._elapsed_timer = QElapsedTimer()
        self._elapsed_update = QTimer(self)
        self._elapsed_update.setInterval(1000)
        self._elapsed_update.timeout.connect(self._update_elapsed_status)
        self.provenance_label.setVisible(False)
        self.source_smiles_label = QLabel("Molécula original (SMILES)", self)
        self.source_smiles_label.setVisible(False)
        self.source_smiles_preview = QPlainTextEdit(self)
        self.source_smiles_preview.setReadOnly(True)
        self.source_smiles_preview.setMaximumHeight(72)
        self.source_smiles_preview.setVisible(False)
        self.proposal_smiles_label = QLabel("Estructura propuesta (SMILES)", self)
        self.proposal_smiles_label.setVisible(False)
        self.smiles_preview = QPlainTextEdit(self)
        self.smiles_preview.setReadOnly(True)
        self.smiles_preview.setMaximumHeight(100)
        self.smiles_preview.setVisible(False)

        self.generate_button = QPushButton("Generar y validar", self)
        self.insert_button = QPushButton(
            "Reemplazar molécula seleccionada"
            if self._transform_mode
            else "Insertar en el documento",
            self,
        )
        self.insert_button.setVisible(False)
        self.insert_variant_button = QPushButton("Insertar variante", self)
        self.insert_variant_button.setVisible(False)
        self.close_button = QPushButton(
            "Cancelar" if self._transform_mode else "Cerrar",
            self,
        )

        buttons = QHBoxLayout()
        buttons.addStretch(1)
        buttons.addWidget(self.generate_button)
        buttons.addWidget(self.insert_variant_button)
        buttons.addWidget(self.insert_button)
        buttons.addWidget(self.close_button)

        layout = QVBoxLayout(self)
        layout.addLayout(form)
        layout.addWidget(advanced_group)
        layout.addWidget(self.status_label)
        layout.addWidget(self.provenance_label)
        layout.addWidget(self.identity_label)
        layout.addWidget(self.source_smiles_label)
        layout.addWidget(self.source_smiles_preview)
        layout.addWidget(self.proposal_smiles_label)
        layout.addWidget(self.smiles_preview)
        layout.addWidget(self.preview_group)
        layout.addLayout(buttons)

        self.provider_combo.currentIndexChanged.connect(self._on_provider_profile_changed)
        self._on_provider_profile_changed(clear_api_key=False)
        self.generate_button.clicked.connect(self._submit)
        self.insert_button.clicked.connect(self.insert_requested.emit)
        self.insert_variant_button.clicked.connect(self.insert_variant_requested.emit)
        self.close_button.clicked.connect(self.reject)

    @property
    def job_id(self) -> int | None:
        """Return the running or completed generation job identifier."""
        return self._job_id

    def set_job_id(self, job_id: int) -> None:
        self._job_id = int(job_id)

    def clear_api_key(self) -> None:
        """Remove the key from the visible form after it is handed to the worker."""
        self.api_key_edit.clear()

    def set_pending(self) -> None:
        """Show bounded background generation without blocking the window."""
        self.generate_button.setEnabled(False)
        self.insert_button.setVisible(False)
        self.insert_variant_button.setVisible(False)
        self._elapsed_timer.start()
        self._update_elapsed_status()
        self._elapsed_update.start()

    def _update_elapsed_status(self) -> None:
        if self._elapsed_timer.isValid():
            seconds = self._elapsed_timer.elapsed() // 1000
            self.status_label.setText(f"Generando y validando… {seconds} s")

    def _stop_elapsed_timer(self) -> None:
        self._elapsed_update.stop()
        self._elapsed_timer.invalidate()

    def show_configuration_error(self) -> None:
        """Report invalid local configuration without exposing provider details."""
        self._stop_elapsed_timer()
        self.generate_button.setEnabled(True)
        self.insert_button.setVisible(False)
        self.status_label.setText("Revisa el endpoint, el modelo y la API key requerida por el proveedor.")

    def show_failure(self, status: str, reason_code: str) -> None:
        """Translate a stable failure code without exposing raw diagnostics."""
        self._stop_elapsed_timer()
        self.generate_button.setEnabled(True)
        self.insert_button.setVisible(False)
        self.insert_variant_button.setVisible(False)
        self.identity_label.setVisible(False)
        self.provenance_label.setVisible(False)
        self.source_smiles_label.setVisible(False)
        self.source_smiles_preview.setVisible(False)
        self.proposal_smiles_label.setVisible(False)
        self.smiles_preview.setVisible(False)
        self.preview_group.setVisible(False)
        messages = {
            "timeout": (
                "El modelo no respondió antes del límite configurado "
                f"({self.timeout_spin.value()} s)."
            ),
            "invalid_json": "El modelo respondió, pero no respetó el formato estructurado requerido.",
            "network_error": "No fue posible contactar el endpoint configurado.",
            "http_error": "El endpoint rechazó la solicitud o devolvió un error HTTP.",
            "source_export_failed": "No fue posible exportar la molécula original para transformarla.",
            "invalid_smiles": "La estructura propuesta no pudo validarse químicamente.",
        }
        self.status_label.setText(
            messages.get(reason_code, "No se pudo generar una estructura válida.")
        )

    def show_preview(
        self,
        *,
        provider_id: str,
        model_id: str | None,
        smiles: str,
        source_smiles: str | None = None,
        identity_status: str = "unverified",
        requested_name: str | None = None,
        reference_identifier: str | None = None,
    ) -> None:
        """Show parser acceptance separately from semantic identity confidence."""
        self._identity_status = identity_status
        self._stop_elapsed_timer()
        self.generate_button.setEnabled(False)
        self.provenance_label.setText(
            f"SMILES válido (ChemIO): ✓ · Proveedor: {provider_id} · "
            f"Modelo: {model_id or 'N/D'}"
        )
        self.provenance_label.setVisible(True)
        identity_text = {
            "not_applicable": "Identidad molecular: no aplicable a una solicitud generativa abierta.",
            "verified": "Identidad solicitada: ✓ verificada mediante referencia química.",
            "mismatch": (
                "Identidad solicitada: ✗ la estructura propuesta es químicamente "
                "interpretable, pero no coincide con la referencia disponible"
                f"{f' para {requested_name}' if requested_name else ''}. "
                "Solo se insertará tras confirmación explícita."
            ),
            "reference_error": "Identidad: no verificada por un error al consultar/canonicalizar la referencia.",
            "unverified": "Identidad: no verificada; no hay una referencia confiable disponible.",
        }
        identity_suffix = (
            f" Referencia: {reference_identifier}." if reference_identifier else ""
        )
        self.identity_label.setText(
            identity_text.get(identity_status, identity_text["unverified"])
            + identity_suffix
        )
        self.identity_label.setVisible(True)
        if self._transform_mode and source_smiles is not None:
            self.source_smiles_preview.setPlainText(source_smiles)
            self.source_smiles_label.setVisible(True)
            self.source_smiles_preview.setVisible(True)
        self.proposal_smiles_label.setVisible(True)
        self.smiles_preview.setPlainText(smiles)
        self.smiles_preview.setVisible(True)
        self.preview_group.setVisible(True)
        self.insert_button.setVisible(True)
        self.insert_button.setEnabled(True)
        self.insert_variant_button.setVisible(self._transform_mode)
        self.insert_variant_button.setEnabled(self._transform_mode)
        if identity_status == "mismatch":
            self.insert_button.setText(
                "Reemplazar de todos modos"
                if self._transform_mode
                else "Insertar de todos modos"
            )
            if self._transform_mode:
                self.insert_variant_button.setText("Insertar variante de todos modos")
        else:
            self.insert_button.setText(
                "Reemplazar molécula seleccionada"
                if self._transform_mode
                else "Insertar en el documento"
            )
            self.insert_variant_button.setText("Insertar variante")
        self.status_label.setText(
            "Revisa la transformación. No se reemplazará la molécula hasta que lo confirmes."
            if self._transform_mode
            else "Revisa la propuesta. No se insertará hasta que lo confirmes."
        )

    def show_insert_notice(self, message: str) -> None:
        """Show a local insertion constraint while keeping the valid preview."""
        self.status_label.setText(message)

    def _on_provider_profile_changed(self, *_args, clear_api_key: bool = True) -> None:
        previous_profile_id = self._active_profile_id
        if previous_profile_id is not None:
            self._capture_current_profile_preferences(previous_profile_id)
        profile_id = str(self.provider_combo.currentData() or "")
        profile = self._profiles_by_id.get(profile_id)
        if profile is None:
            return
        if clear_api_key and previous_profile_id != profile_id:
            self.api_key_edit.clear()
        values = self._profile_preferences.get(profile_id, {})
        self.base_url_edit.setText(
            str(values.get("base_url", profile.default_base_url) or "")
        )
        self.model_edit.setText(str(values.get("model", "") or ""))
        self.timeout_spin.setValue(int(values.get("timeout_s", 60) or 60))
        self.max_tokens_spin.setValue(int(values.get("max_tokens", 4096) or 4096))
        self.json_output_check.setChecked(bool(values.get("supports_json_output", False)))
        self.api_key_label.setText(
            "API key (requerida)" if profile.api_key_required else "API key (opcional)"
        )
        self._active_profile_id = profile_id

    def _capture_current_profile_preferences(self, profile_id: str) -> None:
        self._profile_preferences[profile_id] = {
            "base_url": self.base_url_edit.text(),
            "model": self.model_edit.text(),
            "timeout_s": self.timeout_spin.value(),
            "max_tokens": self.max_tokens_spin.value(),
            "supports_json_output": self.json_output_check.isChecked(),
        }

    def current_profile_preferences(self) -> dict[str, object]:
        """Return the current profile's non-secret preferences for persistence."""
        profile_id = str(self.provider_combo.currentData() or "")
        self._capture_current_profile_preferences(profile_id)
        return dict(self._profile_preferences[profile_id])

    def _submit(self) -> None:
        description = self.description_edit.toPlainText().strip()
        provider_id = self.provider_combo.currentData()
        profile = self._profiles_by_id.get(provider_id)
        base_url = self.base_url_edit.text().strip()
        model = self.model_edit.text().strip()
        if not description:
            self.status_label.setText("Escribe una descripción molecular.")
            return
        if not base_url or not model:
            self.status_label.setText("Indica explícitamente el endpoint base y el modelo.")
            return
        if profile is None:
            self.show_configuration_error()
            return
        if profile.api_key_required and not self.api_key_edit.text().strip():
            self.status_label.setText("El proveedor seleccionado requiere una API key.")
            return
        self.status_label.setText("Preparando solicitud…")
        self.generation_requested.emit(
            description,
            provider_id,
            base_url,
            model,
            self.api_key_edit.text(),
            self.json_output_check.isChecked(),
            self.timeout_spin.value(),
            self.max_tokens_spin.value(),
        )
