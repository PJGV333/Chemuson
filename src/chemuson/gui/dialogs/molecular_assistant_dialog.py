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
    use_ai_proposal_requested = pyqtSignal()

    def __init__(
        self,
        parent=None,
        *,
        profiles: Sequence[ProviderProfileView],
        transform_mode: bool = False,
        profile_preferences: Mapping[str, Mapping[str, object]] | None = None,
        resolution_method: str = "ai_reference",
        allow_external_reference: bool = False,
        identity_verification_enabled: bool | None = None,
        allow_external_identity_reference: bool | None = None,
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
            else "Asistente molecular"
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

        self.resolution_combo = QComboBox(self)
        self.resolution_combo.addItem("IA + referencia", "ai_reference")
        self.resolution_combo.addItem("Solo IA", "ai")
        self.resolution_combo.addItem("Referencia química", "reference")
        if self._transform_mode:
            self.resolution_combo.setCurrentIndex(self.resolution_combo.findData("ai"))
            self.resolution_combo.setVisible(False)

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

        if allow_external_identity_reference is not None:
            allow_external_reference = bool(allow_external_identity_reference)
        self.external_identity_check = QCheckBox(
            "Permitir referencias químicas externas (PubChem)", self
        )
        self.external_identity_check.setToolTip(
            "Sólo se envía el nombre químico extraído de una petición explícita. "
            "No se envían el prompt completo, documentos, SMILES privados ni API keys. "
            "Sin permiso externo sólo se consultan referencias locales y la caché. "
            "Solo IA nunca consulta referencias ni accede a la red."
        )
        self.external_identity_check.setChecked(bool(allow_external_reference))
        self.reference_group = QGroupBox("Referencias químicas", self)
        reference_layout = QVBoxLayout(self.reference_group)
        reference_layout.addWidget(self.external_identity_check)

        self.advanced_group = QGroupBox("Opciones avanzadas", self)
        advanced_form = QFormLayout(self.advanced_group)
        advanced_form.addRow("Timeout (s)", self.timeout_spin)
        advanced_form.addRow("Máximo de tokens de salida", self.max_tokens_spin)
        advanced_form.addRow("", self.json_output_check)

        self.provider_group = QGroupBox("Proveedor IA", self)
        provider_form = QFormLayout(self.provider_group)
        provider_form.addRow("Proveedor", self.provider_combo)
        provider_form.addRow("Endpoint base", self.base_url_edit)
        provider_form.addRow("Modelo", self.model_edit)
        provider_form.addRow(self.api_key_label, self.api_key_edit)

        form = QFormLayout()
        form.addRow("Descripción", self.description_edit)
        if not self._transform_mode:
            form.addRow("Método de resolución", self.resolution_combo)

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
        self.provenance_label.setTextFormat(Qt.TextFormat.PlainText)
        self.provenance_label.setTextInteractionFlags(Qt.TextInteractionFlag.TextSelectableByMouse)
        self.output_diagnostic_label = QLabel(self)
        self.output_diagnostic_label.setWordWrap(True)
        self.output_diagnostic_label.setVisible(False)
        self.identity_label = QLabel(self)
        self.identity_label.setTextFormat(Qt.TextFormat.PlainText)
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
        self.reference_smiles_label = QLabel("Referencia química (SMILES)", self)
        self.reference_smiles_label.setVisible(False)
        self.reference_smiles_preview = QPlainTextEdit(self)
        self.reference_smiles_preview.setReadOnly(True)
        self.reference_smiles_preview.setMaximumHeight(100)
        self.reference_smiles_preview.setVisible(False)

        self.generate_button = QPushButton("Resolver y validar", self)
        self.insert_button = QPushButton(
            "Reemplazar molécula seleccionada"
            if self._transform_mode
            else "Insertar en el documento",
            self,
        )
        self.insert_button.setVisible(False)
        self.insert_variant_button = QPushButton("Insertar variante", self)
        self.insert_variant_button.setVisible(False)
        self.use_ai_proposal_button = QPushButton("Usar propuesta IA de todos modos", self)
        self.use_ai_proposal_button.setVisible(False)
        self.close_button = QPushButton(
            "Cancelar" if self._transform_mode else "Cerrar",
            self,
        )

        buttons = QHBoxLayout()
        buttons.addStretch(1)
        buttons.addWidget(self.generate_button)
        buttons.addWidget(self.use_ai_proposal_button)
        buttons.addWidget(self.insert_variant_button)
        buttons.addWidget(self.insert_button)
        buttons.addWidget(self.close_button)

        layout = QVBoxLayout(self)
        layout.addLayout(form)
        layout.addWidget(self.provider_group)
        layout.addWidget(self.reference_group)
        layout.addWidget(self.advanced_group)
        layout.addWidget(self.status_label)
        layout.addWidget(self.provenance_label)
        layout.addWidget(self.output_diagnostic_label)
        layout.addWidget(self.identity_label)
        layout.addWidget(self.source_smiles_label)
        layout.addWidget(self.source_smiles_preview)
        layout.addWidget(self.proposal_smiles_label)
        layout.addWidget(self.smiles_preview)
        layout.addWidget(self.reference_smiles_label)
        layout.addWidget(self.reference_smiles_preview)
        layout.addWidget(self.preview_group)
        layout.addLayout(buttons)

        self.provider_combo.currentIndexChanged.connect(self._on_provider_profile_changed)
        self.resolution_combo.currentIndexChanged.connect(self._on_resolution_method_changed)
        self._on_provider_profile_changed(clear_api_key=False)
        if resolution_method in {"ai", "ai_reference", "reference"}:
            self.resolution_combo.setCurrentIndex(self.resolution_combo.findData(resolution_method))
        self._on_resolution_method_changed()
        self.generate_button.clicked.connect(self._submit)
        self.insert_button.clicked.connect(self.insert_requested.emit)
        self.insert_variant_button.clicked.connect(self.insert_variant_requested.emit)
        self.use_ai_proposal_button.clicked.connect(self.use_ai_proposal_requested.emit)
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

    @property
    def resolution_method(self) -> str:
        """Return the selected route; whole-molecule transforms remain AI-only."""
        if self._transform_mode:
            return "ai"
        return str(self.resolution_combo.currentData() or "ai_reference")

    @property
    def identity_verification_enabled(self) -> bool:
        """Compatibility view: reference modes always reconcile identity when possible."""
        return self.resolution_method == "ai_reference"

    @property
    def allow_external_reference(self) -> bool:
        """Return the persisted permission; AI-only routes never invoke the resolver."""
        return not self._transform_mode and self.external_identity_check.isChecked()

    @property
    def allow_external_identity_reference(self) -> bool:
        """Compatibility alias for the external chemical-reference preference."""
        return self.allow_external_reference

    def _on_resolution_method_changed(self, *_args) -> None:
        method = self.resolution_method
        is_ai = method != "reference"
        self.provider_group.setVisible(is_ai)
        self.advanced_group.setVisible(is_ai)
        self.reference_group.setVisible(not self._transform_mode)
        self.external_identity_check.setEnabled(
            not self._transform_mode and method in {"ai_reference", "reference"}
        )
        self.generate_button.setText(
            "Resolver referencia" if method == "reference" else "Generar y validar"
        )

    def _restore_request_controls(self) -> None:
        self.resolution_combo.setEnabled(not self._transform_mode)
        self.provider_combo.setEnabled(True)
        self.base_url_edit.setEnabled(True)
        self.model_edit.setEnabled(True)
        self.api_key_edit.setEnabled(True)
        self._on_resolution_method_changed()

    def set_pending(self) -> None:
        """Show bounded background generation without blocking the window."""
        self.generate_button.setEnabled(False)
        self.insert_button.setVisible(False)
        self.insert_variant_button.setVisible(False)
        self.use_ai_proposal_button.setVisible(False)
        self.resolution_combo.setEnabled(False)
        self.external_identity_check.setEnabled(False)
        self.provider_combo.setEnabled(False)
        self.base_url_edit.setEnabled(False)
        self.model_edit.setEnabled(False)
        self.api_key_edit.setEnabled(False)
        self._elapsed_timer.start()
        self._update_elapsed_status()
        self._elapsed_update.start()

    def _update_elapsed_status(self) -> None:
        if self._elapsed_timer.isValid():
            seconds = self._elapsed_timer.elapsed() // 1000
            action = (
                "Resolviendo referencia"
                if self.resolution_method == "reference"
                else "Generando y validando"
            )
            self.status_label.setText(f"{action}… {seconds} s")

    def _stop_elapsed_timer(self) -> None:
        self._elapsed_update.stop()
        self._elapsed_timer.invalidate()

    def show_configuration_error(self) -> None:
        """Report invalid local configuration without exposing provider details."""
        self._stop_elapsed_timer()
        self.generate_button.setEnabled(True)
        self._restore_request_controls()
        self.insert_button.setVisible(False)
        self.status_label.setText(
            "Revisa endpoint, modelo y API key. Las credenciales sólo pueden enviarse "
            "a un endpoint HTTPS o a loopback local."
        )

    def show_failure(
        self,
        status: str,
        reason_code: str,
        *,
        structured_output_requested: bool = False,
        structured_output_native: bool | None = None,
        structured_output_fallback_used: bool = False,
        format_repair_used: bool = False,
        format_repair_succeeded: bool | None = None,
        finish_reason: str | None = None,
        completion_tokens: int | None = None,
        reasoning_tokens: int | None = None,
        reference_failure: str | None = None,
    ) -> None:
        """Translate stable codes and bounded diagnostics without raw provider data."""
        self._stop_elapsed_timer()
        self.generate_button.setEnabled(True)
        self._restore_request_controls()
        self.insert_button.setVisible(False)
        self.insert_variant_button.setVisible(False)
        self.use_ai_proposal_button.setVisible(False)
        self.identity_label.setVisible(False)
        self.provenance_label.setVisible(False)
        self._show_output_diagnostics(
            structured_output_requested=structured_output_requested,
            structured_output_native=structured_output_native,
            structured_output_fallback_used=structured_output_fallback_used,
            format_repair_used=format_repair_used,
            format_repair_succeeded=format_repair_succeeded,
        )
        self.source_smiles_label.setVisible(False)
        self.source_smiles_preview.setVisible(False)
        self.proposal_smiles_label.setVisible(False)
        self.smiles_preview.setVisible(False)
        self.reference_smiles_label.setVisible(False)
        self.reference_smiles_preview.setVisible(False)
        self.preview_group.setVisible(False)
        messages = {
            "timeout": (
                "El modelo no respondió antes del límite configurado "
                f"({self.timeout_spin.value()} s)."
            ),
            "generation_exhausted": (
                "El modelo agotó el límite de generación antes de producir una respuesta final."
            ),
            "reference_name_required": "Indica una petición explícita de una molécula por nombre.",
            "reference_not_found": "No se encontró una referencia química utilizable.",
            "reference_invalid": "La estructura de referencia no superó la validación ChemIO.",
            "reference_unavailable": "No fue posible resolver una referencia química.",
            "invalid_json": "El modelo respondió, pero no respetó el formato estructurado requerido.",
            "network_error": "No fue posible contactar el endpoint configurado.",
            "http_error": "El endpoint rechazó la solicitud o devolvió un error HTTP.",
            "structured_output_unsupported": "El endpoint no admite el modo JSON estructurado.",
            "source_export_failed": "No fue posible exportar la molécula original para transformarla.",
            "invalid_smiles": "La estructura propuesta no pudo validarse químicamente.",
        }
        message = messages.get(reason_code, "No se pudo generar una estructura válida.")
        reference_messages = {
            "reference_not_found": "No se encontró una referencia química utilizable.",
            "reference_error": "No se pudo consultar o validar una referencia química.",
            "invalid_reference": "La referencia disponible no superó ChemIO.",
        }
        if reference_failure in reference_messages:
            message += f" {reference_messages[reference_failure]}"
        self.status_label.setText(message)
        if finish_reason == "length" or completion_tokens is not None or reasoning_tokens is not None:
            self._show_output_diagnostics(
                structured_output_requested=structured_output_requested,
                structured_output_native=structured_output_native,
                structured_output_fallback_used=structured_output_fallback_used,
                format_repair_used=format_repair_used,
                format_repair_succeeded=format_repair_succeeded,
                finish_reason=finish_reason,
                completion_tokens=completion_tokens,
                reasoning_tokens=reasoning_tokens,
            )

    def show_preview(
        self,
        *,
        provider_id: str,
        model_id: str | None,
        smiles: str,
        source_smiles: str | None = None,
        identity_status: str = "unverified",
        identity_reason_code: str | None = None,
        structured_output_requested: bool = False,
        structured_output_native: bool | None = None,
        structured_output_fallback_used: bool = False,
        format_repair_used: bool = False,
        format_repair_succeeded: bool | None = None,
        requested_name: str | None = None,
        reference_identifier: str | None = None,
        structure_origin: str = "ai",
        ai_smiles: str | None = None,
        reference_smiles: str | None = None,
        reference_source: str | None = None,
        reference_resolved_name: str | None = None,
        reference_from_cache: bool = False,
        ai_failure_reason: str | None = None,
        completion_tokens: int | None = None,
        reasoning_tokens: int | None = None,
        ai_provider_id: str | None = None,
        ai_model_id: str | None = None,
    ) -> None:
        """Show AI/reference origin, identity, and both candidates without mutation."""
        self._identity_status = identity_status
        self._insert_candidate = "reference" if structure_origin == "ai_mismatch_reference" else (
            "reference" if structure_origin == "reference" else "ai"
        )
        self._stop_elapsed_timer()
        self.generate_button.setEnabled(False)
        self._restore_request_controls()
        self.generate_button.setEnabled(False)
        self.external_identity_check.setEnabled(False)
        self.resolution_combo.setEnabled(False)
        self.provenance_label.setText(self._provenance_text(
            provider_id,
            model_id,
            structure_origin,
            reference_source,
            reference_resolved_name,
            reference_from_cache,
            ai_provider_id,
            ai_model_id,
        ))
        self.provenance_label.setVisible(True)
        self._show_output_diagnostics(
            structured_output_requested=structured_output_requested,
            structured_output_native=structured_output_native,
            structured_output_fallback_used=structured_output_fallback_used,
            format_repair_used=format_repair_used,
            format_repair_succeeded=format_repair_succeeded,
            finish_reason="length" if ai_failure_reason == "generation_exhausted" else None,
            completion_tokens=completion_tokens,
            reasoning_tokens=reasoning_tokens,
        )
        if identity_status == "verified":
            identity_text = "Identidad: ✓ verificada por comparación InChI aislada."
        elif identity_status == "mismatch":
            identity_text = (
                f"La propuesta IA no coincide con la referencia"
                f"{f' para {requested_name}' if requested_name else ''}. "
                "La referencia es la opción recomendada; la IA requiere override explícito."
            )
        elif identity_status == "reference_error":
            identity_text = "Identidad no verificada: la referencia no pudo compararse de forma segura."
        elif identity_status == "not_applicable" and structure_origin == "reference":
            identity_text = "Estructura de referencia seleccionada; no es una propuesta generada por IA."
        elif identity_status == "not_applicable":
            identity_text = "Identidad: no aplicable a una solicitud generativa abierta."
        elif identity_reason_code == "reference_not_requested":
            identity_text = "Identidad no verificada: método Solo IA, sin consulta de referencia."
        elif identity_reason_code == "reference_not_found_offline":
            identity_text = "Identidad no verificada: no hubo referencia local; no se consultó la red."
        else:
            identity_text = "Identidad no verificada: no hay una referencia química confiable disponible."
        if reference_identifier:
            identity_text += f" Referencia: {reference_identifier}."
        self.identity_label.setText(identity_text)
        self.identity_label.setVisible(True)

        if self._transform_mode and source_smiles is not None:
            self.source_smiles_preview.setPlainText(source_smiles)
            self.source_smiles_label.setVisible(True)
            self.source_smiles_preview.setVisible(True)

        ai_text = ai_smiles if ai_smiles is not None else smiles
        self.smiles_preview.setPlainText(ai_text if structure_origin != "reference" else smiles)
        self.proposal_smiles_label.setText(
            "Estructura propuesta (SMILES)"
            if self._transform_mode
            else "Propuesta IA (SMILES)"
            if structure_origin != "reference"
            else "Estructura de referencia (SMILES)"
        )
        self.proposal_smiles_label.setVisible(True)
        self.smiles_preview.setVisible(True)
        has_reference = bool(reference_smiles)
        self.reference_smiles_preview.setPlainText(reference_smiles or "")
        self.reference_smiles_label.setVisible(has_reference)
        self.reference_smiles_preview.setVisible(has_reference)
        self.preview_group.setVisible(True)

        self.use_ai_proposal_button.setVisible(identity_status == "mismatch" and not self._transform_mode)
        self.insert_variant_button.setVisible(self._transform_mode)
        self.insert_variant_button.setEnabled(self._transform_mode)
        self.insert_button.setVisible(True)
        self.insert_button.setEnabled(True)
        if self._transform_mode:
            self.insert_button.setText("Reemplazar molécula seleccionada")
            self.insert_variant_button.setText("Insertar variante")
            self.status_label.setText("Revisa la transformación; no se aplicará hasta que lo confirmes.")
        elif identity_status == "mismatch":
            self.close_button.setText("Cancelar")
            source_label = {
                "pubchem": "referencia PubChem",
                "offline-common": "referencia local",
            }.get(reference_source or "", "referencia química")
            self.insert_button.setText(f"Usar {source_label}")
            self.use_ai_proposal_button.setText("Usar propuesta IA de todos modos")
            self.status_label.setText("Compara ambas estructuras. La referencia es la opción recomendada.")
        elif structure_origin == "reference":
            self.close_button.setText("Cancelar")
            self.insert_button.setText("Insertar referencia química")
            self.status_label.setText(
                self._ai_failure_message(ai_failure_reason)
                if ai_failure_reason
                else "Revisa la referencia química antes de insertarla."
            )
        else:
            self.close_button.setText("Cerrar")
            self.insert_button.setText("Insertar propuesta IA")
            self.status_label.setText("Revisa la propuesta IA antes de insertarla.")

    @staticmethod
    def _ai_failure_message(reason_code: str) -> str:
        messages = {
            "timeout": "La IA no produjo una estructura utilizable (timeout). Se encontró una referencia química.",
            "generation_exhausted": "La IA agotó el límite de generación antes de producir una respuesta final. Se encontró una referencia química.",
            "invalid_json": "La IA no produjo una respuesta estructurada utilizable. Se encontró una referencia química.",
            "provider_error": "La IA no produjo una estructura utilizable. Se encontró una referencia química.",
            "network_error": "La IA no pudo contactar el endpoint. Se encontró una referencia química.",
            "http_error": "El endpoint de IA devolvió un error. Se encontró una referencia química.",
            "invalid_smiles": "La estructura de IA no superó ChemIO. Se encontró una referencia química.",
        }
        return messages.get(reason_code, "La IA no produjo una estructura utilizable. Se encontró una referencia química.")

    @staticmethod
    def _provenance_text(
        provider_id: str,
        model_id: str | None,
        structure_origin: str,
        reference_source: str | None,
        reference_resolved_name: str | None,
        reference_from_cache: bool,
        ai_provider_id: str | None,
        ai_model_id: str | None,
    ) -> str:
        source_labels = {"pubchem": "PubChem", "offline-common": "Referencia local"}
        source = source_labels.get(reference_source or "", "Referencia química")
        ref_name = reference_resolved_name or "nombre solicitado"
        cache = " · caché local" if reference_from_cache else ""
        if structure_origin == "reference":
            ai_source = (
                f" · IA: {ai_provider_id}/{ai_model_id or 'N/D'}"
                if ai_provider_id
                else ""
            )
            return f"Fuente principal: referencia química · {source} · {ref_name}{cache} · ChemIO: ✓{ai_source}"
        if structure_origin == "ai_verified_by_reference":
            return f"Fuente principal: IA · {provider_id}/{model_id or 'N/D'} · identidad verificada con {source}: {ref_name}"
        if structure_origin == "ai_mismatch_reference":
            ai_source = f"{ai_provider_id or provider_id}/{ai_model_id or model_id or 'N/D'}"
            return f"IA ({ai_source}) y referencia no coinciden · referencia recomendada: {source} · {ref_name}{cache}"
        return f"Fuente principal: IA · {provider_id}/{model_id or 'N/D'} · ChemIO: ✓ · identidad no verificada"

    def _show_output_diagnostics(
        self,
        *,
        structured_output_requested: bool,
        structured_output_native: bool | None,
        structured_output_fallback_used: bool,
        format_repair_used: bool,
        format_repair_succeeded: bool | None,
        finish_reason: str | None = None,
        completion_tokens: int | None = None,
        reasoning_tokens: int | None = None,
    ) -> None:
        messages: list[str] = []
        if structured_output_fallback_used:
            messages.append("JSON nativo no soportado; se usó el contrato textual")
        elif structured_output_native is True:
            messages.append("Salida JSON nativa: usada")
        elif structured_output_requested:
            messages.append("Salida JSON nativa: no confirmada")
        elif structured_output_native is False:
            messages.append("Contrato JSON textual: usado")
        if format_repair_used:
            messages.append(
                "Reintento de formato: exitoso"
                if format_repair_succeeded is True
                else "Reintento de formato: fallido"
            )
        if finish_reason == "length":
            counts = []
            if completion_tokens is not None:
                counts.append(f"completion tokens: {completion_tokens}")
            if reasoning_tokens is not None:
                counts.append(f"reasoning tokens: {reasoning_tokens}")
            messages.append("Generación agotada" + (f" ({'; '.join(counts)})" if counts else ""))
        self.output_diagnostic_label.setText(" · ".join(messages))
        self.output_diagnostic_label.setVisible(bool(messages))

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
        if self.resolution_method != "reference":
            if not base_url or not model:
                self.status_label.setText("Indica explícitamente el endpoint base y el modelo.")
                return
            if profile is None:
                self.show_configuration_error()
                return
            if profile.api_key_required and not self.api_key_edit.text().strip():
                self.status_label.setText("El proveedor seleccionado requiere una API key.")
                return
        else:
            provider_id = str(provider_id or "")
            base_url = ""
            model = ""
        self.status_label.setText("Preparando solicitud…")
        self.generation_requested.emit(
            description,
            provider_id,
            base_url,
            model,
            self.api_key_edit.text() if self.resolution_method != "reference" else "",
            self.json_output_check.isChecked(),
            self.timeout_spin.value(),
            self.max_tokens_spin.value(),
        )
