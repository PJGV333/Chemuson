"""Small modeless UI for requesting and reviewing an AI structure proposal."""

from __future__ import annotations

from PyQt6.QtCore import Qt, pyqtSignal
from PyQt6.QtWidgets import (
    QCheckBox,
    QDialog,
    QFormLayout,
    QHBoxLayout,
    QLabel,
    QLineEdit,
    QPlainTextEdit,
    QPushButton,
    QVBoxLayout,
)


class MolecularAssistantDialog(QDialog):
    """Collect one explicit request/configuration and review before insertion."""

    generation_requested = pyqtSignal(str, str, str, str, bool)
    insert_requested = pyqtSignal()

    def __init__(self, parent=None) -> None:
        super().__init__(parent)
        self.setWindowTitle("Generar estructura con IA")
        self.setModal(False)
        self.setMinimumWidth(520)
        self.setAttribute(Qt.WidgetAttribute.WA_DeleteOnClose, True)
        self._job_id: int | None = None

        self.description_edit = QPlainTextEdit(self)
        self.description_edit.setPlaceholderText("Dibuja cafeína…")
        self.description_edit.document().setMaximumBlockCount(100)
        self.description_edit.setMinimumHeight(76)

        self.base_url_edit = QLineEdit(self)
        self.base_url_edit.setPlaceholderText("https://servidor/v1 o http://localhost:1234/v1")
        self.model_edit = QLineEdit(self)
        self.api_key_edit = QLineEdit(self)
        self.api_key_edit.setEchoMode(QLineEdit.EchoMode.Password)
        self.json_output_check = QCheckBox("Solicitar salida JSON estructurada si el endpoint la admite")

        form = QFormLayout()
        form.addRow("Descripción", self.description_edit)
        form.addRow("Endpoint base", self.base_url_edit)
        form.addRow("Modelo", self.model_edit)
        form.addRow("API key (opcional)", self.api_key_edit)
        form.addRow("", self.json_output_check)

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
        self.provenance_label.setVisible(False)
        self.smiles_preview = QPlainTextEdit(self)
        self.smiles_preview.setReadOnly(True)
        self.smiles_preview.setMaximumHeight(100)
        self.smiles_preview.setVisible(False)

        self.generate_button = QPushButton("Generar y validar", self)
        self.insert_button = QPushButton("Insertar en el documento", self)
        self.insert_button.setVisible(False)
        self.close_button = QPushButton("Cerrar", self)

        buttons = QHBoxLayout()
        buttons.addStretch(1)
        buttons.addWidget(self.generate_button)
        buttons.addWidget(self.insert_button)
        buttons.addWidget(self.close_button)

        layout = QVBoxLayout(self)
        layout.addLayout(form)
        layout.addWidget(self.status_label)
        layout.addWidget(self.provenance_label)
        layout.addWidget(self.smiles_preview)
        layout.addWidget(self.preview_group)
        layout.addLayout(buttons)

        self.generate_button.clicked.connect(self._submit)
        self.insert_button.clicked.connect(self.insert_requested.emit)
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
        self.status_label.setText("Generando y validando en segundo plano…")

    def show_configuration_error(self) -> None:
        """Report invalid local configuration without exposing provider details."""
        self.generate_button.setEnabled(True)
        self.insert_button.setVisible(False)
        self.status_label.setText("Revisa la URL HTTP(S) y el identificador del modelo.")

    def show_failure(self, status: str, reason_code: str) -> None:
        """Present stable result identifiers only; raw diagnostics are never shown."""
        self.generate_button.setEnabled(True)
        self.insert_button.setVisible(False)
        self.provenance_label.setVisible(False)
        self.smiles_preview.setVisible(False)
        self.preview_group.setVisible(False)
        self.status_label.setText(
            f"No se generó una estructura válida. Estado: {status}; "
            f"motivo: {reason_code}."
        )

    def show_preview(
        self,
        *,
        provider_id: str,
        model_id: str | None,
        smiles: str,
    ) -> None:
        """Show the validated proposal, provenance, and semantic limitation."""
        self.generate_button.setEnabled(False)
        self.provenance_label.setText(
            f"Validación ChemIO: aceptada · Proveedor: {provider_id} · "
            f"Modelo: {model_id or 'N/D'}"
        )
        self.provenance_label.setVisible(True)
        self.smiles_preview.setPlainText(smiles)
        self.smiles_preview.setVisible(True)
        self.preview_group.setVisible(True)
        self.insert_button.setVisible(True)
        self.insert_button.setEnabled(True)
        self.status_label.setText("Revisa la propuesta. No se insertará hasta que lo confirmes.")

    def show_insert_notice(self, message: str) -> None:
        """Show a local insertion constraint while keeping the valid preview."""
        self.status_label.setText(message)

    def _submit(self) -> None:
        description = self.description_edit.toPlainText().strip()
        base_url = self.base_url_edit.text().strip()
        model = self.model_edit.text().strip()
        if not description:
            self.status_label.setText("Escribe una descripción molecular.")
            return
        if not base_url or not model:
            self.status_label.setText("Indica explícitamente el endpoint base y el modelo.")
            return
        self.status_label.setText("Preparando solicitud…")
        self.generation_requested.emit(
            description,
            base_url,
            model,
            self.api_key_edit.text(),
            self.json_output_check.isChecked(),
        )
