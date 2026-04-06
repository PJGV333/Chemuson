from __future__ import annotations

from datetime import datetime
import json
import os
from typing import Callable, Optional

from PyQt6.QtWidgets import (
    QAbstractItemView,
    QDialog,
    QDialogButtonBox,
    QHeaderView,
    QHBoxLayout,
    QLabel,
    QMessageBox,
    QPushButton,
    QSizePolicy,
    QTableWidget,
    QTableWidgetItem,
    QTabWidget,
    QVBoxLayout,
    QWidget,
)

from chemuson.chemio.persistence import PersistenceManager
from chemuson.utils.autosave import AutosaveManager


class RecoveryController:
    @staticmethod
    def read_autosave_metadata(filepath: str) -> Optional[dict]:
        try:
            with open(filepath, "r", encoding="utf-8") as f:
                payload = json.load(f)
        except Exception:
            return None
        if payload.get("application") != "Chemuson":
            return None
        metadata = payload.get("autosave_metadata")
        if not isinstance(metadata, dict):
            metadata = {}
        raw_path = metadata.get("original_path")
        raw_timestamp = metadata.get("timestamp")
        return {
            "autosave_path": filepath,
            "original_path": str(raw_path) if raw_path else None,
            "timestamp": str(raw_timestamp) if raw_timestamp else "Desconocida",
        }

    @staticmethod
    def list_autosave_entries(directory: str) -> list[dict]:
        if not os.path.isdir(directory):
            return []
        entries: list[dict] = []
        for name in sorted(os.listdir(directory)):
            if not name.endswith(".json"):
                continue
            filepath = os.path.join(directory, name)
            if not os.path.isfile(filepath):
                continue
            metadata = RecoveryController.read_autosave_metadata(filepath)
            if metadata is None:
                continue
            metadata["filename"] = name
            entries.append(metadata)
        entries.sort(key=lambda entry: os.path.getmtime(entry["autosave_path"]), reverse=True)
        return entries

    @staticmethod
    def archive_autosave(path: str, autosave_dir: str) -> str:
        old_dir = os.path.join(autosave_dir, "old")
        os.makedirs(old_dir, exist_ok=True)
        basename = os.path.basename(path)
        target = os.path.join(old_dir, basename)
        if os.path.exists(target):
            root, ext = os.path.splitext(basename)
            target = os.path.join(old_dir, f"{root}_{datetime.now().strftime('%Y%m%d_%H%M%S_%f')}{ext}")
        os.replace(path, target)
        return target

    def open_autosave_document(self, window, autosave_path: str, original_path: Optional[str] = None) -> bool:
        canvas = window._create_document_tab(make_current=True)
        window._apply_toolbar_defaults_to_canvas(canvas)
        try:
            PersistenceManager.load_from_file(autosave_path, canvas)
            canvas.undo_stack.setClean()
            window._set_canvas_file_path(canvas, original_path)
            window.tabs.setCurrentWidget(canvas)
            window._set_active_canvas(canvas)
            window._update_total_charge_indicator()
            return True
        except Exception as exc:
            if hasattr(window, "_before_canvas_discard"):
                window._before_canvas_discard(canvas)
            if hasattr(window, "_tab_manager"):
                window._tab_manager.discard_canvas(canvas)
            if window.tabs.count() == 0:
                replacement = window._create_document_tab(make_current=True)
                window._apply_toolbar_defaults_to_canvas(replacement)
                window._set_active_canvas(replacement)
            QMessageBox.warning(window, "Error al abrir autosave", f"No se pudo abrir el autosave:\n{exc}")
            return False

    def on_open_recovery_center(self, window) -> None:
        self.show_recovery_center(window, show_only_if_pending=False)

    def show_recovery_center(self, window, show_only_if_pending: bool = False) -> None:
        autosave_dir = AutosaveManager.default_autosave_dir()
        pending_entries = self.list_autosave_entries(autosave_dir)
        recovered_entries = self.list_autosave_entries(os.path.join(autosave_dir, "old"))
        recent_entries = [p for p in window._recent_files if os.path.exists(p)]
        if recent_entries != window._recent_files:
            window._recent_files = recent_entries
            window._save_recent_files()
            window._update_recent_menu()

        if show_only_if_pending and not pending_entries:
            return

        dialog = QDialog(window)
        dialog.setWindowTitle("Centro de recuperación")
        dialog.resize(980, 560)
        layout = QVBoxLayout(dialog)
        layout.addWidget(QLabel("Gestiona autosaves pendientes, sesiones recuperadas y tus archivos recientes."))

        tabs = QTabWidget(dialog)
        layout.addWidget(tabs)

        def _setup_table(table: QTableWidget, action_column: int, action_width: int) -> None:
            table.setEditTriggers(QAbstractItemView.EditTrigger.NoEditTriggers)
            table.setSelectionMode(QAbstractItemView.SelectionMode.NoSelection)
            table.verticalHeader().setVisible(False)
            table.verticalHeader().setDefaultSectionSize(48)
            table.horizontalHeader().setSectionResizeMode(0, QHeaderView.ResizeMode.Stretch)
            if table.columnCount() > 2:
                table.horizontalHeader().setSectionResizeMode(1, QHeaderView.ResizeMode.ResizeToContents)
            table.horizontalHeader().setSectionResizeMode(action_column, QHeaderView.ResizeMode.Fixed)
            table.setColumnWidth(action_column, action_width)

        def _set_actions_cell(table: QTableWidget, row: int, buttons: list[tuple[str, Callable]]) -> None:
            action_widget = QWidget(table)
            action_layout = QHBoxLayout(action_widget)
            action_layout.setContentsMargins(8, 6, 8, 6)
            action_layout.setSpacing(8)
            button_widgets: list[QPushButton] = []
            for text, callback in buttons:
                button = QPushButton(text, action_widget)
                button.setMinimumHeight(38)
                button.setMinimumWidth(136)
                button.setSizePolicy(QSizePolicy.Policy.Expanding, QSizePolicy.Policy.Fixed)
                button.clicked.connect(callback)
                action_layout.addWidget(button)
                button_widgets.append(button)
            table.setCellWidget(row, table.columnCount() - 1, action_widget)

        pending_tab = QWidget()
        pending_layout = QVBoxLayout(pending_tab)
        if not pending_entries:
            pending_layout.addWidget(QLabel("No hay autosaves pendientes."))
        else:
            pending_table = QTableWidget(len(pending_entries), 3, pending_tab)
            pending_table.setHorizontalHeaderLabels(["Archivo original", "Timestamp", "Acciones"])
            _setup_table(pending_table, action_column=2, action_width=360)

            def _mark_pending_done(row: int, status: str) -> None:
                pending_table.setCellWidget(row, 2, QLabel(status))

            for row, entry in enumerate(pending_entries):
                original_path_value = entry.get("original_path")
                original_path = original_path_value or "Sin nombre"
                timestamp = entry.get("timestamp") or "Desconocida"
                autosave_path = entry["autosave_path"]

                pending_table.setItem(row, 0, QTableWidgetItem(original_path))
                pending_table.setItem(row, 1, QTableWidgetItem(timestamp))

                def _recover_handler(_checked: bool = False, r: int = row, p: str = autosave_path, original: Optional[str] = original_path_value) -> None:
                    if not self.open_autosave_document(window, p, original_path=original):
                        return
                    try:
                        self.archive_autosave(p, autosave_dir)
                        _mark_pending_done(r, "Recuperado")
                        window.statusBar().showMessage("Autosave recuperado")
                    except Exception as exc:
                        QMessageBox.warning(window, "Error al mover autosave", f"No se pudo archivar el autosave:\n{exc}")

                def _discard_handler(_checked: bool = False, r: int = row, p: str = autosave_path) -> None:
                    try:
                        self.archive_autosave(p, autosave_dir)
                        _mark_pending_done(r, "Descartado")
                    except Exception as exc:
                        QMessageBox.warning(window, "Error al descartar", f"No se pudo mover el autosave:\n{exc}")

                _set_actions_cell(pending_table, row, [("Recuperar", _recover_handler), ("Descartar", _discard_handler)])
            pending_layout.addWidget(pending_table)
        tabs.addTab(pending_tab, f"Pendientes ({len(pending_entries)})")

        recent_tab = QWidget()
        recent_layout = QVBoxLayout(recent_tab)
        if not recent_entries:
            recent_layout.addWidget(QLabel("No hay archivos recientes disponibles."))
        else:
            recent_table = QTableWidget(len(recent_entries), 2, recent_tab)
            recent_table.setHorizontalHeaderLabels(["Archivo reciente", "Acciones"])
            _setup_table(recent_table, action_column=1, action_width=170)
            for row, filepath in enumerate(recent_entries):
                recent_table.setItem(row, 0, QTableWidgetItem(filepath))

                def _open_recent_handler(_checked: bool = False, p: str = filepath) -> None:
                    if not os.path.exists(p):
                        QMessageBox.warning(window, "Archivo no encontrado", "El archivo no existe.")
                        window._update_recent_menu()
                        return
                    window._open_file_path(p)

                _set_actions_cell(recent_table, row, [("Abrir", _open_recent_handler)])
            recent_layout.addWidget(recent_table)
        tabs.addTab(recent_tab, f"Recientes ({len(recent_entries)})")

        recovered_tab = QWidget()
        recovered_layout = QVBoxLayout(recovered_tab)
        if not recovered_entries:
            recovered_layout.addWidget(QLabel("No hay autosaves recuperados/archivados."))
        else:
            recovered_table = QTableWidget(len(recovered_entries), 3, recovered_tab)
            recovered_table.setHorizontalHeaderLabels(["Archivo original", "Timestamp", "Acciones"])
            _setup_table(recovered_table, action_column=2, action_width=170)
            for row, entry in enumerate(recovered_entries):
                original_path_value = entry.get("original_path")
                original_path = original_path_value or "Sin nombre"
                timestamp = entry.get("timestamp") or "Desconocida"
                autosave_path = entry["autosave_path"]
                recovered_table.setItem(row, 0, QTableWidgetItem(original_path))
                recovered_table.setItem(row, 1, QTableWidgetItem(timestamp))

                def _open_recovered_handler(_checked: bool = False, p: str = autosave_path, original: Optional[str] = original_path_value) -> None:
                    self.open_autosave_document(window, p, original_path=original)

                _set_actions_cell(recovered_table, row, [("Abrir", _open_recovered_handler)])
            recovered_layout.addWidget(recovered_table)
        tabs.addTab(recovered_tab, f"Recuperados ({len(recovered_entries)})")

        buttons = QDialogButtonBox(QDialogButtonBox.StandardButton.Close)
        buttons.rejected.connect(dialog.reject)
        buttons.accepted.connect(dialog.accept)
        layout.addWidget(buttons)
        dialog.exec()

    def check_autosaves(self, window) -> None:
        self.show_recovery_center(window, show_only_if_pending=True)
