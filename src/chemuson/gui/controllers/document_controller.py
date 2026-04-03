from __future__ import annotations

import os
from typing import Optional


class DocumentController:
    """Operaciones de documento separadas de ChemusonWindow para test no-GUI."""

    @staticmethod
    def load_recent_files(window) -> list[str]:
        value = window._settings.value("recent_files", [])
        if isinstance(value, str):
            return [value]
        if isinstance(value, list):
            return [str(p) for p in value if p]
        return []

    @staticmethod
    def save_recent_files(window) -> None:
        window._settings.setValue("recent_files", window._recent_files)

    @staticmethod
    def set_canvas_file_path(window, canvas, filepath: Optional[str]) -> None:
        clean_path = os.path.abspath(filepath) if filepath else None
        window._canvas_file_paths[canvas] = clean_path
        autosave_manager = window._canvas_autosave_managers.get(canvas)
        if autosave_manager is not None:
            autosave_manager.set_original_path(clean_path)
        if clean_path:
            window._canvas_tab_titles[canvas] = os.path.basename(clean_path)
        if canvas is window.canvas:
            window._current_file_path = clean_path
        window._update_tab_title(canvas)

    @staticmethod
    def update_tab_title(window, canvas) -> None:
        if not window._tabs_alive():
            return
        try:
            index = window.tabs.indexOf(canvas)
            if index < 0:
                return
            base_title = window._canvas_tab_titles.get(canvas, "Sin título")
            dirty_suffix = " *" if not canvas.undo_stack.isClean() else ""
            window.tabs.setTabText(index, f"{base_title}{dirty_suffix}")
            path = window._canvas_file_paths.get(canvas)
            window.tabs.setTabToolTip(index, path or "Documento sin guardar")
        except RuntimeError:
            return

    @staticmethod
    def add_recent_file(window, filepath: str) -> None:
        if not filepath:
            return
        path = os.path.abspath(filepath)
        window._recent_files = [p for p in window._recent_files if p != path]
        window._recent_files.insert(0, path)
        window._recent_files = window._recent_files[:10]
        window._save_recent_files()
        window._update_recent_menu()

    @staticmethod
    def confirm_discard_changes(window, canvas=None) -> bool:
        from PyQt6.QtWidgets import QMessageBox

        target = canvas or window.canvas
        if target.undo_stack.isClean():
            return True
        index = window.tabs.indexOf(target)
        title = window._canvas_tab_titles.get(target, "Documento")
        if index >= 0:
            title = window.tabs.tabText(index).replace(" *", "")
        reply = QMessageBox.question(
            window,
            "Cambios sin guardar",
            f"'{title}' tiene cambios sin guardar. ¿Deseas guardar antes de cerrar?",
            QMessageBox.StandardButton.Save
            | QMessageBox.StandardButton.Discard
            | QMessageBox.StandardButton.Cancel,
            QMessageBox.StandardButton.Save,
        )
        if reply == QMessageBox.StandardButton.Save:
            previous_index = window.tabs.currentIndex()
            if index >= 0 and previous_index != index:
                window.tabs.setCurrentIndex(index)
                window._set_active_canvas(target)
            window._on_file_save()
            saved = target.undo_stack.isClean()
            if (
                previous_index >= 0
                and previous_index < window.tabs.count()
                and window.tabs.currentIndex() != previous_index
            ):
                window.tabs.setCurrentIndex(previous_index)
            return saved
        if reply == QMessageBox.StandardButton.Discard:
            return True
        return False
