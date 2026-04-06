from __future__ import annotations

import os
from dataclasses import dataclass
from typing import Callable, Optional


@dataclass(slots=True)
class RecentFilesContext:
    """Dependencias mínimas para manejar archivos recientes."""

    settings: object
    recent_files: list[str]
    persist_recent_files: Callable[[], None]
    refresh_recent_menu: Callable[[], None]


@dataclass(slots=True)
class DocumentDiscardContext:
    """Dependencias mínimas para confirmar descarte de cambios."""

    parent: object
    canvas: object
    tab_manager: object
    save_canvas: Callable[[], None]
    activate_canvas: Callable[[object], None]


@dataclass(slots=True)
class DocumentTabsContext:
    """Compatibilidad temporal para mapear canvas a pestañas y rutas."""

    tab_manager: object
    current_canvas: object | None = None
    current_file_path_setter: Callable[[Optional[str]], None] | None = None


class DocumentController:
    """Operaciones de documento separadas de ChemusonWindow para test no-GUI."""

    @staticmethod
    def _settings_from(target) -> object:
        return getattr(target, "_settings", target)

    @classmethod
    def load_recent_files(cls, target) -> list[str]:
        settings = cls._settings_from(target)
        value = settings.value("recent_files", [])
        if isinstance(value, str):
            return [value]
        if isinstance(value, list):
            return [str(p) for p in value if p]
        return []

    @classmethod
    def save_recent_files(cls, context: RecentFilesContext) -> None:
        context.settings.setValue("recent_files", context.recent_files)

    @staticmethod
    def set_canvas_file_path(target, canvas, filepath: Optional[str]) -> None:
        clean_path = os.path.abspath(filepath) if filepath else None
        if isinstance(target, DocumentTabsContext):
            target.tab_manager.set_canvas_file_path(canvas, clean_path)
            if canvas is target.current_canvas and target.current_file_path_setter is not None:
                target.current_file_path_setter(clean_path)
            return

        manager = getattr(target, "_tab_manager", None)
        if manager is None:
            raise TypeError("DocumentController.set_canvas_file_path requiere DocumentTabsContext")
        manager.set_canvas_file_path(canvas, clean_path)
        if canvas is getattr(target, "canvas", None):
            target._current_file_path = clean_path
        updater = getattr(target, "_update_tab_title", None)
        if callable(updater):
            updater(canvas)

    @staticmethod
    def update_tab_title(target, canvas) -> None:
        if isinstance(target, DocumentTabsContext):
            target.tab_manager.update_tab_title(canvas)
            return

        manager = getattr(target, "_tab_manager", None)
        if manager is None:
            raise TypeError("DocumentController.update_tab_title requiere DocumentTabsContext")
        manager.update_tab_title(canvas)

    @classmethod
    def add_recent_file(cls, context: RecentFilesContext, filepath: str) -> None:
        if not filepath:
            return
        path = os.path.abspath(filepath)
        context.recent_files[:] = [p for p in context.recent_files if p != path]
        context.recent_files.insert(0, path)
        del context.recent_files[10:]
        context.persist_recent_files()
        context.refresh_recent_menu()

    @classmethod
    def confirm_discard_changes(
        cls,
        context: DocumentDiscardContext,
        canvas=None,
    ) -> bool:
        from PyQt6.QtWidgets import QMessageBox

        target_canvas = canvas or context.canvas
        if target_canvas.undo_stack.isClean():
            return True
        tabs = context.tab_manager.tabs
        index = tabs.indexOf(target_canvas)
        title = context.tab_manager.tab_titles.get(target_canvas, "Documento")
        if index >= 0:
            title = tabs.tabText(index).replace(" *", "")
        reply = QMessageBox.question(
            context.parent,
            "Cambios sin guardar",
            f"'{title}' tiene cambios sin guardar. ¿Deseas guardar antes de cerrar?",
            QMessageBox.StandardButton.Save
            | QMessageBox.StandardButton.Discard
            | QMessageBox.StandardButton.Cancel,
            QMessageBox.StandardButton.Save,
        )
        if reply == QMessageBox.StandardButton.Save:
            previous_index = tabs.currentIndex()
            if index >= 0 and previous_index != index:
                tabs.setCurrentIndex(index)
                context.activate_canvas(target_canvas)
            context.save_canvas()
            saved = target_canvas.undo_stack.isClean()
            if (
                previous_index >= 0
                and previous_index < tabs.count()
                and tabs.currentIndex() != previous_index
            ):
                tabs.setCurrentIndex(previous_index)
                previous_canvas = context.tab_manager.canvas_from_tab_index(previous_index)
                if previous_canvas is not None:
                    context.activate_canvas(previous_canvas)
            return saved
        if reply == QMessageBox.StandardButton.Discard:
            return True
        return False
