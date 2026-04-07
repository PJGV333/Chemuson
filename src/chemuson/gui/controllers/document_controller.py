from __future__ import annotations

import os
from dataclasses import dataclass
from typing import Callable, Optional

from chemuson.gui.tab_manager import CanvasTabManager


@dataclass(slots=True)
class RecentFilesContext:
    """Dependencias mínimas para manejar archivos recientes."""

    settings: object
    recent_files: list[str]
    persist_recent_files: Callable[[], None]
    refresh_recent_menu: Callable[[], None]

    def __post_init__(self) -> None:
        if not isinstance(self.recent_files, list):
            raise TypeError("RecentFilesContext.recent_files must be a list[str]")
        if not callable(self.persist_recent_files):
            raise TypeError("RecentFilesContext.persist_recent_files must be callable")
        if not callable(self.refresh_recent_menu):
            raise TypeError("RecentFilesContext.refresh_recent_menu must be callable")


@dataclass(slots=True)
class DocumentDiscardContext:
    """Dependencias mínimas para confirmar descarte de cambios."""

    parent: object
    canvas: object
    tab_manager: object
    save_canvas: Callable[[], None]
    activate_canvas: Callable[[object], None]

    def __post_init__(self) -> None:
        if not isinstance(self.tab_manager, CanvasTabManager):
            raise TypeError("DocumentDiscardContext.tab_manager must be a CanvasTabManager")
        if not callable(self.save_canvas):
            raise TypeError("DocumentDiscardContext.save_canvas must be callable")
        if not callable(self.activate_canvas):
            raise TypeError("DocumentDiscardContext.activate_canvas must be callable")


@dataclass(slots=True)
class DocumentTabsContext:
    """Contrato explícito para mapear canvas a pestañas y rutas."""

    tab_manager: CanvasTabManager
    current_canvas: object | None = None
    current_file_path_setter: Callable[[Optional[str]], None] | None = None

    def __post_init__(self) -> None:
        if not isinstance(self.tab_manager, CanvasTabManager):
            raise TypeError("DocumentTabsContext.tab_manager must be a CanvasTabManager")
        if (
            self.current_file_path_setter is not None
            and not callable(self.current_file_path_setter)
        ):
            raise TypeError(
                "DocumentTabsContext.current_file_path_setter must be callable or None"
            )


class DocumentController:
    """Operaciones de documento separadas de ChemusonWindow para test no-GUI."""

    @classmethod
    def load_recent_files(cls, context: RecentFilesContext) -> list[str]:
        value = context.settings.value("recent_files", [])
        if isinstance(value, str):
            return [value]
        if isinstance(value, list):
            return [str(p) for p in value if p]
        return []

    @classmethod
    def save_recent_files(cls, context: RecentFilesContext) -> None:
        context.settings.setValue("recent_files", context.recent_files)

    @staticmethod
    def set_canvas_file_path(context: DocumentTabsContext, canvas, filepath: Optional[str]) -> None:
        clean_path = os.path.abspath(filepath) if filepath else None
        context.tab_manager.set_canvas_file_path(canvas, clean_path)
        if canvas is context.current_canvas and context.current_file_path_setter is not None:
            context.current_file_path_setter(clean_path)

    @staticmethod
    def update_tab_title(context: DocumentTabsContext, canvas) -> None:
        context.tab_manager.update_tab_title(canvas)

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
        canvas,
    ) -> bool:
        from PyQt6.QtWidgets import QMessageBox

        if canvas is None:
            raise TypeError("confirm_discard_changes requires an explicit canvas")
        if canvas.undo_stack.isClean():
            return True
        tabs = context.tab_manager.tabs
        index = tabs.indexOf(canvas)
        title = context.tab_manager.tab_titles.get(canvas, "Documento")
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
                context.activate_canvas(canvas)
            context.save_canvas()
            saved = canvas.undo_stack.isClean()
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
