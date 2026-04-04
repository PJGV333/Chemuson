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


class DocumentController:
    """Operaciones de documento separadas de ChemusonWindow para test no-GUI."""

    @staticmethod
    def _settings_from(target) -> object:
        return getattr(target, "_settings", target)

    @staticmethod
    def _recent_context_from(target) -> RecentFilesContext:
        if isinstance(target, RecentFilesContext):
            return target
        return RecentFilesContext(
            settings=target._settings,
            recent_files=target._recent_files,
            persist_recent_files=target._save_recent_files,
            refresh_recent_menu=target._update_recent_menu,
        )

    @staticmethod
    def _discard_context_from(target) -> DocumentDiscardContext:
        if isinstance(target, DocumentDiscardContext):
            return target
        return DocumentDiscardContext(
            parent=target,
            canvas=target.canvas,
            tab_manager=getattr(target, "_tab_manager", target),
            save_canvas=target._on_file_save,
            activate_canvas=target._set_active_canvas,
        )

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
    def save_recent_files(cls, target) -> None:
        context = cls._recent_context_from(target)
        context.settings.setValue("recent_files", context.recent_files)

    @staticmethod
    def set_canvas_file_path(target, canvas, filepath: Optional[str]) -> None:
        clean_path = os.path.abspath(filepath) if filepath else None
        manager = getattr(target, "_tab_manager", None)
        if manager is not None:
            manager.set_canvas_file_path(canvas, clean_path)
        else:
            target._canvas_file_paths[canvas] = clean_path
            autosave_manager = target._canvas_autosave_managers.get(canvas)
            if autosave_manager is not None:
                autosave_manager.set_original_path(clean_path)
            if clean_path:
                target._canvas_tab_titles[canvas] = os.path.basename(clean_path)
        if canvas is getattr(target, "canvas", None):
            target._current_file_path = clean_path
        updater = getattr(target, "_update_tab_title", None)
        if callable(updater):
            updater(canvas)

    @staticmethod
    def update_tab_title(target, canvas) -> None:
        manager = getattr(target, "_tab_manager", None)
        if manager is not None:
            manager.update_tab_title(canvas)
            return
        if not target._tabs_alive():
            return
        try:
            index = target.tabs.indexOf(canvas)
            if index < 0:
                return
            base_title = target._canvas_tab_titles.get(canvas, "Sin título")
            dirty_suffix = " *" if not canvas.undo_stack.isClean() else ""
            target.tabs.setTabText(index, f"{base_title}{dirty_suffix}")
            path = target._canvas_file_paths.get(canvas)
            target.tabs.setTabToolTip(index, path or "Documento sin guardar")
        except RuntimeError:
            return

    @classmethod
    def add_recent_file(cls, target, filepath: str) -> None:
        if not filepath:
            return
        context = cls._recent_context_from(target)
        path = os.path.abspath(filepath)
        context.recent_files[:] = [p for p in context.recent_files if p != path]
        context.recent_files.insert(0, path)
        del context.recent_files[10:]
        context.persist_recent_files()
        context.refresh_recent_menu()

    @classmethod
    def confirm_discard_changes(cls, target, canvas=None) -> bool:
        from PyQt6.QtWidgets import QMessageBox

        context = cls._discard_context_from(target)
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
