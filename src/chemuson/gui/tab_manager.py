from __future__ import annotations

import os
from typing import Callable, Iterable, Optional

from PyQt6.QtCore import QObject, QTimer
from PyQt6.QtWidgets import QTabWidget

from chemuson.gui.canvas import ChemusonCanvas
from chemuson.chemio.persistence import PersistenceManager
from chemuson.utils.autosave import (
    AutosaveController,
    AutosaveManager,
    AutosaveSerializer,
    AutosaveUndoStack,
)


AutosaveFactory = Callable[
    [QObject, ChemusonCanvas, AutosaveUndoStack, AutosaveSerializer[ChemusonCanvas]],
    AutosaveController,
]


class _QtAutosaveManager(QObject):
    """Restore Qt ownership around the dependency-free autosave core."""

    def __init__(
        self,
        parent: QObject,
        document: ChemusonCanvas,
        undo_stack: AutosaveUndoStack,
        serializer: AutosaveSerializer[ChemusonCanvas],
        *,
        backup_limit: int = AutosaveManager.MAX_BACKUPS,
    ) -> None:
        super().__init__(parent)

        def create_timer(
            interval_ms: int,
            single_shot: bool,
            callback: Callable[[], None],
        ) -> QTimer:
            timer = QTimer(self)
            timer.setInterval(interval_ms)
            timer.setSingleShot(single_shot)
            timer.timeout.connect(callback)
            return timer

        self._core = AutosaveManager[ChemusonCanvas](
            document,
            undo_stack,
            serializer,
            create_timer,
            backup_limit=backup_limit,
        )

    @property
    def core(self) -> AutosaveManager[ChemusonCanvas]:
        """Return the dependency-free core for focused composition tests."""
        return self._core

    def start(self) -> None:
        self._core.start()

    def stop(self) -> None:
        self._core.stop()

    def set_original_path(self, filepath: Optional[str]) -> None:
        self._core.set_original_path(filepath)

    def restart_debounce(self) -> None:
        self._core.restart_debounce()

    def cancel_debounce(self) -> None:
        self._core.cancel_debounce()

    def cleanup_after_save(self) -> None:
        self._core.cleanup_after_save()


class CanvasTabManager:
    """Gestiona pestañas de documento, metadatos y autosave asociados."""

    def __init__(
        self,
        tabs: QTabWidget,
        *,
        autosave_parent: QObject,
        canvas_factory: Callable[[], ChemusonCanvas] = ChemusonCanvas,
        autosave_factory: AutosaveFactory = _QtAutosaveManager,
    ) -> None:
        self.tabs = tabs
        self._autosave_parent = autosave_parent
        self._canvas_factory = canvas_factory
        self._autosave_factory = autosave_factory
        self.file_paths: dict[ChemusonCanvas, Optional[str]] = {}
        self.tab_titles: dict[ChemusonCanvas, str] = {}
        self.autosave_managers: dict[ChemusonCanvas, AutosaveController] = {}
        self._untitled_counter = 1

    def tabs_alive(self) -> bool:
        """Indica si el widget de pestañas sigue disponible."""
        try:
            self.tabs.count()
        except RuntimeError:
            return False
        return True

    def canvas_from_tab_index(self, index: int) -> Optional[ChemusonCanvas]:
        """Obtiene el canvas asociado al índice de pestaña."""
        if not self.tabs_alive() or index < 0:
            return None
        widget = self.tabs.widget(index)
        return widget if isinstance(widget, ChemusonCanvas) else None

    def current_canvas(self) -> Optional[ChemusonCanvas]:
        """Devuelve el canvas activo actualmente en el tab widget."""
        return self.canvas_from_tab_index(self.tabs.currentIndex())

    def iter_canvases(self) -> Iterable[ChemusonCanvas]:
        """Itera los canvases actualmente montados en pestañas."""
        if not self.tabs_alive():
            return []
        canvases: list[ChemusonCanvas] = []
        for index in range(self.tabs.count()):
            canvas = self.canvas_from_tab_index(index)
            if canvas is not None:
                canvases.append(canvas)
        return canvases

    def file_path_for(self, canvas: ChemusonCanvas) -> Optional[str]:
        """Ruta asociada a un canvas, si existe."""
        return self.file_paths.get(canvas)

    def autosave_manager_for(self, canvas: ChemusonCanvas) -> Optional[AutosaveController]:
        """Autosave manager asociado a un canvas, si existe."""
        return self.autosave_managers.get(canvas)

    def create_document_tab(
        self,
        *,
        make_current: bool = False,
        prepare_canvas: Optional[Callable[[ChemusonCanvas], None]] = None,
        clean_state_changed: Optional[Callable[[ChemusonCanvas, Optional[bool]], None]] = None,
    ) -> ChemusonCanvas:
        """Crea una nueva pestaña con su canvas y autosave listos."""
        canvas = self._canvas_factory()
        tab_title = (
            "Sin título"
            if self._untitled_counter == 1
            else f"Sin título {self._untitled_counter}"
        )
        self._untitled_counter += 1
        self.file_paths[canvas] = None
        self.tab_titles[canvas] = tab_title
        if prepare_canvas is not None:
            prepare_canvas(canvas)
        index = self.tabs.addTab(canvas, tab_title)
        self.tabs.setTabToolTip(index, "Documento sin guardar")
        autosave_manager = self._autosave_factory(
            self._autosave_parent,
            canvas,
            canvas.undo_stack,
            PersistenceManager.save_to_dict,
        )
        autosave_manager.start()
        self.autosave_managers[canvas] = autosave_manager
        callback = clean_state_changed or self.on_canvas_clean_state_changed
        canvas.undo_stack.cleanChanged.connect(
            lambda clean, c=canvas: callback(c, bool(clean))
        )
        self.update_tab_title(canvas)
        if make_current:
            self.tabs.setCurrentIndex(index)
        return canvas

    def on_canvas_clean_state_changed(
        self,
        canvas: ChemusonCanvas,
        clean: Optional[bool] = None,
    ) -> None:
        """Sincroniza autosave y título cuando cambia el estado clean/dirty."""
        autosave_manager = self.autosave_managers.get(canvas)
        if autosave_manager is not None:
            is_clean = canvas.undo_stack.isClean() if clean is None else bool(clean)
            if is_clean:
                autosave_manager.cancel_debounce()
            else:
                autosave_manager.restart_debounce()
        self.update_tab_title(canvas)

    def set_canvas_file_path(
        self,
        canvas: ChemusonCanvas,
        filepath: Optional[str],
    ) -> Optional[str]:
        """Asigna ruta a un canvas y refresca su presentación."""
        clean_path = os.path.abspath(filepath) if filepath else None
        self.file_paths[canvas] = clean_path
        autosave_manager = self.autosave_managers.get(canvas)
        if autosave_manager is not None:
            autosave_manager.set_original_path(clean_path)
        if clean_path:
            self.tab_titles[canvas] = os.path.basename(clean_path)
        self.update_tab_title(canvas)
        return clean_path

    def update_tab_title(self, canvas: ChemusonCanvas) -> None:
        """Actualiza título y tooltip del canvas dentro del tab widget."""
        if not self.tabs_alive():
            return
        try:
            index = self.tabs.indexOf(canvas)
            if index < 0:
                return
            base_title = self.tab_titles.get(canvas, "Sin título")
            dirty_suffix = " *" if not canvas.undo_stack.isClean() else ""
            self.tabs.setTabText(index, f"{base_title}{dirty_suffix}")
            path = self.file_paths.get(canvas)
            self.tabs.setTabToolTip(index, path or "Documento sin guardar")
        except RuntimeError:
            return

    def discard_canvas(self, canvas: ChemusonCanvas) -> bool:
        """Elimina una pestaña y limpia su estado asociado sin pedir confirmación."""
        index = self.tabs.indexOf(canvas)
        if index < 0:
            return False
        self.tabs.removeTab(index)
        autosave_manager = self.autosave_managers.pop(canvas, None)
        if autosave_manager is not None:
            autosave_manager.stop()
        self.file_paths.pop(canvas, None)
        self.tab_titles.pop(canvas, None)
        canvas.deleteLater()
        return True

    def close_canvas_tab(
        self,
        canvas: ChemusonCanvas,
        *,
        confirm_discard: Callable[[ChemusonCanvas], bool],
        before_discard: Optional[Callable[[ChemusonCanvas], None]] = None,
        create_replacement: Optional[Callable[[], ChemusonCanvas]] = None,
        activate_canvas: Optional[Callable[[ChemusonCanvas], None]] = None,
    ) -> bool:
        """Cierra una pestaña, opcionalmente creando un reemplazo si queda vacía."""
        if not confirm_discard(canvas):
            return False
        was_active = canvas is self.current_canvas()
        if before_discard is not None:
            before_discard(canvas)
        if not self.discard_canvas(canvas):
            return False
        if self.tabs.count() == 0 and create_replacement is not None:
            replacement = create_replacement()
            if activate_canvas is not None:
                activate_canvas(replacement)
            return True
        if was_active and activate_canvas is not None:
            current = self.current_canvas()
            if current is not None:
                activate_canvas(current)
        return True
