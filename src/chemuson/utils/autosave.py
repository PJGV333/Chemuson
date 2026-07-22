"""Gestión de autosave y backups rotativos para Chemuson."""

from __future__ import annotations

import hashlib
import json
import os
import uuid
from datetime import datetime, timezone
from typing import Callable, Optional, Protocol


class AutosaveDocument(Protocol):
    """Opaque document accepted by an injected autosave serializer."""


class AutosaveUndoStack(Protocol):
    """Minimal undo state needed to decide whether a snapshot is due."""

    def isClean(self) -> bool: ...


class AutosaveTimer(Protocol):
    """Minimal timer lifecycle owned by the composition root."""

    def start(self) -> None: ...

    def stop(self) -> None: ...


AutosaveSerializer = Callable[[AutosaveDocument], dict[str, object]]
AutosaveTimerFactory = Callable[[int, bool, Callable[[], None]], AutosaveTimer]


class AutosaveManager:
    """Encapsula la lógica de autosave de un documento."""

    AUTOSAVE_INTERVAL_MS = 2 * 60 * 1000
    DEBOUNCE_INTERVAL_MS = 8 * 1000
    MAX_BACKUPS = 5

    def __init__(
        self,
        document: AutosaveDocument,
        undo_stack: AutosaveUndoStack,
        serializer: AutosaveSerializer,
        timer_factory: AutosaveTimerFactory,
        *,
        backup_limit: int = MAX_BACKUPS,
    ) -> None:
        self._document = document
        self._undo_stack = undo_stack
        self._serializer = serializer
        self._backup_limit = max(1, int(backup_limit))
        self.autosave_dir = self.default_autosave_dir()
        os.makedirs(self.autosave_dir, exist_ok=True)

        # UUID estable por pestaña para documentos sin ruta asignada.
        self._session_doc_id = str(uuid.uuid4())
        self._doc_id_aliases: set[str] = {self._session_doc_id}
        self._original_path: Optional[str] = None
        self._running = False

        self._periodic_timer = timer_factory(
            self.AUTOSAVE_INTERVAL_MS,
            False,
            self._maybe_autosave,
        )
        self._debounce_timer = timer_factory(
            self.DEBOUNCE_INTERVAL_MS,
            True,
            self._maybe_autosave,
        )

    @classmethod
    def default_base_dir(cls) -> str:
        """Resuelve el directorio base de Chemuson en entorno local."""
        xdg_home = os.environ.get("XDG_CONFIG_HOME")
        if xdg_home:
            return os.path.join(xdg_home, "chemuson")
        return os.path.expanduser("~/.chemuson")

    @classmethod
    def default_autosave_dir(cls) -> str:
        """Devuelve la ruta base de autosaves de Chemuson."""
        return os.path.join(cls.default_base_dir(), "autosave")

    def start(self) -> None:
        """Inicia timers de autosave periódico y debounce."""
        if self._running:
            return
        self._running = True
        self._periodic_timer.start()
        if not self._undo_stack.isClean():
            self._debounce_timer.start()

    def stop(self) -> None:
        """Detiene completamente el mecanismo de autosave."""
        self._running = False
        self._periodic_timer.stop()
        self._debounce_timer.stop()

    def set_original_path(self, filepath: Optional[str]) -> None:
        """Actualiza la ruta de origen del documento para metadatos y hash."""
        clean_path = os.path.abspath(filepath) if filepath else None
        self._original_path = clean_path
        if clean_path:
            self._doc_id_aliases.add(clean_path)

    def restart_debounce(self) -> None:
        """Reinicia el debounce para guardar tras unos segundos sin actividad."""
        if not self._running:
            return
        if self._undo_stack.isClean():
            self._debounce_timer.stop()
            return
        self._debounce_timer.start()

    def cancel_debounce(self) -> None:
        """Cancela un autosave diferido pendiente."""
        self._debounce_timer.stop()

    def cleanup_after_save(self) -> None:
        """Limpia autosaves obsoletos del documento y guarda backup rotativo."""
        for stale_path in self._list_document_autosaves():
            try:
                os.remove(stale_path)
            except OSError:
                continue
        # Conservamos un snapshot post-guardado como respaldo de seguridad.
        self._write_autosave(force=True)

    def _doc_id(self) -> str:
        return self._original_path or self._session_doc_id

    def _doc_hash(self, doc_id: Optional[str] = None) -> str:
        value = (doc_id or self._doc_id()).encode("utf-8", errors="ignore")
        return hashlib.sha256(value).hexdigest()[:16]

    def _autosave_name(self) -> str:
        now = datetime.now(timezone.utc)
        return f"{self._doc_hash()}_{now.strftime('%Y%m%d_%H%M%S_%f')}.json"

    def _list_document_autosaves(self) -> list[str]:
        prefixes = {self._doc_hash(alias) for alias in self._doc_id_aliases}
        if self._original_path:
            prefixes.add(self._doc_hash(self._original_path))
        paths: list[str] = []
        try:
            names = os.listdir(self.autosave_dir)
        except OSError:
            return paths
        for name in names:
            if not name.endswith(".json"):
                continue
            if not any(name.startswith(f"{prefix}_") for prefix in prefixes):
                continue
            full_path = os.path.join(self.autosave_dir, name)
            if os.path.isfile(full_path):
                paths.append(full_path)
        paths.sort(key=lambda p: os.path.getmtime(p), reverse=True)
        return paths

    def _rotate_backups(self) -> None:
        """Elimina backups antiguos cuando supera el límite configurado."""
        stale_files = self._list_document_autosaves()[self._backup_limit :]
        for filepath in stale_files:
            try:
                os.remove(filepath)
            except OSError:
                continue

    def _maybe_autosave(self) -> None:
        """Ejecuta autosave solo cuando hay cambios pendientes."""
        if self._undo_stack.isClean():
            return
        self._write_autosave(force=False)

    def _write_autosave(self, *, force: bool = False) -> Optional[str]:
        """Serializa el estado actual y lo persiste en JSON con metadatos."""
        if not force and self._undo_stack.isClean():
            return None

        os.makedirs(self.autosave_dir, exist_ok=True)
        payload = self._serializer(self._document)
        now = datetime.now(timezone.utc)
        payload["autosave_metadata"] = {
            "original_path": self._original_path,
            "timestamp": now.isoformat(),
            "doc_id": self._doc_id(),
            "doc_hash": self._doc_hash(),
        }

        autosave_path = os.path.join(self.autosave_dir, self._autosave_name())
        with open(autosave_path, "w", encoding="utf-8") as f:
            json.dump(payload, f, indent=2, ensure_ascii=False)

        self._rotate_backups()
        return autosave_path
