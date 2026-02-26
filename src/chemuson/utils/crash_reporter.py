"""Registro y notificación de excepciones no controladas en Chemuson."""

from __future__ import annotations

import os
import platform
import sys
import traceback
from datetime import datetime
from importlib import metadata
from typing import Optional

from PyQt6.QtWidgets import QApplication, QMessageBox

_INSTALLED = False


def _base_dir() -> str:
    """Resuelve directorio base de Chemuson para archivos locales."""
    xdg_home = os.environ.get("XDG_CONFIG_HOME")
    if xdg_home:
        return os.path.join(xdg_home, "chemuson")
    return os.path.expanduser("~/.chemuson")


def _chemuson_version() -> str:
    """Obtiene la versión instalada; usa fallback si no está disponible."""
    try:
        return metadata.version("chemuson")
    except Exception:
        return "desconocida"


def write_crash_log(exc: BaseException) -> str:
    """Escribe un reporte de excepción en ~/.chemuson/crash_logs/."""
    logs_dir = os.path.join(_base_dir(), "crash_logs")
    os.makedirs(logs_dir, exist_ok=True)

    now = datetime.now()
    filepath = os.path.join(logs_dir, f"crash_{now.strftime('%Y%m%d_%H%M%S')}.txt")

    trace_text = "".join(traceback.format_exception(type(exc), exc, exc.__traceback__))
    content = [
        "Chemuson Crash Report",
        "=" * 80,
        f"Timestamp: {now.isoformat()}",
        f"Chemuson version: {_chemuson_version()}",
        f"Python: {platform.python_version()}",
        f"Platform: {platform.platform()}",
        "",
        "Traceback:",
        trace_text,
    ]

    with open(filepath, "w", encoding="utf-8") as f:
        f.write("\n".join(content))

    return filepath


def install() -> None:
    """Instala un excepthook global para registrar crashes y notificar al usuario."""
    global _INSTALLED
    if _INSTALLED:
        return

    previous_hook = sys.excepthook

    def _hook(exc_type, exc_value, exc_tb) -> None:
        if issubclass(exc_type, KeyboardInterrupt):
            previous_hook(exc_type, exc_value, exc_tb)
            return

        exception = exc_value if isinstance(exc_value, BaseException) else Exception(str(exc_value))
        log_path: Optional[str] = None
        try:
            log_path = write_crash_log(exception)
        except Exception:
            pass

        message = "Ocurrió un error inesperado."
        if log_path:
            message += f"\nSe guardó un reporte en:\n{log_path}"

        app = QApplication.instance()
        if app is not None:
            QMessageBox.critical(None, "Chemuson - Error crítico", message)
        else:
            sys.stderr.write(message + "\n")

        # Conserva el comportamiento estándar de impresión en stderr.
        if previous_hook is not None:
            previous_hook(exc_type, exc_value, exc_tb)

    sys.excepthook = _hook
    _INSTALLED = True
