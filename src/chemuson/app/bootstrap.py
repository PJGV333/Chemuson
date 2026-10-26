"""Composition root de la aplicación gráfica Chemuson."""

from __future__ import annotations

import sys
from pathlib import Path

from PyQt6.QtGui import QIcon
from PyQt6.QtWidgets import QApplication, QMessageBox

from chemuson.gui.main_window import ChemusonWindow
from chemuson.resilience import crash_reporter
from chemuson.version import get_app_version


def _application_icon_path() -> Path:
    """Return the project logo bundled beside the ChemUSON package."""
    return (
        Path(__file__).resolve().parents[1]
        / "assets"
        / "branding"
        / "chemuson-icon-512.png"
    )


def run_app() -> None:
    """Punto de entrada de la GUI de Chemuson."""
    crash_reporter.install()
    try:
        app = QApplication(sys.argv)
        app.setApplicationName("Chemuson")
        app.setApplicationDisplayName("ChemUSON")
        app.setApplicationVersion(get_app_version())
        icon_path = _application_icon_path()
        if not icon_path.is_file():
            raise FileNotFoundError(f"ChemUSON application icon is missing: {icon_path}")
        app_icon = QIcon(str(icon_path))
        if app_icon.isNull():
            raise RuntimeError(f"ChemUSON application icon could not be loaded: {icon_path}")
        app.setWindowIcon(app_icon)
        window = ChemusonWindow()
        window.setWindowIcon(app_icon)
        ChemusonWindow.check_autosaves(window)
        window.show()
        exit_code = app.exec()
    except Exception as exc:
        log_path = crash_reporter.write_crash_log(exc)
        if QApplication.instance() is not None:
            QMessageBox.critical(
                None,
                "ChemUSON - Error crítico",
                "No se pudo iniciar la aplicación.\n"
                f"Se guardó un reporte en:\n{log_path}",
            )
        else:
            sys.stderr.write(
                "No se pudo iniciar la aplicación.\n"
                f"Se guardó un reporte en: {log_path}\n"
            )
        return
    sys.exit(exit_code)
