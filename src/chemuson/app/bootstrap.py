"""Composition root de la aplicación gráfica Chemuson."""

from __future__ import annotations

import sys

from PyQt6.QtWidgets import QApplication, QMessageBox

from chemuson.gui.main_window import ChemusonWindow
from chemuson.utils import crash_reporter
from chemuson.version import get_app_version


def run_app() -> None:
    """Punto de entrada de la GUI de Chemuson."""
    crash_reporter.install()
    try:
        app = QApplication(sys.argv)
        app.setApplicationName("Chemuson")
        app.setApplicationVersion(get_app_version())
        window = ChemusonWindow()
        ChemusonWindow.check_autosaves(window)
        window.show()
        exit_code = app.exec()
    except Exception as exc:
        log_path = crash_reporter.write_crash_log(exc)
        if QApplication.instance() is not None:
            QMessageBox.critical(
                None,
                "Chemuson - Error crítico",
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
