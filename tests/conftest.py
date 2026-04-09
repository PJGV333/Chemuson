"""Configuración compartida para pruebas Qt."""

from __future__ import annotations

import gc
import os

import pytest
from PyQt6.QtWidgets import QApplication


os.environ.setdefault("QT_QPA_PLATFORM", "offscreen")


@pytest.fixture(scope="session", autouse=True)
def _session_qapp():
    """Mantiene un único QApplication vivo durante toda la sesión."""
    app = QApplication.instance() or QApplication([])
    app.setQuitOnLastWindowClosed(False)
    yield app
    QApplication.processEvents()


@pytest.fixture(autouse=True)
def _cleanup_qt_widgets(_session_qapp):
    """Cierra widgets residuales entre pruebas para evitar abortos en GC/Qt."""
    yield
    for widget in list(QApplication.topLevelWidgets()):
        widget.close()
        widget.deleteLater()
    QApplication.processEvents()
    gc.collect()
    QApplication.processEvents()
