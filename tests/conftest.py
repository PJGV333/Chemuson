"""Configuración compartida para pruebas Qt."""

from __future__ import annotations

import gc
import os

import pytest


os.environ.setdefault("QT_QPA_PLATFORM", "offscreen")

try:
    from PyQt6.QtWidgets import QApplication
except Exception:  # pragma: no cover - jobs sin Qt (p. ej. windows-smoke)
    QApplication = None


@pytest.fixture(scope="session", autouse=True)
def _session_qapp():
    """Mantiene un único QApplication vivo durante toda la sesión cuando exista Qt."""
    if QApplication is None:
        yield None
        return
    app = QApplication.instance() or QApplication([])
    app.setQuitOnLastWindowClosed(False)
    yield app
    QApplication.processEvents()


@pytest.fixture(autouse=True)
def _cleanup_qt_widgets(_session_qapp):
    """Cierra widgets residuales entre pruebas para evitar abortos en GC/Qt."""
    yield
    if QApplication is None:
        return
    for widget in list(QApplication.topLevelWidgets()):
        widget.close()
        widget.deleteLater()
    QApplication.processEvents()
    gc.collect()
    QApplication.processEvents()


@pytest.fixture
def isolated_shell_settings(tmp_path, monkeypatch):
    """Aísla ``application_settings`` en un INI temporal con defaults fijos.

    Necesario para contratos de *defaults* del shell (barra de menús y barra
    de acciones visibles): sin esto la ventana leería el INI real del equipo
    y el resultado dependería de la configuración local.
    Deja ``ui/onboarding/completed`` en ``True`` para no abrir el overlay.
    """
    from PyQt6.QtCore import QSettings

    settings_path = tmp_path / "Chemuson-shell-settings.ini"

    def settings_factory():
        store = QSettings(str(settings_path), QSettings.Format.IniFormat)
        store.clear()
        store.setValue("ui/onboarding/completed", True)
        store.sync()
        return store

    monkeypatch.setattr(
        "chemuson.platform.settings.application_settings", settings_factory
    )
    monkeypatch.setattr(
        "chemuson.gui.shell.assembly.application_settings", settings_factory
    )
    return settings_factory
