"""Pruebas de texto para mensajes del updater."""

from __future__ import annotations

import os
import sys

sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), "..", "src")))

from chemuson.gui.dialogs import format_update_behavior_summary
from chemuson.gui.main_window import (
    FLATPAK_APP_ID,
    format_no_update_message,
    format_update_disabled_message,
)


def test_format_no_update_message_explains_published_releases() -> None:
    message = format_no_update_message(
        current_version="0.2.1",
        channel="beta",
        reason="La version remota no es elegible para el canal o no es mas nueva.",
        latest_version="0.2.2-beta.2",
        source="remote",
    )

    assert "version publicada mas nueva" in message.lower()
    assert "releases publicadas" in message.lower()
    assert "version instalada: 0.2.1" in message.lower()
    assert "canal: beta" in message.lower()
    assert "ultima version consultada: 0.2.2-beta.2" in message.lower()
    assert "origen de datos: github" in message.lower()
    assert "detalle:" in message.lower()


def test_format_no_update_message_warns_when_using_cache() -> None:
    message = format_no_update_message(
        current_version="0.2.2-beta.1",
        channel="beta",
        reason="La version remota no es elegible para el canal o no es mas nueva.",
        latest_version="0.2.2-beta.1",
        source="cache",
    )

    assert "origen de datos: caché local" in message.lower()
    assert "aviso:" in message.lower()
    assert "caché local" in message.lower()


def test_format_update_disabled_message_for_flatpak() -> None:
    message = format_update_disabled_message(flatpak=True, app_id=FLATPAK_APP_ID)

    assert "flatpak" in message.lower()
    assert "no usa auto-update real" in message.lower()
    assert f"flatpak update {FLATPAK_APP_ID}".lower() in message.lower()


def test_format_update_disabled_message_generic() -> None:
    message = format_update_disabled_message(flatpak=False)

    assert "deshabilitado por entorno" in message.lower()


def test_format_update_behavior_summary_clarifies_delivery_modes() -> None:
    message = format_update_behavior_summary(FLATPAK_APP_ID)

    assert "auto-update real" in message.lower()
    assert "appimage" in message.lower()
    assert "instalación silenciosa al cerrar" in message.lower()
    assert "windows instalado con setup" in message.lower()
    assert f"flatpak update {FLATPAK_APP_ID}".lower() in message.lower()
