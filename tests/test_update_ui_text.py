"""Pruebas de texto para mensajes del updater."""

from __future__ import annotations

import os
import sys

sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), "..", "src")))

from chemuson.gui.main_window import format_no_update_message


def test_format_no_update_message_explains_published_releases() -> None:
    message = format_no_update_message(
        current_version="0.2.1",
        channel="beta",
        reason="La version remota no es elegible para el canal o no es mas nueva.",
    )

    assert "version publicada mas nueva" in message.lower()
    assert "releases publicadas" in message.lower()
    assert "version instalada: 0.2.1" in message.lower()
    assert "canal: beta" in message.lower()
    assert "detalle:" in message.lower()
