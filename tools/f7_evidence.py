"""Genera evidencia visual de la Fase 7 (offscreen).

Genera screenshots del tema claro/oscuro en dos tamaños, el onboarding,
la página de plantillas del panel lateral y la paleta de comandos, además
de un smoke HiDPI (QT_SCALE_FACTOR=2).
"""
from __future__ import annotations

import os
import sys

_REPO_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
OUT = os.path.join(
    _REPO_ROOT,
    "openspec", "changes", "2026-10-01-modernize-ui-polish", "evidence",
)
os.makedirs(OUT, exist_ok=True)

from PyQt6.QtCore import QSettings  # noqa: E402
from PyQt6.QtWidgets import QApplication  # noqa: E402

app = QApplication(sys.argv)

_ONBOARD_KEY = "ui/onboarding/completed"


def _set_onboarding(completed: bool) -> None:
    s = QSettings("Chemuson", "Chemuson")
    s.setValue(_ONBOARD_KEY, "true" if completed else "false")
    s.sync()


def build_window(width: int, height: int, theme: str, onboarding: bool = False):
    from chemuson.gui.main_window import ChemusonWindow

    _set_onboarding(not onboarding)  # las vistas base se muestran limpias
    win = ChemusonWindow()
    win.toggle_theme(theme == "dark")
    win.resize(width, height)
    win.show()
    app.processEvents()
    return win


def save(win, name: str) -> None:
    win.grab().save(os.path.join(OUT, name))
    win.close()
    app.processEvents()
    print(f"OK {name}", flush=True)


# --- Ventana base light / dark en dos tamaños ----------------------------
for theme in ("light", "dark"):
    for w, h in ((1440, 900), (980, 600)):
        win = build_window(w, h, theme)
        save(win, f"main_{theme}_{w}x{h}.png")

# --- Onboarding (primer uso) ---------------------------------------------
# Con el flag en False, ``assemble_application_shell`` muestra el overlay
# de onboarding automáticamente (flujo real de primera ejecución).
win = build_window(1440, 900, "light", onboarding=True)
app.processEvents()
win.grab().save(os.path.join(OUT, "onboarding_step1.png"))
win.close()
app.processEvents()
print("OK onboarding_step1.png", flush=True)

# --- Plantillas (panel lateral) ------------------------------------------
win = build_window(1440, 900, "light")
win.side_panel.show_page("templates")
app.processEvents()
win.grab().save(os.path.join(OUT, "side_panel_templates.png"))
win.close()
app.processEvents()
print("OK side_panel_templates.png", flush=True)

# --- Paleta de comandos (Ctrl+P) ------------------------------------------
win = build_window(1440, 900, "dark")
win._open_command_palette()
app.processEvents()
win.grab().save(os.path.join(OUT, "command_palette.png"))
win.close()
app.processEvents()
print("OK command_palette.png", flush=True)

print("DONE", flush=True)
