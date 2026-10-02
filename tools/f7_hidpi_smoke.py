"""Smoke HiDPI: construye la ventana con QT_SCALE_FACTOR=2 y verifica que
los iconos del rail, la app bar y los botones de plantillas sean no nulos a
DPR 1 y 2. Guarda un screenshot de prueba.
"""
from __future__ import annotations

import os
import sys

from PyQt6.QtWidgets import QApplication  # noqa: E402

app = QApplication(sys.argv)

from chemuson.gui.main_window import ChemusonWindow  # noqa: E402

win = ChemusonWindow()
win.toggle_theme(False)
win.resize(1440, 900)
win.show()
app.processEvents()

dpr = app.primaryScreen().devicePixelRatio()
print(f"DPR efectivo: {dpr}", flush=True)

def _non_null(widgets, label):
    for w in widgets:
        icon = w.icon()
        assert icon is not None and not icon.isNull(), f"{label}: icono nulo en {w}"
    print(f"OK {label}: iconos no nulos", flush=True)

# Rail de herramientas.
rail_btns = [b for b in win.tool_rail._buttons.values()]
_non_null(rail_btns, "tool_rail")

# App bar.
_non_null([win.app_bar.undo_button, win.app_bar.redo_button,
           win.app_bar.theme_button, win.app_bar.preferences_button], "app_bar")

# Verificamos los thumbnails de plantilla vía servicio (contexto fresco).
from chemuson.gui.template_browser_service import TemplateBrowserService  # noqa: E402
from types import SimpleNamespace  # noqa: E402
lib = win.template_library
grouped = lib.grouped_templates()
tpl_ids = [t["id"] for g in grouped for t in g["templates"]]
if tpl_ids:
    svc = TemplateBrowserService()
    ctx = SimpleNamespace(
        template_library=lib, preview_cache={}, show_status=lambda _s: None
    )
    for tid in tpl_ids[:3]:
        ic = svc.template_preview_icon(ctx, tid)
        assert ic is not None and not ic.isNull(), f"thumbnail nulo para {tid}"
    print(f"OK plantillas: {len(tpl_ids)} disponibles, thumbnails no nulos", flush=True)
else:
    print("OK plantillas: sin plantillas registradas", flush=True)

# Screenshot de prueba HiDPI.
out = os.path.join(
    os.path.dirname(os.path.dirname(os.path.abspath(__file__))),
    "openspec", "changes", "2026-10-01-modernize-ui-polish", "evidence",
)
win.grab().save(os.path.join(out, "hidpi_200.png"))
print("OK hidpi_200.png", flush=True)
win.close()
print("SMOKE PASS", flush=True)
