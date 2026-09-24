"""Capturas de referencia de la Fase 2 (iconos SVG): light/dark, ventana real.

Uso (desde la raíz del repo):

    QT_QPA_PLATFORM=offscreen PYTHONPATH=src \
        python docs/ui-modernization/icon-phase-shots/make_shots.py <outdir>

Crea ``ChemusonWindow``, aplica ``light``→``dark`` vía ``_apply_theme()`` y
captura la ventana completa y cada toolbar con ``grab()`` (sin tocar
producción; solo lectura del estado de la UI).
"""
import os
import sys
from pathlib import Path

os.environ["QT_QPA_PLATFORM"] = "offscreen"
sys.path.insert(0, str(Path(__file__).resolve().parents[3] / "src"))

from PyQt6.QtWidgets import QApplication
from PyQt6.QtCore import QEventLoop, QTimer

app = QApplication([])
from chemuson.gui.main_window import ChemusonWindow

win = ChemusonWindow()
win.resize(1280, 860)
win.inspector_dock.show()
win.show()
outdir = Path(sys.argv[1])
outdir.mkdir(parents=True, exist_ok=True)


def capture(theme_name: str) -> None:
    loop = QEventLoop()

    def _snap() -> None:
        win.grab().save(str(outdir / f"icon-phase-{theme_name}.png"))
        for attr, fname in (
            ("main_toolbar", "maintoolbar"),
            ("toolbar", "toolbar-draw"),
            ("symbols_toolbar", "toolbar-symbols"),
            ("text_toolbar", "toolbar-text"),
        ):
            widget = getattr(win, attr, None)
            if widget is not None:
                widget.grab().save(str(outdir / f"{fname}-{theme_name}.png"))
        loop.quit()

    win.current_theme = theme_name
    win._apply_theme()
    app.processEvents()
    QTimer.singleShot(500, _snap)
    loop.exec()
    print(f"capture {theme_name} ok")


capture("light")
capture("dark")
win.close()
print("SHOTS: OK")
