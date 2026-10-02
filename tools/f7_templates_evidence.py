"""Evidencia visual de los thumbnails de plantillas (Fase 7, HiDPI).

Regenera, en offscreen y con el mismo ``.venv``:

- ``templates_after_100.png`` / ``templates_after_200.png``: captura del
  SidePanel (página *Plantillas*) a DPR 1 y DPR 2.
- ``templates_montage_100.png`` / ``templates_montage_200.png``: montaje
  comparativo antes (``templates_before_*.png``, estado sin ``iconSize``)
  vs. después (la captura regenerada por este script).

El DPR se fuerza con ``QT_SCALE_FACTOR`` en el entorno del hijo; cada escala
se captura en un subprocess porque el ratio de dispositivo se fija al crear
la ``QApplication``. Los ``templates_before_*.png`` NO se regeneran: son el
baseline visual comprometido en el OpenSpec.

Uso: ``python tools/f7_templates_evidence.py``
"""
from __future__ import annotations

import os
import subprocess
import sys

_REPO_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
OUT = os.path.join(
    _REPO_ROOT,
    "openspec",
    "changes",
    "2026-10-01-modernize-ui-polish",
    "evidence",
)
os.makedirs(OUT, exist_ok=True)

_ONBOARD_KEY = "ui/onboarding/completed"
_SCALERS = (1, 2)


def _child_capture() -> None:
    """Captura el SidePanel (página Plantillas) con el DPR del entorno."""
    sys.path.insert(0, os.path.join(_REPO_ROOT, "src"))

    from PyQt6.QtCore import QSettings
    from PyQt6.QtWidgets import QApplication

    app = QApplication(sys.argv)
    scale = int(os.environ.get("QT_SCALE_FACTOR", "1"))

    settings = QSettings("Chemuson", "Chemuson")
    settings.setValue(_ONBOARD_KEY, "true")
    settings.sync()

    from chemuson.gui.main_window import ChemusonWindow

    win = ChemusonWindow()
    # Tema claro: coincide con la evidencia antes/después comprometida, de modo
    # que el montaje comparativo sea comparable píxel a píxel.
    win.toggle_theme(False)
    win.resize(1440, 900)
    win.show()
    win.side_panel.show_page("templates")
    app.processEvents()

    panel = win.side_panel.grab()
    path = os.path.join(OUT, f"templates_after_{scale * 100}.png")
    panel.save(path)
    print(
        f"OK templates_after_{scale * 100}.png "
        f"({panel.width()}x{panel.height()}, DPR {app.primaryScreen().devicePixelRatio()})",
        flush=True,
    )
    win.close()


def _montage(scale: int) -> None:
    """Montaje comparativo antes/después de la misma captura."""
    from PIL import Image, ImageDraw

    before = Image.open(os.path.join(OUT, f"templates_before_{scale * 100}.png")).convert("RGBA")
    after = Image.open(os.path.join(OUT, f"templates_after_{scale * 100}.png")).convert("RGBA")

    margin = 32 * scale
    label_h = 42 * scale
    gap = 32 * scale
    width = margin * 2 + before.width + gap + after.width
    height = margin * 2 + label_h * 2 + max(before.height, after.height)

    sheet = Image.new("RGBA", (width, height), (248, 250, 252, 255))
    draw = ImageDraw.Draw(sheet)
    draw.text((margin, margin), f"SidePanel · Plantillas · {scale * 100} %", fill=(26, 26, 26))

    y = margin + label_h
    draw.text((margin, y), "Antes", fill=(26, 26, 26))
    draw.text((margin + before.width + gap, y), "Después", fill=(26, 26, 26))
    y += label_h
    sheet.paste(before, (margin, y))
    sheet.paste(after, (margin + before.width + gap, y))

    path = os.path.join(OUT, f"templates_montage_{scale * 100}.png")
    sheet.save(path)
    print(f"OK templates_montage_{scale * 100}.png ({width}x{height})", flush=True)


def main() -> None:
    if os.environ.get("_F7_TEMPLATES_EVIDENCE_CHILD") == "1":
        _child_capture()
        return

    env = dict(os.environ, QT_QPA_PLATFORM="offscreen", _F7_TEMPLATES_EVIDENCE_CHILD="1")
    for scale in _SCALERS:
        subprocess.run(
            [sys.executable, os.path.abspath(__file__)],
            env=dict(env, QT_SCALE_FACTOR=str(scale)),
            check=True,
        )
    for scale in _SCALERS:
        _montage(scale)
    print("DONE", flush=True)


if __name__ == "__main__":
    main()
