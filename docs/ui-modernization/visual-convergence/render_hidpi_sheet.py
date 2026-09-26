"""Renderiza una hoja local de iconos a DPR 1/1.25/1.5/2 sin ventana grande."""
from __future__ import annotations

import os
from pathlib import Path

os.environ.setdefault("QT_QPA_PLATFORM", "wayland")

from PyQt6.QtGui import QPainter
from PyQt6.QtWidgets import QApplication
from PIL import Image, ImageDraw

from chemuson.gui.theme.icon_provider import IconProvider

ROOT = Path(__file__).resolve().parent
OUT = ROOT / "real-kde" / "hidpi-icon-sheet.png"
OUT.parent.mkdir(parents=True, exist_ok=True)
app = QApplication([])

# Los estáticos de la barra/rail y builders reales para ring/atom.
icons = [
    ("pointer", "pointer"),
    ("bond", "bond-single"),
    ("ring", "dynamic:ring"),
    ("atom", "dynamic:atom"),
    ("flask", "flask"),
    ("plus", "plus"),
    ("undo", "undo"),
    ("redo", "redo"),
    ("moon", "moon"),
    ("sun", "sun"),
    ("settings", "sliders"),
]
sizes = (16, 18, 20, 21, 24, 32)
dprs = (1.0, 1.25, 1.5, 2.0)
cell_w, cell_h, pad = 112, 78, 12
width = pad + len(sizes) * cell_w
height = pad + len(icons) * len(dprs) * cell_h + 30
sheet = Image.new("RGB", (width, height), (245, 247, 250))
draw = ImageDraw.Draw(sheet)
draw.text((pad, 5), "ChemUSON SVG HiDPI — real Qt backend", fill=(15, 23, 42))

for row, (label, icon_name) in enumerate(icons):
    for di, dpr in enumerate(dprs):
        y = pad + (row * len(dprs) + di) * cell_h + 24
        draw.text((pad, y + 5), f"{label} @ {dpr:g}", fill=(30, 41, 59))
        provider = IconProvider(dpr=dpr)
        for si, size in enumerate(sizes):
            x = pad + 112 + si * cell_w
            if icon_name.startswith("dynamic:"):
                key = icon_name.split(":", 1)[1]
                if key == "ring":
                    pm = provider.pixmap_dynamic("ring", "#334155", size, variant="benzene")
                else:
                    pm = provider.pixmap_dynamic("atom", "#334155", size, element="N")
            else:
                pm = provider.pixmap(icon_name, "#334155", size)
            if not pm.isNull():
                q = QPainter(pm)
                q.end()
                data = pm.toImage().convertToFormat(pm.toImage().Format.Format_RGBA8888)
                ptr = data.bits(); ptr.setsize(data.sizeInBytes())
                im = Image.frombytes("RGBA", (data.width(), data.height()), bytes(ptr))
                # Aplanar sobre blanco conservando DPR visual/lógico.
                logical_w = max(1, round(pm.width() / pm.devicePixelRatioF()))
                logical_h = max(1, round(pm.height() / pm.devicePixelRatioF()))
                im = im.resize((logical_w, logical_h), Image.Resampling.LANCZOS)
                bg = Image.new("RGBA", im.size, (255, 255, 255, 255)); bg.alpha_composite(im)
                sheet.paste(bg.convert("RGB"), (x + (cell_w - logical_w)//2, y + 12))
            draw.text((x + 4, y + 48), str(size), fill=(71, 85, 105))

sheet.save(OUT)
print(OUT)
