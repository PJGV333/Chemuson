"""Harness de regresión visual: producción vs. spike PyQt6 aprobado.

El spike (``docs/ui-modernization/pyqt6-spike/``) es el **contrato visual**.
Este script produce, en offscreen, con el mismo ``.venv``:

1. **Chequeos numéricos** a 1440×900 (appbar 54, rail 58, status 34,
   menubar no visible, text toolbar no visible por defecto, rail = QWidget
   no-QToolBar, botones 42 px, icono de hamburguesa) + fit a 980×600.
2. **Capturas de producción** (tema claro y oscuro): ventana completa
   1440×900, 980×600, close-ups de rail y app bar, flyout de enlaces.
3. **Capturas del spike** (claro y oscuro) 1440×900, tal cual.
4. **Contact sheets** (montaje PIL) y ``checks.json`` con el resultado
   numérico de cada chequeo (para la matriz de aceptación del README).

Salida: ``docs/ui-modernization/visual-convergence/captures/``.
"""
from __future__ import annotations

import importlib.util
import json
import os
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
REPO = HERE.parents[2]
OUT = HERE / "captures"
OUT.mkdir(parents=True, exist_ok=True)

os.environ.setdefault("QT_QPA_PLATFORM", "offscreen")
sys.path.insert(0, str(REPO / "src"))
sys.path.insert(0, str(HERE.parent / "pyqt6-spike"))

from PyQt6.QtWidgets import QApplication  # noqa: E402

app = QApplication(sys.argv)

# ---------------------------------------------------------------------------
# Carga de la ventana del spike (script standalone)
# ---------------------------------------------------------------------------
_spike_path = REPO / "docs" / "ui-modernization" / "pyqt6-spike" / "app.py"
_spec = importlib.util.spec_from_file_location("spike_app", _spike_path)
spike_app = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(spike_app)


def _grab(widget, name: str) -> None:
    png = widget.grab()
    path = OUT / name
    png.save(str(path))
    print(f"captured: {name} ({png.width()}x{png.height()})")


def _set_theme_production(win, theme: str) -> None:
    win.current_theme = theme
    win._apply_theme()
    app.processEvents()
    app.processEvents()


def _close_all_flyouts(win) -> None:
    for flyout in win.tool_rail._flyouts.values():
        if flyout is not None and flyout.isVisible():
            flyout.close_with(None)
    app.processEvents()


# ---------------------------------------------------------------------------
# 1. Chequeos numéricos (producción, 1440×900)
# ---------------------------------------------------------------------------
checks: list[dict] = []


def check(name: str, expected, actual) -> bool:
    ok = expected == actual
    checks.append({"name": name, "expected": expected, "actual": actual, "ok": bool(ok)})
    print(f"{'OK  ' if ok else 'FAIL'} {name}: {actual!r} (esperado {expected!r})")
    return ok


def _measure_production(win, size_w: int, size_h: int) -> dict:
    win.resize(size_w, size_h)
    app.processEvents()
    app.processEvents()
    rail = win.tool_rail
    return {
        "appbar_h": win.app_bar.height(),
        "rail_w": rail.width(),
        "rail_btn_w": rail._buttons["select"].width(),
        "status_h": win.statusBar().height(),
        "menubar_visible": win.menuBar().isVisible(),
        "text_toolbar_visible": win.text_toolbar.isVisible(),
        "n_rail_buttons": len(rail._buttons),
        "window_w": win.width(),
        "window_h": win.height(),
        "hamburger_icon": not win.app_bar.menu_button.icon().isNull(),
        "min_w": win.minimumWidth(),
        "min_h": win.minimumHeight(),
    }


from chemuson.gui.main_window import ChemusonWindow  # noqa: E402

win = ChemusonWindow()
win.resize(1440, 900)
win.show()
app.processEvents()
app.processEvents()
m = _measure_production(win, 1440, 900)

check("appbar_h", 54, m["appbar_h"])
check("rail_w", 58, m["rail_w"])
check("rail_btn_w", 42, m["rail_btn_w"])
check("status_h", 34, m["status_h"])
check("menubar_visible", False, m["menubar_visible"])
check("text_toolbar_visible", False, m["text_toolbar_visible"])
check("n_rail_buttons", 15, m["n_rail_buttons"])
check("hamburger_icon", True, m["hamburger_icon"])
check("min_size", (900, 560), (m["min_w"], m["min_h"]))
from PyQt6.QtWidgets import QToolBar  # noqa: E402

check("rail_is_qtoolbar", False, isinstance(win.tool_rail, QToolBar))
visible_tb = [t.objectName() for t in win.findChildren(QToolBar) if t.isVisible()]
check("visible_qtoolbars", [], visible_tb)

# Fit 980×600
m600 = _measure_production(win, 980, 600)
check("fits_980x600", (True, True), (m600["window_w"] <= 980, m600["window_h"] <= 600))
check("rail_btn_42_at_980x600", 42, m600["rail_btn_w"])

# ---------------------------------------------------------------------------
# 2. Capturas de producción (claro)
# ---------------------------------------------------------------------------
_set_theme_production(win, "light")
win.resize(1440, 900)
app.processEvents()
app.processEvents()
_grab(win, "prod-full-light.png")

# Rail + app bar (close-ups)
rail = win.tool_rail
_rl = rail.mapToGlobal(rail.rect().topLeft())
_rb = rail.mapToGlobal(rail.rect().bottomRight())
win.grab().save("/tmp/_win_tmp.png")
from PIL import Image  # noqa: E402

_im = Image.open("/tmp/_win_tmp.png")
_im.crop((_rl.x(), _rl.y(), _rb.x() + 1, _rb.y() + 1)).save(str(OUT / "prod-rail-light.png"))
ab = win.app_bar
_al = ab.mapToGlobal(ab.rect().topLeft())
_ab = ab.mapToGlobal(ab.rect().bottomRight())
_im.crop((_al.x(), _al.y(), _ab.x() + 1, _ab.y() + 1)).save(str(OUT / "prod-appbar-light.png"))

# Flyout de enlaces (primer + segundo clic sobre la categoría activa)
win.toolbar.tool_changed.emit("tool_select")
app.processEvents()
win.tool_rail._buttons["bond"].clicked.emit()
app.processEvents()
win.tool_rail._buttons["bond"].clicked.emit()
app.processEvents()
fly = win.tool_rail._flyouts["bond"]
assert fly is not None and fly.isVisible()
_grab(fly, "prod-bond-flyout-light.png")
fly.close_with(None)
app.processEvents()

# 980×600
win.resize(980, 600)
app.processEvents()
app.processEvents()
_grab(win, "prod-980x600-light.png")
_rail = win.tool_rail
_rl = _rail.mapToGlobal(_rail.rect().topLeft())
_rb = _rail.mapToGlobal(_rail.rect().bottomRight())
win.grab().save("/tmp/_win_tmp.png")
_im = Image.open("/tmp/_win_tmp.png")
_im.crop((_rl.x(), _rl.y(), _rb.x() + 1, _rb.y() + 1)).save(str(OUT / "prod-rail-980x600-light.png"))

# Producción (oscuro)
_close_all_flyouts(win)
_set_theme_production(win, "dark")
win.resize(1440, 900)
app.processEvents()
app.processEvents()
_grab(win, "prod-full-dark.png")
_rail = win.tool_rail
_rl = _rail.mapToGlobal(_rail.rect().topLeft())
_rb = _rail.mapToGlobal(_rail.rect().bottomRight())
win.grab().save("/tmp/_win_tmp.png")
_im = Image.open("/tmp/_win_tmp.png")
_im.crop((_rl.x(), _rl.y(), _rb.x() + 1, _rb.y() + 1)).save(str(OUT / "prod-rail-dark.png"))

# ---------------------------------------------------------------------------
# 3. Capturas del spike (claro / oscuro, 1440×900)
# ---------------------------------------------------------------------------
spike = spike_app.SpikeWindow(app)
spike.set_theme("light")
spike.resize(1440, 900)
spike.show()
app.processEvents()
app.processEvents()
_grab(spike, "spike-full-light.png")
spike.set_theme("dark")
app.processEvents()
app.processEvents()
_grab(spike, "spike-full-dark.png")
spike.close()

# ---------------------------------------------------------------------------
# 4. Contact sheets
# ---------------------------------------------------------------------------


def _contact(names: list[str], out_name: str, thumb: int = 900) -> None:
    imgs: list[Image.Image] = []
    for name in names:
        im = Image.open(OUT / name).convert("RGB")
        if im.size[0] > thumb:
            im = im.resize((thumb, int(im.size[1] * thumb / im.size[0])))
        imgs.append(im)
    # rejilla 2x2
    hmax = max(im.size[1] for im in imgs)
    grid = Image.new("RGB", (thumb * 2 + 10, hmax * 2 + 30), (24, 26, 32))
    from PIL import ImageDraw  # noqa: E402

    draw = ImageDraw.Draw(grid)
    for i, (im, name) in enumerate(zip(imgs, names)):
        r, c = divmod(i, 2)
        x, y = c * (thumb + 10), r * (hmax + 30)
        grid.paste(im, (x, y + 20))
        draw.text((x + 4, y + 4), name, fill=(226, 232, 240))
    grid.save(OUT / out_name)
    print(f"contact: {out_name}")


_contact(
    ["prod-full-light.png", "prod-full-dark.png", "spike-full-light.png", "spike-full-dark.png"],
    "contact-window.png",
)
_contact(
    ["prod-rail-light.png", "prod-rail-dark.png"],
    "contact-rail.png",
)
_contact(
    ["prod-appbar-light.png", "prod-bond-flyout-light.png"],
    "contact-appbar-flyout.png",
)
_contact(
    ["prod-980x600-light.png", "prod-rail-980x600-light.png"],
    "contact-980x600.png",
)

# ---------------------------------------------------------------------------
# 5. checks.json + resumen
# ---------------------------------------------------------------------------
with open(OUT / "checks.json", "w") as fh:
    json.dump(checks, fh, indent=2)

n_ok = sum(1 for c in checks if c["ok"])
print(f"\nCHECKS: {n_ok}/{len(checks)} passed")
if n_ok != len(checks):
    for c in checks:
        if not c["ok"]:
            print("  FAIL:", c["name"], c["actual"], "!= expected", c["expected"])
# Evitamos el teardown de Qt (segfault conocido en offscreen al cerrar).
os._exit(0 if n_ok == len(checks) else 1)
