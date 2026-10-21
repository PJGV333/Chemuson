from __future__ import annotations

import argparse
import json
from pathlib import Path

from PyQt6.QtCore import QPoint, Qt
from PyQt6.QtGui import QActionGroup, QRawFont
from PyQt6.QtTest import QTest
from PyQt6.QtWidgets import QApplication, QMainWindow

from chemuson.gui.styles import get_tool_palette_stylesheet
from chemuson.gui.theme import apply_theme
from chemuson.gui.toolbar import SymbolPaletteToolbar


def _submenu(parent_menu, title: str):
    for action in parent_menu.actions():
        if action.text() == title and action.menu() is not None:
            return action.menu()
    raise RuntimeError(f"Submenu not found: {title}")


def _capture_menu(app: QApplication, menu, out_dir: Path, stem: str) -> dict[str, object]:
    menu.ensurePolished()
    menu.adjustSize()
    screen = app.primaryScreen()
    geometry = screen.availableGeometry() if screen is not None else None
    point = QPoint(24, 36)
    if geometry is not None:
        point += geometry.topLeft()
    menu.popup(point)
    QTest.qWait(120)
    app.processEvents()
    menu.setActiveAction(None)
    app.processEvents()

    # QWidget.grab captures the Qt-painted menu content. QScreen.grabWindow
    # is intentionally avoided because Wayland may return a black surface for
    # another top-level popup even when it is visibly rendered.
    pixmap = menu.grab()
    image_path = out_dir / f"{stem}.png"
    if not pixmap.save(str(image_path)):
        raise RuntimeError(f"Could not save menu capture: {image_path}")

    actions = []
    minimum_width = 0.0
    for action in menu.actions():
        text = action.text()
        if not text:
            continue
        rect = menu.actionGeometry(action)
        metrics = menu.fontMetrics()
        advance = float(metrics.horizontalAdvance(text))
        minimum_width = max(minimum_width, advance + 32.0)
        actions.append(
            {
                "text": text,
                "rect": [rect.x(), rect.y(), rect.width(), rect.height()],
                "text_advance_px": advance,
                "width_margin_px": float(rect.width()) - advance,
                "text_fits_with_margin": float(rect.width()) >= advance + 20.0,
                "has_submenu": action.menu() is not None,
            }
        )
    font = menu.font()
    raw_font = QRawFont.fromFont(font)
    codepoints = "₂⁻²³⁴⁵⁶⁸"
    result = {
        "image": image_path.name,
        "menu_size": [menu.width(), menu.height()],
        "screen_available": [geometry.width(), geometry.height()] if geometry else None,
        "screen_dpr": float(screen.devicePixelRatio()) if screen else None,
        "platform": app.platformName(),
        "font_family": font.family(),
        "font_pixel_size": font.pixelSize(),
        "unicode_glyph_coverage": {
            char: bool(raw_font.supportsCharacter(char)) for char in codepoints
        },
        "minimum_text_width_px": minimum_width,
        "actions": actions,
    }
    menu.hide()
    app.processEvents()
    return result


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--out-dir", type=Path, required=True)
    parser.add_argument("--theme", choices=("light", "dark"), required=True)
    args = parser.parse_args()
    args.out_dir.mkdir(parents=True, exist_ok=True)

    app = QApplication.instance() or QApplication([])
    window = QMainWindow()
    apply_theme(window, args.theme)
    toolbar = SymbolPaletteToolbar(QActionGroup(window), window)
    toolbar.setStyleSheet(get_tool_palette_stylesheet(args.theme))
    window.addToolBar(Qt.ToolBarArea.LeftToolBarArea, toolbar)
    window.resize(900, 600)
    window.show()
    app.processEvents()

    root = toolbar.energy_diagram_button.menu()
    presets = _submenu(root, "Electronic Diagram Presets")
    mo_menu = _submenu(presets, "Diatomic MO Presets")
    ligand_menu = _submenu(presets, "Ligand Field Presets")
    captures = {
        "energy_palette": _capture_menu(app, root, args.out_dir, f"energy-palette-{args.theme}"),
        "molecular_orbitals": _capture_menu(app, mo_menu, args.out_dir, f"mo-presets-{args.theme}"),
        "ligand_field": _capture_menu(app, ligand_menu, args.out_dir, f"ligand-presets-{args.theme}"),
    }
    (args.out_dir / f"menu-metrics-{args.theme}.json").write_text(
        json.dumps(captures, ensure_ascii=False, indent=2) + "\n",
        encoding="utf-8",
    )
    window.close()
    app.processEvents()


if __name__ == "__main__":
    main()
