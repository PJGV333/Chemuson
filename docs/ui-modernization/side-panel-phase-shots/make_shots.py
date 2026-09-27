"""Reproducible Phase 5 screenshots for the production side panel.

Usage (from the repository root):

    QT_QPA_PLATFORM=offscreen PYTHONPATH=src \
        python docs/ui-modernization/side-panel-phase-shots/make_shots.py

The script captures the real ChemusonWindow in light/dark themes at 1440x900,
with Validation, Properties, Appearance, the overflow menu, and at 980x600.
Configuration is isolated under TMPDIR; user preferences are not changed.
These offscreen captures are evidence, not visual approval.
"""
from __future__ import annotations

import os
import sys
import tempfile
from pathlib import Path

os.environ["QT_QPA_PLATFORM"] = "offscreen"
HERE = Path(__file__).resolve().parent
REPO = HERE.parents[2]
sys.path.insert(0, str(REPO / "src"))
settings_home = Path(tempfile.mkdtemp(prefix="chemuson-phase5-settings-"))
os.environ["XDG_CONFIG_HOME"] = str(settings_home)

from PyQt6.QtCore import QEventLoop, QPoint, QTimer  # noqa: E402
from PyQt6.QtGui import QPainter  # noqa: E402
from PyQt6.QtWidgets import QApplication  # noqa: E402

app = QApplication([])
from chemuson.gui.main_window import ChemusonWindow  # noqa: E402


outdir = Path(sys.argv[1]).resolve() if len(sys.argv) > 1 else HERE
outdir.mkdir(parents=True, exist_ok=True)
window = ChemusonWindow()
window.resize(1440, 900)
window.show()
app.processEvents()


def settle(milliseconds: int = 350) -> None:
    loop = QEventLoop()
    QTimer.singleShot(milliseconds, loop.quit)
    loop.exec()
    app.processEvents()


def capture(
    filename: str,
    *,
    theme: str = "light",
    page: str = "inspector",
    size: tuple[int, int] = (1440, 900),
    overflow_open: bool = False,
) -> None:
    window.resize(*size)
    window.current_theme = theme
    window._apply_theme()
    panel = getattr(window, "side_panel")
    panel.show_page(page)
    app.processEvents()
    if overflow_open:
        panel._show_overflow_menu()
    settle()

    image = window.grab()
    if overflow_open:
        menu = panel.overflow_menu
        if not menu.isVisible():
            raise RuntimeError("Overflow menu did not become visible")
        menu_image = menu.grab()
        # The offscreen virtual screen is smaller than the 1440x900 window and
        # clamps QMenu.popup(); composite the real menu widget at its intended
        # anchor in the captured window for a reproducible screenshot.
        button = panel.overflow_button
        button_origin = button.mapTo(window, QPoint(0, 0))
        menu_x = button_origin.x() + button.width() - menu_image.width()
        menu_x = max(0, min(menu_x, window.width() - menu_image.width()))
        menu_point = QPoint(menu_x, button_origin.y() + button.height())
        painter = QPainter(image)
        painter.drawPixmap(menu_point, menu_image)
        painter.end()
        menu.hide()
        app.processEvents()

    path = outdir / filename
    if not image.save(str(path), "PNG"):
        raise OSError(f"Could not save screenshot: {path}")
    print(f"captured {path} ({image.width()}x{image.height()})")


try:
    capture("side-panel-light-1440x900-inspector.png")
    capture("side-panel-dark-1440x900-inspector.png", theme="dark")
    capture("side-panel-light-validation.png", page="validation")
    capture("side-panel-dark-properties.png", theme="dark", page="properties")
    capture(
        "side-panel-overflow-open.png",
        page="spectroscopy",
        overflow_open=True,
    )
    capture("side-panel-appearance.png", page="appearance")
    capture(
        "side-panel-980x600.png",
        size=(980, 600),
    )
finally:
    window.close()
    app.processEvents()
    for thread, _worker, _base_rows in list(
        getattr(window, "_descriptor_jobs", {}).values()
    ):
        if thread.isRunning():
            thread.wait(6000)
    app.processEvents()

print("SHOTS: OK", flush=True)
# PyQt6's offscreen backend can segfault during interpreter teardown after a
# successful render; other UI capture harnesses use the same safe exit pattern.
os._exit(0)
