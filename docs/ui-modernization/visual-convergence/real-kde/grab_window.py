"""Brief QWidget.grab capture using the real current Qt session backend."""
from __future__ import annotations
import os, sys
from pathlib import Path
if os.environ.get("QT_QPA_PLATFORM", "").lower() == "offscreen":
    raise SystemExit("Refusing offscreen capture")
from PyQt6.QtCore import QTimer
from PyQt6.QtWidgets import QApplication
ROOT=Path(__file__).resolve().parents[5]
sys.path.insert(0,str(ROOT/"src"))
from chemuson.gui.main_window import ChemusonWindow
mode=sys.argv[1] if len(sys.argv)>1 else "light"
if mode not in {"light","dark"}: raise SystemExit("usage: grab_window.py light|dark")
app=QApplication(sys.argv[:1]); win=ChemusonWindow()
win.resize(1200,760); win.show()
out=Path(__file__).resolve().parent/f"real-qt-wayland-{mode}.png"
def capture():
    win.current_theme=mode; win._apply_theme(); app.processEvents()
    ok=win.grab().save(str(out)); handle=win.windowHandle()
    dpr=handle.devicePixelRatio() if handle else "unknown"
    print(f"platform={app.platformName()} dpr={dpr} saved={out} ok={ok}",flush=True)
    win.close(); app.quit()
QTimer.singleShot(500,capture); QTimer.singleShot(1800,app.quit); app.exec()
