"""Capturas de referencia de la Fase 3 (app bar + pestañas de documento).

Uso (desde la raíz del repo):

    QT_QPA_PLATFORM=offscreen PYTHONPATH=src \
        python docs/ui-modernization/app-bar-phase-shots/make_shots.py <outdir>

Crea ``ChemusonWindow``, monta 3 pestañas (una nombrada, una modificada
vía ``QUndoStack`` real) y captura: ventana completa light/dark a 1440×900,
ventana completa light a 980×600 (comportamiento estrecho) y close-ups de
la ``AppBar`` light/dark (con y sin historial undo/redo).

Solo lectura del estado de la UI; no modifica producción.
"""
import os
import sys
from pathlib import Path

os.environ["QT_QPA_PLATFORM"] = "offscreen"
sys.path.insert(0, str(Path(__file__).resolve().parents[3] / "src"))

from PyQt6.QtCore import QEventLoop, QTimer
from PyQt6.QtGui import QUndoCommand
from PyQt6.QtWidgets import QApplication

app = QApplication([])
from chemuson.gui.main_window import ChemusonWindow


class _NoopCommand(QUndoCommand):
    def redo(self) -> None:
        return

    def undo(self) -> None:
        return


win = ChemusonWindow()
outdir = Path(sys.argv[1])
outdir.mkdir(parents=True, exist_ok=True)

# El tema persistido (M21) puede ser cualquiera: la captura termina
# restaurando el valor original para no alterar la preferencia del usuario.
_original_theme = win.current_theme

# Escenario: 3 pestañas — una con nombre de archivo (solo título), una
# modificada (comando undo real -> suciedad real), una limpia.
win._on_file_new()
win._on_file_new()
app.processEvents()
win._tab_manager.set_canvas_file_path(win._canvas_from_tab_index(0), "/ficticio/Proyecto A.cmsn")
win._tab_manager.update_tab_title(win._canvas_from_tab_index(0))
win._canvas_from_tab_index(1).undo_stack.push(_NoopCommand())
win._tab_manager.update_tab_title(win._canvas_from_tab_index(1))
app.processEvents()


def _snap(closures: list) -> None:
    for label, widget in closures:
        widget.grab().save(str(outdir / f"{label}"))
    print("saved:", ", ".join(label for label, _ in closures))


def capture(theme_name: str, size: tuple, label_prefix: str, with_history: bool = False) -> None:
    loop = QEventLoop()

    def _do() -> None:
        win.current_theme = theme_name
        win._apply_theme()
        app.processEvents()
        win.resize(*size)
        app.processEvents()
        if with_history:
            win.canvas.undo_stack.push(_NoopCommand())
            win._update_tab_title(win.canvas)
        else:
            win.canvas.undo_stack.setClean()
            win._update_tab_title(win.canvas)
        app.processEvents()

        def _take() -> None:
            _snap(
                [
                    (f"{label_prefix}-full-{theme_name}.png", win),
                    (f"{label_prefix}-appbar-{theme_name}.png", win.app_bar),
                ]
            )
            loop.quit()

        QTimer.singleShot(600, _take)

    _do()
    loop.exec()
    print(f"capture {theme_name} {size} ok")


# Ventana completa 1440x900 light/dark (con historial: undo habilitado).
capture("light", (1440, 900), "appbar-phase", with_history=True)
capture("dark", (1440, 900), "appbar-phase", with_history=True)
# Comportamiento estrecho 980x600 light.
capture("light", (980, 600), "appbar-phase-narrow")

# Limpia todo antes de cerrar (evita el diálogo de cambios sin guardar)
y restaura el tema original persistido.
for canvas in list(win._tab_manager.iter_canvases()):
    canvas.undo_stack.setClean()
win.current_theme = _original_theme
win._apply_theme()
win.close()
print("SHOTS: OK")
