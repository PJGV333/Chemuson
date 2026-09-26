"""Capturas de referencia de la Fase 4 (rail de herramientas + flyouts).

Uso (desde la raíz del repo):

    QT_QPA_PLATFORM=offscreen PYTHONPATH=src \
        python docs/ui-modernization/tool-rail-phase-shots/make_shots.py <outdir>

Crea ``ChemusonWindow`` (con el rail montado y las toolbars históricas
ocultas) y captura:

- Ventana completa light/dark a 1440×900 (estado inicial: herramienta de
  selección activa);
- Close-up del rail (``win.tool_rail``) light/dark;
- Flyout de enlaces abierto (celda "Enlace doble" resaltada tras el clic)
  y flyout de diagramas de energía abierto (con pie de preset), light/dark.

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
_original_theme = win.current_theme


def _snap(items: list) -> None:
    for label, widget in items:
        widget.grab().save(str(outdir / f"{label}"))
    print("saved:", ", ".join(label for label, _ in items))


def capture(theme_name: str, label_prefix: str, stage) -> None:
    loop = QEventLoop()

    def _do() -> None:
        win.current_theme = theme_name
        win._apply_theme()
        app.processEvents()
        win.resize(1440, 900)
        app.processEvents()
        stage_widget = stage()  # aplica el escenario y devuelve el widget a capturar (None = ventana)
        app.processEvents()

        def _take() -> None:
            if stage_widget is None:
                _snap(
                    [
                        (f"{label_prefix}-full-{theme_name}.png", win),
                        (f"{label_prefix}-rail-{theme_name}.png", win.tool_rail),
                    ]
                )
            else:
                _snap([(f"{label_prefix}-{theme_name}.png", stage_widget)])
            loop.quit()

        QTimer.singleShot(600, _take)

    _do()
    loop.exec()
    print(f"capture {theme_name} {label_prefix} ok")


def _stage_full() -> None:
    # Estado inicial: herramienta de selección activa (devuelve None: se
    # captura la ventana completa + close-up del rail).
    win.toolbar.select_action.trigger()
    app.processEvents()
    for flyout in win.tool_rail._flyouts.values():
        flyout.close_with(None)
    return None


def _stage_bond() -> "QWidget":
    # Flyout de enlaces abierto con la celda «Enlace doble» resaltada.
    win.tool_rail._buttons["bond"].set_active(True)
    flyout = win.tool_rail.open_flyout("bond")
    labels = [c._label.text().replace("\n", " ").lower() for c in flyout._cells]
    index = next(i for i, text in enumerate(labels) if "doble" in text)
    flyout._cells[index].set_active(True)
    return flyout


def _stage_energy() -> "QWidget":
    # Flyout de diagramas de energía abierto con su pie de preset.
    win.tool_rail._buttons["energy"].set_active(True)
    return win.tool_rail.open_flyout("energy")


# Ventana completa (light/dark) con el rail.
for theme in ("light", "dark"):
    capture(theme, "rail-phase-full", _stage_full)
# Close-ups de flyouts (light/dark).
for theme in ("light", "dark"):
    capture(theme, "rail-phase-bond", _stage_bond)
    capture(theme, "rail-phase-energy", _stage_energy)

# Limpia todo antes de cerrar (evita el diálogo de cambios sin guardar)
# y restaura el tema original persistido.
for canvas in list(win._tab_manager.iter_canvases()):
    canvas.undo_stack.setClean()
for flyout in win.tool_rail._flyouts.values():
    flyout.close_with(None)
win.current_theme = _original_theme
win._apply_theme()
win.close()
print("SHOTS: OK")
