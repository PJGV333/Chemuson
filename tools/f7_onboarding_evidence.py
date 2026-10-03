"""Evidencia visual del onboarding y del clic simple en Plantillas (Fase 7).

Genera, en offscreen (KDE-like: el equipo de la fase corre KDE/Wayland a DPR 2,
aquí se simula con ``QT_SCALE_FACTOR``):

- ``onboarding_light_step{1,2,3}.png`` y ``onboarding_dark_step{1,2,3}.png``:
  ventana completa con el overlay (máscara + agujero + tarjeta) en cada paso.
  ``onboarding_step1.png`` (el nombre ya comprometido en el OpenSpec) se
  regenera con la captura light paso 1.
- ``onboarding_light_step{1,2,3}_dpr2.png`` y ``onboarding_dark_step{1,2,3}_dpr2.png``:
  las mismas capturas a DPR 2 (el backend donde aparecían las franjas negras).
- ``templates_click_simple.png``: SidePanel (página Plantillas) con un item de
  plantilla activado por un clic simple (QTest), mostrando la selección.

Cada escala se captura en un subprocess porque el devicePixelRatio se fija al
crear la ``QApplication``.

Uso: ``python tools/f7_onboarding_evidence.py``
"""
from __future__ import annotations

import os
import subprocess
import sys
import tempfile

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
_THEMES = ("light", "dark")
_SCALERS = (1, 2)


def _build_window(theme: str):
    from PyQt6.QtCore import QSettings

    from chemuson.gui.main_window import ChemusonWindow

    # QSettings aislado en un directorio temporal: la evidencia no toca la
    # configuración real del usuario (la clave ``ui/onboarding/completed``).
    scratch = os.environ.get("_F7_EVIDENCE_SETTINGS_DIR")
    if scratch:
        QSettings.setPath(QSettings.Format.NativeFormat, QSettings.Scope.UserScope, scratch)

    settings = QSettings("Chemuson", "Chemuson")
    settings.remove(_ONBOARD_KEY)
    settings.sync()

    win = ChemusonWindow()
    win.toggle_theme(theme == "dark")
    win.resize(1440, 900)
    win.show()
    return win


def _capture(app, name: str, widget) -> None:
    pixmap = widget.grab()
    path = os.path.join(OUT, name)
    pixmap.save(path)
    print(
        f"OK {name} ({pixmap.width()}x{pixmap.height()}, "
        f"DPR {app.primaryScreen().devicePixelRatio()})",
        flush=True,
    )


def _child_capture() -> None:
    sys.path.insert(0, os.path.join(_REPO_ROOT, "src"))

    from PyQt6.QtWidgets import QApplication

    app = QApplication(sys.argv)
    scale = int(os.environ.get("QT_SCALE_FACTOR", "1"))
    suffix = "" if scale == 1 else "_dpr2"

    from chemuson.gui.onboarding import OnboardingOverlay

    wins = []  # se mantienen vivas: destruir una ventana con su overlay hijo
    # activo durante el bucle provoca un cierre inestable en offscreen.
    for theme in _THEMES:
        win = _build_window(theme)
        app.processEvents()
        overlay = OnboardingOverlay(win, [win.tool_rail, win.canvas, win.side_panel])
        overlay.show()
        app.processEvents()
        for step in range(3):
            # Se navega con la API pública (``advance``), que actualiza la
            # tarjeta y el layout: así la captura muestra el paso real.
            if step > 0:
                overlay.advance()
            app.processEvents()
            _capture(app, f"onboarding_{theme}_step{step + 1}{suffix}.png", win)
            if theme == "light" and step == 0 and scale == 1:
                # Nombre ya comprometido en el OpenSpec (evidence/onboarding_step1.png).
                _capture(app, "onboarding_step1.png", win)
        wins.append((win, overlay))

    if scale != 1:
        os._exit(0)

    # Clic simple sobre una plantilla (SidePanel, página Plantillas).
    from PyQt6.QtCore import Qt
    from PyQt6.QtTest import QTest

    win = _build_window("light")
    win.side_panel.show_page("templates")
    app.processEvents()

    dock = win.templates_dock
    # No se muta la biblioteca: se usa la primera plantilla ya existente.
    target = None
    for i in range(dock.tree.topLevelItemCount()):
        group = dock.tree.topLevelItem(i)
        for j in range(group.childCount()):
            child = group.child(j)
            payload = child.data(0, Qt.ItemDataRole.UserRole)
            if isinstance(payload, dict) and payload.get("kind") == "template":
                target = child
                break
        if target is not None:
            break

    assert target is not None, "el árbol no tiene plantillas que seleccionar"
    dock.tree.setCurrentItem(target)
    QTest.mouseClick(
        dock.tree.viewport(),
        Qt.MouseButton.LeftButton,
        Qt.KeyboardModifier.NoModifier,
        dock.tree.visualItemRect(target).center(),
    )
    app.processEvents()
    _capture(app, "templates_click_simple.png", win.side_panel)

    # Los widgets se mantienen vivos hasta aquí. Destruir la ventana con el
    # overlay translúcido activo en offscreen a DPR 2 produce un cierre
    # inestable (SIGSEGV en la limpieza de Qt), así que el proceso termina sin
    # limpieza explícita: las capturas ya están escritas y volcadas.
    os._exit(0)


def main() -> None:
    if os.environ.get("_F7_ONBOARDING_EVIDENCE_CHILD") == "1":
        _child_capture()
        return

    env = dict(os.environ, QT_QPA_PLATFORM="offscreen", _F7_ONBOARDING_EVIDENCE_CHILD="1")
    scratch = tempfile.mkdtemp(prefix="f7_onboarding_evidence_")
    env["_F7_EVIDENCE_SETTINGS_DIR"] = scratch
    for scale in _SCALERS:
        subprocess.run(
            [sys.executable, os.path.abspath(__file__)],
            env=dict(env, QT_SCALE_FACTOR=str(scale)),
            check=True,
        )
    print("DONE", flush=True)


if __name__ == "__main__":
    main()
