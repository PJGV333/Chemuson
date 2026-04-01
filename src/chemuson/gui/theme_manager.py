"""Gestión de temas internos (light/dark/auto) para Chemuson."""

from __future__ import annotations

import os
import platform
import subprocess

from PyQt6.QtCore import QSettings
from PyQt6.QtGui import QColor, QPalette
from PyQt6.QtWidgets import QApplication, QStyleFactory

from chemuson.gui.styles import build_main_stylesheet, theme_tokens

THEME_MODE_KEY = "appearance/theme_mode"
VALID_THEME_MODES = {"light", "dark", "auto"}
VALID_THEMES = {"light", "dark"}
_CURRENT_THEME = "light"


def load_theme_preference(settings: QSettings) -> str:
    mode = str(settings.value(THEME_MODE_KEY, "auto") or "auto").strip().lower()
    return mode if mode in VALID_THEME_MODES else "auto"


def save_theme_preference(settings: QSettings, mode: str) -> str:
    normalized = normalize_mode(mode)
    settings.setValue(THEME_MODE_KEY, normalized)
    return normalized


def normalize_mode(mode: str) -> str:
    value = str(mode or "").strip().lower()
    return value if value in VALID_THEME_MODES else "auto"


def resolve_effective_theme(mode_or_theme: str) -> str:
    value = str(mode_or_theme or "").strip().lower()
    if value in VALID_THEMES:
        return value
    return "dark" if system_prefers_dark() else "light"


def current_theme() -> str:
    return _CURRENT_THEME


def _build_palette(theme: str) -> QPalette:
    tokens = theme_tokens(theme)
    p = QPalette()
    p.setColor(QPalette.ColorRole.Window, QColor(tokens.window))
    p.setColor(QPalette.ColorRole.WindowText, QColor(tokens.text))
    p.setColor(QPalette.ColorRole.Base, QColor(tokens.surface))
    p.setColor(QPalette.ColorRole.AlternateBase, QColor(tokens.surface_alt))
    p.setColor(QPalette.ColorRole.ToolTipBase, QColor(tokens.tooltip_bg))
    p.setColor(QPalette.ColorRole.ToolTipText, QColor(tokens.tooltip_text))
    p.setColor(QPalette.ColorRole.Text, QColor(tokens.text))
    p.setColor(QPalette.ColorRole.Button, QColor(tokens.panel))
    p.setColor(QPalette.ColorRole.ButtonText, QColor(tokens.text))
    p.setColor(QPalette.ColorRole.BrightText, QColor(tokens.text_inverse))
    p.setColor(QPalette.ColorRole.Highlight, QColor(tokens.accent))
    p.setColor(QPalette.ColorRole.HighlightedText, QColor(tokens.text_inverse))
    p.setColor(QPalette.ColorRole.Link, QColor(tokens.accent))
    return p


def apply_theme(app: QApplication, mode_or_theme: str) -> str:
    global _CURRENT_THEME
    theme = resolve_effective_theme(mode_or_theme)
    fusion = QStyleFactory.create("Fusion")
    if fusion is not None:
        app.setStyle(fusion)
    app.setPalette(_build_palette(theme))
    app.setStyleSheet(build_main_stylesheet(theme_tokens(theme)))
    _CURRENT_THEME = theme
    return theme


def system_prefers_dark() -> bool:
    qt_app = QApplication.instance()
    if qt_app is not None:
        try:
            base = qt_app.styleHints().colorScheme().name.lower()
            if "dark" in base:
                return True
            if "light" in base:
                return False
        except Exception:
            pass
    if platform.system().lower() == "linux":
        plasma = _plasma_prefers_dark()
        if plasma is not None:
            return plasma
    return False


def _plasma_prefers_dark() -> bool | None:
    config = os.path.expanduser("~/.config/kdeglobals")
    if os.path.exists(config):
        try:
            with open(config, "r", encoding="utf-8") as fh:
                in_section = False
                for raw in fh:
                    line = raw.strip()
                    if line.startswith("["):
                        in_section = line.lower() == "[general]"
                        continue
                    if in_section and line.lower().startswith("colorscheme="):
                        scheme = line.split("=", 1)[1].strip().lower()
                        if scheme:
                            return "dark" in scheme or scheme.endswith("night")
        except Exception:
            return None
    try:
        result = subprocess.run(
            ["gsettings", "get", "org.gnome.desktop.interface", "color-scheme"],
            check=False,
            capture_output=True,
            text=True,
            timeout=1.2,
        )
        out = (result.stdout or "").strip().lower()
        if "dark" in out:
            return True
        if "light" in out or "default" in out:
            return False
    except Exception:
        return None
    return None
