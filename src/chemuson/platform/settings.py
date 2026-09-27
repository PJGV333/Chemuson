"""Application configuration backed by Qt's portable settings store."""

from __future__ import annotations

from dataclasses import dataclass
from typing import Protocol

from PyQt6.QtCore import QSettings


class SettingsStore(Protocol):
    """Minimal key/value contract used by application configuration policy."""

    def value(self, key: str, default=None): ...

    def setValue(self, key: str, value: object) -> None: ...

    def remove(self, key: str) -> None: ...


@dataclass(frozen=True, slots=True)
class NamingPreferences:
    advanced_enabled: bool = True
    rdkit_isolated: bool = True


@dataclass(frozen=True, slots=True)
class NumberingPreferences:
    mode: str = "atoms"
    include_export: bool = True


@dataclass(frozen=True, slots=True)
class UiPreferences:
    """Preferencias de interfaz (tema) persistidas.

    ``theme`` admite ``"light"``, ``"dark"`` o ``"system"`` (valor lógico
    preparado para la opción "seguir sistema"; el resolve a light/dark lo
    hace la capa de temas, no este módulo).
    """

    theme: str = "light"


@dataclass(frozen=True, slots=True)
class SidePanelPreferences:
    """Persisted page and visibility for the modern right-side panel."""

    active_tab: str = "inspector"
    visible: bool = True


#: Stable page keys accepted by ``ui/side_panel/active_tab``.
SIDE_PANEL_TAB_KEYS: tuple[str, ...] = (
    "inspector",
    "validation",
    "properties",
    "templates",
    "appearance",
    "spectroscopy",
    "compchem",
)


#: Valores válidos para ``ui/theme`` (incluye el valor lógico ``system``).
UI_THEME_CHOICES: tuple[str, ...] = ("light", "dark", "system")


def application_settings() -> QSettings:
    """Create the application's persistent settings store."""
    return QSettings("Chemuson", "Chemuson")


def setting_bool(value: object, default: bool) -> bool:
    """Normalize values returned by QSettings and legacy stores."""
    if value is None:
        return bool(default)
    if isinstance(value, bool):
        return value
    if isinstance(value, str):
        lowered = value.strip().lower()
        if lowered in {"1", "true", "yes", "on", "si", "sí"}:
            return True
        if lowered in {"0", "false", "no", "off"}:
            return False
    try:
        return bool(int(str(value)))
    except (TypeError, ValueError):
        return bool(value)


def load_naming_preferences(settings: SettingsStore) -> NamingPreferences:
    """Read persistent naming preferences without GUI knowledge."""
    return NamingPreferences(
        advanced_enabled=setting_bool(settings.value("naming/advanced_enabled", True), True),
        rdkit_isolated=setting_bool(settings.value("naming/rdkit_isolated", True), True),
    )


def save_naming_preferences(settings: SettingsStore, preferences: NamingPreferences) -> None:
    """Persist naming preferences using the historical keys."""
    settings.setValue("naming/advanced_enabled", bool(preferences.advanced_enabled))
    settings.setValue("naming/rdkit_isolated", bool(preferences.rdkit_isolated))


def load_numbering_preferences(settings: SettingsStore) -> NumberingPreferences:
    """Read and normalize global numbering preferences."""
    mode = str(settings.value("numbering/mode", "atoms") or "atoms").strip().lower()
    if mode not in {"atoms", "structures", "both"}:
        mode = "atoms"
    include_export = setting_bool(settings.value("numbering/include_export", True), True)
    return NumberingPreferences(mode=mode, include_export=include_export)


def save_numbering_preferences(settings: SettingsStore, preferences: NumberingPreferences) -> None:
    """Persist numbering preferences using the historical keys."""
    settings.remove("numbering/enabled")
    settings.setValue("numbering/mode", str(preferences.mode))
    settings.setValue("numbering/include_export", bool(preferences.include_export))


def load_ui_preferences(settings: SettingsStore) -> UiPreferences:
    """Read and normalize persistent UI preferences (theme).

    Valores ausentes o inválidos resuelven a ``"light"`` (comportamiento
    histórico de la aplicación).
    """
    raw = str(settings.value("ui/theme", "light") or "light").strip().lower()
    if raw not in UI_THEME_CHOICES:
        raw = "light"
    return UiPreferences(theme=raw)


def save_ui_preferences(settings: SettingsStore, preferences: UiPreferences) -> None:
    """Persist UI preferences under the ``ui/*`` keys."""
    theme = str(preferences.theme).strip().lower()
    if theme not in UI_THEME_CHOICES:
        theme = "light"
    settings.setValue("ui/theme", theme)


def _side_panel_visible(value: object) -> bool:
    """Normalize side-panel visibility, defaulting unknown values to visible."""
    if isinstance(value, bool):
        return value
    if isinstance(value, int) and value in (0, 1):
        return bool(value)
    if isinstance(value, str):
        lowered = value.strip().lower()
        if lowered in {"1", "true", "yes", "on", "si", "sí"}:
            return True
        if lowered in {"0", "false", "no", "off"}:
            return False
    return True


def load_side_panel_preferences(settings: SettingsStore) -> SidePanelPreferences:
    """Read and normalize the right-side panel's persistent state."""
    active_tab = str(
        settings.value("ui/side_panel/active_tab", "inspector") or "inspector"
    ).strip().lower()
    if active_tab not in SIDE_PANEL_TAB_KEYS:
        active_tab = "inspector"
    visible = _side_panel_visible(settings.value("ui/side_panel/visible", True))
    return SidePanelPreferences(active_tab=active_tab, visible=visible)


def save_side_panel_preferences(
    settings: SettingsStore,
    preferences: SidePanelPreferences,
) -> None:
    """Persist the active page and visibility under ``ui/side_panel/*``."""
    active_tab = str(preferences.active_tab or "inspector").strip().lower()
    if active_tab not in SIDE_PANEL_TAB_KEYS:
        active_tab = "inspector"
    settings.setValue("ui/side_panel/active_tab", active_tab)
    settings.setValue("ui/side_panel/visible", bool(preferences.visible))
