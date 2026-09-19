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
