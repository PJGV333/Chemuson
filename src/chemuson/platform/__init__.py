"""Platform-neutral application settings and preference policies."""

from .settings import (
    NamingPreferences,
    NumberingPreferences,
    application_settings,
    load_naming_preferences,
    load_numbering_preferences,
    save_naming_preferences,
    save_numbering_preferences,
    setting_bool,
)

__all__ = [
    "NamingPreferences",
    "NumberingPreferences",
    "application_settings",
    "load_naming_preferences",
    "load_numbering_preferences",
    "save_naming_preferences",
    "save_numbering_preferences",
    "setting_bool",
]
