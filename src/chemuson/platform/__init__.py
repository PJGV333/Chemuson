"""Platform-neutral application settings and preference policies."""

from .settings import (
    NamingPreferences,
    NumberingPreferences,
    UI_THEME_CHOICES,
    UiPreferences,
    application_settings,
    load_naming_preferences,
    load_numbering_preferences,
    load_ui_preferences,
    save_naming_preferences,
    save_numbering_preferences,
    save_ui_preferences,
    setting_bool,
)

__all__ = [
    "NamingPreferences",
    "NumberingPreferences",
    "UI_THEME_CHOICES",
    "UiPreferences",
    "application_settings",
    "load_naming_preferences",
    "load_numbering_preferences",
    "load_ui_preferences",
    "save_naming_preferences",
    "save_numbering_preferences",
    "save_ui_preferences",
    "setting_bool",
]
