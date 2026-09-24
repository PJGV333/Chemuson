from __future__ import annotations

from dataclasses import dataclass

from chemuson.platform.settings import (
    NamingPreferences,
    NumberingPreferences,
    UI_THEME_CHOICES,
    UiPreferences,
    load_naming_preferences,
    load_numbering_preferences,
    load_ui_preferences,
    save_naming_preferences,
    save_numbering_preferences,
    save_ui_preferences,
    setting_bool,
)


@dataclass
class FakeSettings:
    values: dict[str, object]

    def value(self, key: str, default=None):
        return self.values.get(key, default)

    def setValue(self, key: str, value: object) -> None:
        self.values[key] = value

    def remove(self, key: str) -> None:
        self.values.pop(key, None)


def test_setting_bool_preserves_legacy_qsettings_values() -> None:
    assert setting_bool("sí", False) is True
    assert setting_bool("0", True) is False
    assert setting_bool(None, True) is True


def test_naming_preferences_load_and_save_roundtrip() -> None:
    settings = FakeSettings({"naming/advanced_enabled": "0", "naming/rdkit_isolated": True})

    assert load_naming_preferences(settings) == NamingPreferences(False, True)

    save_naming_preferences(settings, NamingPreferences(True, False))
    assert settings.values["naming/advanced_enabled"] is True
    assert settings.values["naming/rdkit_isolated"] is False


def test_numbering_preferences_normalize_and_save() -> None:
    settings = FakeSettings({"numbering/mode": "invalid", "numbering/include_export": "no"})

    assert load_numbering_preferences(settings) == NumberingPreferences("atoms", False)

    save_numbering_preferences(settings, NumberingPreferences("both", True))
    assert settings.values["numbering/mode"] == "both"
    assert settings.values["numbering/include_export"] is True
    assert "numbering/enabled" not in settings.values


def test_ui_preferences_roundtrip_and_normalize() -> None:
    assert set(UI_THEME_CHOICES) == {"light", "dark", "system"}

    assert load_ui_preferences(FakeSettings({"ui/theme": "dark"})) == UiPreferences(theme="dark")
    assert load_ui_preferences(FakeSettings({"ui/theme": "system"})) == UiPreferences(theme="system")

    # Valores inválidos o ausentes resuelven a light (comportamiento histórico).
    assert load_ui_preferences(FakeSettings({"ui/theme": "neon"})) == UiPreferences(theme="light")
    assert load_ui_preferences(FakeSettings({})) == UiPreferences(theme="light")

    settings = FakeSettings({})
    save_ui_preferences(settings, UiPreferences(theme="system"))
    assert settings.values["ui/theme"] == "system"

    save_ui_preferences(settings, UiPreferences(theme="invalid"))
    assert settings.values["ui/theme"] == "light"

    # El valor persistido se relee normalizado.
    settings = FakeSettings({"ui/theme": "DARK"})
    assert load_ui_preferences(settings).theme == "dark"
