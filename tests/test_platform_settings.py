from __future__ import annotations

from dataclasses import dataclass

from chemuson.platform.settings import (
    NamingPreferences,
    NumberingPreferences,
    load_naming_preferences,
    load_numbering_preferences,
    save_naming_preferences,
    save_numbering_preferences,
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
