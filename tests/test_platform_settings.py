from __future__ import annotations

from dataclasses import dataclass

import pytest
from PyQt6.QtCore import QSettings

from chemuson.platform.settings import (
    AIProviderPreferences,
    IdentityVerificationPreferences,
    MolecularAssistantPreferences,
    NamingPreferences,
    NumberingPreferences,
    UI_THEME_CHOICES,
    UiPreferences,
    SidePanelPreferences,
    load_ai_provider_preferences,
    load_identity_verification_preferences,
    load_molecular_assistant_preferences,
    load_naming_preferences,
    load_numbering_preferences,
    load_side_panel_preferences,
    load_ui_preferences,
    save_ai_provider_preferences,
    save_identity_verification_preferences,
    save_molecular_assistant_preferences,
    save_naming_preferences,
    save_numbering_preferences,
    save_side_panel_preferences,
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


def test_ai_provider_preferences_roundtrip_per_profile_without_secrets() -> None:
    settings = FakeSettings({
        "ai/providers/llama-cpp/api_key": "legacy-secret",
    })
    defaults = load_ai_provider_preferences(
        settings,
        "llama-cpp",
        default_base_url="http://127.0.0.1:8080/v1",
    )
    assert defaults == AIProviderPreferences(
        "llama-cpp", "http://127.0.0.1:8080/v1", "", 60, False, 4096
    )

    saved = AIProviderPreferences(
        "llama-cpp", "http://127.0.0.1:8081/v1", "local-model", 180, True, 2048
    )
    save_ai_provider_preferences(settings, saved)
    assert load_ai_provider_preferences(settings, "llama-cpp") == saved
    assert "ai/providers/llama-cpp/api_key" not in settings.values
    assert "ai/providers/openai/base_url" not in settings.values


def test_identity_verification_preferences_default_offline_and_roundtrip() -> None:
    settings = FakeSettings({})
    assert load_identity_verification_preferences(settings) == IdentityVerificationPreferences()

    preferences = IdentityVerificationPreferences(
        enabled=True,
        allow_external_reference=True,
    )
    save_identity_verification_preferences(settings, preferences)
    assert load_identity_verification_preferences(settings) == preferences
    assert settings.values == {
        "ai/identity_verification/enabled": True,
        "ai/identity_verification/allow_external_reference": True,
    }

    save_identity_verification_preferences(
        settings,
        IdentityVerificationPreferences(enabled=False, allow_external_reference=True),
    )
    assert load_identity_verification_preferences(settings) == IdentityVerificationPreferences(
        enabled=False,
        allow_external_reference=False,
    )


def test_qsettings_provider_preferences_never_persist_api_keys(tmp_path) -> None:
    path = tmp_path / "settings.ini"
    settings = QSettings(str(path), QSettings.Format.IniFormat)
    settings.setValue("ai/providers/openai/api_key", "legacy-sensitive-key")
    save_ai_provider_preferences(
        settings,
        AIProviderPreferences("openai", "https://api.openai.com/v1", "model"),
    )
    save_identity_verification_preferences(
        settings,
        IdentityVerificationPreferences(True, False),
    )
    settings.sync()

    assert not settings.contains("ai/providers/openai/api_key")
    assert "legacy-sensitive-key" not in path.read_text(encoding="utf-8")
    assert load_identity_verification_preferences(settings) == IdentityVerificationPreferences()


def test_molecular_assistant_route_preferences_are_offline_by_default_and_nonsecret():
    settings = FakeSettings({})
    assert load_molecular_assistant_preferences(settings) == MolecularAssistantPreferences()
    assert load_molecular_assistant_preferences(
        FakeSettings({"ai/identity_verification/enabled": False})
    ).resolution_method == "ai"

    preferences = MolecularAssistantPreferences("reference", True)
    save_molecular_assistant_preferences(settings, preferences)
    assert load_molecular_assistant_preferences(settings) == preferences
    assert settings.values == {
        "ai/molecular_assistant/resolution_method": "reference",
        "ai/identity_verification/allow_external_reference": True,
    }

    settings.values["ai/molecular_assistant/resolution_method"] = "arbitrary"
    assert load_molecular_assistant_preferences(settings).resolution_method == "ai_reference"
    with pytest.raises(ValueError):
        save_molecular_assistant_preferences(
            settings,
            MolecularAssistantPreferences("ai", "yes"),
        )


def test_ai_provider_preferences_normalize_out_of_range_settings() -> None:
    settings = FakeSettings({
        "ai/providers/lm-studio/timeout_s": 601,
        "ai/providers/lm-studio/max_tokens": 1,
        "ai/providers/lm-studio/supports_json_output": "yes",
    })
    preferences = load_ai_provider_preferences(settings, "lm-studio")
    assert preferences.timeout_s == 60
    assert preferences.max_tokens == 4096
    assert preferences.supports_json_output is True


def test_ai_provider_preferences_reject_invalid_save() -> None:
    settings = FakeSettings({})
    with pytest.raises(ValueError):
        save_ai_provider_preferences(
            settings,
            AIProviderPreferences("openai", timeout_s=9),
        )


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


def test_side_panel_preferences_default_and_roundtrip() -> None:
    assert load_side_panel_preferences(FakeSettings({})) == SidePanelPreferences(
        active_tab="inspector",
        visible=True,
    )

    settings = FakeSettings(
        {
            "ui/side_panel/active_tab": "compchem",
            "ui/side_panel/visible": "0",
        }
    )
    assert load_side_panel_preferences(settings) == SidePanelPreferences(
        active_tab="compchem",
        visible=False,
    )

    save_side_panel_preferences(
        settings,
        SidePanelPreferences(active_tab="appearance", visible=True),
    )
    assert settings.values["ui/side_panel/active_tab"] == "appearance"
    assert settings.values["ui/side_panel/visible"] is True


def test_side_panel_preferences_invalid_values_use_safe_defaults() -> None:
    settings = FakeSettings(
        {
            "ui/side_panel/active_tab": "unknown",
            "ui/side_panel/visible": "",
        }
    )
    assert load_side_panel_preferences(settings) == SidePanelPreferences(
        active_tab="inspector",
        visible=True,
    )

    save_side_panel_preferences(
        settings,
        SidePanelPreferences(active_tab="unknown", visible=False),
    )
    assert settings.values["ui/side_panel/active_tab"] == "inspector"
    assert settings.values["ui/side_panel/visible"] is False
