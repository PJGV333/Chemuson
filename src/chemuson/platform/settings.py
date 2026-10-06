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


@dataclass(frozen=True, slots=True)
class AIProviderPreferences:
    """Non-secret, per-profile settings for an OpenAI-compatible endpoint."""

    profile_id: str
    base_url: str = ""
    model: str = ""
    timeout_s: int = 60
    supports_json_output: bool = False
    max_tokens: int = 4096


@dataclass(frozen=True, slots=True)
class IdentityVerificationPreferences:
    """Legacy non-secret identity lookup policy; external access defaults off."""

    enabled: bool = True
    allow_external_reference: bool = False


@dataclass(frozen=True, slots=True)
class MolecularAssistantPreferences:
    """Non-secret reference-resolution method and network permission."""

    resolution_method: str = "ai_reference"
    allow_external_reference: bool = False


MOLECULAR_RESOLUTION_METHODS: tuple[str, ...] = (
    "ai",
    "ai_reference",
    "reference",
)


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


def load_ai_provider_preferences(
    settings: SettingsStore,
    profile_id: str,
    *,
    default_base_url: str = "",
) -> AIProviderPreferences:
    """Load bounded provider settings; credentials and prompts are never read."""
    _validate_ai_profile_id(profile_id)
    prefix = f"ai/providers/{profile_id}"
    timeout = _bounded_int(settings.value(f"{prefix}/timeout_s", 60), 60, 10, 600)
    max_tokens = _bounded_int(settings.value(f"{prefix}/max_tokens", 4096), 4096, 64, 8192)
    base_url = settings.value(f"{prefix}/base_url", default_base_url)
    model = settings.value(f"{prefix}/model", "")
    return AIProviderPreferences(
        profile_id=profile_id,
        base_url=str(base_url or "")[:2048],
        model=str(model or "")[:256],
        timeout_s=timeout,
        supports_json_output=setting_bool(
            settings.value(f"{prefix}/supports_json_output", False), False
        ),
        max_tokens=max_tokens,
    )


def save_ai_provider_preferences(
    settings: SettingsStore,
    preferences: AIProviderPreferences,
) -> None:
    """Persist only endpoint/runtime choices, explicitly removing legacy key storage."""
    _validate_ai_profile_id(preferences.profile_id)
    if not 10 <= int(preferences.timeout_s) <= 600:
        raise ValueError("timeout_s must be between 10 and 600")
    if not 64 <= int(preferences.max_tokens) <= 8192:
        raise ValueError("max_tokens must be between 64 and 8192")
    prefix = f"ai/providers/{preferences.profile_id}"
    settings.setValue(f"{prefix}/base_url", str(preferences.base_url)[:2048])
    settings.setValue(f"{prefix}/model", str(preferences.model)[:256])
    settings.setValue(f"{prefix}/timeout_s", int(preferences.timeout_s))
    settings.setValue(
        f"{prefix}/supports_json_output", bool(preferences.supports_json_output)
    )
    settings.setValue(f"{prefix}/max_tokens", int(preferences.max_tokens))
    settings.remove(f"{prefix}/api_key")


def _validate_ai_profile_id(profile_id: str) -> None:
    if (
        not isinstance(profile_id, str)
        or not profile_id
        or len(profile_id) > 64
        or any(not (char.isascii() and (char.isalnum() or char in "-_")) for char in profile_id)
    ):
        raise ValueError("profile_id is invalid")


def _bounded_int(value: object, default: int, minimum: int, maximum: int) -> int:
    if isinstance(value, bool):
        return default
    try:
        normalized = int(value)
    except (TypeError, ValueError, OverflowError):
        return default
    return normalized if minimum <= normalized <= maximum else default


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


def load_identity_verification_preferences(
    settings: SettingsStore,
) -> IdentityVerificationPreferences:
    """Load offline-by-default molecular identity policy."""
    enabled = setting_bool(settings.value("ai/identity_verification/enabled", True), True)
    allow_external = setting_bool(
        settings.value("ai/identity_verification/allow_external_reference", False),
        False,
    )
    return IdentityVerificationPreferences(
        enabled=enabled,
        allow_external_reference=enabled and allow_external,
    )


def save_identity_verification_preferences(
    settings: SettingsStore,
    preferences: IdentityVerificationPreferences,
) -> None:
    """Persist identity policy only; it contains no provider credential."""
    settings.setValue("ai/identity_verification/enabled", bool(preferences.enabled))
    settings.setValue(
        "ai/identity_verification/allow_external_reference",
        bool(preferences.enabled and preferences.allow_external_reference),
    )


def load_molecular_assistant_preferences(
    settings: SettingsStore,
) -> MolecularAssistantPreferences:
    """Load the offline-by-default resolution route without reading secrets."""
    stored_method = settings.value("ai/molecular_assistant/resolution_method", None)
    if stored_method is None:
        legacy_enabled = setting_bool(
            settings.value("ai/identity_verification/enabled", True),
            True,
        )
        method = "ai_reference" if legacy_enabled else "ai"
    else:
        method = str(stored_method or "ai_reference").strip()
    if method not in MOLECULAR_RESOLUTION_METHODS:
        method = "ai_reference"
    allow_external = setting_bool(
        settings.value("ai/identity_verification/allow_external_reference", False),
        False,
    )
    return MolecularAssistantPreferences(method, allow_external)


def save_molecular_assistant_preferences(
    settings: SettingsStore,
    preferences: MolecularAssistantPreferences,
) -> None:
    """Persist only the selected route and explicit PubChem network permission."""
    method = str(preferences.resolution_method or "ai_reference").strip()
    if method not in MOLECULAR_RESOLUTION_METHODS:
        raise ValueError("resolution_method is invalid")
    if not isinstance(preferences.allow_external_reference, bool):
        raise ValueError("allow_external_reference must be boolean")
    settings.setValue("ai/molecular_assistant/resolution_method", method)
    settings.setValue(
        "ai/identity_verification/allow_external_reference",
        bool(preferences.allow_external_reference),
    )


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
