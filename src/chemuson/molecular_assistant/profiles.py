"""Explicit endpoint presets for OpenAI-compatible Chat Completions APIs."""

from __future__ import annotations

from dataclasses import dataclass


@dataclass(frozen=True)
class OpenAICompatibleProfile:
    """Provider identity and editable default endpoint; model IDs stay endpoint-owned."""

    profile_id: str
    display_name: str
    default_base_url: str
    api_key_required: bool = False


OPENAI_COMPATIBLE_PROFILES: tuple[OpenAICompatibleProfile, ...] = (
    OpenAICompatibleProfile(
        profile_id="openai-compatible",
        display_name="OpenAI-compatible (personalizado)",
        default_base_url="",
    ),
    OpenAICompatibleProfile(
        profile_id="openai",
        display_name="OpenAI",
        default_base_url="https://api.openai.com/v1",
        api_key_required=True,
    ),
    OpenAICompatibleProfile(
        profile_id="lm-studio",
        display_name="LM Studio",
        default_base_url="http://127.0.0.1:1234/v1",
    ),
    OpenAICompatibleProfile(
        profile_id="llama-cpp",
        display_name="llama.cpp server",
        default_base_url="http://127.0.0.1:8080/v1",
    ),
)

_PROFILES_BY_ID = {profile.profile_id: profile for profile in OPENAI_COMPATIBLE_PROFILES}


def get_openai_compatible_profile(profile_id: str) -> OpenAICompatibleProfile | None:
    """Return a known profile, or ``None`` for an unsupported identifier."""
    if not isinstance(profile_id, str):
        return None
    return _PROFILES_BY_ID.get(profile_id)
