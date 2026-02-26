"""Pruebas de utilidades SemVer para el módulo de update."""

import os
import sys

import pytest

sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), "..", "src")))

from chemuson.update.semver import (
    channel_accepts_version,
    compare_versions,
    is_newer_version,
    is_prerelease,
    parse_semver,
)
from chemuson.update.types import UpdateChannel


def test_compare_versions_respects_prerelease_order() -> None:
    assert compare_versions("1.2.0", "1.2.0") == 0
    assert compare_versions("1.2.1", "1.2.0") > 0
    assert compare_versions("1.2.0-beta.1", "1.2.0-alpha.9") > 0
    assert compare_versions("1.2.0", "1.2.0-rc.1") > 0


def test_is_newer_version_and_prerelease_detection() -> None:
    assert is_newer_version("2.0.0", "1.9.9") is True
    assert is_newer_version("1.0.0-beta.1", "1.0.0-alpha.5") is True
    assert is_prerelease("1.0.0-beta.1") is True
    assert is_prerelease("v1.0.0") is False


def test_channel_filtering_prerelease() -> None:
    assert channel_accepts_version(UpdateChannel.STABLE, "1.3.0") is True
    assert channel_accepts_version(UpdateChannel.STABLE, "1.3.0-beta.1") is False
    assert channel_accepts_version(UpdateChannel.BETA, "1.3.0-beta.1") is True


def test_parse_semver_rejects_invalid_values() -> None:
    with pytest.raises(ValueError):
        parse_semver("1.2")
    with pytest.raises(ValueError):
        parse_semver("abc")

