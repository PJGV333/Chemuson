"""Smoke tests para distribución y updater en Windows."""

import os
import sys

sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), "..", "src")))

from chemuson.update.core import AutoUpdateCore
from chemuson.update.provider import GitHubReleasesProvider
from chemuson.update.types import (
    ReleaseAsset,
    ReleaseInfo,
    UpdateChannel,
    UpdateMode,
    UpdateSettings,
)
from chemuson.update.windows import (
    INSTALL_MARKER_FILENAME,
    build_inno_installer_command,
    choose_windows_asset_flavor,
    detect_windows_install_context,
)


def _windows_release(version: str = "1.5.0") -> ReleaseInfo:
    return ReleaseInfo(
        tag_name=f"v{version}",
        version=version,
        prerelease=False,
        published_at="2026-02-26T00:00:00Z",
        html_url="https://example.invalid/release",
        assets=[
            ReleaseAsset(
                name=f"Chemuson-v{version}-windows-x86_64-portable.exe",
                url="https://example.invalid/portable.exe",
            ),
            ReleaseAsset(
                name=f"Chemuson-v{version}-windows-x86_64-setup.exe",
                url="https://example.invalid/setup.exe",
            ),
        ],
    )


def test_detect_windows_install_context_with_marker(tmp_path) -> None:
    exe_path = tmp_path / "Chemuson.exe"
    exe_path.write_text("stub", encoding="utf-8")
    marker = tmp_path / INSTALL_MARKER_FILENAME
    marker.write_text("installed=true", encoding="utf-8")

    context = detect_windows_install_context(
        executable_path=str(exe_path),
        platform_name="windows",
    )
    assert context.is_windows is True
    assert context.installed is True
    assert choose_windows_asset_flavor(context) == "installer"


def test_detect_windows_install_context_portable(tmp_path) -> None:
    exe_path = tmp_path / "Chemuson.exe"
    exe_path.write_text("stub", encoding="utf-8")

    context = detect_windows_install_context(
        executable_path=str(exe_path),
        platform_name="windows",
    )
    assert context.is_windows is True
    assert context.installed is False
    assert choose_windows_asset_flavor(context) == "portable"


def test_build_inno_installer_command_flags() -> None:
    silent_command = build_inno_installer_command(r"C:\tmp\Chemuson-setup.exe", silent=True)
    normal_command = build_inno_installer_command(r"C:\tmp\Chemuson-setup.exe", silent=False)

    assert "/VERYSILENT" in silent_command
    assert "/SUPPRESSMSGBOXES" in silent_command
    assert "/SILENT" in normal_command
    assert "/CLOSEAPPLICATIONS" in normal_command


def test_provider_prefers_requested_windows_flavor() -> None:
    provider = GitHubReleasesProvider("PJGV333", "Chemuson")
    release = _windows_release()

    installer = provider.find_platform_asset(
        release,
        platform_tag="windows-x86_64",
        preferred_flavor="installer",
    )
    portable = provider.find_platform_asset(
        release,
        platform_tag="windows-x86_64",
        preferred_flavor="portable",
    )

    assert installer is not None
    assert portable is not None
    assert "setup" in installer.name.lower()
    assert "portable" in portable.name.lower()


class _ProviderCapturingFlavor:
    def __init__(self, release: ReleaseInfo):
        self.release = release
        self.last_flavor = None

    def latest_for_channel(self, _channel):
        return self.release

    def find_platform_asset(self, release, platform_tag=None, preferred_flavor=None):
        _ = platform_tag
        self.last_flavor = preferred_flavor
        return release.assets[0]


def test_core_passes_preferred_flavor_to_provider() -> None:
    release = _windows_release("1.6.0")
    provider = _ProviderCapturingFlavor(release)
    settings = UpdateSettings(
        enabled=True,
        channel=UpdateChannel.STABLE,
        mode=UpdateMode.NOTIFY,
        check_interval_hours=24,
    )
    core = AutoUpdateCore(provider=provider, settings=settings)

    result = core.check_for_updates(
        current_version="1.5.0",
        platform_tag="windows-x86_64",
        preferred_asset_flavor="installer",
    )

    assert result.available is True
    assert provider.last_flavor == "installer"
