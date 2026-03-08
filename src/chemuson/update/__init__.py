"""API pública del subsistema de actualización."""

from chemuson.update.core import AutoUpdateCore, VerificationResult
from chemuson.update.policy import (
    can_offer_update,
    is_downgrade_suspected,
    mark_checked,
    should_check_now,
)
from chemuson.update.portable import (
    APPIMAGE_ENV_VAR,
    PortableUpdateContext,
    build_portable_update_launcher_command,
    detect_portable_update_context,
    is_portable_target_writable,
    launch_portable_update_script,
    write_portable_update_script,
)
from chemuson.update.provider import GitHubReleasesProvider, detect_platform_tag
from chemuson.update.rollback import RollbackManager
from chemuson.update.security import SignatureVerifier
from chemuson.update.telemetry import UpdateTelemetryLogger
from chemuson.update.semver import (
    channel_accepts_version,
    compare_versions,
    is_newer_version,
    is_prerelease,
    parse_semver,
)
from chemuson.update.types import (
    DownloadedUpdate,
    ReleaseAsset,
    ReleaseInfo,
    UpdateCandidate,
    UpdateChannel,
    UpdateCheckResult,
    UpdateMode,
    UpdateSettings,
    coerce_update_channel,
    coerce_update_mode,
)
from chemuson.update.windows import (
    INSTALL_MARKER_FILENAME,
    WindowsInstallContext,
    build_inno_installer_command,
    choose_windows_asset_flavor,
    detect_windows_install_context,
    is_windows_installer_asset,
    launch_inno_installer,
)

__all__ = [
    "AutoUpdateCore",
    "VerificationResult",
    "GitHubReleasesProvider",
    "detect_platform_tag",
    "RollbackManager",
    "SignatureVerifier",
    "UpdateTelemetryLogger",
    "parse_semver",
    "compare_versions",
    "is_newer_version",
    "is_prerelease",
    "channel_accepts_version",
    "should_check_now",
    "mark_checked",
    "can_offer_update",
    "is_downgrade_suspected",
    "APPIMAGE_ENV_VAR",
    "PortableUpdateContext",
    "detect_portable_update_context",
    "is_portable_target_writable",
    "build_portable_update_launcher_command",
    "write_portable_update_script",
    "launch_portable_update_script",
    "DownloadedUpdate",
    "ReleaseAsset",
    "ReleaseInfo",
    "UpdateCandidate",
    "UpdateChannel",
    "UpdateCheckResult",
    "UpdateMode",
    "UpdateSettings",
    "coerce_update_channel",
    "coerce_update_mode",
    "INSTALL_MARKER_FILENAME",
    "WindowsInstallContext",
    "detect_windows_install_context",
    "choose_windows_asset_flavor",
    "is_windows_installer_asset",
    "build_inno_installer_command",
    "launch_inno_installer",
]
