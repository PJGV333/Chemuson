"""Pruebas de helpers para actualización portable/AppImage."""

from __future__ import annotations

import os
import stat
import sys

sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), "..", "src")))

from chemuson.update.portable import (
    APPIMAGE_ENV_VAR,
    build_portable_update_launcher_command,
    detect_portable_update_context,
    is_portable_target_writable,
    write_portable_update_script,
)
from chemuson.update.windows import detect_windows_install_context


def test_detect_portable_update_context_prefers_appimage_env(tmp_path) -> None:
    mounted_exe = tmp_path / "squashfs-root" / "AppRun"
    mounted_exe.parent.mkdir(parents=True, exist_ok=True)
    mounted_exe.write_text("stub", encoding="utf-8")
    appimage = tmp_path / "Chemuson.AppImage"
    appimage.write_text("appimage", encoding="utf-8")

    context = detect_portable_update_context(
        executable_path=str(mounted_exe),
        platform_name="linux",
        env={APPIMAGE_ENV_VAR: str(appimage)},
        is_frozen=True,
    )

    assert context.can_self_update is True
    assert context.is_appimage is True
    assert context.target_path == str(appimage.resolve())
    assert context.display_name == "AppImage actual"


def test_detect_portable_update_context_for_windows_portable(tmp_path) -> None:
    exe_path = tmp_path / "Chemuson.exe"
    exe_path.write_text("stub", encoding="utf-8")
    windows_context = detect_windows_install_context(
        executable_path=str(exe_path),
        platform_name="windows",
    )

    context = detect_portable_update_context(
        executable_path=str(exe_path),
        platform_name="windows",
        windows_context=windows_context,
        is_frozen=True,
    )

    assert context.can_self_update is True
    assert context.is_windows is True
    assert context.is_appimage is False
    assert context.target_path == str(exe_path.resolve())


def test_is_portable_target_writable_checks_parent_directory(tmp_path) -> None:
    target = tmp_path / "Chemuson.AppImage"
    target.write_text("stub", encoding="utf-8")

    assert is_portable_target_writable(str(target)) is True


def test_write_portable_update_script_for_unix(tmp_path) -> None:
    script_path = tmp_path / "apply-portable-update.sh"
    source_path = tmp_path / "Chemuson-v1.2.3.AppImage"
    target_path = tmp_path / "Chemuson.AppImage"
    backup_path = tmp_path / "rollback" / "Chemuson.AppImage.bak"
    log_path = tmp_path / "portable-update.log"

    write_portable_update_script(
        str(script_path),
        source_path=str(source_path),
        target_path=str(target_path),
        backup_path=str(backup_path),
        log_path=str(log_path),
        pid=1234,
        relaunch=True,
        platform_name="linux",
    )

    content = script_path.read_text(encoding="utf-8")
    assert 'nohup "$TARGET"' in content
    assert 'mv -f "$SOURCE" "$TARGET"' in content
    assert build_portable_update_launcher_command(str(script_path), platform_name="linux") == [
        "sh",
        str(script_path.resolve()),
    ]
    assert bool(script_path.stat().st_mode & stat.S_IXUSR)


def test_write_portable_update_script_for_windows_escapes_batch_values(tmp_path) -> None:
    script_path = tmp_path / "apply-portable-update.cmd"
    source_path = tmp_path / "Chemuson 100%.exe"
    target_path = tmp_path / "Current Chemuson.exe"
    backup_path = tmp_path / "rollback" / "Chemuson.exe.bak"
    log_path = tmp_path / "portable-update.log"

    write_portable_update_script(
        str(script_path),
        source_path=str(source_path),
        target_path=str(target_path),
        backup_path=str(backup_path),
        log_path=str(log_path),
        pid=4321,
        relaunch=False,
        platform_name="windows",
    )

    content = script_path.read_text(encoding="utf-8")
    assert "Chemuson 100%%.exe" in content
    assert 'move /Y "%SOURCE%" "%TARGET%"' in content
    assert build_portable_update_launcher_command(str(script_path), platform_name="windows") == [
        "cmd",
        "/c",
        str(script_path.resolve()),
    ]
