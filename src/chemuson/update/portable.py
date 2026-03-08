"""Helpers de actualización para AppImage y binarios portables."""

from __future__ import annotations

import os
import platform
import shlex
import stat
import subprocess
import sys
from dataclasses import dataclass
from pathlib import Path
from typing import Mapping

from chemuson.update.windows import WindowsInstallContext, detect_windows_install_context

APPIMAGE_ENV_VAR = "APPIMAGE"


@dataclass(slots=True)
class PortableUpdateContext:
    """Describe si la app actual puede auto-reemplazarse como binario portable."""

    is_portable: bool
    is_windows: bool
    is_appimage: bool
    target_path: str
    executable_path: str
    display_name: str

    @property
    def can_self_update(self) -> bool:
        """Indica si hay una ruta objetivo concreta para auto-update."""
        return bool(self.is_portable and self.target_path)


def detect_portable_update_context(
    executable_path: str | None = None,
    platform_name: str | None = None,
    env: Mapping[str, str] | None = None,
    *,
    is_frozen: bool | None = None,
    windows_context: WindowsInstallContext | None = None,
) -> PortableUpdateContext:
    """Detecta si la app corre como AppImage o binario portable auto-actualizable."""
    sys_name = str(platform_name or platform.system() or "").strip().lower()
    is_windows = sys_name.startswith("win")
    env_map = env if env is not None else os.environ

    resolved_exe = Path(executable_path or sys.executable).resolve()
    appimage_raw = str(env_map.get(APPIMAGE_ENV_VAR, "") or "").strip()
    if appimage_raw:
        target = Path(appimage_raw).expanduser()
        if target.exists():
            target = target.resolve()
        return PortableUpdateContext(
            is_portable=True,
            is_windows=is_windows,
            is_appimage=True,
            target_path=str(target),
            executable_path=str(resolved_exe),
            display_name="AppImage actual",
        )

    if str(resolved_exe).lower().endswith(".appimage"):
        return PortableUpdateContext(
            is_portable=True,
            is_windows=is_windows,
            is_appimage=True,
            target_path=str(resolved_exe),
            executable_path=str(resolved_exe),
            display_name="AppImage actual",
        )

    if windows_context is None:
        windows_context = detect_windows_install_context(
            executable_path=str(resolved_exe),
            platform_name=sys_name,
        )
    if windows_context.is_windows and not windows_context.installed:
        return PortableUpdateContext(
            is_portable=True,
            is_windows=True,
            is_appimage=False,
            target_path=str(windows_context.executable_path),
            executable_path=str(windows_context.executable_path),
            display_name="ejecutable portable actual",
        )

    frozen = bool(getattr(sys, "frozen", False) if is_frozen is None else is_frozen)
    if frozen and not is_windows:
        return PortableUpdateContext(
            is_portable=True,
            is_windows=False,
            is_appimage=False,
            target_path=str(resolved_exe),
            executable_path=str(resolved_exe),
            display_name="ejecutable portable actual",
        )

    return PortableUpdateContext(
        is_portable=False,
        is_windows=is_windows,
        is_appimage=False,
        target_path="",
        executable_path=str(resolved_exe),
        display_name="",
    )


def is_portable_target_writable(target_path: str) -> bool:
    """Valida si la ruta del binario actual parece escribible."""
    target = Path(str(target_path or "").strip()).expanduser()
    if not str(target):
        return False
    parent = target.parent
    if not parent.exists() or not os.access(parent, os.W_OK):
        return False
    if target.exists() and not os.access(target, os.W_OK):
        return False
    return True


def _escape_windows_batch_value(value: str) -> str:
    return str(value or "").replace("%", "%%")


def _build_windows_replace_script(
    target_path: str,
    source_path: str,
    backup_path: str,
    log_path: str,
    pid: int,
    relaunch: bool,
) -> str:
    target = _escape_windows_batch_value(target_path)
    source = _escape_windows_batch_value(source_path)
    backup = _escape_windows_batch_value(backup_path)
    log = _escape_windows_batch_value(log_path)
    relaunch_flag = "1" if relaunch else "0"
    return "\r\n".join(
        [
            "@echo off",
            "setlocal enableextensions",
            f'set "TARGET={target}"',
            f'set "SOURCE={source}"',
            f'set "BACKUP={backup}"',
            f'set "LOG={log}"',
            f'set "PID={int(pid)}"',
            f'set "RELAUNCH={relaunch_flag}"',
            '>>"%LOG%" echo [%DATE% %TIME%] Starting portable update helper.',
            "for /L %%I in (1,1,300) do (",
            '  tasklist /FI "PID eq %PID%" 2>NUL | find "%PID%" >NUL',
            "  if errorlevel 1 goto after_wait",
            "  timeout /T 1 /NOBREAK >NUL",
            ")",
            '>>"%LOG%" echo Timed out waiting for Chemuson to exit.',
            "exit /b 1",
            ":after_wait",
            'if exist "%TARGET%" copy /Y "%TARGET%" "%BACKUP%" >NUL',
            'move /Y "%SOURCE%" "%TARGET%" >NUL',
            "if errorlevel 1 (",
            '  >>"%LOG%" echo Failed to replace target binary.',
            '  if exist "%BACKUP%" copy /Y "%BACKUP%" "%TARGET%" >NUL',
            "  exit /b 1",
            ")",
            'if "%RELAUNCH%"=="1" start "" "%TARGET%"',
            'del "%~f0" >NUL 2>NUL',
            "exit /b 0",
            "",
        ]
    )


def _build_unix_replace_script(
    target_path: str,
    source_path: str,
    backup_path: str,
    log_path: str,
    pid: int,
    relaunch: bool,
) -> str:
    target = shlex.quote(str(target_path))
    source = shlex.quote(str(source_path))
    backup = shlex.quote(str(backup_path))
    log = shlex.quote(str(log_path))
    relaunch_flag = "1" if relaunch else "0"
    return "\n".join(
        [
            "#!/usr/bin/env sh",
            "set -eu",
            f"TARGET={target}",
            f"SOURCE={source}",
            f"BACKUP={backup}",
            f"LOG={log}",
            f"PID={int(pid)}",
            f"RELAUNCH={relaunch_flag}",
            'mkdir -p "$(dirname "$BACKUP")"',
            'touch "$LOG"',
            'exec >>"$LOG" 2>&1',
            'echo "[$(date -u +%Y-%m-%dT%H:%M:%SZ)] Starting portable update helper."',
            "COUNT=0",
            'while kill -0 "$PID" 2>/dev/null; do',
            '  COUNT=$((COUNT + 1))',
            '  if [ "$COUNT" -ge 300 ]; then',
            '    echo "Timed out waiting for Chemuson to exit."',
            "    exit 1",
            "  fi",
            "  sleep 1",
            "done",
            'if [ -f "$TARGET" ]; then',
            '  cp -p "$TARGET" "$BACKUP"',
            "fi",
            'chmod +x "$SOURCE" || true',
            'if ! mv -f "$SOURCE" "$TARGET"; then',
            '  echo "Failed to replace target binary."',
            '  if [ -f "$BACKUP" ]; then',
            '    cp -f "$BACKUP" "$TARGET" || true',
            "  fi",
            "  exit 1",
            "fi",
            'chmod +x "$TARGET" || true',
            'if [ "$RELAUNCH" = "1" ]; then',
            '  nohup "$TARGET" >/dev/null 2>&1 &',
            "fi",
            'rm -f -- "$0"',
            "",
        ]
    )


def write_portable_update_script(
    script_path: str,
    *,
    source_path: str,
    target_path: str,
    backup_path: str,
    log_path: str,
    pid: int,
    relaunch: bool = False,
    platform_name: str | None = None,
) -> str:
    """Escribe script helper que reemplaza el binario luego de cerrar la app."""
    sys_name = str(platform_name or platform.system() or "").strip().lower()
    is_windows = sys_name.startswith("win")
    path = Path(script_path).resolve()
    path.parent.mkdir(parents=True, exist_ok=True)
    script_text = (
        _build_windows_replace_script(
            target_path=target_path,
            source_path=source_path,
            backup_path=backup_path,
            log_path=log_path,
            pid=pid,
            relaunch=relaunch,
        )
        if is_windows
        else _build_unix_replace_script(
            target_path=target_path,
            source_path=source_path,
            backup_path=backup_path,
            log_path=log_path,
            pid=pid,
            relaunch=relaunch,
        )
    )
    path.write_text(script_text, encoding="utf-8", newline="" if is_windows else "\n")
    if not is_windows:
        current_mode = path.stat().st_mode
        path.chmod(current_mode | stat.S_IXUSR | stat.S_IWUSR | stat.S_IRUSR)
    return str(path)


def build_portable_update_launcher_command(
    script_path: str,
    platform_name: str | None = None,
) -> list[str]:
    """Construye comando para lanzar helper de reemplazo desacoplado."""
    sys_name = str(platform_name or platform.system() or "").strip().lower()
    resolved = str(Path(script_path).resolve())
    if sys_name.startswith("win"):
        return ["cmd", "/c", resolved]
    return ["sh", resolved]


def launch_portable_update_script(
    script_path: str,
    platform_name: str | None = None,
) -> list[str]:
    """Lanza helper desacoplado y devuelve el comando usado."""
    command = build_portable_update_launcher_command(script_path, platform_name=platform_name)
    sys_name = str(platform_name or platform.system() or "").strip().lower()
    if sys_name.startswith("win"):
        creationflags = 0
        if hasattr(subprocess, "DETACHED_PROCESS"):
            creationflags |= int(subprocess.DETACHED_PROCESS)
        if hasattr(subprocess, "CREATE_NEW_PROCESS_GROUP"):
            creationflags |= int(subprocess.CREATE_NEW_PROCESS_GROUP)
        subprocess.Popen(command, close_fds=True, creationflags=creationflags)
    else:
        subprocess.Popen(command, close_fds=True, start_new_session=True)
    return command
