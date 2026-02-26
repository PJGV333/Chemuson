"""Utilidades de actualización para instalaciones Windows."""

from __future__ import annotations

import os
import platform
import subprocess
import sys
from dataclasses import dataclass
from pathlib import Path

INSTALL_MARKER_FILENAME = ".chemuson-installed"


@dataclass(slots=True)
class WindowsInstallContext:
    """Describe si la app se ejecuta como instalación Windows o portable."""

    is_windows: bool
    installed: bool
    executable_path: str
    install_dir: str
    marker_path: str


def _is_truthy(value: str) -> bool:
    lowered = str(value or "").strip().lower()
    return lowered in {"1", "true", "yes", "on", "si", "sí"}


def detect_windows_install_context(
    executable_path: str | None = None,
    platform_name: str | None = None,
) -> WindowsInstallContext:
    """Detecta si Chemuson corre como instalación formal en Windows.

    Se considera "instalado" cuando existe marcador de instalación,
    cuando una variable de entorno fuerza el modo o cuando el ejecutable
    está ubicado bajo Program Files.
    """

    sys_name = str(platform_name or platform.system() or "").strip().lower()
    is_windows = sys_name.startswith("win")

    resolved_exe = Path(executable_path or sys.executable).resolve()
    install_dir = resolved_exe.parent
    marker_path = install_dir / INSTALL_MARKER_FILENAME

    forced = _is_truthy(os.getenv("CHEMUSON_FORCE_INSTALLED", ""))
    marker_exists = marker_path.exists()
    exe_lower = str(resolved_exe).replace("/", "\\").lower()
    in_program_files = "\\program files\\" in exe_lower or "\\program files (x86)\\" in exe_lower

    installed = bool(is_windows and (forced or marker_exists or in_program_files))
    return WindowsInstallContext(
        is_windows=is_windows,
        installed=installed,
        executable_path=str(resolved_exe),
        install_dir=str(install_dir),
        marker_path=str(marker_path),
    )


def choose_windows_asset_flavor(context: WindowsInstallContext) -> str:
    """Devuelve flavor de asset preferido para Windows.

    - `installer` cuando la app está instalada.
    - `portable` cuando se ejecuta en modo portable.
    """

    if not context.is_windows:
        return ""
    return "installer" if context.installed else "portable"


def is_windows_installer_asset(asset_name: str) -> bool:
    """Indica si un asset corresponde a instalador Windows."""

    lowered = str(asset_name or "").strip().lower()
    return lowered.endswith(".msi") or "setup" in lowered or "installer" in lowered


def build_inno_installer_command(installer_path: str, silent: bool = True) -> list[str]:
    """Construye comando seguro para ejecutar update vía Inno Setup."""

    command = [
        str(installer_path),
        "/SP-",
        "/NORESTART",
        "/CLOSEAPPLICATIONS",
        "/FORCECLOSEAPPLICATIONS",
    ]
    if silent:
        command.extend(["/VERYSILENT", "/SUPPRESSMSGBOXES"])
    else:
        command.append("/SILENT")
    return command


def launch_inno_installer(installer_path: str, silent: bool = True) -> list[str]:
    """Lanza instalador Inno en proceso desacoplado y devuelve comando usado."""

    path = Path(installer_path).resolve()
    if not path.exists():
        raise FileNotFoundError(f"No existe instalador: {path}")

    command = build_inno_installer_command(str(path), silent=silent)
    creationflags = 0
    if hasattr(subprocess, "DETACHED_PROCESS"):
        creationflags |= int(subprocess.DETACHED_PROCESS)
    if hasattr(subprocess, "CREATE_NEW_PROCESS_GROUP"):
        creationflags |= int(subprocess.CREATE_NEW_PROCESS_GROUP)

    # Se usa proceso desacoplado para permitir que la app actual cierre sin bloquear.
    subprocess.Popen(command, close_fds=True, creationflags=creationflags)
    return command
