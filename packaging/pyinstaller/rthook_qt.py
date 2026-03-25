"""Runtime hook para asegurar descubrimiento de plugins Qt en binarios congelados."""

from __future__ import annotations

import os
import sys
from pathlib import Path


def _is_frozen_app() -> bool:
    return bool(getattr(sys, "frozen", False))


def _first_existing(paths: list[Path]) -> Path | None:
    for candidate in paths:
        if candidate.exists():
            return candidate
    return None


def _has_wayland_plugin(platforms_dir: Path) -> bool:
    return any(p.name.startswith(("libqwayland", "qwayland")) for p in platforms_dir.glob("*"))


def _has_xcb_plugin(platforms_dir: Path) -> bool:
    return any("qxcb" in p.name for p in platforms_dir.glob("*"))


def _configure_ssl_cert_file() -> None:
    try:
        import certifi
    except Exception:
        return
    try:
        os.environ.setdefault("SSL_CERT_FILE", certifi.where())
    except Exception:
        pass


def _configure_qt_runtime() -> None:
    if not _is_frozen_app():
        return

    _configure_ssl_cert_file()

    base_dir = Path(getattr(sys, "_MEIPASS", Path(sys.executable).resolve().parent))
    plugin_root = _first_existing(
        [
            base_dir / "PyQt6" / "Qt6" / "plugins",
            base_dir / "PyQt6" / "Qt" / "plugins",
            base_dir / "qt6_plugins",
            base_dir / "qt_plugins",
            base_dir / "plugins",
        ]
    )
    if plugin_root is None:
        return

    os.environ.setdefault("QT_PLUGIN_PATH", str(plugin_root))
    platforms_dir = plugin_root / "platforms"
    if platforms_dir.exists():
        os.environ.setdefault("QT_QPA_PLATFORM_PLUGIN_PATH", str(platforms_dir))

        # En Linux, fuerza xcb cuando falta plugin Wayland para evitar fallo de arranque.
        if sys.platform.startswith("linux"):
            if not _has_wayland_plugin(platforms_dir) and _has_xcb_plugin(platforms_dir):
                os.environ.setdefault("QT_QPA_PLATFORM", "xcb")


_configure_qt_runtime()
