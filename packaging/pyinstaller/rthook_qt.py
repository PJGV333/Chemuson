import os
import sys
from pathlib import Path


def _set_env_if_empty(key: str, value: str) -> None:
    if not os.environ.get(key):
        os.environ[key] = value


def _first_existing(paths: list[Path]) -> Path | None:
    for path in paths:
        if path.exists():
            return path
    return None


if getattr(sys, "frozen", False):
    meipass = Path(getattr(sys, "_MEIPASS", Path(sys.executable).resolve().parent))
    base_candidates = [
        meipass,
        meipass / "_internal",
        Path(sys.executable).resolve().parent,
    ]

    plugin_root = None
    for base_dir in base_candidates:
        plugin_root = _first_existing(
            [
                base_dir / "PyQt6" / "Qt6" / "plugins",
                base_dir / "PyQt6" / "Qt" / "plugins",
                base_dir / "qt6_plugins",
                base_dir / "qt_plugins",
                base_dir / "plugins",
            ]
        )
        if plugin_root is not None:
            break

    if plugin_root is not None:
        platforms = plugin_root / "platforms"

        _set_env_if_empty("QT_PLUGIN_PATH", str(plugin_root))
        if platforms.exists():
            _set_env_if_empty("QT_QPA_PLATFORM_PLUGIN_PATH", str(platforms))

            wayland_exists = any(
                p.name.startswith(("libqwayland", "qwayland")) for p in platforms.glob("*")
            )
            xcb_exists = any("qxcb" in p.name for p in platforms.glob("*"))
            if xcb_exists and not wayland_exists:
                _set_env_if_empty("QT_QPA_PLATFORM", "xcb")
