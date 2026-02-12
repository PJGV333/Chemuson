import os
import sys
from pathlib import Path

def _set_env_if_empty(key: str, value: str) -> None:
    if not os.environ.get(key):
        os.environ[key] = value

if getattr(sys, "frozen", False) and hasattr(sys, "_MEIPASS"):
    base = Path(sys._MEIPASS)

    plugins = base / "PyQt6" / "Qt6" / "plugins"
    platforms = plugins / "platforms"

    if plugins.exists():
        _set_env_if_empty("QT_PLUGIN_PATH", str(plugins))
    if platforms.exists():
        _set_env_if_empty("QT_QPA_PLATFORM_PLUGIN_PATH", str(platforms))

    wayland_plugin = platforms / "libqwayland-egl.so"
    xcb_plugin = platforms / "libqxcb.so"
    if xcb_plugin.exists() and not wayland_plugin.exists():
        _set_env_if_empty("QT_QPA_PLATFORM", "xcb")
