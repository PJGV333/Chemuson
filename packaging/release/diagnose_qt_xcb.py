"""Fail early with Qt xcb linker diagnostics before packaged GUI smoke tests."""

from __future__ import annotations

import importlib.util
import os
import re
import shutil
import subprocess
import sys
import tempfile
from pathlib import Path


def _run_and_print(args: list[str], *, env: dict[str, str] | None = None) -> subprocess.CompletedProcess[str]:
    print("$ " + " ".join(args), flush=True)
    result = subprocess.run(
        args,
        check=False,
        capture_output=True,
        text=True,
        env=env,
        timeout=60,
    )
    if result.stdout:
        print(result.stdout, end="", flush=True)
    if result.stderr:
        print(result.stderr, end="", file=sys.stderr, flush=True)
    if result.returncode:
        raise RuntimeError(f"Command exited {result.returncode}: {args[0]}")
    return result


def diagnose_qt_xcb() -> None:
    package = importlib.util.find_spec("PyQt6")
    if package is None or package.origin is None:
        raise RuntimeError("PyQt6 is not installed in the current Python environment.")
    qt_root = Path(package.origin).resolve().parent / "Qt6"
    libraries = (
        qt_root / "lib/libQt6XcbQpa.so.6",
        qt_root / "plugins/platforms/libqxcb.so",
    )
    readelf = shutil.which("readelf")
    ldd = shutil.which("ldd")
    xvfb_run = shutil.which("xvfb-run")
    if not readelf or not ldd or not xvfb_run:
        raise RuntimeError("readelf, ldd, and xvfb-run are required for the Qt xcb diagnostic.")

    for library in libraries:
        if not library.is_file():
            raise RuntimeError(f"Qt xcb runtime library is missing: {library}")
        print(f"\n=== Qt xcb dependency inspection: {library} ===", flush=True)
        _run_and_print([readelf, "-d", str(library)])
        resolved = _run_and_print([ldd, "-r", str(library)])
        if re.search(r"=>\s+not found|undefined symbol", resolved.stdout + resolved.stderr):
            raise RuntimeError(f"Qt xcb dependencies are unresolved for {library}")

    smoke_env = os.environ.copy()
    for key in (
        "APPIMAGE",
        "APPDIR",
        "LD_PRELOAD",
        "LD_LIBRARY_PATH",
        "QT_PLUGIN_PATH",
        "QT_QPA_PLATFORM_PLUGIN_PATH",
    ):
        smoke_env.pop(key, None)
    smoke_env.update(
        {
            "QT_QPA_PLATFORM": "xcb",
            "QT_DEBUG_PLUGINS": "1",
        }
    )
    print("\n=== Qt xcb plugin initialization under Xvfb ===", flush=True)
    with tempfile.TemporaryDirectory(prefix="chemuson-xcb-loader-") as temp:
        trace_prefix = Path(temp) / "ld-debug"
        result = _run_and_print(
            [
                xvfb_run,
                "-a",
                "-s",
                "-screen 0 1280x900x24 +extension GLX +render -noreset",
                "env",
                "LD_DEBUG=libs",
                f"LD_DEBUG_OUTPUT={trace_prefix}",
                sys.executable,
                "-c",
                (
                    "from PyQt6.QtWidgets import QApplication; "
                    "app = QApplication([]); "
                    "print('platform=' + app.platformName()); "
                    "app.quit()"
                ),
            ],
            env=smoke_env,
        )
        if "platform=xcb" not in result.stdout:
            raise RuntimeError("Qt did not initialize the xcb platform plugin under Xvfb.")
        loader_lines = [
            line.rstrip()
            for trace in Path(temp).glob("ld-debug.*")
            for line in trace.read_text(encoding="utf-8", errors="replace").splitlines()
            if any(
                library in line.lower()
                for library in (
                    "libqxcb.so",
                    "libqt6xcbqpa.so",
                    "libxcb-keysyms.so",
                    "libxcb-icccm.so",
                    "libxcb-shape.so",
                )
            )
        ]
        if not loader_lines:
            raise RuntimeError("LD_DEBUG=libs produced no trace for the Qt xcb libraries.")
        print("\n=== LD_DEBUG=libs Qt xcb resolution trace ===", flush=True)
        print("\n".join(loader_lines), flush=True)
    print("Qt xcb loader diagnostic passed; packaged graphics smoke remains required.")


def main() -> int:
    try:
        diagnose_qt_xcb()
    except (OSError, RuntimeError, subprocess.SubprocessError) as exc:
        print(f"Qt xcb loader diagnostic failed: {exc}", file=sys.stderr)
        return 1
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
