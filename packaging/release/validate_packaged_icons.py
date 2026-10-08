"""Run the opt-in SVG icon diagnostic inside a built PyInstaller executable."""

from __future__ import annotations

import argparse
import json
import os
import subprocess
import tempfile
from pathlib import Path
from typing import Any

EXPECTED_SVG_COUNT = 69
EXPECTED_ESSENTIAL_ICONS = {
    "pointer/select",
    "single bond",
    "aromatic bond",
    "aromatic ring",
    "search",
    "undo",
    "redo",
    "new document",
    "clean",
    "carbon atom",
    "coordination sphere",
}


def validate_executable(executable: Path, *, timeout: int = 120) -> dict[str, Any]:
    """Require the actual frozen executable to find and rasterize its SVGs."""
    binary = executable.resolve()
    if not binary.is_file() or binary.stat().st_size == 0:
        raise ValueError(f"Frozen executable is missing or empty: {binary}")

    with tempfile.TemporaryDirectory(prefix="chemuson-icons-smoke-") as temporary:
        isolated_cwd = Path(temporary).resolve()
        environment = os.environ.copy()
        environment.update(
            {
                "CHEMUSON_ICON_SMOKE_TEST": "1",
                "QT_QPA_PLATFORM": "offscreen",
                "HOME": str(isolated_cwd / "home"),
                "XDG_CONFIG_HOME": str(isolated_cwd / "config"),
                "XDG_DATA_HOME": str(isolated_cwd / "data"),
                "XDG_CACHE_HOME": str(isolated_cwd / "cache"),
            }
        )
        for variable in ("HOME", "XDG_CONFIG_HOME", "XDG_DATA_HOME", "XDG_CACHE_HOME"):
            Path(environment[variable]).mkdir(parents=True, exist_ok=True)
        try:
            result = subprocess.run(
                [str(binary), "--icon-smoke-test"],
                cwd=isolated_cwd,
                env=environment,
                check=False,
                capture_output=True,
                text=True,
                timeout=timeout,
            )
        except (OSError, subprocess.TimeoutExpired) as exc:
            raise ValueError(f"Could not run frozen icon smoke for {binary}: {exc}") from exc
        if result.returncode != 0:
            raise ValueError(
                f"Frozen icon smoke failed ({result.returncode}) for {binary}:\n"
                f"stdout:\n{result.stdout[-4000:]}\nstderr:\n{result.stderr[-4000:]}"
            )
        try:
            report = json.loads(result.stdout.strip())
        except json.JSONDecodeError as exc:
            raise ValueError(
                f"Frozen icon smoke did not return a JSON report: {result.stdout[-4000:]}"
            ) from exc

        if report.get("status") != "success" or report.get("frozen") != "true":
            raise ValueError("Icon diagnostic did not report a successful frozen process.")
        if report.get("static_svg_count") != EXPECTED_SVG_COUNT:
            raise ValueError(
                f"Frozen process reported {report.get('static_svg_count')} SVGs; "
                f"expected {EXPECTED_SVG_COUNT}."
            )
        if set(report.get("essential_icons", ())) != EXPECTED_ESSENTIAL_ICONS:
            raise ValueError("Frozen process did not verify the complete essential icon set.")
        if Path(report.get("cwd", "")).resolve() != isolated_cwd:
            raise ValueError("Frozen icon smoke did not run from its isolated working directory.")
        if Path(report.get("icons_dir", "")).resolve() == isolated_cwd:
            raise ValueError("Icon resource lookup unexpectedly depends on the working directory.")
        meipass = Path(report.get("meipass", "")).resolve()
        expected_icons_dir = meipass / "chemuson/gui/theme/icons"
        if Path(report.get("icons_dir", "")).resolve() != expected_icons_dir:
            raise ValueError("Frozen SVG path does not match the recorded sys._MEIPASS directory.")
        if not Path(report.get("icon_provider_file", "")).resolve().is_relative_to(meipass):
            raise ValueError("Frozen IconProvider.__file__ is outside sys._MEIPASS.")
        if not Path(report.get("qt_svg_module", "")).resolve().is_relative_to(meipass):
            raise ValueError("Frozen PyQt6.QtSvg module is outside sys._MEIPASS.")

        themes = report.get("themes", {})
        if set(themes) != {"light", "dark"}:
            raise ValueError("Frozen icon smoke did not validate both light and dark themes.")
        for theme, values in themes.items():
            if (
                values.get("static_svg_count") != EXPECTED_SVG_COUNT
                or values.get("static_icons_with_visible_pixels") != EXPECTED_SVG_COUNT
            ):
                raise ValueError(f"Not all 69 static SVGs rendered visible pixels in {theme} theme.")
            for name in EXPECTED_ESSENTIAL_ICONS:
                if values.get(f"essential:{name}", 0) <= 0:
                    raise ValueError(f"Essential {name!r} icon is blank in {theme} theme.")
                if values.get(f"toolbar:{name}", 0) <= 0:
                    raise ValueError(f"QToolButton displayed no {name!r} icon pixels in {theme} theme.")
        if report.get("hidpi_dpr") != 2.0 or report.get("hidpi_visible_pixels", 0) <= 0:
            raise ValueError("Frozen executable failed the DPR 2 HiDPI raster check.")

    return report


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--executable", required=True, type=Path)
    parser.add_argument("--timeout", type=int, default=120)
    args = parser.parse_args()
    report = validate_executable(args.executable, timeout=args.timeout)
    print(json.dumps(report, indent=2, sort_keys=True))


if __name__ == "__main__":
    try:
        main()
    except Exception as exc:
        raise SystemExit(f"Packaged icon validation failed: {exc}") from exc
