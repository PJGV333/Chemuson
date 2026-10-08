"""Opt-in packaged-icon smoke used by the Windows/Linux build pipelines.

This diagnostic is intentionally inert during normal application startup. It
runs only when the hidden ``--icon-smoke-test`` CLI option is combined with
``CHEMUSON_ICON_SMOKE_TEST=1``.
"""

from __future__ import annotations

import sys
from pathlib import Path
from typing import Any

from PyQt6 import QtSvg
from PyQt6.QtCore import QSize, Qt
from PyQt6.QtGui import QIcon, QImage, QPainter
from PyQt6.QtWidgets import QApplication, QToolButton

from chemuson.gui import icons as icon_facade
from chemuson.gui.theme.icon_provider import DEFAULT_ICONS_DIR, IconProvider

EXPECTED_STATIC_SVG_COUNT = 69
SMOKE_ICON_NAMES = (
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
)


def _raster_signature(icon: Any, size: int = 32) -> tuple[frozenset[int], int]:
    pixmap = icon.pixmap(QSize(size, size))
    if pixmap.isNull():
        return frozenset(), 0
    image = pixmap.toImage()
    visible: list[int] = []
    for y in range(image.height()):
        for x in range(image.width()):
            pixel = image.pixelColor(x, y)
            if pixel.alpha() > 0:
                visible.append(pixel.rgba())
    return frozenset(visible), len(visible)


def _tool_button_image(icon: QIcon, app: QApplication) -> QImage:
    button = QToolButton()
    button.setAutoRaise(True)
    button.setToolButtonStyle(Qt.ToolButtonStyle.ToolButtonIconOnly)
    button.setIconSize(QSize(32, 32))
    button.resize(48, 48)
    button.setIcon(icon)
    button.show()
    app.processEvents()

    image = QImage(button.size(), QImage.Format.Format_ARGB32_Premultiplied)
    image.fill(Qt.GlobalColor.transparent)
    painter = QPainter(image)
    button.render(painter)
    painter.end()
    button.close()
    button.deleteLater()
    app.processEvents()
    return image


def _tool_button_ink_difference(image: QImage, baseline: QImage) -> int:
    """Count changed visible pixels in the icon area of a QToolButton."""
    changed = 0
    for y in range(8, min(40, image.height(), baseline.height())):
        for x in range(8, min(40, image.width(), baseline.width())):
            pixel = image.pixelColor(x, y)
            if pixel.alpha() > 0 and pixel.rgba() != baseline.pixelColor(x, y).rgba():
                changed += 1
    return changed


def _runtime_paths() -> dict[str, str | None]:
    icon_directory = DEFAULT_ICONS_DIR.resolve()
    provider_file = Path(
        sys.modules["chemuson.gui.theme.icon_provider"].__file__
    ).resolve()
    frozen = bool(getattr(sys, "frozen", False))
    meipass: Path | None = None
    if frozen:
        raw_meipass = getattr(sys, "_MEIPASS", None)
        if not raw_meipass:
            raise RuntimeError("Frozen process has no sys._MEIPASS directory.")
        meipass = Path(raw_meipass).resolve()
        expected_directory = meipass / "chemuson/gui/theme/icons"
        if icon_directory != expected_directory:
            raise RuntimeError(
                "IconProvider path does not resolve under sys._MEIPASS: "
                f"expected {expected_directory}, got {icon_directory}."
            )
        if not provider_file.is_relative_to(meipass):
            raise RuntimeError(
                f"IconProvider.__file__ is outside sys._MEIPASS: {provider_file}."
            )
        qt_svg_file = Path(QtSvg.__file__).resolve()
        if not qt_svg_file.is_relative_to(meipass):
            raise RuntimeError(f"PyQt6.QtSvg is not loaded from sys._MEIPASS: {qt_svg_file}.")

    return {
        "icons_dir": str(icon_directory),
        "icon_provider_file": str(provider_file),
        "meipass": str(meipass) if meipass is not None else None,
        "qt_svg_module": str(Path(QtSvg.__file__).resolve()),
        "frozen": str(frozen).lower(),
    }


def run_icon_smoke_test() -> dict[str, object]:
    """Check packaged paths, SVG inventory, QtSvg raster output and themes."""
    app = QApplication.instance()
    if app is None:
        app = QApplication(["chemuson-icon-smoke"])

    paths = _runtime_paths()
    icon_directory = Path(paths["icons_dir"] or "")
    svg_files = sorted(icon_directory.glob("i-*.svg"))
    if len(svg_files) != EXPECTED_STATIC_SVG_COUNT:
        raise RuntimeError(
            f"Expected {EXPECTED_STATIC_SVG_COUNT} packaged SVG icons, found {len(svg_files)} "
            f"in {icon_directory}."
        )
    if any(not path.is_file() or path.stat().st_size == 0 for path in svg_files):
        raise RuntimeError("At least one packaged SVG icon is missing or empty.")

    provider = IconProvider(dpr=1.0, icons_dir=icon_directory)
    theme_results: dict[str, dict[str, int]] = {}
    essential_signatures: dict[str, dict[str, frozenset[int]]] = {}

    try:
        for theme in ("light", "dark"):
            icon_facade.set_icon_theme(theme)
            tint = icon_facade.icon_foreground_color()
            static_counts: dict[str, int] = {}
            for path in svg_files:
                icon_name = path.stem.removeprefix("i-")
                signature, pixel_count = _raster_signature(
                    provider.icon(icon_name, tint, 32)
                )
                if not signature or pixel_count == 0:
                    raise RuntimeError(
                        f"QtSvg rendered no visible pixels for {path.name} in {theme} theme."
                    )
                static_counts[icon_name] = pixel_count

            essential_icons = {
                "pointer/select": icon_facade.draw_generic_icon("pointer"),
                "single bond": icon_facade.draw_bond_icon("single"),
                "aromatic bond": icon_facade.draw_bond_icon("aromatic"),
                "aromatic ring": icon_facade.draw_ring_icon(6, aromatic=True),
                "search": icon_facade._static("search"),
                "undo": icon_facade.draw_generic_icon("undo"),
                "redo": icon_facade.draw_generic_icon("redo"),
                "new document": icon_facade.draw_generic_icon("document_new"),
                "clean": icon_facade.draw_generic_icon("clean"),
                "carbon atom": icon_facade.draw_atom_icon("C"),
                "coordination sphere": icon_facade.draw_coordination_sphere_icon(),
            }
            essential_counts: dict[str, int] = {}
            toolbar_counts: dict[str, int] = {}
            essential_signatures[theme] = {}
            blank_button = _tool_button_image(QIcon(), app)
            for name, icon in essential_icons.items():
                signature, pixel_count = _raster_signature(icon)
                if icon.isNull() or not signature or pixel_count == 0:
                    raise RuntimeError(
                        f"Essential icon {name!r} rendered blank in {theme} theme."
                    )
                essential_counts[name] = pixel_count
                essential_signatures[theme][name] = signature
                button_image = _tool_button_image(icon, app)
                changed_pixels = _tool_button_ink_difference(button_image, blank_button)
                if changed_pixels == 0:
                    raise RuntimeError(
                        f"QToolButton displayed no visible icon pixels for {name!r} in {theme} theme."
                    )
                toolbar_counts[name] = changed_pixels
            theme_results[theme] = {
                "static_svg_count": len(static_counts),
                "static_icons_with_visible_pixels": sum(
                    count > 0 for count in static_counts.values()
                ),
                **{f"essential:{name}": count for name, count in essential_counts.items()},
                **{f"toolbar:{name}": count for name, count in toolbar_counts.items()},
            }

        for name in (
            "pointer/select",
            "single bond",
            "aromatic ring",
            "search",
            "undo",
            "redo",
            "new document",
            "clean",
        ):
            if essential_signatures["light"][name] == essential_signatures["dark"][name]:
                raise RuntimeError(f"Theme tint did not change the {name!r} icon raster.")
        for name in ("carbon atom", "coordination sphere"):
            if essential_signatures["light"][name] != essential_signatures["dark"][name]:
                raise RuntimeError(f"Domain-colored {name!r} changed unexpectedly with theme.")
    finally:
        icon_facade.set_icon_theme("light")

    hidpi_provider = IconProvider(dpr=2.0, icons_dir=icon_directory)
    hidpi_pixmap = hidpi_provider.pixmap("pointer", IconProvider.theme_color("light"), 32)
    _hidpi_signature, hidpi_pixel_count = _raster_signature(
        hidpi_provider.icon("pointer", IconProvider.theme_color("light"), 32), 64
    )
    if (
        hidpi_pixmap.isNull()
        or hidpi_pixmap.width() != 64
        or hidpi_pixmap.height() != 64
        or hidpi_pixmap.devicePixelRatio() != 2.0
        or hidpi_pixel_count == 0
    ):
        raise RuntimeError("The HiDPI SVG raster smoke failed at DPR 2.")

    return {
        "status": "success",
        **paths,
        "static_svg_count": len(svg_files),
        "themes": theme_results,
        "essential_icons": list(SMOKE_ICON_NAMES),
        "hidpi_dpr": 2.0,
        "hidpi_visible_pixels": hidpi_pixel_count,
        "cwd": str(Path.cwd().resolve()),
    }
