"""HiDPI regressions for SVG -> QPixmap -> QIcon rendering."""
from __future__ import annotations

import pytest
from PyQt6.QtCore import QSize
from PyQt6.QtWidgets import QApplication

from chemuson.gui.theme.icon_provider import IconProvider


DPRS = (1.0, 1.25, 1.5, 2.0)
REPRESENTATIVE = ("pointer", "bond-single", "ring", "flask", "undo", "sliders")


def _ink_bounds(pixmap):
    image = pixmap.toImage()
    points = [
        (x, y)
        for y in range(image.height())
        for x in range(image.width())
        if image.pixelColor(x, y).alpha() > 16
    ]
    assert points, "icon pixmap must contain visible ink"
    return (
        min(x for x, _ in points), min(y for _, y in points),
        max(x for x, _ in points), max(y for _, y in points),
    )


@pytest.mark.parametrize("dpr", DPRS)
@pytest.mark.parametrize("name", REPRESENTATIVE)
def test_hidpi_svg_backing_store_is_physical_and_icon_is_not_clipped(
    _session_qapp: QApplication, dpr: float, name: str
) -> None:
    provider = IconProvider(dpr=dpr)
    logical_size = 21
    if name == "ring":
        pixmap = provider.pixmap_dynamic(
            "ring", "#334155", logical_size, variant="benzene"
        )
        icon = provider.icon_dynamic(
            "ring", "#334155", logical_size, variant="benzene"
        )
    else:
        pixmap = provider.pixmap(name, "#334155", logical_size)
        icon = provider.icon(name, "#334155", logical_size)

    physical_size = round(logical_size * dpr)
    assert not pixmap.isNull()
    assert (pixmap.width(), pixmap.height()) == (physical_size, physical_size)
    assert pixmap.devicePixelRatioF() == pytest.approx(dpr)
    assert pixmap.deviceIndependentSize().width() == pytest.approx(
        physical_size / dpr
    )
    assert pixmap.deviceIndependentSize().width() == pytest.approx(
        logical_size, abs=0.4
    )

    left, top, right, bottom = _ink_bounds(pixmap)
    # Los SVGs representativos ocupan un área significativa, pero el trazo
    # conserva margen interior (al menos 1 physical pixel) en todos los DPR.
    ink_w = right - left + 1
    ink_h = bottom - top + 1
    assert ink_w >= physical_size * 0.35
    assert ink_h >= physical_size * 0.35
    assert left >= 1 and top >= 1
    assert right <= physical_size - 2 and bottom <= physical_size - 2
    assert abs((left + right) / 2 - (physical_size - 1) / 2) <= physical_size * 0.2
    assert abs((top + bottom) / 2 - (physical_size - 1) / 2) <= physical_size * 0.2

    assert not icon.isNull()
    assert icon.actualSize(QSize(logical_size, logical_size)) == QSize(
        logical_size, logical_size
    )


def test_cache_identity_includes_dpr_for_pixmaps_and_icons(
    _session_qapp: QApplication,
) -> None:
    low = IconProvider(dpr=1.0)
    high = IconProvider(dpr=2.0)

    low_pm = low.pixmap("pointer", "#334155", 21)
    high_pm = high.pixmap("pointer", "#334155", 21)
    assert low_pm.size() != high_pm.size()
    assert low_pm.width() == 21
    assert high_pm.width() == 42
    assert low.icon("pointer", "#334155", 21) is not high.icon(
        "pointer", "#334155", 21
    )


def test_same_provider_can_change_dpr_without_stale_cache(
    _session_qapp: QApplication,
) -> None:
    provider = IconProvider(dpr=1.0)
    first = provider.pixmap("pointer", "#334155", 21)
    provider.set_device_pixel_ratio(2.0)
    second = provider.pixmap("pointer", "#334155", 21)
    assert first.width() == 21
    assert second.width() == 42
    assert second is not first
    assert second.devicePixelRatioF() == pytest.approx(2.0)
