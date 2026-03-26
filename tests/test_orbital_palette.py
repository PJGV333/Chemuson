"""Pruebas funcionales y visuales para la paleta de orbitales."""

from __future__ import annotations

import os
import sys
from pathlib import Path

import pytest
from PyQt6.QtCore import QPointF
from PyQt6.QtGui import QImage
from PyQt6.QtWidgets import QApplication

sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), "..", "src")))

from chemuson.gui.canvas import ChemusonCanvas
from chemuson.gui.main_window import ChemusonWindow
from chemuson.gui.orbitals import (
    ORBITAL_PALETTE_MODEL,
    orbital_display_name,
    orbital_tool_id,
    render_orbital_palette_image,
)


_GOLDEN_PATH = Path(__file__).with_name("data") / "orbitals" / "palette_golden.png"


@pytest.fixture(scope="module", autouse=True)
def _qapp():
    return QApplication.instance() or QApplication([])


@pytest.fixture(autouse=True)
def _isolated_config_home(tmp_path, monkeypatch):
    monkeypatch.setenv("XDG_CONFIG_HOME", str(tmp_path))


def _pixel_delta(actual: QImage, expected: QImage) -> tuple[int, int]:
    bad_pixels = 0
    worst = 0
    for y in range(actual.height()):
        for x in range(actual.width()):
            a = actual.pixelColor(x, y)
            b = expected.pixelColor(x, y)
            diff = max(
                abs(a.red() - b.red()),
                abs(a.green() - b.green()),
                abs(a.blue() - b.blue()),
                abs(a.alpha() - b.alpha()),
            )
            worst = max(worst, diff)
            if diff > 6:
                bad_pixels += 1
    return bad_pixels, worst


def test_orbital_tool_id_maps_to_generic_canvas_tool() -> None:
    canvas = ChemusonCanvas()
    try:
        canvas.set_current_tool(orbital_tool_id("dz2_shaded"))

        assert canvas.current_tool == "tool_orbital"
        assert canvas.state.active_tool == "tool_orbital"
        assert canvas.state.active_orbital_kind == "dz2_shaded"
    finally:
        canvas.close()


def test_inserted_orbital_is_persistent_and_undoable() -> None:
    canvas = ChemusonCanvas()
    restored = None
    try:
        canvas.state.active_orbital_kind = "pi_bonding_shaded"
        center = QPointF(160.0, 140.0)

        item = canvas._insert_orbital_item(center)

        assert item is not None
        assert len(canvas.orbital_items) == 1
        assert canvas.orbital_items[0].kind() == "pi_bonding_shaded"
        assert canvas.orbital_items[0].anchor0() == center

        data = canvas.get_persistence_data()
        restored = ChemusonCanvas()
        restored.load_persistence_data(data)

        assert len(restored.orbital_items) == 1
        assert restored.orbital_items[0].kind() == "pi_bonding_shaded"
        assert restored.orbital_items[0].anchor0() == center
        assert restored.orbital_items[0].anchor1().y() < center.y()

        canvas.undo_stack.undo()
        assert len(canvas.orbital_items) == 0

        canvas.undo_stack.redo()
        assert len(canvas.orbital_items) == 1
        assert canvas.orbital_items[0].kind() == "pi_bonding_shaded"
    finally:
        if restored is not None:
            restored.close()
        canvas.close()


def test_orbital_toolbar_selection_updates_active_canvas() -> None:
    window = ChemusonWindow()
    try:
        tool_id = orbital_tool_id("sigma_bonding_solid")
        window.symbols_toolbar._select_orbital_tool(tool_id)

        assert window._current_tool_id == tool_id
        assert window.canvas.current_tool == "tool_orbital"
        assert window.canvas.state.active_orbital_kind == "sigma_bonding_solid"
        assert window.symbols_toolbar.orbital_action.isChecked() is True
        assert window.symbols_toolbar.orbital_action.toolTip() == orbital_display_name("sigma_bonding_solid")
    finally:
        window.canvas.undo_stack.clear()
        window.canvas.undo_stack.setClean()
        window.close()


def test_orbital_palette_snapshot_matches_golden() -> None:
    image = render_orbital_palette_image(ORBITAL_PALETTE_MODEL).convertToFormat(
        QImage.Format.Format_ARGB32
    )
    expected = QImage(str(_GOLDEN_PATH)).convertToFormat(QImage.Format.Format_ARGB32)

    assert expected.isNull() is False, f"Golden faltante: {_GOLDEN_PATH}"
    assert image.size() == expected.size()

    bad_pixels, worst = _pixel_delta(image, expected)
    total_pixels = image.width() * image.height()

    assert worst <= 24
    assert bad_pixels / total_pixels <= 0.004
