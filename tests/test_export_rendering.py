"""Regresiones para la ruta de exportación completa."""

import os
import sys

import pytest
from PyQt6.QtCore import Qt
from PyQt6.QtGui import QColor, QBrush, QPen
from PyQt6.QtWidgets import QApplication, QGraphicsDropShadowEffect, QGraphicsRectItem

sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), "..", "src")))

from chemuson.gui.canvas import ChemusonCanvas


@pytest.fixture(scope="module", autouse=True)
def _qapp():
    return QApplication.instance() or QApplication([])


def _add_orphan_paper(canvas: ChemusonCanvas) -> QGraphicsRectItem:
    """Inserta una hoja legacy huérfana para simular escenas antiguas."""
    orphan = QGraphicsRectItem(0, 0, canvas.paper_width, canvas.paper_height)
    orphan.setBrush(QBrush(Qt.GlobalColor.white))
    orphan.setPen(QPen(QColor("#CCCCCC"), 1))
    orphan.setZValue(-10)
    shadow = QGraphicsDropShadowEffect()
    shadow.setBlurRadius(20)
    shadow.setColor(QColor(0, 0, 0, 80))
    shadow.setOffset(5, 5)
    orphan.setGraphicsEffect(shadow)
    canvas.scene.addItem(orphan)
    return orphan


def test_create_paper_removes_orphan_legacy_pages() -> None:
    """Recrear la hoja no debe acumular varias páginas en la escena."""
    canvas = ChemusonCanvas()
    _add_orphan_paper(canvas)

    assert len(canvas._paper_scene_items()) == 2

    canvas._create_paper()

    papers = canvas._paper_scene_items()
    assert len(papers) == 1
    assert papers[0] is canvas.paper


def test_full_export_render_hides_orphan_paper_shadow_band() -> None:
    """Exportar la escena completa no debe arrastrar bandas de hojas huérfanas."""
    canvas = ChemusonCanvas()
    _add_orphan_paper(canvas)

    previous = None
    for x in [100.0, 350.0, 600.0, 850.0, 1100.0]:
        atom = canvas.model.add_atom("C", x, 200.0, is_explicit=True)
        if previous is not None:
            canvas.model.add_bond(previous, atom.id, order=1)
        previous = atom.id
    canvas._rebuild_items_from_model()

    image = canvas._render_scene_image(background=QColor("white"))
    assert image is not None

    tall_gray_columns = 0
    for x in range(image.width()):
        grayish = 0
        for y in range(image.height()):
            color = image.pixelColor(x, y)
            if color == QColor("white"):
                continue
            if abs(color.red() - color.green()) < 10 and abs(color.green() - color.blue()) < 10:
                grayish += 1
        if grayish >= max(1, int(image.height() * 0.8)):
            tall_gray_columns += 1

    assert tall_gray_columns == 0
