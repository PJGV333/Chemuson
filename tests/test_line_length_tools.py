"""Pruebas para homogeneización de longitud en líneas de anotación."""

import math
import os

import pytest
from PyQt6.QtCore import QPointF
from PyQt6.QtWidgets import QApplication


from chemuson.gui.canvas import ChemusonCanvas


@pytest.fixture(scope="module", autouse=True)
def _qapp():
    os.environ.setdefault("QT_QPA_PLATFORM", "offscreen")
    return QApplication.instance() or QApplication([])


def _length(item) -> float:
    start = item.start_point()
    end = item.end_point()
    return math.hypot(end.x() - start.x(), end.y() - start.y())


def _center(item) -> QPointF:
    start = item.start_point()
    end = item.end_point()
    return QPointF((start.x() + end.x()) * 0.5, (start.y() + end.y()) * 0.5)


def test_set_line_arrow_length_preserves_center_and_undo() -> None:
    canvas = ChemusonCanvas()
    try:
        line = canvas.add_arrow_item(QPointF(10.0, 20.0), QPointF(50.0, 20.0), kind="line")
        line.setSelected(True)

        initial_center = _center(line)
        initial_length = _length(line)
        assert initial_length == pytest.approx(40.0)

        assert canvas._set_line_arrow_length(80.0) is True
        assert _length(line) == pytest.approx(80.0)
        assert _center(line).x() == pytest.approx(initial_center.x())
        assert _center(line).y() == pytest.approx(initial_center.y())

        canvas.undo_stack.undo()
        assert _length(line) == pytest.approx(initial_length)
        assert _center(line).x() == pytest.approx(initial_center.x())
        assert _center(line).y() == pytest.approx(initial_center.y())
    finally:
        canvas.close()


def test_equalize_selected_line_arrow_lengths_ignores_non_line_arrows() -> None:
    canvas = ChemusonCanvas()
    try:
        reference = canvas.add_arrow_item(
            QPointF(10.0, 20.0),
            QPointF(70.0, 20.0),
            kind="line",
        )
        dashed = canvas.add_arrow_item(
            QPointF(100.0, 50.0),
            QPointF(130.0, 80.0),
            kind="line_dashed",
        )
        reaction = canvas.add_arrow_item(
            QPointF(180.0, 40.0),
            QPointF(250.0, 40.0),
            kind="forward",
        )

        reference.setSelected(True)
        dashed.setSelected(True)
        reaction.setSelected(True)

        dashed_center = _center(dashed)
        dashed_before = _length(dashed)
        reaction_before = _length(reaction)

        assert canvas._equalize_selected_line_arrow_lengths(reference_item=reference) is True
        assert _length(reference) == pytest.approx(60.0)
        assert _length(dashed) == pytest.approx(60.0)
        assert _center(dashed).x() == pytest.approx(dashed_center.x())
        assert _center(dashed).y() == pytest.approx(dashed_center.y())
        assert _length(reaction) == pytest.approx(reaction_before)

        canvas.undo_stack.undo()
        assert _length(dashed) == pytest.approx(dashed_before)
        assert _length(reaction) == pytest.approx(reaction_before)
    finally:
        canvas.close()
