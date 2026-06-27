"""Regresiones para selección y ciclo de vida de edición de texto."""

from __future__ import annotations


import pytest
from PyQt6.QtCore import QPoint, QPointF, Qt
from PyQt6.QtGui import QFont, QKeySequence
from PyQt6.QtTest import QTest
from PyQt6.QtWidgets import QApplication


from chemuson.gui.canvas import ChemusonCanvas
from chemuson.gui.items import TextAnnotationItem
from chemuson.gui.main_window import ChemusonWindow
from chemuson.gui.text_toolbar import TextFormatToolbar


@pytest.fixture(scope="module", autouse=True)
def _qapp():
    return QApplication.instance() or QApplication([])


@pytest.fixture(autouse=True)
def _isolated_config_home(tmp_path, monkeypatch):
    monkeypatch.setenv("XDG_CONFIG_HOME", str(tmp_path))


def _text_items(canvas: ChemusonCanvas) -> list[TextAnnotationItem]:
    return [item for item in canvas.scene.items() if isinstance(item, TextAnnotationItem)]


def _format_flag(fmt, property_name: str) -> bool:
    if property_name == "bold":
        return int(fmt.fontWeight()) >= int(QFont.Weight.Bold)
    if property_name == "italic":
        return bool(fmt.fontItalic())
    if property_name == "underline":
        return bool(fmt.fontUnderline())
    raise AssertionError(f"Unsupported property: {property_name}")


def test_programmatic_text_format_does_not_leave_document_selected() -> None:
    canvas = ChemusonCanvas()
    try:
        item = TextAnnotationItem("Condiciones", 120.0, 120.0)
        canvas.add_text_item(item)
        item.setSelected(True)
        canvas._sync_selection_from_scene()

        canvas.update_text_format("Arial", 14, True, False, False, False, False, "bold")

        assert item in canvas.scene.selectedItems()
        assert not item.textCursor().hasSelection()
    finally:
        canvas.close()


def test_switching_tool_finishes_text_editing_and_clears_selection() -> None:
    canvas = ChemusonCanvas()
    try:
        item = TextAnnotationItem("Intermedio", 140.0, 140.0)
        canvas.add_text_item(item)
        item.setSelected(True)
        item.setTextInteractionFlags(Qt.TextInteractionFlag.TextEditorInteraction)
        item.setFocus()
        cursor = item.textCursor()
        cursor.select(cursor.SelectionType.Document)
        item.setTextCursor(cursor)
        canvas.remember_text_edit_item(item)
        canvas.set_current_tool("tool_text")

        canvas.set_current_tool("tool_bond")

        assert canvas.current_tool == "tool_bond"
        assert item.scene() is canvas.scene
        assert item.textInteractionFlags() == Qt.TextInteractionFlag.NoTextInteraction
        assert not item.textCursor().hasSelection()
        assert not canvas.scene.selectedItems()
    finally:
        canvas.close()


def test_clicking_outside_active_text_editor_finishes_edit_without_creating_new_box() -> None:
    canvas = ChemusonCanvas()
    try:
        canvas.resize(900, 700)
        canvas.show()
        QApplication.processEvents()

        item = TextAnnotationItem("Producto", 150.0, 150.0)
        canvas.add_text_item(item)
        item.setSelected(True)
        item.setTextInteractionFlags(Qt.TextInteractionFlag.TextEditorInteraction)
        item.setFocus()
        canvas.remember_text_edit_item(item)
        canvas.set_current_tool("tool_text")

        outside = canvas.mapFromScene(QPointF(320.0, 260.0))
        QTest.mouseClick(
            canvas.viewport(),
            Qt.MouseButton.LeftButton,
            Qt.KeyboardModifier.NoModifier,
            outside,
        )
        QApplication.processEvents()

        assert len(_text_items(canvas)) == 1
        assert item.textInteractionFlags() == Qt.TextInteractionFlag.NoTextInteraction
        assert not canvas.scene.selectedItems()
    finally:
        canvas.close()


def test_text_toolbar_uses_conventional_shortcuts() -> None:
    toolbar = TextFormatToolbar()

    assert toolbar.action_bold.shortcut() == QKeySequence("Ctrl+B")
    assert toolbar.action_italic.shortcut() == QKeySequence("Ctrl+I")
    assert toolbar.action_underline.shortcut() == QKeySequence("Ctrl+U")


@pytest.mark.parametrize("property_name", ["bold", "italic", "underline"])
def test_text_format_undo_redo_restores_selected_range(property_name: str) -> None:
    canvas = ChemusonCanvas()
    try:
        item = TextAnnotationItem("abc", 120.0, 120.0)
        canvas.add_text_item(item)
        item.setSelected(True)
        item.setTextInteractionFlags(Qt.TextInteractionFlag.TextEditorInteraction)
        item.setFocus()
        canvas.remember_text_edit_item(item)

        cursor = item.textCursor()
        cursor.setPosition(0)
        cursor.setPosition(1, cursor.MoveMode.KeepAnchor)
        item.setTextCursor(cursor)

        before_count = canvas.undo_stack.count()
        settings = dict(canvas._current_text_settings)
        settings[property_name] = True
        canvas.update_text_format(
            settings["family"],
            settings["size"],
            settings["bold"],
            settings["italic"],
            settings["underline"],
            settings["sub"],
            settings["sup"],
            property_name,
        )

        probe = item.textCursor()
        probe.setPosition(0)
        probe.setPosition(1, probe.MoveMode.KeepAnchor)
        fmt_after = probe.charFormat()
        assert canvas.undo_stack.count() == before_count + 1
        assert _format_flag(fmt_after, property_name) is True

        canvas.undo_stack.undo()

        probe = item.textCursor()
        probe.setPosition(0)
        probe.setPosition(1, probe.MoveMode.KeepAnchor)
        fmt_undo = probe.charFormat()
        assert _format_flag(fmt_undo, property_name) is False

        canvas.undo_stack.redo()

        probe = item.textCursor()
        probe.setPosition(0)
        probe.setPosition(1, probe.MoveMode.KeepAnchor)
        fmt_redo = probe.charFormat()
        assert _format_flag(fmt_redo, property_name) is True
    finally:
        canvas.close()


@pytest.mark.parametrize(
    ("property_name", "action_name"),
    [
        ("bold", "action_bold"),
        ("italic", "action_italic"),
        ("underline", "action_underline"),
    ],
)
def test_text_toolbar_button_resyncs_on_new_plain_selection(
    property_name: str,
    action_name: str,
) -> None:
    window = ChemusonWindow()
    try:
        canvas = window.canvas
        item = TextAnnotationItem("abc", 120.0, 120.0)
        canvas.add_text_item(item)
        item.setSelected(True)
        item.setTextInteractionFlags(Qt.TextInteractionFlag.TextEditorInteraction)
        item.setFocus()
        canvas.remember_text_edit_item(item)
        canvas.sync_text_selection_state()

        cursor = item.textCursor()
        cursor.setPosition(0)
        cursor.setPosition(1, cursor.MoveMode.KeepAnchor)
        item.setTextCursor(cursor)
        canvas.sync_text_selection_state()

        action = getattr(window.text_toolbar, action_name)
        action.trigger()

        assert action.isChecked() is True

        cursor = item.textCursor()
        cursor.setPosition(1)
        cursor.setPosition(2, cursor.MoveMode.KeepAnchor)
        item.setTextCursor(cursor)
        canvas.sync_text_selection_state()

        probe = item.textCursor()
        probe.setPosition(1)
        probe.setPosition(2, probe.MoveMode.KeepAnchor)
        fmt = probe.charFormat()
        assert _format_flag(fmt, property_name) is False
        assert action.isChecked() is False
    finally:
        window.canvas.undo_stack.clear()
        window.canvas.undo_stack.setClean()
        window.close()


def test_dragging_text_scale_handle_reflows_text_inside_new_box() -> None:
    canvas = ChemusonCanvas()
    try:
        canvas.resize(960, 720)
        canvas.show()
        QApplication.processEvents()

        item = TextAnnotationItem(
            (
                "Solid-phase heteroditopic receptors incorporating urea/thiourea "
                "anion-binding sites and polyether cation-binding spacers were "
                "designed for ion-pair recognition."
            ),
            100.0,
            120.0,
        )
        font = item.font()
        font.setPointSizeF(12.0)
        item.setFont(font)
        canvas.add_text_item(item)
        item.setSelected(True)
        canvas._sync_selection_from_scene()
        QApplication.processEvents()

        before_font_size = item.font().pointSizeF()
        before_rect = item.boundingRect()
        assert item.textWidth() < 0.0
        assert canvas._selection_scale_handle is not None

        handle_view = canvas.mapFromScene(canvas._selection_scale_handle.scenePos())
        shrunk_view = QPoint(handle_view.x() - 120, handle_view.y() - 10)
        QTest.mousePress(
            canvas.viewport(),
            Qt.MouseButton.LeftButton,
            Qt.KeyboardModifier.NoModifier,
            handle_view,
        )
        QTest.mouseMove(canvas.viewport(), shrunk_view)
        QTest.mouseRelease(
            canvas.viewport(),
            Qt.MouseButton.LeftButton,
            Qt.KeyboardModifier.NoModifier,
            shrunk_view,
        )
        QApplication.processEvents()

        assert item.font().pointSizeF() == pytest.approx(before_font_size)
        assert item.textWidth() > 0.0
        assert item.textWidth() < before_rect.width()
        assert item.boundingRect().height() > before_rect.height()

        canvas.undo_stack.undo()

        assert item.font().pointSizeF() == pytest.approx(before_font_size)
        assert item.textWidth() < 0.0
    finally:
        canvas.close()


def test_justify_alignment_preserves_text_box_width_and_updates_blocks() -> None:
    canvas = ChemusonCanvas()
    try:
        item = TextAnnotationItem(
            (
                "Solid-phase heteroditopic receptors incorporating urea/thiourea "
                "anion-binding sites and polyether cation-binding spacers were "
                "designed for ion-pair recognition."
            ),
            80.0,
            90.0,
        )
        item.setTextWidth(220.0)
        canvas.add_text_item(item)
        item.setSelected(True)
        canvas._sync_selection_from_scene()

        before_width = item.textWidth()
        canvas.update_text_alignment(Qt.AlignmentFlag.AlignJustify)

        assert item.textWidth() == pytest.approx(before_width)
        block = item.document().begin()
        assert block.isValid()
        assert bool(block.blockFormat().alignment() & Qt.AlignmentFlag.AlignJustify)
        assert item.document().defaultTextOption().alignment() == Qt.AlignmentFlag.AlignJustify
    finally:
        canvas.close()
