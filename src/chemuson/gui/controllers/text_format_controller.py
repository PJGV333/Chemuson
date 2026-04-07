from __future__ import annotations

from PyQt6.QtCore import Qt
from PyQt6.QtGui import QBrush, QColor, QFont, QTextCharFormat, QTextCursor
from PyQt6.QtWidgets import QFontDialog, QInputDialog, QTextEdit


class TextFormatController:
    """Maneja formato de texto externo sin acoplar al hub principal."""

    def on_text_format_changed(
        self,
        window,
        family: str,
        size: int,
        bold: bool,
        italic: bool,
        underline: bool,
        sub: bool,
        sup: bool,
        property_name: str,
    ) -> None:
        if self.apply_text_format_to_external_editor(
            window,
            family,
            size,
            bold,
            italic,
            underline,
            sub,
            sup,
            property_name,
        ):
            return
        window.canvas.update_text_format(
            family,
            size,
            bold,
            italic,
            underline,
            sub,
            sup,
            property_name,
        )

    def on_text_color_changed(self, window, color: QColor) -> None:
        if self.apply_text_color_to_external_editor(window, color):
            return
        window.canvas.update_text_color(color)

    def on_text_alignment_changed(self, window, alignment: Qt.AlignmentFlag) -> None:
        if self.apply_text_alignment_to_external_editor(window, alignment):
            return
        window.canvas.update_text_alignment(alignment)

    def on_external_text_copy_available(self, window, _available: bool) -> None:
        self.remember_external_text_cursor(window)

    def set_external_text_editor(self, window, editor: QTextEdit | None) -> None:
        previous = getattr(window, "_external_text_editor", None)
        if previous is editor:
            if editor is not None:
                self.sync_text_toolbar_from_external_editor(window)
            return
        if previous is not None:
            for signal, slot in (
                (previous.cursorPositionChanged, window._sync_text_toolbar_from_external_editor),
                (previous.textChanged, window._sync_text_toolbar_from_external_editor),
                (previous.cursorPositionChanged, window._remember_external_text_cursor),
                (previous.textChanged, window._remember_external_text_cursor),
                (previous.copyAvailable, window._on_external_text_copy_available),
            ):
                try:
                    signal.disconnect(slot)
                except (TypeError, RuntimeError):
                    pass
        window._external_text_editor = editor
        window._external_text_cursor_state = None
        window._external_text_selected_range = None
        window.text_toolbar.opacity_spin.setEnabled(editor is None)
        if editor is not None:
            editor.cursorPositionChanged.connect(window._sync_text_toolbar_from_external_editor)
            editor.textChanged.connect(window._sync_text_toolbar_from_external_editor)
            editor.cursorPositionChanged.connect(window._remember_external_text_cursor)
            editor.textChanged.connect(window._remember_external_text_cursor)
            editor.copyAvailable.connect(window._on_external_text_copy_available)
            self.remember_external_text_cursor(window)
            self.sync_text_toolbar_from_external_editor(window)
            return
        window.text_toolbar.set_opacity_percent(window.canvas.current_opacity_percent())
        window.canvas.sync_text_selection_state()

    def remember_external_text_cursor(self, window) -> None:
        editor = getattr(window, "_external_text_editor", None)
        if editor is None:
            window._external_text_cursor_state = None
            window._external_text_selected_range = None
            return
        cursor = editor.textCursor()
        anchor = int(cursor.anchor())
        position = int(cursor.position())
        window._external_text_cursor_state = (anchor, position)
        if cursor.hasSelection():
            window._external_text_selected_range = (anchor, position)
        elif editor.hasFocus():
            window._external_text_selected_range = None

    def external_text_cursor_for_formatting(self, window, editor: QTextEdit) -> QTextCursor:
        cursor = editor.textCursor()
        if cursor.hasSelection():
            start = cursor.selectionStart()
            end = cursor.selectionEnd()
            cursor.setPosition(start)
            cursor.setPosition(end, QTextCursor.MoveMode.KeepAnchor)
            return cursor

        stored_selected_range = getattr(window, "_external_text_selected_range", None)
        if stored_selected_range:
            anchor, position = stored_selected_range
            text_len = editor.document().characterCount() - 1
            if text_len < 0:
                text_len = 0
            anchor = max(0, min(int(anchor), text_len))
            position = max(0, min(int(position), text_len))
            cursor.setPosition(anchor)
            cursor.setPosition(position, QTextCursor.MoveMode.KeepAnchor)
            if cursor.hasSelection():
                return cursor

        stored = getattr(window, "_external_text_cursor_state", None)
        if stored:
            anchor, position = stored
            text_len = editor.document().characterCount() - 1
            if text_len < 0:
                text_len = 0
            anchor = max(0, min(int(anchor), text_len))
            position = max(0, min(int(position), text_len))
            cursor.setPosition(anchor)
            cursor.setPosition(position, QTextCursor.MoveMode.KeepAnchor)
            if cursor.hasSelection():
                start = cursor.selectionStart()
                end = cursor.selectionEnd()
                cursor.setPosition(start)
                cursor.setPosition(end, QTextCursor.MoveMode.KeepAnchor)

        return cursor

    def sync_text_toolbar_from_external_editor(self, window) -> None:
        editor = getattr(window, "_external_text_editor", None)
        if editor is None:
            return
        cursor = self.external_text_cursor_for_formatting(window, editor)
        fmt = cursor.charFormat()
        font = QFont(editor.currentFont())
        families = fmt.fontFamilies()
        if families:
            font.setFamily(families[0])
        point_size = float(fmt.fontPointSize() or 0.0)
        if point_size <= 0.0:
            point_size = float(font.pointSizeF() or font.pointSize() or 12.0)
        font.setPointSizeF(point_size)
        weight = int(fmt.fontWeight() or font.weight())
        font.setWeight(weight)
        font.setItalic(bool(fmt.fontItalic()))
        font.setUnderline(bool(fmt.fontUnderline()))
        color = fmt.foreground().color()
        if not color.isValid():
            color = QColor(Qt.GlobalColor.black)
        window.text_toolbar.set_state(
            font,
            {
                "color": color,
                "sub": fmt.verticalAlignment() == QTextCharFormat.VerticalAlignment.AlignSubScript,
                "sup": fmt.verticalAlignment() == QTextCharFormat.VerticalAlignment.AlignSuperScript,
            },
        )

    def apply_text_format_to_external_editor(
        self,
        window,
        family: str,
        size: int,
        bold: bool,
        italic: bool,
        underline: bool,
        sub: bool,
        sup: bool,
        property_name: str,
    ) -> bool:
        editor = getattr(window, "_external_text_editor", None)
        if editor is None:
            return False

        cursor = self.external_text_cursor_for_formatting(window, editor)
        fmt = QTextCharFormat()

        if property_name in ("all", "family"):
            fmt.setFontFamilies([family])
        if property_name in ("all", "size"):
            fmt.setFontPointSize(float(size))
        if property_name in ("all", "bold"):
            fmt.setFontWeight(QFont.Weight.Bold if bold else QFont.Weight.Normal)
        if property_name in ("all", "italic"):
            fmt.setFontItalic(bool(italic))
        if property_name in ("all", "underline"):
            fmt.setFontUnderline(bool(underline))
        if property_name in ("all", "sub", "sup"):
            if sub:
                fmt.setVerticalAlignment(QTextCharFormat.VerticalAlignment.AlignSubScript)
            elif sup:
                fmt.setVerticalAlignment(QTextCharFormat.VerticalAlignment.AlignSuperScript)
            else:
                fmt.setVerticalAlignment(QTextCharFormat.VerticalAlignment.AlignNormal)

        if cursor.hasSelection():
            cursor.mergeCharFormat(fmt)
            editor.setTextCursor(cursor)
        else:
            editor.mergeCurrentCharFormat(fmt)

        if not editor.hasFocus():
            editor.activateWindow()
            editor.raise_()
            editor.setFocus()

        self.remember_external_text_cursor(window)
        self.sync_text_toolbar_from_external_editor(window)
        return True

    def apply_text_color_to_external_editor(self, window, color: QColor) -> bool:
        editor = getattr(window, "_external_text_editor", None)
        if editor is None:
            return False
        editor.activateWindow()
        editor.raise_()
        editor.setFocus()
        fmt = QTextCharFormat()
        fmt.setForeground(QBrush(color))
        cursor = self.external_text_cursor_for_formatting(window, editor)
        if cursor.hasSelection():
            cursor.mergeCharFormat(fmt)
            editor.setTextCursor(cursor)
        editor.mergeCurrentCharFormat(fmt)
        self.remember_external_text_cursor(window)
        self.sync_text_toolbar_from_external_editor(window)
        return True

    def apply_text_alignment_to_external_editor(self, window, alignment: Qt.AlignmentFlag) -> bool:
        editor = getattr(window, "_external_text_editor", None)
        if editor is None:
            return False
        editor.activateWindow()
        editor.raise_()
        editor.setFocus()
        editor.setAlignment(alignment)
        self.remember_external_text_cursor(window)
        self.sync_text_toolbar_from_external_editor(window)
        return True

    def sync_label_menu_state(self, window) -> None:
        window.action_label_bold.setChecked(window.canvas.state.label_font_bold)
        window.action_label_italic.setChecked(window.canvas.state.label_font_italic)
        window.action_label_underline.setChecked(window.canvas.state.label_font_underline)
        window.action_label_color_element.setChecked(window.canvas.state.use_element_colors)
        window.action_label_color_black.setChecked(not window.canvas.state.use_element_colors)

    def on_label_font(self, window) -> None:
        font, ok = QFontDialog.getFont(
            window.canvas.label_font(),
            window,
            "Fuente de etiquetas",
        )
        if not ok:
            return
        window.canvas.apply_label_font(font)
        self.sync_label_menu_state(window)

    def on_label_font_size_dialog(self, window) -> None:
        size, ok = QInputDialog.getDouble(
            window,
            "Tamaño de etiquetas",
            "Tamaño (pt):",
            float(window.canvas.current_label_size_value()),
            6.0,
            72.0,
            1,
        )
        if not ok:
            return
        window.canvas.apply_label_font_size(float(size))
        self.sync_label_menu_state(window)

    def on_label_bold(self, window, checked: bool) -> None:
        if self.apply_text_format_to_external_editor(
            window,
            "",
            0,
            checked,
            False,
            False,
            False,
            False,
            "bold",
        ):
            return
        font = window.canvas.label_font()
        font.setBold(checked)
        window.canvas.apply_label_font(font)
        self.sync_label_menu_state(window)

    def on_label_italic(self, window, checked: bool) -> None:
        if self.apply_text_format_to_external_editor(
            window,
            "",
            0,
            False,
            checked,
            False,
            False,
            False,
            "italic",
        ):
            return
        font = window.canvas.label_font()
        font.setItalic(checked)
        window.canvas.apply_label_font(font)
        self.sync_label_menu_state(window)

    def on_label_underline(self, window, checked: bool) -> None:
        if self.apply_text_format_to_external_editor(
            window,
            "",
            0,
            False,
            False,
            checked,
            False,
            False,
            "underline",
        ):
            return
        font = window.canvas.label_font()
        font.setUnderline(checked)
        window.canvas.apply_label_font(font)
        self.sync_label_menu_state(window)

    def on_label_font_size(self, window, delta: float) -> None:
        window.canvas.adjust_label_font_size(delta)

    def on_label_subscript(self, window) -> None:
        if self.apply_text_format_to_external_editor(
            window,
            "",
            0,
            False,
            False,
            False,
            True,
            False,
            "sub",
        ):
            return
        window.canvas.update_text_format(sub=True, property_name="sub")

    def on_label_superscript(self, window) -> None:
        if self.apply_text_format_to_external_editor(
            window,
            "",
            0,
            False,
            False,
            False,
            False,
            True,
            "sup",
        ):
            return
        window.canvas.update_text_format(sup=True, property_name="sup")

    def on_label_color_mode(self, window, use_element_colors: bool) -> None:
        window.canvas.set_use_element_colors(use_element_colors)
        self.sync_label_menu_state(window)
        window.statusBar().showMessage(
            "Etiquetas: por elemento" if use_element_colors else "Etiquetas: negro"
        )
