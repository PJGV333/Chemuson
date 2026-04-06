from __future__ import annotations

import re
from typing import Callable

from PyQt6.QtCore import QEventLoop, Qt
from PyQt6.QtGui import QFont, QKeySequence, QShortcut, QTextCharFormat, QTextCursor
from PyQt6.QtWidgets import (
    QDialog,
    QDialogButtonBox,
    QLabel,
    QTextEdit,
    QVBoxLayout,
    QWidget,
)


def rich_text_editor_value(editor: QTextEdit) -> str:
    """Devuelve texto plano cuando no hay formato relevante, u HTML si sí lo hay."""
    html_value = str(editor.toHtml() or "")
    plain_value = str(editor.toPlainText() or "")
    rich_markers = (
        "<span style=",
        "vertical-align:",
        "font-weight:700",
        "font-weight:600",
        "font-style:italic",
        "text-decoration:",
        "text-align:center",
        "text-align:right",
        "text-align:justify",
        "<ul",
        "<ol",
    )
    if not any(marker in html_value for marker in rich_markers):
        return plain_value
    html_value = html_value.replace("<!--StartFragment-->", "").replace("<!--EndFragment-->", "")
    html_value = re.sub(r"\s*font-family:'[^']*';", " ", html_value)
    html_value = re.sub(r"\s*font-size:[0-9.]+pt;", " ", html_value)
    html_value = re.sub(r"\s*font-weight:400;", " ", html_value)
    html_value = re.sub(r"\s*font-style:normal;", " ", html_value)
    html_value = re.sub(r"\s{2,}", " ", html_value)
    return html_value


def open_rich_text_value_dialog(
    parent: QWidget,
    *,
    title: str,
    label: str,
    initial_text: str = "",
    set_external_text_editor: Callable[[QTextEdit | None], None],
    sync_text_toolbar: Callable[[], None],
) -> tuple[str, bool]:
    """Abre un editor enriquecido temporal conectado al toolbar principal."""
    dialog = QDialog(parent)
    dialog.setWindowTitle(title)
    dialog.setModal(False)
    dialog.setWindowModality(Qt.WindowModality.NonModal)
    dialog.resize(420, 240)
    layout = QVBoxLayout(dialog)
    layout.addWidget(QLabel(label, dialog))

    editor = QTextEdit(dialog)
    editor.setAcceptRichText(True)
    editor.setMinimumHeight(96)
    editor.setFont(QFont("Arial", 10))
    if Qt.mightBeRichText(str(initial_text or "")):
        editor.setHtml(str(initial_text or ""))
    else:
        editor.setPlainText(str(initial_text or ""))
    layout.addWidget(editor)

    def handle_format(prop_name: str) -> None:
        editor.setFocus()
        cursor = editor.textCursor()
        fmt = QTextCharFormat()
        current_fmt = editor.currentCharFormat()

        if prop_name == "bold":
            is_bold = current_fmt.fontWeight() != QFont.Weight.Bold
            fmt.setFontWeight(QFont.Weight.Bold if is_bold else QFont.Weight.Normal)
        elif prop_name == "italic":
            fmt.setFontItalic(not current_fmt.fontItalic())
        elif prop_name == "underline":
            fmt.setFontUnderline(not current_fmt.fontUnderline())
        elif prop_name == "sub":
            current_align = current_fmt.verticalAlignment()
            if current_align == QTextCharFormat.VerticalAlignment.AlignSubScript:
                fmt.setVerticalAlignment(QTextCharFormat.VerticalAlignment.AlignNormal)
            else:
                fmt.setVerticalAlignment(QTextCharFormat.VerticalAlignment.AlignSubScript)
        elif prop_name == "sup":
            current_align = current_fmt.verticalAlignment()
            if current_align == QTextCharFormat.VerticalAlignment.AlignSuperScript:
                fmt.setVerticalAlignment(QTextCharFormat.VerticalAlignment.AlignNormal)
            else:
                fmt.setVerticalAlignment(QTextCharFormat.VerticalAlignment.AlignSuperScript)

        if cursor.hasSelection():
            start = cursor.selectionStart()
            end = cursor.selectionEnd()
            cursor.setPosition(start)
            cursor.setPosition(end, QTextCursor.MoveMode.KeepAnchor)
            cursor.mergeCharFormat(fmt)
            editor.setTextCursor(cursor)
        else:
            editor.mergeCurrentCharFormat(fmt)
        sync_text_toolbar()

    shortcuts = (
        (QKeySequence.StandardKey.Bold, lambda: handle_format("bold")),
        (QKeySequence.StandardKey.Italic, lambda: handle_format("italic")),
        (QKeySequence.StandardKey.Underline, lambda: handle_format("underline")),
        (QKeySequence("Ctrl+="), lambda: handle_format("sub")),
        (QKeySequence("Ctrl+,"), lambda: handle_format("sub")),
        (QKeySequence("Ctrl+Shift+B"), lambda: handle_format("sub")),
        (QKeySequence("Ctrl+Shift+="), lambda: handle_format("sup")),
        (QKeySequence("Ctrl++"), lambda: handle_format("sup")),
        (QKeySequence("Ctrl+."), lambda: handle_format("sup")),
        (QKeySequence("Ctrl+Shift+P"), lambda: handle_format("sup")),
        (QKeySequence("Alt+-"), lambda: editor.insertPlainText("—")),
        (QKeySequence("Alt+Shift+-"), lambda: editor.insertPlainText("—")),
    )
    for sequence, handler in shortcuts:
        shortcut = QShortcut(sequence, dialog)
        shortcut.setContext(Qt.ShortcutContext.WidgetWithChildrenShortcut)
        shortcut.activated.connect(handler)

    buttons = QDialogButtonBox(
        QDialogButtonBox.StandardButton.Ok | QDialogButtonBox.StandardButton.Cancel,
        parent=dialog,
    )
    buttons.accepted.connect(dialog.accept)
    buttons.rejected.connect(dialog.reject)
    layout.addWidget(buttons)

    loop = QEventLoop(dialog)
    dialog.finished.connect(loop.quit)
    set_external_text_editor(editor)
    dialog.show()
    dialog.raise_()
    dialog.activateWindow()
    editor.setFocus()
    loop.exec()

    value = rich_text_editor_value(editor)
    accepted = dialog.result() == QDialog.DialogCode.Accepted
    set_external_text_editor(None)
    dialog.deleteLater()
    return value, accepted
