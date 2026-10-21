from __future__ import annotations

import json
import os
from pathlib import Path

from PyQt6.QtCore import QPointF, Qt
from PyQt6.QtGui import QBrush, QColor, QPainterPath
from PyQt6.QtTest import QTest
from PyQt6.QtWidgets import QApplication

from chemuson.gui.canvas import ChemusonCanvas
from chemuson.gui.geom import endpoint_from_angle_len
from chemuson.gui.theme import apply_theme

OUT = Path(__file__).parent
app = QApplication.instance() or QApplication([])


def new_canvas() -> ChemusonCanvas:
    canvas = ChemusonCanvas()
    canvas.resize(1000, 760)
    canvas.show()
    app.processEvents()
    canvas.center_on_paper()
    app.processEvents()
    canvas.state.fixed_angles = True
    canvas.state.fixed_lengths = True
    canvas.state.angle_step_deg = 30
    canvas.state.bond_length = 40.0
    return canvas


def scene_pos(canvas: ChemusonCanvas, atom_id: int) -> QPointF:
    atom = canvas.model.get_atom(atom_id)
    return QPointF(atom.x, atom.y)


def atom_click(canvas: ChemusonCanvas, point: QPointF, element: str) -> int:
    canvas.set_current_tool("tool_atom")
    canvas.state.default_element = element
    pos = canvas.mapFromScene(point)
    before = set(canvas.model.atoms)
    QTest.mousePress(canvas.viewport(), Qt.MouseButton.LeftButton, pos=pos)
    QTest.mouseRelease(canvas.viewport(), Qt.MouseButton.LeftButton, pos=pos)
    app.processEvents()
    created = set(canvas.model.atoms) - before
    if len(created) != 1:
        raise RuntimeError(f"atom click expected one atom, got {created}")
    return created.pop()


def bond_click(canvas: ChemusonCanvas, anchor_id: int, element: str, order: int = 1) -> int:
    canvas.set_current_tool("tool_bond")
    canvas.state.default_element = element
    canvas.state.active_bond_order = order
    canvas.state.active_bond_aromatic = False
    canvas.state.active_bond_stereo = 0
    pos = canvas.mapFromScene(scene_pos(canvas, anchor_id))
    before = set(canvas.model.atoms)
    QTest.mousePress(canvas.viewport(), Qt.MouseButton.LeftButton, pos=pos)
    app.processEvents()
    preview = canvas._preview_bond_item.path()
    preview_end = QPointF(preview.elementAt(1).x, preview.elementAt(1).y)
    QTest.mouseRelease(canvas.viewport(), Qt.MouseButton.LeftButton, pos=pos)
    app.processEvents()
    created = set(canvas.model.atoms) - before
    if len(created) != 1:
        raise RuntimeError(f"bond click expected one atom, got {created}")
    atom_id = created.pop()
    atom = canvas.model.get_atom(atom_id)
    if (atom.x, atom.y) != (preview_end.x(), preview_end.y()):
        raise RuntimeError("preview/commit endpoint changed")
    return atom_id


def bond_drag(canvas: ChemusonCanvas, anchor_id: int, element: str, angle: float, order: int = 1) -> int:
    canvas.set_current_tool("tool_bond")
    canvas.state.default_element = element
    canvas.state.active_bond_order = order
    canvas.state.active_bond_aromatic = False
    start = scene_pos(canvas, anchor_id)
    end = endpoint_from_angle_len(start, angle, canvas.state.bond_length)
    p0, p1 = canvas.mapFromScene(start), canvas.mapFromScene(end)
    before = set(canvas.model.atoms)
    QTest.mousePress(canvas.viewport(), Qt.MouseButton.LeftButton, pos=p0)
    QTest.mouseMove(canvas.viewport(), p1, 10)
    app.processEvents()
    preview = canvas._preview_bond_item.path()
    preview_end = QPointF(preview.elementAt(1).x, preview.elementAt(1).y)
    QTest.mouseRelease(canvas.viewport(), Qt.MouseButton.LeftButton, pos=p1)
    app.processEvents()
    created = set(canvas.model.atoms) - before
    if len(created) != 1:
        raise RuntimeError(f"bond drag expected one atom, got {created}")
    atom_id = created.pop()
    atom = canvas.model.get_atom(atom_id)
    if (atom.x, atom.y) != (preview_end.x(), preview_end.y()):
        raise RuntimeError("preview/commit endpoint changed")
    return atom_id


def draw_ring(canvas: ChemusonCanvas) -> int:
    canvas.set_current_tool("tool_ring")
    canvas.state.active_ring_size = 6
    canvas.state.active_ring_aromatic = True
    center = canvas.mapToScene(canvas.viewport().rect().center())
    end = endpoint_from_angle_len(center, 0.0, canvas.state.bond_length)
    p0, p1 = canvas.mapFromScene(center), canvas.mapFromScene(end)
    QTest.mousePress(canvas.viewport(), Qt.MouseButton.LeftButton, pos=p0)
    QTest.mouseMove(canvas.viewport(), p1, 10)
    QTest.mouseRelease(canvas.viewport(), Qt.MouseButton.LeftButton, pos=p1)
    app.processEvents()
    return max((atom for atom in canvas.model.atoms.values()), key=lambda atom: atom.x).id


def draw_case(name: str) -> tuple[ChemusonCanvas, dict[str, int]]:
    canvas = new_canvas()
    ids: dict[str, int] = {}
    if name in {"urea", "acetanilide"}:
        aryl_carbon = draw_ring(canvas)
        ids["aryl_carbon"] = aryl_carbon
        nitrogen = bond_click(canvas, aryl_carbon, "N")
        ids["nitrogen"] = nitrogen
        carbonyl = bond_click(canvas, nitrogen, "C")
        ids["carbonyl"] = carbonyl
        oxygen = bond_click(canvas, carbonyl, "O", 2)
        ids["oxygen"] = oxygen
        if name == "urea":
            second_n = bond_click(canvas, carbonyl, "N")
            methyl = bond_click(canvas, second_n, "C")
            ids.update(second_nitrogen=second_n, methyl=methyl)
        else:
            methyl = bond_click(canvas, carbonyl, "C")
            ids["methyl"] = methyl
    else:
        center = canvas.mapToScene(canvas.viewport().rect().center())
        nitrogen = atom_click(canvas, center, "N")
        ids["nitrogen"] = nitrogen
        first = bond_click(canvas, nitrogen, "C")
        ids["substituent1"] = first
        if name == "secondary-amine":
            second = bond_click(canvas, nitrogen, "C")
            ids["substituent2"] = second
    canvas.state.default_element = "C"
    canvas._update_hover(QPointF(-100000.0, -100000.0))
    app.processEvents()
    return canvas, ids


def label_metrics(canvas: ChemusonCanvas, atom_id: int) -> dict:
    atom = canvas.model.get_atom(atom_id)
    item = canvas.atom_items[atom_id]
    text_item = item.label
    block = text_item.document().begin()
    text = text_item.toPlainText()
    layout = block.layout()
    document_layout = text_item.document().documentLayout()
    block_rect = document_layout.blockBoundingRect(block)
    glyphs = {}
    for index, char in enumerate(text):
        if not char.isalpha():
            continue
        line = layout.lineForTextPosition(index)
        if not line.isValid():
            continue
        x = line.cursorToX(index)
        if isinstance(x, tuple):
            x = x[0]
        cursor = text_item.textCursor()
        cursor.setPosition(index)
        cursor.setPosition(index + 1, cursor.MoveMode.KeepAnchor)
        glyph_font = cursor.charFormat().font()
        path = QPainterPath()
        path.addText(block_rect.x() + line.position().x() + float(x), block_rect.y() + line.position().y() + line.ascent(), glyph_font, char)
        rect = text_item.mapToScene(path).boundingRect()
        glyphs.setdefault(char, []).append([rect.x(), rect.y(), rect.width(), rect.height()])
    rect = text_item.mapToScene(text_item.shape()).boundingRect()
    return {
        "text": text,
        "center": [atom.x, atom.y],
        "label_rect": [rect.x(), rect.y(), rect.width(), rect.height()],
        "glyph_rects": glyphs,
    }


cases = {}
for name in ("urea", "acetanilide", "primary-amine", "secondary-amine"):
    canvas, ids = draw_case(name)
    nitrogens = [atom_id for key, atom_id in ids.items() if "nitrogen" in key]
    cases[name] = {
        "ids": ids,
        "atoms": {
            str(atom_id): {
                "element": canvas.model.get_atom(atom_id).element,
                "coords": [canvas.model.get_atom(atom_id).x, canvas.model.get_atom(atom_id).y],
                "implicit_h": canvas.model.implicit_h_count(atom_id),
            }
            for atom_id in ids.values()
        },
        "labels": {str(atom_id): label_metrics(canvas, atom_id) for atom_id in nitrogens},
    }
    for theme, background in (("light", "#FFFFFF"), ("dark", "#111827")):
        apply_theme(canvas, theme)
        canvas.scene.setBackgroundBrush(QBrush(QColor(background)))
        app.processEvents()
        image = canvas.viewport().grab().toImage()
        path = OUT / f"{name}-before-{theme}.png"
        if not image.save(str(path)):
            raise RuntimeError(f"failed to save {path}")
    canvas.close()
    app.processEvents()

(OUT / "layout-before-metrics.json").write_text(json.dumps(cases, indent=2) + "\n", encoding="utf-8")
print(json.dumps({name: data["labels"] for name, data in cases.items()}, indent=2))
