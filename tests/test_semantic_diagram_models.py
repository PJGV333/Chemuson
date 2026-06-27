"""Pruebas para modelos semánticos de diagramas electrónicos."""

from __future__ import annotations


import pytest
from PyQt6.QtCore import QPoint, QPointF, Qt
from PyQt6.QtTest import QTest
from PyQt6.QtWidgets import QGraphicsTextItem
from PyQt6.QtWidgets import QApplication


from chemuson.gui.canvas import ChemusonCanvas
from chemuson.gui.composite_diagram_item import CompositeDiagramItem
from chemuson.gui.diagram_layout import (
    SEMANTIC_TEXT_KIND_ROLE,
    build_items_from_semantic_diagram,
)
from chemuson.gui.diagram_models import (
    DiagramConnector,
    DiagramLane,
    DiagramLevel,
    SemanticDiagram,
)
from chemuson.gui.energy_diagrams import build_atomic_subshell_diagram
from chemuson.gui.items import EnergyDiagramItem


@pytest.fixture(scope="module", autouse=True)
def _qapp():
    return QApplication.instance() or QApplication([])


def test_semantic_diagram_json_roundtrip() -> None:
    diagram = SemanticDiagram(
        kind="molecular_orbital",
        title="Sample",
        lanes=[
            DiagramLane(id="left", title="Left", x=-120.0),
            DiagramLane(id="center", title="Center", x=0.0),
        ],
        levels=[
            DiagramLevel(
                id="sigma_2s",
                lane_id="center",
                energy=1.5,
                label="σ2s",
                representation="boxes",
                degeneracy=1,
                occupancies=[2],
                metadata={"bonding": True},
            ),
        ],
        connectors=[DiagramConnector("sigma_2s", "sigma_2s", "solid")],
        metadata={"total_electrons": 2},
    )

    restored = SemanticDiagram.from_json_dict(diagram.to_json_dict())

    assert restored.to_json_dict() == diagram.to_json_dict()
    assert restored.levels[0].occupancies == [2]


def test_diagram_level_normalizes_box_occupancies() -> None:
    level = DiagramLevel(
        id="2p",
        lane_id="p_lane",
        energy=2.0,
        label="2p",
        representation="boxes",
        degeneracy=3,
        occupancies=["pair", "up", 9, "ignored"],
    )

    assert level.occupancies == [2, 1, 2]


def test_semantic_diagram_canvas_persistence_roundtrip() -> None:
    canvas = ChemusonCanvas()
    restored = None
    try:
        diagram = build_atomic_subshell_diagram(7)
        inserted = canvas.insert_semantic_diagram(diagram, QPointF(180.0, 160.0))

        assert len(inserted) == 1
        assert len(canvas.semantic_diagram_items) == 1

        data = canvas.get_persistence_data()
        assert len(data["annotations"]["semantic_diagrams"]) == 1

        restored = ChemusonCanvas()
        restored.load_persistence_data(data)

        assert len(restored.semantic_diagram_items) == 1
        payload = restored.semantic_diagram_items[0].to_json()
        assert payload["semantic_diagram"]["kind"] == "atomic"
        assert payload["semantic_diagram"]["metadata"]["electron_count"] == 7
    finally:
        if restored is not None:
            restored.close()
        canvas.close()


def test_semantic_diagram_click_resolves_to_composite_and_copies() -> None:
    source = ChemusonCanvas()
    target = ChemusonCanvas()
    try:
        inserted = source.insert_semantic_diagram(
            build_atomic_subshell_diagram(8),
            QPointF(220.0, 180.0),
        )

        assert len(inserted) == 1
        composite = inserted[0]
        assert isinstance(composite, CompositeDiagramItem)

        level_item = next(
            item
            for item in composite.level_items()
            if item.metadata().get("semantic_level_id") == "2p"
        )
        child_center = level_item.mapToScene(level_item.boundingRect().center())

        assert source._get_item_at(child_center) is composite
        assert source._resolve_click_item(child_center) is composite

        source.scene.clearSelection()
        composite.setSelected(True)
        source._sync_selection_from_scene()
        assert source.has_copyable_selection() is True
        payload = source._build_selection_payload()

        assert payload is not None
        assert payload["energy_diagrams"] == []
        assert len(payload["semantic_diagrams"]) == 1
        source.copy_to_clipboard()
        assert target.can_paste_from_clipboard() is True

        target._last_scene_pos = QPointF(360.0, 260.0)
        target._paste_selection_payload(payload)

        assert len(target.semantic_diagram_items) == 1
        pasted = target.semantic_diagram_items[0]
        assert pasted.isSelected() is True
        assert pasted.to_json()["semantic_diagram"]["metadata"]["electron_count"] == 8
    finally:
        source.close()
        target.close()


def test_semantic_layout_keeps_titles_above_content_and_compacts_single_levels() -> None:
    diagram = build_atomic_subshell_diagram(8)
    items = build_items_from_semantic_diagram(diagram)
    levels = [item for item in items if isinstance(item, EnergyDiagramItem)]
    texts = [item for item in items if isinstance(item, QGraphicsTextItem)]

    assert levels
    assert texts

    content_top = min(item.display_rect().top() for item in levels)
    header_texts = [
        item
        for item in texts
        if str(item.data(SEMANTIC_TEXT_KIND_ROLE) or "") != "diagram_summary"
    ]
    assert all(
        item.pos().y() + item.boundingRect().height() <= content_top
        for item in header_texts
    )

    widths = {
        str(item.metadata().get("semantic_level_id")): item.display_rect().width()
        for item in levels
    }
    assert widths["1s"] < 70.0
    assert widths["2p"] > widths["1s"]


def test_semantic_diagram_labels_are_editable_and_persisted() -> None:
    composite = CompositeDiagramItem(build_atomic_subshell_diagram(8))

    assert composite.set_diagram_title("Diagrama actualizado") is True
    assert composite.set_lane_title("s_lane", "Subniveles s") is True
    assert composite.set_level_label("2p", "2p valencia") is True

    payload = composite.to_json()["semantic_diagram"]

    assert payload["title"] == "Diagrama actualizado"
    assert payload["lanes"][0]["title"] == "Subniveles s"
    level_payload = next(level for level in payload["levels"] if level["id"] == "2p")
    assert level_payload["label"] == "2p valencia"


def test_semantic_layout_renders_rich_text_titles() -> None:
    diagram = build_atomic_subshell_diagram(8)
    diagram.title = "O<sub>2</sub> Diagram"

    items = build_items_from_semantic_diagram(diagram)
    title_item = next(
        item
        for item in items
        if isinstance(item, QGraphicsTextItem)
        and str(item.data(SEMANTIC_TEXT_KIND_ROLE) or "") == "diagram_title"
    )

    assert "sub" in title_item.toHtml().lower()


def test_semantic_diagram_copy_exports_png_for_external_apps() -> None:
    canvas = ChemusonCanvas()
    try:
        inserted = canvas.insert_semantic_diagram(
            build_atomic_subshell_diagram(8),
            QPointF(220.0, 180.0),
        )
        assert inserted
        canvas.copy_to_clipboard()

        mime = QApplication.clipboard().mimeData()
        assert mime is not None
        assert mime.hasFormat("application/x-chemuson-selection")
        assert mime.hasFormat("image/png")
        assert mime.hasImage()
        assert bytes(mime.data("image/png"))
    finally:
        canvas.close()


def test_pasted_semantic_diagram_can_drag_and_delete_with_keyboard() -> None:
    source = ChemusonCanvas()
    target = ChemusonCanvas()
    try:
        source.insert_semantic_diagram(
            build_atomic_subshell_diagram(8),
            QPointF(220.0, 180.0),
        )
        source.copy_to_clipboard()

        target.resize(900, 700)
        target.show()
        QApplication.processEvents()
        target._last_scene_pos = QPointF(320.0, 240.0)
        target.paste_from_clipboard()
        QApplication.processEvents()

        assert len(target.semantic_diagram_items) == 1
        item = target.semantic_diagram_items[0]
        before_rect = item.display_rect()

        start_view = target.mapFromScene(before_rect.center())
        end_view = QPoint(start_view.x() + 32, start_view.y() + 18)
        QTest.mousePress(
            target.viewport(),
            Qt.MouseButton.LeftButton,
            Qt.KeyboardModifier.NoModifier,
            start_view,
        )
        QTest.mouseMove(target.viewport(), end_view)
        QTest.mouseRelease(
            target.viewport(),
            Qt.MouseButton.LeftButton,
            Qt.KeyboardModifier.NoModifier,
            end_view,
        )
        QApplication.processEvents()

        moved_rect = item.display_rect()
        assert moved_rect.x() == pytest.approx(before_rect.x() + 32.0, abs=1.0)
        assert moved_rect.y() == pytest.approx(before_rect.y() + 18.0, abs=1.0)

        target.setFocus()
        target.viewport().setFocus()
        QTest.keyClick(target.viewport(), Qt.Key.Key_Delete)
        QApplication.processEvents()

        assert len(target.semantic_diagram_items) == 0
    finally:
        source.close()
        target.close()
