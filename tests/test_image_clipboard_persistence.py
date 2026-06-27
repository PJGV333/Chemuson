"""Regresiones para imágenes pegadas desde portapapeles y su persistencia."""

from __future__ import annotations


import pytest
from PyQt6.QtCore import QBuffer, QMimeData, QPointF, QRectF, QUrl
from PyQt6.QtGui import QImage
from PyQt6.QtCore import QPoint, Qt
from PyQt6.QtTest import QTest
from PyQt6.QtWidgets import QApplication


from chemuson.chemio.persistence import PersistenceManager
from chemuson.gui.commands import AddAtomCommand
from chemuson.gui.canvas import ChemusonCanvas
from chemuson.gui.items import ImageAnnotationItem


@pytest.fixture(scope="module", autouse=True)
def _qapp():
    return QApplication.instance() or QApplication([])


def _png_bytes(width: int = 24, height: int = 16) -> bytes:
    image = QImage(width, height, QImage.Format.Format_ARGB32)
    image.fill(0xFF3366CC)
    buffer = QBuffer()
    buffer.open(QBuffer.OpenModeFlag.WriteOnly)
    image.save(buffer, "PNG")
    return bytes(buffer.data())


def test_paste_local_png_file_url_selects_image_and_supports_undo(tmp_path):
    """Pegar un archivo PNG copiado desde el explorador debe crear una imagen seleccionada."""
    image_path = tmp_path / "clipboard-image.png"
    image_path.write_bytes(_png_bytes())

    mime = QMimeData()
    mime.setUrls([QUrl.fromLocalFile(str(image_path))])
    QApplication.clipboard().setMimeData(mime)

    canvas = ChemusonCanvas()
    canvas._last_scene_pos = QPointF(220.0, 180.0)
    canvas.paste_from_clipboard()

    assert len(canvas.image_items) == 1
    item = canvas.image_items[0]
    assert item in canvas.scene.selectedItems()
    assert item.mime_type() == "image/png"

    rect = item.display_rect()
    assert rect.center().x() == pytest.approx(220.0, abs=0.6)
    assert rect.center().y() == pytest.approx(180.0, abs=0.6)

    canvas.undo_stack.undo()
    assert len(canvas.image_items) == 0

    canvas.undo_stack.redo()
    assert len(canvas.image_items) == 1
    assert canvas.image_items[0] in canvas.scene.selectedItems()


def test_image_copy_payload_and_persistence_roundtrip():
    """Las imágenes deben copiarse internamente y persistir en archivos .cmsn."""
    canvas = ChemusonCanvas()
    item = ImageAnnotationItem(_png_bytes(), "image/png", width=48.0, height=24.0, source_name="test.png")
    item.set_display_rect(QRectF(80.0, 120.0, 48.0, 24.0))
    item.setRotation(18.0)
    canvas.readd_image_item(item)

    item.setSelected(True)
    canvas._sync_selection_from_scene()
    payload = canvas._build_selection_payload()
    assert payload is not None
    assert len(payload["images"]) == 1

    target = ChemusonCanvas()
    target._last_scene_pos = QPointF(260.0, 260.0)
    target._paste_selection_payload(payload)
    assert len(target.image_items) == 1
    pasted = target.image_items[0]
    assert pasted in target.scene.selectedItems()
    assert pasted.mime_type() == "image/png"

    data = PersistenceManager.save_to_dict(canvas)
    assert len(data["canvas"]["annotations"]["images"]) == 1

    restored = ChemusonCanvas()
    PersistenceManager.load_from_dict(data, restored)

    assert len(restored.image_items) == 1
    restored_item = restored.image_items[0]
    restored_rect = restored_item.display_rect()
    assert restored_item.mime_type() == "image/png"
    assert restored_item.source_name() == "test.png"
    assert restored_rect.x() == pytest.approx(80.0, abs=0.2)
    assert restored_rect.y() == pytest.approx(120.0, abs=0.2)
    assert restored_rect.width() == pytest.approx(48.0, abs=0.2)
    assert restored_rect.height() == pytest.approx(24.0, abs=0.2)
    assert restored_item.rotation() == pytest.approx(18.0, abs=0.2)


def test_paste_local_svg_file_url_when_supported(tmp_path):
    """Si el runtime soporta SVG, el pegado desde archivo debe preservarlo como tal."""
    svg_bytes = (
        b'<svg xmlns="http://www.w3.org/2000/svg" width="40" height="20" viewBox="0 0 40 20">'
        b'<rect x="1" y="1" width="38" height="18" fill="#88ccff" stroke="#003366" />'
        b"</svg>"
    )
    probe = ImageAnnotationItem(svg_bytes, "image/svg+xml", width=40.0, height=20.0)
    if not probe.is_valid_image():
        pytest.skip("SVG no soportado por el runtime Qt actual")

    image_path = tmp_path / "clipboard-image.svg"
    image_path.write_bytes(svg_bytes)

    mime = QMimeData()
    mime.setUrls([QUrl.fromLocalFile(str(image_path))])
    QApplication.clipboard().setMimeData(mime)

    canvas = ChemusonCanvas()
    canvas.paste_from_clipboard()

    assert len(canvas.image_items) == 1
    assert canvas.image_items[0].mime_type() == "image/svg+xml"


def test_mouse_drag_and_scale_on_selected_image(tmp_path):
    """La imagen pegada debe poder moverse y escalarse con el mouse."""
    image_path = tmp_path / "interactive-image.png"
    image_path.write_bytes(_png_bytes(120, 60))

    mime = QMimeData()
    mime.setUrls([QUrl.fromLocalFile(str(image_path))])
    QApplication.clipboard().setMimeData(mime)

    canvas = ChemusonCanvas()
    canvas.resize(900, 700)
    canvas.show()
    QApplication.processEvents()
    canvas._last_scene_pos = QPointF(260.0, 220.0)
    canvas.paste_from_clipboard()
    QApplication.processEvents()

    item = canvas.image_items[0]
    before_rect = item.display_rect()

    start_scene = before_rect.center()
    start_view = canvas.mapFromScene(start_scene)
    end_view = QPoint(start_view.x() + 36, start_view.y() + 18)
    QTest.mousePress(
        canvas.viewport(),
        Qt.MouseButton.LeftButton,
        Qt.KeyboardModifier.NoModifier,
        start_view,
    )
    QTest.mouseMove(canvas.viewport(), end_view)
    QTest.mouseRelease(
        canvas.viewport(),
        Qt.MouseButton.LeftButton,
        Qt.KeyboardModifier.NoModifier,
        end_view,
    )
    QApplication.processEvents()

    moved_rect = item.display_rect()
    assert moved_rect.x() == pytest.approx(before_rect.x() + 36.0, abs=1.0)
    assert moved_rect.y() == pytest.approx(before_rect.y() + 18.0, abs=1.0)

    handle_scene = canvas._selection_scale_handle.scenePos()
    handle_view = canvas.mapFromScene(handle_scene)
    scaled_view = QPoint(handle_view.x() + 28, handle_view.y() + 14)
    QTest.mousePress(
        canvas.viewport(),
        Qt.MouseButton.LeftButton,
        Qt.KeyboardModifier.NoModifier,
        handle_view,
    )
    QTest.mouseMove(canvas.viewport(), scaled_view)
    QTest.mouseRelease(
        canvas.viewport(),
        Qt.MouseButton.LeftButton,
        Qt.KeyboardModifier.NoModifier,
        scaled_view,
    )
    QApplication.processEvents()

    scaled_rect = item.display_rect()
    assert scaled_rect.width() > moved_rect.width()
    assert scaled_rect.height() > moved_rect.height()


def test_rotation_handle_and_delete_key_work_for_selected_image(tmp_path):
    """La imagen seleccionada debe rotar desde el handle superior y borrarse con Supr."""
    image_path = tmp_path / "interactive-rotate-image.png"
    image_path.write_bytes(_png_bytes(96, 48))

    mime = QMimeData()
    mime.setUrls([QUrl.fromLocalFile(str(image_path))])
    QApplication.clipboard().setMimeData(mime)

    canvas = ChemusonCanvas()
    canvas.resize(900, 700)
    canvas.show()
    QApplication.processEvents()
    canvas._last_scene_pos = QPointF(280.0, 210.0)
    canvas.paste_from_clipboard()
    QApplication.processEvents()

    item = canvas.image_items[0]
    initial_rotation = item.rotation()

    rotate_scene = canvas._selection_handle.scenePos()
    rotate_view = canvas.mapFromScene(rotate_scene)
    rotated_view = QPoint(rotate_view.x() + 34, rotate_view.y() + 22)
    QTest.mousePress(
        canvas.viewport(),
        Qt.MouseButton.LeftButton,
        Qt.KeyboardModifier.NoModifier,
        rotate_view,
    )
    QTest.mouseMove(canvas.viewport(), rotated_view)
    QTest.mouseRelease(
        canvas.viewport(),
        Qt.MouseButton.LeftButton,
        Qt.KeyboardModifier.NoModifier,
        rotated_view,
    )
    QApplication.processEvents()

    assert item.rotation() != pytest.approx(initial_rotation, abs=0.2)

    QTest.keyClick(canvas.viewport(), Qt.Key.Key_Delete)
    QApplication.processEvents()

    assert len(canvas.image_items) == 0


def test_handles_work_for_image_with_zoom_and_scroll(tmp_path):
    """Los handles deben seguir funcionando con zoom y scroll activo."""
    image_path = tmp_path / "zoom-scroll-image.png"
    image_path.write_bytes(_png_bytes(800, 400))

    mime = QMimeData()
    mime.setUrls([QUrl.fromLocalFile(str(image_path))])
    QApplication.clipboard().setMimeData(mime)

    canvas = ChemusonCanvas()
    canvas.resize(1200, 800)
    canvas.show()
    QApplication.processEvents()
    canvas.zoom_in()
    canvas.zoom_in()
    canvas.horizontalScrollBar().setValue(700)
    canvas.verticalScrollBar().setValue(600)
    QApplication.processEvents()

    canvas._last_scene_pos = QPointF(850.0, 900.0)
    canvas.paste_from_clipboard()
    QApplication.processEvents()

    item = canvas.image_items[0]
    moved_before = item.display_rect()

    move_scene = canvas._selection_move_handle.scenePos()
    move_view = canvas.mapFromScene(move_scene)
    moved_view = QPoint(move_view.x() + 48, move_view.y() + 26)
    QTest.mousePress(
        canvas.viewport(),
        Qt.MouseButton.LeftButton,
        Qt.KeyboardModifier.NoModifier,
        move_view,
    )
    QTest.mouseMove(canvas.viewport(), moved_view)
    QTest.mouseRelease(
        canvas.viewport(),
        Qt.MouseButton.LeftButton,
        Qt.KeyboardModifier.NoModifier,
        moved_view,
    )
    QApplication.processEvents()

    moved_after = item.display_rect()
    assert moved_after.x() > moved_before.x()
    assert moved_after.y() > moved_before.y()

    scale_scene = canvas._selection_scale_handle.scenePos()
    scale_view = canvas.mapFromScene(scale_scene)
    scaled_view = QPoint(scale_view.x() + 36, scale_view.y() + 18)
    QTest.mousePress(
        canvas.viewport(),
        Qt.MouseButton.LeftButton,
        Qt.KeyboardModifier.NoModifier,
        scale_view,
    )
    QTest.mouseMove(canvas.viewport(), scaled_view)
    QTest.mouseRelease(
        canvas.viewport(),
        Qt.MouseButton.LeftButton,
        Qt.KeyboardModifier.NoModifier,
        scaled_view,
    )
    QApplication.processEvents()

    scaled_after = item.display_rect()
    assert scaled_after.width() > moved_after.width()
    assert scaled_after.height() > moved_after.height()


def test_selected_image_remains_transformable_over_busy_chemistry(tmp_path):
    """Una imagen seleccionada debe seguir mandando aunque haya química superpuesta."""
    image_path = tmp_path / "busy-canvas-image.png"
    image_path.write_bytes(_png_bytes(720, 360))

    mime = QMimeData()
    mime.setUrls([QUrl.fromLocalFile(str(image_path))])
    QApplication.clipboard().setMimeData(mime)

    canvas = ChemusonCanvas()
    canvas.resize(1200, 900)
    canvas.show()
    QApplication.processEvents()

    atom_a = canvas.model.add_atom("C", 500.0, 500.0)
    atom_b = canvas.model.add_atom("C", 560.0, 500.0)
    canvas.model.add_bond(atom_a.id, atom_b.id)
    canvas._rebuild_items_from_model()
    canvas._insert_wavy_anchor(QPointF(530.0, 500.0), atom_a.id)
    QApplication.processEvents()

    canvas._last_scene_pos = QPointF(500.0, 500.0)
    canvas.paste_from_clipboard()
    QApplication.processEvents()

    item = canvas.image_items[0]
    before_rect = item.display_rect()

    center_view = canvas.mapFromScene(before_rect.center())
    moved_view = QPoint(center_view.x() + 42, center_view.y() + 24)
    QTest.mousePress(
        canvas.viewport(),
        Qt.MouseButton.LeftButton,
        Qt.KeyboardModifier.NoModifier,
        center_view,
    )
    QTest.mouseMove(canvas.viewport(), moved_view)
    QTest.mouseRelease(
        canvas.viewport(),
        Qt.MouseButton.LeftButton,
        Qt.KeyboardModifier.NoModifier,
        moved_view,
    )
    QApplication.processEvents()

    moved_rect = item.display_rect()
    assert moved_rect.x() == pytest.approx(before_rect.x() + 42.0, abs=1.0)
    assert moved_rect.y() == pytest.approx(before_rect.y() + 24.0, abs=1.0)
    assert item in canvas.scene.selectedItems()

    scale_scene = canvas._selection_scale_handle.scenePos()
    scale_view = canvas.mapFromScene(scale_scene)
    scaled_view = QPoint(scale_view.x() + 30, scale_view.y() + 16)
    QTest.mousePress(
        canvas.viewport(),
        Qt.MouseButton.LeftButton,
        Qt.KeyboardModifier.NoModifier,
        scale_view,
    )
    QTest.mouseMove(canvas.viewport(), scaled_view)
    QTest.mouseRelease(
        canvas.viewport(),
        Qt.MouseButton.LeftButton,
        Qt.KeyboardModifier.NoModifier,
        scaled_view,
    )
    QApplication.processEvents()

    scaled_rect = item.display_rect()
    assert scaled_rect.width() > moved_rect.width()
    assert scaled_rect.height() > moved_rect.height()

    rotate_scene = canvas._selection_handle.scenePos()
    rotate_view = canvas.mapFromScene(rotate_scene)
    rotated_view = QPoint(rotate_view.x() + 28, rotate_view.y() + 18)
    initial_rotation = item.rotation()
    QTest.mousePress(
        canvas.viewport(),
        Qt.MouseButton.LeftButton,
        Qt.KeyboardModifier.NoModifier,
        rotate_view,
    )
    QTest.mouseMove(canvas.viewport(), rotated_view)
    QTest.mouseRelease(
        canvas.viewport(),
        Qt.MouseButton.LeftButton,
        Qt.KeyboardModifier.NoModifier,
        rotated_view,
    )
    QApplication.processEvents()

    assert item.rotation() != pytest.approx(initial_rotation, abs=0.2)


def test_selected_image_survives_structure_edit_refresh(tmp_path):
    """Editar química existente no debe romper ni deseleccionar la imagen activa."""
    image_path = tmp_path / "refresh-after-edit-image.png"
    image_path.write_bytes(_png_bytes(720, 360))

    mime = QMimeData()
    mime.setUrls([QUrl.fromLocalFile(str(image_path))])
    QApplication.clipboard().setMimeData(mime)

    canvas = ChemusonCanvas()
    canvas.resize(1200, 900)
    canvas.show()
    QApplication.processEvents()

    atom_a = canvas.model.add_atom("C", 420.0, 420.0)
    atom_b = canvas.model.add_atom("N", 480.0, 420.0)
    canvas.model.add_bond(atom_a.id, atom_b.id)
    canvas._rebuild_items_from_model()
    QApplication.processEvents()

    canvas._last_scene_pos = QPointF(520.0, 520.0)
    canvas.paste_from_clipboard()
    QApplication.processEvents()

    item = canvas.image_items[0]
    assert item in canvas.scene.selectedItems()

    canvas.undo_stack.push(
        AddAtomCommand(canvas.model, canvas, "O", 760.0, 360.0, is_explicit=True)
    )
    QApplication.processEvents()

    assert item in canvas.scene.selectedItems()

    before_rect = item.display_rect()
    center_view = canvas.mapFromScene(before_rect.center())
    moved_view = QPoint(center_view.x() + 32, center_view.y() + 20)
    QTest.mousePress(
        canvas.viewport(),
        Qt.MouseButton.LeftButton,
        Qt.KeyboardModifier.NoModifier,
        center_view,
    )
    QTest.mouseMove(canvas.viewport(), moved_view)
    QTest.mouseRelease(
        canvas.viewport(),
        Qt.MouseButton.LeftButton,
        Qt.KeyboardModifier.NoModifier,
        moved_view,
    )
    QApplication.processEvents()

    moved_rect = item.display_rect()
    assert moved_rect.x() == pytest.approx(before_rect.x() + 32.0, abs=1.0)
    assert moved_rect.y() == pytest.approx(before_rect.y() + 20.0, abs=1.0)


def test_selected_image_can_transform_even_with_non_select_tool(tmp_path):
    """Pegar una imagen con otra herramienta activa no debe bloquear sus handles."""
    image_path = tmp_path / "non-select-tool-image.png"
    image_path.write_bytes(_png_bytes(720, 360))

    mime = QMimeData()
    mime.setUrls([QUrl.fromLocalFile(str(image_path))])
    QApplication.clipboard().setMimeData(mime)

    canvas = ChemusonCanvas()
    canvas.resize(1200, 900)
    canvas.show()
    QApplication.processEvents()

    atom_a = canvas.model.add_atom("C", 420.0, 420.0)
    atom_b = canvas.model.add_atom("N", 480.0, 420.0)
    canvas.model.add_bond(atom_a.id, atom_b.id)
    canvas._rebuild_items_from_model()
    canvas.set_current_tool("tool_bond")
    QApplication.processEvents()

    canvas._last_scene_pos = QPointF(540.0, 540.0)
    canvas.paste_from_clipboard()
    QApplication.processEvents()

    item = canvas.image_items[0]
    before_rect = item.display_rect()

    center_view = canvas.mapFromScene(before_rect.center())
    moved_view = QPoint(center_view.x() + 36, center_view.y() + 18)
    QTest.mousePress(
        canvas.viewport(),
        Qt.MouseButton.LeftButton,
        Qt.KeyboardModifier.NoModifier,
        center_view,
    )
    QTest.mouseMove(canvas.viewport(), moved_view)
    QTest.mouseRelease(
        canvas.viewport(),
        Qt.MouseButton.LeftButton,
        Qt.KeyboardModifier.NoModifier,
        moved_view,
    )
    QApplication.processEvents()

    moved_rect = item.display_rect()
    assert moved_rect.x() == pytest.approx(before_rect.x() + 36.0, abs=1.0)
    assert moved_rect.y() == pytest.approx(before_rect.y() + 18.0, abs=1.0)

    scale_scene = canvas._selection_scale_handle.scenePos()
    scale_view = canvas.mapFromScene(scale_scene)
    scaled_view = QPoint(scale_view.x() + 26, scale_view.y() + 14)
    QTest.mousePress(
        canvas.viewport(),
        Qt.MouseButton.LeftButton,
        Qt.KeyboardModifier.NoModifier,
        scale_view,
    )
    QTest.mouseMove(canvas.viewport(), scaled_view)
    QTest.mouseRelease(
        canvas.viewport(),
        Qt.MouseButton.LeftButton,
        Qt.KeyboardModifier.NoModifier,
        scaled_view,
    )
    QApplication.processEvents()

    scaled_rect = item.display_rect()
    assert scaled_rect.width() > moved_rect.width()
    assert scaled_rect.height() > moved_rect.height()
