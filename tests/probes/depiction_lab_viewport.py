from __future__ import annotations

import json
import os
from pathlib import Path
import sys

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "src"))
sys.path.insert(0, str(ROOT))
os.environ.setdefault("QT_QPA_PLATFORM", "offscreen")

from PyQt6.QtCore import QPointF
from PyQt6.QtGui import QColor, QFont, QImage, QPainter
from PyQt6.QtTest import QTest
from PyQt6.QtWidgets import QApplication

from tools.visual.depiction_lab_app.ui import DepictionLabWindow


def _events(app: QApplication) -> None:
    for _ in range(4):
        app.processEvents()
        QTest.qWait(15)


def _inspect(canvas, scene_rect) -> dict[str, int | float]:
    viewport = canvas.viewport()
    view_rect = viewport.rect()
    assert not view_rect.isEmpty()
    atom_centers = [
        canvas.mapFromScene(QPointF(atom.x, atom.y))
        for atom in canvas.model.atoms.values()
    ]
    visible_atoms = sum(view_rect.contains(point) for point in atom_centers)
    bond_centers = [
        canvas.mapFromScene(item.mapToScene(item.path().boundingRect().center()))
        for item in canvas.bond_items.values()
    ]
    visible_bonds = sum(view_rect.contains(point) for point in bond_centers)
    target = canvas.mapFromScene(scene_rect.center())
    center = view_rect.center()
    dx, dy = abs(target.x() - center.x()), abs(target.y() - center.y())

    image = viewport.grab().toImage().convertToFormat(QImage.Format.Format_ARGB32)
    dark_pixels = 0
    for y in range(image.height()):
        for x in range(image.width()):
            color = image.pixelColor(x, y)
            if max(color.red(), color.green(), color.blue()) < 170:
                dark_pixels += 1
    scale = (canvas.transform().m11() ** 2 + canvas.transform().m12() ** 2) ** 0.5
    visible_fraction = 1.0 if scale <= 1.0 else 0.3
    assert visible_atoms >= len(canvas.model.atoms) * visible_fraction, (
        visible_atoms, len(canvas.model.atoms), scale
    )
    assert visible_bonds >= len(canvas.model.bonds) * visible_fraction, (
        visible_bonds, len(canvas.model.bonds), scale
    )
    assert dark_pixels > 20, dark_pixels
    assert dx <= 10 and dy <= 10, (
        dx, dy, canvas.transform().m11(), canvas.transform().m12(),
        canvas.transform().m21(), canvas.transform().m22(),
        canvas.horizontalScrollBar().value(), canvas.horizontalScrollBar().minimum(),
        canvas.horizontalScrollBar().maximum(), canvas.scene.sceneRect(),
    )
    return {
        "atoms": visible_atoms,
        "bonds": visible_bonds,
        "dark_pixels": dark_pixels,
        "center_error_x": dx,
        "center_error_y": dy,
    }


def _save_synced_viewports(path: Path, left, right) -> None:
    left_image = left.viewport().grab().toImage()
    right_image = right.viewport().grab().toImage()
    image = QImage(
        left_image.width() + right_image.width(),
        max(left_image.height(), right_image.height()) + 38,
        QImage.Format.Format_ARGB32_Premultiplied,
    )
    image.fill(QColor("#E2E8F0"))
    painter = QPainter(image)
    painter.setFont(QFont("DejaVu Sans", 12, QFont.Weight.DemiBold))
    painter.setPen(QColor("#0F172A"))
    painter.drawText(12, 25, "CURRENT · C3 Haworth · 100% · 0°")
    painter.drawText(left_image.width() + 12, 25, "WORKING · C3 Haworth · 100% · 0°")
    painter.drawImage(0, 38, left_image)
    painter.drawImage(left_image.width(), 38, right_image)
    painter.end()
    image.save(str(path), "PNG")


def _panel_image(canvas) -> bytes:
    image = canvas.viewport().grab().toImage().convertToFormat(QImage.Format.Format_ARGB32)
    bits = image.bits()
    bits.setsize(image.sizeInBytes())
    return bytes(bits)


app = QApplication([])
window = DepictionLabWindow()
window.resize(1500, 1000)
window.show()
_events(app)

observations = []
for case_id in ("A1", "C3"):
    window.case_combo.setCurrentIndex(window.case_combo.findData(case_id))
    _events(app)
    for size in ((1650, 1080), (1120, 850)):
        window.resize(*size)
        _events(app)
        for zoom_index, angle_index in ((2, 0), (0, 0), (4, 4)):
            window.zoom_combo.setCurrentIndex(zoom_index)
            window.angle_combo.setCurrentIndex(angle_index)
            _events(app)
            current = window.current_panel.adapter
            working = window.working_panel.adapter
            assert current is not None and working is not None
            left = _inspect(current.canvas, current.scene_rect)
            right = _inspect(working.canvas, working.scene_rect)
            assert current.canvas.transform() == working.canvas.transform()
            assert left["atoms"] == right["atoms"] and left["bonds"] == right["bonds"]
            if case_id == "A1" and size == (1650, 1080) and zoom_index == 2 and angle_index == 0:
                assert _panel_image(current.canvas) == _panel_image(working.canvas)
            evidence_path = os.environ.get("CHEMUSON_VIEWPORT_EVIDENCE_PATH")
            if (
                evidence_path
                and case_id == "C3"
                and size == (1650, 1080)
                and zoom_index == 2
                and angle_index == 0
            ):
                _save_synced_viewports(Path(evidence_path), current.canvas, working.canvas)
            observations.append({"case": case_id, "size": size, "zoom_index": zoom_index, "angle_index": angle_index})

window.close()
_events(app)
print(json.dumps({"observations": len(observations), "closed": not window.isVisible()}))
