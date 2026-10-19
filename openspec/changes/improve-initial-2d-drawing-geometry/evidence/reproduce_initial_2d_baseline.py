from __future__ import annotations

import json
import math
import os
from pathlib import Path

from PyQt6.QtCore import QPointF, QRectF, Qt
from PyQt6.QtGui import QImage, QPainter
from PyQt6.QtTest import QTest
from PyQt6.QtWidgets import QApplication

from chemuson.gui.canvas import ChemusonCanvas
from chemuson.gui.geom import angle_deg, angle_distance_deg, endpoint_from_angle_len

OUT = Path(os.environ.get("INITIAL_2D_OUT", "/tmp/chemuson-initial-2d-baseline"))
OUT.mkdir(parents=True, exist_ok=True)
app = QApplication.instance() or QApplication([])
canvas = ChemusonCanvas()
canvas.resize(1000, 760)
canvas.show()
app.processEvents()
canvas.center_on_paper()
app.processEvents()
canvas.state.fixed_angles = True
canvas.state.fixed_lengths = True
canvas.state.angle_step_deg = 30
canvas.state.active_ring_size = 6
canvas.state.active_ring_aromatic = True


def drag(start: QPointF, end: QPointF) -> None:
    p0 = canvas.mapFromScene(start)
    p1 = canvas.mapFromScene(end)
    if not canvas.viewport().rect().contains(p0) or not canvas.viewport().rect().contains(p1):
        raise RuntimeError(f"gesture outside viewport: {p0}, {p1}")
    QTest.mousePress(canvas.viewport(), Qt.MouseButton.LeftButton, pos=p0)
    QTest.mouseMove(canvas.viewport(), p1, 10)
    QTest.mouseRelease(canvas.viewport(), Qt.MouseButton.LeftButton, pos=p1)
    app.processEvents()


def directed_bond(atom_id: int, theta_deg: float, order: int = 1) -> int:
    canvas.state.active_bond_order = order
    canvas.state.active_bond_aromatic = False
    canvas.state.default_element = "O" if order == 2 else "C"
    canvas.set_current_tool("tool_bond")
    atom = canvas.model.get_atom(atom_id)
    p0 = QPointF(atom.x, atom.y)
    p1 = endpoint_from_angle_len(p0, theta_deg, canvas.state.bond_length)
    before = set(canvas.model.atoms)
    drag(p0, p1)
    created = set(canvas.model.atoms) - before
    if len(created) != 1:
        raise RuntimeError(f"expected one atom after directed bond, got {created}")
    return created.pop()

# Isolated benceno by the actual ring-tool gesture.
canvas.set_current_tool("tool_ring")
center = QPointF(400.0, 300.0)
drag(center, QPointF(440.0, 300.0))
ring_ids = sorted(canvas.model.atoms)
ring_atoms = [canvas.model.get_atom(atom_id) for atom_id in ring_ids]


def save_capture(name: str, scene_source: QRectF) -> None:
    canvas.state.selected_atoms.clear()
    canvas.state.selected_bonds.clear()
    for item in canvas.atom_items.values():
        item.setSelected(False)
    for item in canvas.bond_items.values():
        item.setSelected(False)
    canvas._update_hover(QPointF(-10000.0, -10000.0))
    for scale in (1.0, 1.25):
        width = round(scene_source.width() * scale)
        height = round(scene_source.height() * scale)
        image = QImage(width, height, QImage.Format.Format_ARGB32_Premultiplied)
        image.fill(Qt.GlobalColor.white)
        painter = QPainter(image)
        canvas.scene.render(painter, QRectF(0, 0, width, height), scene_source)
        painter.end()
        image.save(str(OUT / f"{name}-{scale:g}x.png"))


save_capture("before-benzene-isolated", QRectF(340.0, 240.0, 120.0, 120.0))

# Ph-CH2-C(=O)-CH2-CH3, each bond created by a directed canvas gesture;
# no coordinate values are read from or inferred from the raster output.
ring_anchor = max(ring_atoms, key=lambda atom: atom.x)
chain_ids = [directed_bond(ring_anchor.id, 0.0)]
save_capture("before-benzene-substituent", QRectF(340.0, 240.0, 180.0, 120.0))
chain_ids.append(directed_bond(chain_ids[-1], 60.0))  # carbonyl C
chain_ids.append(directed_bond(chain_ids[-1], 120.0, order=2))  # carbonyl O
chain_ids.append(directed_bond(chain_ids[-2], 0.0))  # alpha CH2
chain_ids.append(directed_bond(chain_ids[-1], 60.0))  # terminal CH3

all_atoms = [canvas.model.get_atom(atom_id) for atom_id in sorted(canvas.model.atoms)]
ring_positions = [(canvas.model.get_atom(i).x, canvas.model.get_atom(i).y) for i in ring_ids]
ring_lengths = [
    math.dist(ring_positions[i], ring_positions[(i + 1) % len(ring_positions)])
    for i in range(len(ring_positions))
]
ring_angles = []
for i, (x, y) in enumerate(ring_positions):
    p = QPointF(x, y)
    prev = QPointF(*ring_positions[(i - 1) % len(ring_positions)])
    nxt = QPointF(*ring_positions[(i + 1) % len(ring_positions)])
    a = angle_deg(p, prev)
    b = angle_deg(p, nxt)
    ring_angles.append(min(angle_distance_deg(a, b), 360 - angle_distance_deg(a, b)))


def atom_angle(atom_id: int) -> list[float]:
    atom = canvas.model.get_atom(atom_id)
    neighbors = []
    for bond in canvas.model.bonds.values():
        if bond.a1_id == atom_id:
            other = canvas.model.get_atom(bond.a2_id)
        elif bond.a2_id == atom_id:
            other = canvas.model.get_atom(bond.a1_id)
        else:
            continue
        neighbors.append(angle_deg(QPointF(atom.x, atom.y), QPointF(other.x, other.y)))
    return sorted(
        round(angle_distance_deg(a, b), 8)
        for idx, a in enumerate(neighbors)
        for b in neighbors[idx + 1 :]
    )

scene_source = QRectF(280.0, 180.0, 360.0, 250.0)
save_capture("before-1-phenylbutan-2-one", scene_source)

snapshot = {
    "atoms": [
        {"id": atom.id, "element": atom.element, "x": atom.x, "y": atom.y}
        for atom in all_atoms
    ],
    "bonds": [
        {
            "id": bond.id,
            "a1": bond.a1_id,
            "a2": bond.a2_id,
            "order": bond.order,
            "aromatic": bond.is_aromatic,
            "ring_id": bond.ring_id,
        }
        for bond in canvas.model.bonds.values()
    ],
    "metrics": {
        "ring_side_lengths": ring_lengths,
        "ring_internal_angles_deg": ring_angles,
        "local_pair_angles_deg": {str(atom_id): atom_angle(atom_id) for atom_id in chain_ids},
        "bond_lengths": {
            str(bond.id): math.hypot(
                canvas.model.get_atom(bond.a1_id).x - canvas.model.get_atom(bond.a2_id).x,
                canvas.model.get_atom(bond.a1_id).y - canvas.model.get_atom(bond.a2_id).y,
            )
            for bond in canvas.model.bonds.values()
        },
    },
}
(OUT / "before-1-phenylbutan-2-one.json").write_text(json.dumps(snapshot, indent=2) + "\n")
print(json.dumps({"atoms": len(canvas.model.atoms), "bonds": len(canvas.model.bonds), **snapshot["metrics"]}, indent=2))
canvas.close()
app.processEvents()
