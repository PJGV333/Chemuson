from __future__ import annotations

import json
import math
import os
from pathlib import Path

from PyQt6.QtCore import QPointF, QRectF, Qt
from PyQt6.QtGui import QBrush, QColor, QImage, QPainter
from PyQt6.QtTest import QTest
from PyQt6.QtWidgets import QApplication

from chemuson.chemio.rdkit_io import molgraph_to_smiles_isolated_or_error
from chemuson.gui.canvas import ChemusonCanvas
from chemuson.gui.geom import angle_deg, angle_distance_deg, endpoint_from_angle_len

OUT = Path(os.environ.get("INITIAL_2D_OUT", Path(__file__).parent / "after"))
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


def _cross(p: QPointF, q: QPointF, r: QPointF) -> float:
    return (q.x() - p.x()) * (r.y() - p.y()) - (q.y() - p.y()) * (r.x() - p.x())


def _non_adjacent_crossings() -> int:
    crossings = 0
    bonds = list(canvas.model.bonds.values())
    for index, first in enumerate(bonds):
        first_ids = {first.a1_id, first.a2_id}
        p0 = canvas.model.get_atom(first.a1_id)
        p1 = canvas.model.get_atom(first.a2_id)
        a, b = QPointF(p0.x, p0.y), QPointF(p1.x, p1.y)
        for second in bonds[index + 1 :]:
            if first_ids & {second.a1_id, second.a2_id}:
                continue
            q0 = canvas.model.get_atom(second.a1_id)
            q1 = canvas.model.get_atom(second.a2_id)
            c, d = QPointF(q0.x, q0.y), QPointF(q1.x, q1.y)
            if (
                _cross(a, b, c) * _cross(a, b, d) < -1e-8
                and _cross(c, d, a) * _cross(c, d, b) < -1e-8
            ):
                crossings += 1
    return crossings


def drag(start: QPointF, end: QPointF) -> None:
    p0 = canvas.mapFromScene(start)
    p1 = canvas.mapFromScene(end)
    if not canvas.viewport().rect().contains(p0) or not canvas.viewport().rect().contains(p1):
        raise RuntimeError(f"gesture outside viewport: {p0}, {p1}")
    QTest.mousePress(canvas.viewport(), Qt.MouseButton.LeftButton, pos=p0)
    QTest.mouseMove(canvas.viewport(), p1, 10)
    QTest.mouseRelease(canvas.viewport(), Qt.MouseButton.LeftButton, pos=p1)
    app.processEvents()


def short_bond(anchor_id: int, *, element: str = "C", order: int = 1) -> int:
    canvas.set_current_tool("tool_bond")
    canvas.state.default_element = element
    canvas.state.active_bond_order = order
    canvas.state.active_bond_aromatic = False
    atom = canvas.model.get_atom(anchor_id)
    click = canvas.mapFromScene(QPointF(atom.x, atom.y))
    before = set(canvas.model.atoms)
    QTest.mouseClick(canvas.viewport(), Qt.MouseButton.LeftButton, pos=click)
    app.processEvents()
    created = set(canvas.model.atoms) - before
    if len(created) != 1:
        raise RuntimeError(f"short click at {anchor_id} created {created}")
    return created.pop()


def clear_transient_items() -> None:
    canvas.state.selected_atoms.clear()
    canvas.state.selected_bonds.clear()
    for item in canvas.atom_items.values():
        item.setSelected(False)
    for item in canvas.bond_items.values():
        item.setSelected(False)
    canvas._update_hover(QPointF(-10000.0, -10000.0))


def capture(name: str, source: QRectF) -> None:
    clear_transient_items()
    backgrounds = {"light": "#F1F5F9", "dark": "#111827"}
    for theme, background in backgrounds.items():
        canvas.scene.setBackgroundBrush(QBrush(QColor(background)))
        for scale in (1.0, 1.25):
            width = round(source.width() * scale)
            height = round(source.height() * scale)
            image = QImage(width, height, QImage.Format.Format_ARGB32_Premultiplied)
            image.fill(QColor(background))
            painter = QPainter(image)
            painter.setRenderHint(QPainter.RenderHint.Antialiasing, True)
            painter.setRenderHint(QPainter.RenderHint.TextAntialiasing, True)
            canvas.scene.render(
                painter,
                QRectF(0.0, 0.0, width, height),
                source,
                Qt.AspectRatioMode.IgnoreAspectRatio,
            )
            painter.end()
            if not image.save(str(OUT / f"{name}-{theme}-{scale:g}x.png")):
                raise RuntimeError(f"could not save {name} {theme} {scale}x image")


# Draw the regular aromatic ring with a deliberately rotated scene orientation.
canvas.set_current_tool("tool_ring")
center = canvas.mapToScene(canvas.viewport().rect().center())
cursor = endpoint_from_angle_len(center, 17.0, canvas.state.bond_length)
drag(center, cursor)
ring_ids = sorted(canvas.model.atoms)
ring_anchor = max((canvas.model.get_atom(atom_id) for atom_id in ring_ids), key=lambda atom: atom.x)
ring_center = QPointF(
    sum(canvas.model.get_atom(atom_id).x for atom_id in ring_ids) / 6.0,
    sum(canvas.model.get_atom(atom_id).y for atom_id in ring_ids) / 6.0,
)
ring_points = [
    QPointF(canvas.model.get_atom(atom_id).x, canvas.model.get_atom(atom_id).y)
    for atom_id in ring_ids
]
ring_lengths = [
    math.hypot(
        ring_points[(index + 1) % 6].x() - point.x(),
        ring_points[(index + 1) % 6].y() - point.y(),
    )
    for index, point in enumerate(ring_points)
]
ring_angles = []
for index, point in enumerate(ring_points):
    previous = ring_points[(index - 1) % 6]
    following = ring_points[(index + 1) % 6]
    ring_angles.append(
        angle_distance_deg(angle_deg(point, previous), angle_deg(point, following))
    )

source = QRectF(250.0, center.y() - 150.0, 600.0, 300.0)
capture("benzene-isolated", source)

# Continue with the actual short-click auto-assist, no Clean2D invocation.
phenyl_methylene = short_bond(ring_anchor.id)
capture("benzene-substituent", source)
carbonyl_carbon = short_bond(phenyl_methylene)
oxygen = short_bond(carbonyl_carbon, element="O", order=2)
alpha_methylene = short_bond(carbonyl_carbon)
terminal_methyl = short_bond(alpha_methylene)
canvas.state.default_element = "C"

atoms = [canvas.model.get_atom(atom_id) for atom_id in sorted(canvas.model.atoms)]
bonds = sorted(canvas.model.bonds.values(), key=lambda bond: bond.id)
local_angles: dict[str, list[float]] = {}
for atom_id in (ring_anchor.id, phenyl_methylene, carbonyl_carbon, alpha_methylene):
    atom = canvas.model.get_atom(atom_id)
    directions = []
    for bond in bonds:
        if atom_id not in (bond.a1_id, bond.a2_id):
            continue
        other_id = bond.a2_id if bond.a1_id == atom_id else bond.a1_id
        other = canvas.model.get_atom(other_id)
        directions.append(angle_deg(QPointF(atom.x, atom.y), QPointF(other.x, other.y)))
    local_angles[str(atom_id)] = sorted(
        angle_distance_deg(first, second)
        for index, first in enumerate(directions)
        for second in directions[index + 1 :]
    )

smiles = molgraph_to_smiles_isolated_or_error(canvas.model, timeout_s=2.0)
if smiles != "CCC(=O)Cc1ccccc1":
    raise RuntimeError(f"unexpected structure identity: {smiles}")
if len(atoms) != 11 or len(bonds) != 11:
    raise RuntimeError(f"unexpected graph size: {len(atoms)} atoms, {len(bonds)} bonds")
crossings = _non_adjacent_crossings()
if crossings != 0:
    raise RuntimeError(f"unexpected non-adjacent bond crossings: {crossings}")

snapshot = {
    "smiles": smiles,
    "atoms": [
        {"id": atom.id, "element": atom.element, "x": atom.x, "y": atom.y}
        for atom in atoms
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
        for bond in bonds
    ],
    "metrics": {
        "non_adjacent_crossings": crossings,
        "ring_side_lengths": ring_lengths,
        "ring_internal_angles_deg": ring_angles,
        "local_pair_angles_deg": local_angles,
        "bond_lengths": {
            str(bond.id): math.hypot(
                canvas.model.get_atom(bond.a1_id).x - canvas.model.get_atom(bond.a2_id).x,
                canvas.model.get_atom(bond.a1_id).y - canvas.model.get_atom(bond.a2_id).y,
            )
            for bond in bonds
        },
    },
}
(OUT / "1-phenylbutan-2-one-model.json").write_text(
    json.dumps(snapshot, indent=2) + "\n",
    encoding="utf-8",
)
(OUT / "1-phenylbutan-2-one-metrics.txt").write_text(
    json.dumps({"atoms": len(atoms), "bonds": len(bonds), **snapshot["metrics"]}, indent=2)
    + "\n",
    encoding="utf-8",
)
capture("1-phenylbutan-2-one", source)
print(json.dumps({"smiles": smiles, "atoms": len(atoms), "bonds": len(bonds), **snapshot["metrics"]}, indent=2))
canvas.close()
app.processEvents()
