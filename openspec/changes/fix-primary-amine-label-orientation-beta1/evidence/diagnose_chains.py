"""Diagnóstico (fuera del repo): aminas secundarias en cadena y percepción de
longitud de enlace. Imprime posición del glifo H, dx/dy, y longitud dibujada de
cada enlace (path recortado por etiquetas) para detectar acortamientos."""
from __future__ import annotations

import math
import sys

from PyQt6.QtWidgets import QApplication

from chemuson.gui.canvas import ChemusonCanvas

app = QApplication.instance() or QApplication(sys.argv)

BL = 40.0


def build(atoms, bonds) -> ChemusonCanvas:
    canvas = ChemusonCanvas()
    canvas.resize(900, 700)
    canvas.show()
    app.processEvents()
    ids = {}
    for name, element, x, y in atoms:
        ids[name] = canvas.model.add_atom(element, x, y, is_explicit=False).id
    for a, b, order, arom in bonds:
        canvas.model.add_bond(ids[a], ids[b], order=order, is_aromatic=arom)
    canvas._rebuild_items_from_model()
    app.processEvents()
    return canvas


def bond_path_length(item) -> float:
    lengths = []
    for polygon in item.path().toSubpathPolygons():
        lengths.append(
            sum(
                math.hypot(b.x() - a.x(), b.y() - a.y())
                for a, b in zip(polygon, polygon[1:])
            )
        )
    return max(lengths, default=0.0)


def describe(canvas: ChemusonCanvas, label: str) -> None:
    print(f"--- {label} ---")
    for atom_id, atom in canvas.model.atoms.items():
        if atom.element != "N":
            continue
        item = canvas.atom_items[atom_id]
        pos = item.implicit_hydrogen_position
        n_item, h_item = item.label, item.implicit_hydrogen_label
        n_path = n_item.mapToScene(n_item.shape()).boundingRect()
        if h_item.isVisible():
            h_rect = h_item.mapToScene(h_item.boundingRect()).boundingRect()
            dcx = h_rect.center().x() - n_path.center().x()
            dcy = h_rect.center().y() - n_path.center().y()
            stack = "VERTICAL" if abs(dcy) > abs(dcx) else "HORIZONTAL"
            print(f"  N#{atom_id} label={item._display_label!r} hpos={pos} "
                  f"dx={dcx:+.1f} dy={dcy:+.1f} {stack}")
        else:
            print(f"  N#{atom_id} label={item._display_label!r} sin glifo H")
    for bond in canvas.model.bonds.values():
        item = canvas.bond_items[bond.id]
        e1 = canvas.model.get_atom(bond.a1_id).element
        e2 = canvas.model.get_atom(bond.a2_id).element
        print(f"  bond {bond.a1_id}({e1})-{bond.a2_id}({e2}) "
              f"drawn={bond_path_length(item):.1f}px "
              f"shrink={item._label_shrink_start + item._label_shrink_end:.1f}")


def snapshot(canvas: ChemusonCanvas, path: str) -> None:
    from PyQt6.QtCore import QRectF, Qt
    from PyQt6.QtGui import QImage, QPainter

    bounds = canvas.scene.itemsBoundingRect().adjusted(-20, -20, 20, 20)
    image = QImage(800, 800, QImage.Format.Format_ARGB32)
    image.fill(0xFFFFFFFF)
    painter = QPainter(image)
    painter.setRenderHint(QPainter.RenderHint.Antialiasing)
    canvas.scene.render(painter, QRectF(0, 0, 800, 800), bounds, Qt.AspectRatioMode.KeepAspectRatio)
    painter.end()
    image.save(path)


CASES = {}

# zigzag con N en vértice superior (dietilamina en cadena)
chain = [(0, 0), (34.6, 20), (69.3, 0), (103.9, 20), (138.6, 0)]
atoms = [(f"c{i}", "C", x, y) for i, (x, y) in enumerate(chain)]
atoms[2] = ("n2", "N", *chain[2])
bonds = [(f"c{i}", f"c{i+1}", 1, False) for i in range(4)]
bonds[1] = ("n2", "c1", 1, False)
bonds[2] = ("n2", "c3", 1, False)
CASES["zigzag-sec-amine"] = (atoms, bonds)

# segmento horizontal recto a través de N
atoms = [("c0", "C", 0.0, 0.0), ("n1", "N", 40.0, 0.0), ("c2", "C", 80.0, 0.0)]
bonds = [("c0", "n1", 1, False), ("n1", "c2", 1, False)]
CASES["linear-horizontal"] = (atoms, bonds)

# dietales con N arriba (V invertida)
atoms = [("c1", "C", 60.0, 120.0), ("n", "N", 100.0, 80.0), ("c2", "C", 140.0, 120.0)]
bonds = [("c1", "n", 1, False), ("n", "c2", 1, False)]
CASES["diethylamine-V"] = (atoms, bonds)

# poliamina lineal: H2N-CH2-CH2-NH-CH2-CH2-NH2 (zigzag)
pts = [(i * 34.6, 20.0 if i % 2 else 0.0) for i in range(7)]
atoms = [(f"p{i}", "C", x, y) for i, (x, y) in enumerate(pts)]
atoms[0] = ("n0", "N", *pts[0])
atoms[3] = ("n3", "N", *pts[3])
atoms[6] = ("n6", "N", *pts[6])
bonds = [(f"p{i}", f"p{i+1}", 1, False) for i in range(6)]
bonds[0] = ("n0", "p1", 1, False)
bonds[2] = ("p2", "n3", 1, False)
bonds[3] = ("n3", "p4", 1, False)
bonds[5] = ("p5", "n6", 1, False)
CASES["polyamine-linear"] = (atoms, bonds)

# N-metilanilina: N secundaria sobre anillo
center = (300.0, 300.0)
ring = [
    (center[0] + 40 * math.cos(math.radians(60 * i - 90)),
     center[1] + 40 * math.sin(math.radians(60 * i - 90)))
    for i in range(6)
]
atoms = [(f"r{i}", "C", x, y) for i, (x, y) in enumerate(ring)]
bonds = [(f"r{i}", f"r{(i+1)%6}", 1, True) for i in range(6)]
nx, ny = ring[0]
nx, ny = center[0] + 2 * (nx - center[0]), center[1] + 2 * (ny - center[1])
atoms.append(("n", "N", nx, ny))
bonds.append(("r0", "n", 1, False))
mx, my = nx + 40, ny + 20
atoms.append(("me", "C", mx, my))
bonds.append(("n", "me", 1, False))
CASES["n-methylaniline"] = (atoms, bonds)

tag = sys.argv[1] if len(sys.argv) > 1 else "before"
theme = __import__("os").environ.get("DIAG_THEME", "light")
if __name__ == "__main__":
    if theme != "light":
        from chemuson.gui.theme import apply_theme
    for name, (at, bd) in CASES.items():
        canvas = build(at, bd)
        if theme != "light":
            apply_theme(canvas, theme)
        describe(canvas, name)
        snapshot(canvas, f"/tmp/chemuson-amine2/{tag}_{name}.png")
        canvas.deleteLater()
    print("done")
