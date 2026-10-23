"""Diagnóstico de etiquetas NH2/H2 en aminas vs amidas/ureas (fuera del repo).

Imprime para cada N: etiqueta, posición H, rectángulos de N y H2 en escena,
vector N->H2 (¿apilado vertical?) y si el glifo H choca con enlaces.
"""
from __future__ import annotations

import math
import os
import sys

from PyQt6.QtCore import QRectF
from PyQt6.QtWidgets import QApplication

from chemuson.core.model import MolGraph
from chemuson.gui.canvas import ChemusonCanvas

app = QApplication.instance() or QApplication(sys.argv)

BL = 40.0


def pol(x: float, y: float, angle_deg: float, length: float = BL):
    a = math.radians(angle_deg)
    return x + length * math.cos(a), y + length * math.sin(a)


def build(atoms, bonds) -> ChemusonCanvas:
    canvas = ChemusonCanvas()
    canvas.resize(900, 700)
    canvas.show()
    app.processEvents()
    ids = {}
    for entry in atoms:
        name, element, x, y = entry[:4]
        charge = entry[4] if len(entry) > 4 else 0
        atom = canvas.model.add_atom(element, x, y, formal_charge=charge)
        ids[name] = atom.id
    for a, b, order, arom in bonds:
        canvas.model.add_bond(ids[a], ids[b], order=order, is_aromatic=arom)
    for a in atoms:
        canvas.add_atom_item(canvas.model.get_atom(ids[a[0]]))
    for a, b, _o, _r in bonds:
        canvas.add_bond_item(canvas.model.get_bond(
            next(bid for bid, bd in canvas.model.bonds.items()
                 if {bd.a1_id, bd.a2_id} == {ids[a], ids[b]})))
    canvas.refresh_atom_labels()
    app.processEvents()
    return canvas


def describe(canvas: ChemusonCanvas, label: str) -> None:
    print(f"--- {label} ---")
    for atom_id, atom in canvas.model.atoms.items():
        if atom.element != "N":
            continue
        item = canvas.atom_items[atom_id]
        pos = item.implicit_hydrogen_position
        n_item = item.label
        h_item = item.implicit_hydrogen_label
        n_rect = n_item.mapToScene(n_item.boundingRect()).boundingRect()
        n_path = n_item.mapToScene(n_item.shape()).boundingRect()
        stack = "n/a"
        if h_item.isVisible():
            h_rect = h_item.mapToScene(h_item.boundingRect()).boundingRect()
            dcx = h_rect.center().x() - n_path.center().x()
            dcy = h_rect.center().y() - n_path.center().y()
            stack = ("STACKED-VERTICAL" if abs(dcy) > abs(dcx) else "horizontal")
            overlap = any(
                h_rect.intersects(bd.mapToScene(bd.shape()).boundingRect())
                for bd in canvas.bond_items.values()
            )
            print(f"  N#{atom_id} label={item._display_label!r} hpos={pos} "
                  f"dx={dcx:+.1f} dy={dcy:+.1f} {stack} h_over_bond={overlap}")
        else:
            print(f"  N#{atom_id} label={item._display_label!r} hpos={pos} (sin glifo H)")


def snapshot(canvas: ChemusonCanvas, path: str) -> None:
    canvas.grab().save(path)


CASES = {}

# 1. Etilamina: enlace inclinado hacia arriba y hacia abajo
at = [("c1", "C", 100, 100), ("c2", "C", 140, 100), ("n", "N", *pol(140, 100, -60))]
bd = [("c1", "c2", 1, False), ("c2", "n", 1, False)]
CASES["etilamina"] = (at, bd)

# 2. Anilina: anillo aromático + N
ring = []
for i in range(6):
    a = math.radians(60 * i - 90)
    ring.append((f"r{i}", 200 + 40 * math.cos(a), 200 + 40 * math.sin(a)))
at = [("c", "C", x, y) for (_n, x, y) in ring[:0]]
at = [(f"r{i}", "C", x, y) for i, (_n, x, y) in enumerate(ring)]
at.append(("n", "N", *pol(at[0][2], at[0][3], -90)))
bd = [(f"r{i}", f"r{(i + 1) % 6}", 1, True) for i in range(6)]
bd.append(("r0", "n", 1, False))
CASES["anilina"] = (at, bd)

# 3. Dietilamina (N en vértice, dos enlaces)
at = [("c1", "C", 100, 300), ("n", "N", 140, 300), ("c2", "C", *pol(140, 300, -60))]
bd = [("c1", "n", 1, False), ("n", "c2", 1, False)]
CASES["dietilamina"] = (at, bd)

# 4. 1,1-diamino-3-fenilpropan-2-ona: N arriba y N abajo del mismo C
x1, y1 = 300, 300
n_a = pol(x1, y1, -75)
n_b = pol(x1, y1, 75)
c2 = (x1 + BL, y1)
o2 = pol(*c2, -90)
c3 = pol(*c2, 0)
at = [("c1", "C", x1, y1), ("na", "N", *n_a), ("nb", "N", *n_b),
      ("c2", "C", *c2), ("o2", "O", *o2), ("c3", "C", *c3)]
bd = [("c1", "na", 1, False), ("c1", "nb", 1, False),
      ("c1", "c2", 1, False), ("c2", "o2", 2, False), ("c2", "c3", 1, False)]
ring = []
for i in range(6):
    a = math.radians(60 * i - 90)
    ring.append((c3[0] + 40 * math.cos(a), c3[1] + 40 * math.sin(a)))
at += [(f"r{i}", "C", x, y) for i, (x, y) in enumerate(ring)]
bd += [("c3", "r0", 1, False)] + [(f"r{i}", f"r{(i + 1) % 6}", 1, True) for i in range(6)]
CASES["1,1-diamino-3-fenilpropan-2-ona"] = (at, bd)

# 5. Acetamida
at = [("c1", "C", 100, 500), ("c2", "C", 140, 500), ("o", "O", *pol(140, 500, -90)),
      ("n", "N", *pol(140, 500, 60))]
bd = [("c1", "c2", 1, False), ("c2", "o", 2, False), ("c2", "n", 1, False)]
CASES["acetamida"] = (at, bd)

# 6. N,N'-dietilurea
cu = (400, 500)
o6 = pol(*cu, -90)
n1 = pol(*cu, 150)
n2 = pol(*cu, 30)
at = [("c0", "C", *pol(*n1, -60)), ("n1", "N", *n1), ("cu", "C", *cu),
      ("o", "O", *o6), ("n2", "N", *n2), ("c3", "C", *pol(*n2, -60))]
bd = [("c0", "n1", 1, False), ("n1", "cu", 1, False), ("cu", "o", 2, False),
      ("cu", "n2", 1, False), ("n2", "c3", 1, False)]
CASES["n,n'-dietilurea"] = (at, bd)

# 7. Amoniaco cargado (control de no regresión)
at = [("c1", "C", 600, 300), ("n", "N", 640, 300, 0, 1)]
bd = [("c1", "n", 1, False)]
CASES["amonio"] = (at, bd)

tag = sys.argv[1] if len(sys.argv) > 1 else "after"
theme = os.environ.get("DIAG_THEME", "light")
if theme != "light":
    from chemuson.gui.theme import apply_theme
for name, (at, bd) in CASES.items():
    canvas = build(at, bd)
    if theme != "light":
        apply_theme(canvas, theme)
    if name == "amonio":
        nid = next(aid for aid, a in canvas.model.atoms.items() if a.element == "N")
        canvas.model.get_atom(nid).formal_charge = 1
        canvas.model.get_atom(nid).charge = 1
        canvas.refresh_atom_labels()
    describe(canvas, name)
    snapshot(canvas, f"/tmp/chemuson-amine-label/{tag}_{name}.png")
    canvas.deleteLater()
print("done")
