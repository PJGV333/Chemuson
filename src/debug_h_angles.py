#!/usr/bin/env python3
"""Script de depuración para probar ángulos de hidrógenos."""

from PyQt6.QtWidgets import QApplication
from gui.canvas import ChemusonCanvas

app = QApplication([])
canvas = ChemusonCanvas()

# Simula una cadena zigzag simple.
# Primer átomo
canvas.model.add_atom('C', 100, 100)
# Segundo átomo (30 grados abajo-derecha para el zigzag)
canvas.model.add_atom('C', 140, 120)
# Tercer átomo (30 grados arriba-derecha)
canvas.model.add_atom('C', 180, 100)

canvas.model.add_bond(1, 2)
canvas.model.add_bond(2, 3)

# Añadir ítems al canvas
for atom in canvas.model.atoms.values():
    canvas.add_atom_item(atom)
for bond in canvas.model.bonds.values():
    canvas.add_bond_item(bond)

print("\n=== Testing zigzag chain ===")
print("Atom 1 at (100,100) - terminal C")
print("Atom 2 at (140,120) - middle C with 2 bonds")
print("Atom 3 at (180,100) - terminal C")
print()

# Activar visualización de H implícitos
canvas.state.show_implicit_hydrogens = True
canvas.refresh_atom_visibility()

print("\n=== Done ===")
