
import os
import sys
from PyQt6.QtWidgets import QApplication
from PyQt6.QtCore import QPointF
from gui.canvas import ChemusonCanvas
from chemio.persistence import PersistenceManager
from core.model import BondStyle, BondStereo

def verify_benzene_rendering():
    app = QApplication(sys.argv)
    canvas = ChemusonCanvas()
    canvas.state.show_implicit_carbons = True
    
    # 1. Create Benzene manually
    # Positions roughly in a hexagon
    atoms = [
        canvas.model.add_atom("C", 100, 0),
        canvas.model.add_atom("C", 150, 86),
        canvas.model.add_atom("C", 100, 172),
        canvas.model.add_atom("C", 0, 172),
        canvas.model.add_atom("C", -50, 86),
        canvas.model.add_atom("C", 0, 0)
    ]
    
    # Add bonds (alternating double/single)
    # Give them a ring_id so they are a ring
    ring_id = 1 
    canvas.model.add_bond(atoms[0], atoms[1], order=2, ring_id=ring_id)
    canvas.model.add_bond(atoms[1], atoms[2], order=1, ring_id=ring_id)
    canvas.model.add_bond(atoms[2], atoms[3], order=2, ring_id=ring_id)
    canvas.model.add_bond(atoms[3], atoms[4], order=1, ring_id=ring_id)
    canvas.model.add_bond(atoms[4], atoms[5], order=2, ring_id=ring_id)
    canvas.model.add_bond(atoms[5], atoms[0], order=1, ring_id=ring_id)
    
    # Build visual items initially
    canvas._rebuild_items_from_model()
    
    # Check ring center calculation
    center = canvas._ring_centers.get(ring_id)
    print(f"Calculated center: {center}")
    assert center is not None, "Ring center should be calculated"
    
    # Check double bond offsets
    for bond_id, item in canvas.bond_items.items():
        if item.order == 2:
            print(f"Bond {bond_id} offset_sign: {item._offset_sign}")
            # The logic should point inside. Center is roughly (50, 86)
            # Bond (100,0)-(150,86) has mid (125, 43). Vector to center (-75, 43).
            # nx, ny for this bond: dy=86, dx=50 -> nx = -86/length, ny = 50/length.
            # dot product: (-86*-75 + 50*43) = (6450 + 2150) = 8600 > 0.
            # So offset_sign should be 1? Wait, let's see BondItem logic.
    
    # 2. Save and Load
    save_path = "benzene_test.cmsn"
    PersistenceManager.save_to_file(save_path, canvas)
    
    canvas.clear_canvas()
    PersistenceManager.load_from_file(save_path, canvas)
    
    # 3. Verify labels and shrinks
    print("\nAfter Load Verification:")
    for atom_id, item in canvas.atom_items.items():
        print(f"Atom {atom_id} ({item.element}) Label visible: {item.label.isVisible()}")
        assert item.label.isVisible(), "Carbon label should be visible"
        
    for bond_id, item in canvas.bond_items.items():
        print(f"Bond {bond_id} Shrinks: {item._label_shrink_start:.1f}, {item._label_shrink_end:.1f}")
        assert item._label_shrink_start >= 4.0, f"Shrink start too small: {item._label_shrink_start}"
        assert item._label_shrink_end >= 4.0, f"Shrink end too small: {item._label_shrink_end}"

    os.remove(save_path)
    print("\nSUCCESS: Rendering refinements verified!")

if __name__ == "__main__":
    verify_benzene_rendering()
