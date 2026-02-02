from __future__ import annotations
import json
import os
from typing import Any, Dict, TYPE_CHECKING
from core.model import MolGraph, Atom, Bond, BondStyle, BondStereo

if TYPE_CHECKING:
    from gui.canvas import ChemusonCanvas

class PersistenceManager:
    """Manages saving and loading of Chemuson (.cmsn) files."""
    
    VERSION = "0.1.0"

    @staticmethod
    def save_to_dict(canvas: 'ChemusonCanvas') -> Dict[str, Any]:
        """Collects all data from the canvas and model into a serializable dictionary."""
        graph = canvas.model
        
        # 1. Serialize MolGraph
        atoms_data = []
        for atom in graph.atoms.values():
            atoms_data.append({
                "id": atom.id,
                "element": atom.element,
                "x": atom.x,
                "y": atom.y,
                "charge": atom.charge,
                "isotope": atom.isotope,
                "explicit_h": atom.explicit_h,
                "mapping": atom.mapping,
                "is_query": atom.is_query,
                "is_explicit": atom.is_explicit
            })
            
        bonds_data = []
        for bond in graph.bonds.values():
            bonds_data.append({
                "id": bond.id,
                "a1_id": bond.a1_id,
                "a2_id": bond.a2_id,
                "order": bond.order,
                "style": bond.style.value,
                "stereo": bond.stereo.value,
                "is_aromatic": bond.is_aromatic,
                "display_order": bond.display_order,
                "is_query": bond.is_query,
                "ring_id": bond.ring_id,
                "length_px": bond.length_px
            })
            
        # 2. Serialize Canvas Items (Arrows, Brackets, Text)
        canvas_data = canvas.get_persistence_data()
        
        # 3. Combine everything
        return {
            "application": "Chemuson",
            "version": PersistenceManager.VERSION,
            "model": {
                "atoms": atoms_data,
                "bonds": bonds_data,
                "_next_atom_id": graph._next_atom_id,
                "_next_bond_id": graph._next_bond_id
            },
            "canvas": canvas_data
        }

    @staticmethod
    def load_from_dict(data: Dict[str, Any], canvas: 'ChemusonCanvas') -> None:
        """Restores the canvas and model from a dictionary."""
        if data.get("application") != "Chemuson":
            raise ValueError("Not a valid Chemuson file")
            
        # 1. Restore MolGraph
        canvas.model.clear()
        model_data = data.get("model", {})
        
        for atom_d in model_data.get("atoms", []):
            canvas.model.add_atom(
                element=atom_d["element"],
                x=atom_d["x"],
                y=atom_d["y"],
                atom_id=atom_d["id"],
                charge=atom_d.get("charge", 0),
                isotope=atom_d.get("isotope"),
                explicit_h=atom_d.get("explicit_h"),
                mapping=atom_d.get("mapping"),
                is_query=atom_d.get("is_query", False),
                is_explicit=atom_d.get("is_explicit", False)
            )
            
        for bond_d in model_data.get("bonds", []):
            canvas.model.add_bond(
                a1_id=bond_d["a1_id"],
                a2_id=bond_d["a2_id"],
                order=bond_d.get("order", 1),
                bond_id=bond_d["id"],
                style=BondStyle(bond_d.get("style", "plain")),
                stereo=BondStereo(bond_d.get("stereo", "none")),
                is_aromatic=bond_d.get("is_aromatic", False),
                display_order=bond_d.get("display_order"),
                is_query=bond_d.get("is_query", False),
                ring_id=bond_d.get("ring_id"),
                length_px=bond_d.get("length_px")
            )
            
        canvas.model._next_atom_id = model_data.get("_next_atom_id", canvas.model._next_atom_id)
        canvas.model._next_bond_id = model_data.get("_next_bond_id", canvas.model._next_bond_id)
        
        # 2. Restore Canvas Items and Settings
        canvas.load_persistence_data(data.get("canvas", {}))
        
        # 3. Visual Reconstruction of molecules
        canvas._rebuild_items_from_model()

    @staticmethod
    def save_to_file(filepath: str, canvas: 'ChemusonCanvas') -> None:
        data = PersistenceManager.save_to_dict(canvas)
        with open(filepath, 'w', encoding='utf-8') as f:
            json.dump(data, f, indent=2)

    @staticmethod
    def load_from_file(filepath: str, canvas: 'ChemusonCanvas') -> None:
        with open(filepath, 'r', encoding='utf-8') as f:
            data = json.load(f)
        PersistenceManager.load_from_dict(data, canvas)
