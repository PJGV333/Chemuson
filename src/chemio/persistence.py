"""Persistencia del estado del editor en archivos `.cmsn`.

Este módulo serializa y deserializa el modelo químico y los elementos
de la interfaz para reconstruir el lienzo al abrir un archivo.
"""

from __future__ import annotations
import json
import os
from typing import Any, Dict, TYPE_CHECKING
from core.model import MolGraph, Atom, Bond, BondStyle, BondStereo

if TYPE_CHECKING:
    from gui.canvas import ChemusonCanvas

class PersistenceManager:
    """Gestiona el guardado y carga de archivos `.cmsn` de Chemuson."""
    
    VERSION = "0.1.0"

    @staticmethod
    def save_to_dict(canvas: 'ChemusonCanvas') -> Dict[str, Any]:
        """Serializa el estado del canvas y su modelo en un diccionario.

        Args:
            canvas: Instancia activa del lienzo de Chemuson.

        Returns:
            Diccionario serializable con información de modelo y GUI.

        Side Effects:
            No tiene efectos laterales; solo lee el estado del canvas.
        """
        graph = canvas.model
        
        # 1. Serializar MolGraph (átomos y enlaces)
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
                "length_px": bond.length_px,
                "stroke_px": bond.stroke_px,
                "color": bond.color,
            })
            
        # 2. Serializar elementos del canvas (flechas, brackets, texto)
        canvas_data = canvas.get_persistence_data()
        
        # 3. Combinar todo en un único documento
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
        """Restaura el estado del canvas y el modelo desde un diccionario.

        Args:
            data: Diccionario de estado (resultado de `save_to_dict`).
            canvas: Instancia de lienzo donde se cargará la información.

        Raises:
            ValueError: Si el archivo no corresponde a Chemuson.

        Side Effects:
            Modifica el modelo y reconstruye la vista del canvas.
        """
        if data.get("application") != "Chemuson":
            raise ValueError("Not a valid Chemuson file")
            
        # 1. Restaurar MolGraph
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
                length_px=bond_d.get("length_px"),
                stroke_px=bond_d.get("stroke_px"),
                color=bond_d.get("color"),
            )
            
        canvas.model._next_atom_id = model_data.get("_next_atom_id", canvas.model._next_atom_id)
        canvas.model._next_bond_id = model_data.get("_next_bond_id", canvas.model._next_bond_id)
        
        # 2. Restaurar elementos del canvas y preferencias visuales
        canvas.load_persistence_data(data.get("canvas", {}))
        
        # 3. Reconstrucción visual de moléculas en el lienzo
        canvas._rebuild_items_from_model()

    @staticmethod
    def save_to_file(filepath: str, canvas: 'ChemusonCanvas') -> None:
        """Guarda el estado del canvas en un archivo `.cmsn`.

        Args:
            filepath: Ruta de destino.
            canvas: Instancia del lienzo a serializar.

        Side Effects:
            Escribe en disco el archivo indicado.
        """
        data = PersistenceManager.save_to_dict(canvas)
        with open(filepath, 'w', encoding='utf-8') as f:
            json.dump(data, f, indent=2)

    @staticmethod
    def load_from_file(filepath: str, canvas: 'ChemusonCanvas') -> None:
        """Carga un archivo `.cmsn` y restaura el lienzo.

        Args:
            filepath: Ruta del archivo de entrada.
            canvas: Instancia del lienzo donde se cargará el estado.

        Side Effects:
            Lee desde disco y modifica el canvas.
        """
        with open(filepath, 'r', encoding='utf-8') as f:
            data = json.load(f)
        PersistenceManager.load_from_dict(data, canvas)
