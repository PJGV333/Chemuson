"""Persistencia del estado del editor en archivos `.cmsn`.

Este módulo serializa y deserializa el modelo químico y los elementos
de la interfaz para reconstruir el lienzo al abrir un archivo.
"""

from __future__ import annotations
import json
import os
from typing import Any, Dict, TYPE_CHECKING
from chemuson.core.model import MolGraph, Atom, Bond, BondStyle, BondStereo

if TYPE_CHECKING:
    from chemuson.gui.canvas import ChemusonCanvas

class PersistenceManager:
    """Gestiona el guardado y carga de archivos `.cmsn` de Chemuson."""
    
    VERSION = "0.1.0"

    @staticmethod
    def _parse_bond_style(bond_data: Dict[str, Any]) -> BondStyle:
        """Resuelve el estilo de enlace aceptando llaves nuevas y legadas."""
        raw_style = bond_data.get("style")
        if raw_style is None:
            raw_style = bond_data.get("type", BondStyle.PLAIN.value)
        try:
            return BondStyle(raw_style)
        except Exception:
            return BondStyle.PLAIN

    @staticmethod
    def _parse_bond_stereo(bond_data: Dict[str, Any]) -> BondStereo:
        """Resuelve estereoquímica de enlace con fallback seguro."""
        raw_stereo = bond_data.get("stereo", BondStereo.NONE.value)
        try:
            return BondStereo(raw_stereo)
        except Exception:
            return BondStereo.NONE

    @staticmethod
    def _parse_pi_offset_sign(bond_data: Dict[str, Any]) -> int | None:
        """Normaliza orientación manual de dobles enlaces (+1/-1)."""
        raw_sign = bond_data.get("pi_offset_sign")
        try:
            parsed = int(raw_sign)
        except Exception:
            return None
        return parsed if parsed in {-1, 1} else None

    @staticmethod
    def _infer_coordination_donor(
        graph: MolGraph,
        a1_id: int,
        a2_id: int,
        preferred: Any = None,
    ) -> int:
        """Infere donador para enlaces coordinativos cuando no está serializado."""
        try:
            preferred_id = int(preferred) if preferred is not None else None
        except Exception:
            preferred_id = None
        if preferred_id in {a1_id, a2_id}:
            return preferred_id
        atom1 = graph.atoms.get(a1_id)
        atom2 = graph.atoms.get(a2_id)
        a1_is_center = bool(getattr(atom1, "is_coordination_center", False))
        a2_is_center = bool(getattr(atom2, "is_coordination_center", False))
        if a1_is_center and not a2_is_center:
            return a2_id
        if a2_is_center and not a1_is_center:
            return a1_id
        return a1_id

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
                "formal_charge": int(getattr(atom, "formal_charge", atom.charge)),
                "isotope": atom.isotope,
                "radical_electrons": int(getattr(atom, "radical_electrons", 0) or 0),
                "oxidation_state": getattr(atom, "oxidation_state", None),
                "stereo_cip": getattr(atom, "stereo_cip", None),
                "stereo_axial": getattr(atom, "stereo_axial", None),
                "stereo_helical": getattr(atom, "stereo_helical", None),
                "stereo_si_re": getattr(atom, "stereo_si_re", None),
                "explicit_h": atom.explicit_h,
                "group_h_cap": getattr(atom, "group_h_cap", None),
                "mapping": atom.mapping,
                "is_query": atom.is_query,
                "is_explicit": atom.is_explicit,
                "no_implicit": bool(getattr(atom, "no_implicit", False)),
                "label_scale": getattr(atom, "label_scale", None),
                "is_coordination_center": getattr(atom, "is_coordination_center", False),
                "sphere_radius": getattr(atom, "sphere_radius", None),
                "sphere_color": getattr(atom, "sphere_color", None),
                "sphere_filled": bool(getattr(atom, "sphere_filled", True)),
                "sphere_transparent": bool(getattr(atom, "sphere_transparent", False)),
            })
            
        bonds_data = []
        for bond in graph.bonds.values():
            bonds_data.append({
                "id": bond.id,
                "a1_id": bond.a1_id,
                "a2_id": bond.a2_id,
                "order": bond.order,
                "style": bond.style.value,
                "type": bond.style.value,
                "stereo": bond.stereo.value,
                "stereo_ez": getattr(bond, "stereo_ez", None),
                "stereo_axial": getattr(bond, "stereo_axial", None),
                "stereo_endo_exo": getattr(bond, "stereo_endo_exo", None),
                "stereo_helical": getattr(bond, "stereo_helical", None),
                "is_aromatic": bond.is_aromatic,
                "display_order": bond.display_order,
                "is_query": bond.is_query,
                "ring_id": bond.ring_id,
                "length_px": bond.length_px,
                "stroke_px": bond.stroke_px,
                "color": bond.color,
                "donor_atom_id": getattr(bond, "donor_atom_id", None),
                "flex_curve_1": getattr(bond, "flex_curve_1", None),
                "flex_curve_2": getattr(bond, "flex_curve_2", None),
                "pi_offset_sign": getattr(bond, "pi_offset_sign", None),
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
                charge=atom_d.get("formal_charge", atom_d.get("charge", 0)),
                isotope=atom_d.get("isotope"),
                radical_electrons=int(atom_d.get("radical_electrons", 0) or 0),
                oxidation_state=atom_d.get("oxidation_state"),
                stereo_cip=atom_d.get("stereo_cip"),
                stereo_axial=atom_d.get("stereo_axial"),
                stereo_helical=atom_d.get("stereo_helical"),
                stereo_si_re=atom_d.get("stereo_si_re"),
                explicit_h=atom_d.get("explicit_h"),
                group_h_cap=atom_d.get("group_h_cap"),
                mapping=atom_d.get("mapping"),
                is_query=atom_d.get("is_query", False),
                is_explicit=atom_d.get("is_explicit", False),
                no_implicit=bool(atom_d.get("no_implicit", False)),
                label_scale=atom_d.get("label_scale"),
                is_coordination_center=atom_d.get("is_coordination_center", False),
                sphere_radius=atom_d.get("sphere_radius"),
                sphere_color=atom_d.get("sphere_color"),
                sphere_filled=bool(atom_d.get("sphere_filled", True)),
                sphere_transparent=bool(atom_d.get("sphere_transparent", False)),
            )
            
        for bond_d in model_data.get("bonds", []):
            a1_id = bond_d["a1_id"]
            a2_id = bond_d["a2_id"]
            style = PersistenceManager._parse_bond_style(bond_d)
            stereo = PersistenceManager._parse_bond_stereo(bond_d)
            donor = bond_d.get("donor_atom_id")
            if style == BondStyle.COORDINATION:
                donor = PersistenceManager._infer_coordination_donor(
                    canvas.model,
                    a1_id,
                    a2_id,
                    preferred=donor,
                )
            else:
                donor = None
            canvas.model.add_bond(
                a1_id=a1_id,
                a2_id=a2_id,
                order=bond_d.get("order", 1),
                bond_id=bond_d["id"],
                style=style,
                stereo=stereo,
                stereo_ez=bond_d.get("stereo_ez"),
                stereo_axial=bond_d.get("stereo_axial"),
                stereo_endo_exo=bond_d.get("stereo_endo_exo"),
                stereo_helical=bond_d.get("stereo_helical"),
                is_aromatic=bond_d.get("is_aromatic", False),
                display_order=bond_d.get("display_order"),
                is_query=bond_d.get("is_query", False),
                ring_id=bond_d.get("ring_id"),
                length_px=bond_d.get("length_px"),
                stroke_px=bond_d.get("stroke_px"),
                color=bond_d.get("color"),
                donor_atom_id=donor,
                flex_curve_1=bond_d.get("flex_curve_1"),
                flex_curve_2=bond_d.get("flex_curve_2"),
                pi_offset_sign=PersistenceManager._parse_pi_offset_sign(bond_d),
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
        temp_path = f"{filepath}.tmp"
        with open(temp_path, 'w', encoding='utf-8') as f:
            json.dump(data, f, indent=2)
        os.replace(temp_path, filepath)

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
