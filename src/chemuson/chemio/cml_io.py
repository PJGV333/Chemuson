from __future__ import annotations

"""Importación/exportación CML básica para interoperabilidad semántica."""

import xml.etree.ElementTree as ET

from chemuson.core.model import BondStyle, BondStereo, MolGraph

_CML_NS = "http://www.xml-cml.org/schema"
_CHEMUSON_NS = "https://chemuson.org/schema"
ET.register_namespace("cml", _CML_NS)
ET.register_namespace("chemuson", _CHEMUSON_NS)


def molgraph_to_cml(graph: MolGraph) -> str:
    """Serializa un ``MolGraph`` a CML 2D con extensiones Chemuson."""
    root = ET.Element(_q("cml"))
    molecule = ET.SubElement(root, _q("molecule"), {"id": "m1"})
    atom_array = ET.SubElement(molecule, _q("atomArray"))
    atom_ids: dict[int, str] = {}
    for atom in sorted(graph.atoms.values(), key=lambda item: item.id):
        atom_ref = f"a{int(atom.id)}"
        atom_ids[int(atom.id)] = atom_ref
        attrs = {
            "id": atom_ref,
            "elementType": str(atom.element),
            "x2": _num(atom.x),
            "y2": _num(atom.y),
        }
        if int(getattr(atom, "charge", 0) or 0):
            attrs["formalCharge"] = str(int(atom.charge))
        if getattr(atom, "isotope", None) is not None:
            attrs["isotopeNumber"] = str(int(atom.isotope))
        if int(getattr(atom, "radical_electrons", 0) or 0):
            attrs[_cq("radicalElectrons")] = str(int(atom.radical_electrons))
        if getattr(atom, "explicit_h", None) is not None:
            attrs[_cq("explicitH")] = str(int(atom.explicit_h))
        if bool(getattr(atom, "is_explicit", False)):
            attrs[_cq("isExplicit")] = "true"
        if bool(getattr(atom, "no_implicit", False)):
            attrs[_cq("noImplicit")] = "true"
        r_subs = getattr(atom, "r_group_substituents", ()) or ()
        if r_subs:
            attrs[_cq("rGroupSubstituents")] = ",".join(str(item) for item in r_subs)
        ET.SubElement(atom_array, _q("atom"), attrs)

    bond_array = ET.SubElement(molecule, _q("bondArray"))
    for bond in sorted(graph.bonds.values(), key=lambda item: item.id):
        a1 = atom_ids.get(int(bond.a1_id))
        a2 = atom_ids.get(int(bond.a2_id))
        if not a1 or not a2:
            continue
        attrs = {
            "id": f"b{int(bond.id)}",
            "atomRefs2": f"{a1} {a2}",
            "order": "A" if bool(getattr(bond, "is_aromatic", False)) else str(int(getattr(bond, "order", 1) or 1)),
        }
        style = getattr(bond, "style", BondStyle.PLAIN)
        stereo = getattr(bond, "stereo", BondStereo.NONE)
        attrs[_cq("style")] = str(getattr(style, "value", style))
        attrs[_cq("stereo")] = str(getattr(stereo, "value", stereo))
        if getattr(bond, "display_order", None) is not None:
            attrs[_cq("displayOrder")] = str(int(bond.display_order))
        if getattr(bond, "donor_atom_id", None) is not None:
            attrs[_cq("donorAtomId")] = atom_ids.get(int(bond.donor_atom_id), "")
        ET.SubElement(bond_array, _q("bond"), attrs)
    return ET.tostring(root, encoding="unicode", xml_declaration=True)


def cml_to_molgraph(cml: str) -> MolGraph:
    """Importa CML básico a ``MolGraph``."""
    root = ET.fromstring(str(cml or ""))
    graph = MolGraph()
    atom_id_map: dict[str, int] = {}
    for atom_el in _elements(root, "atom"):
        cml_id = atom_el.attrib.get("id") or f"a{len(atom_id_map) + 1}"
        element = atom_el.attrib.get("elementType", "C")
        atom = graph.add_atom(
            element,
            _float_attr(atom_el, "x2", 0.0),
            _float_attr(atom_el, "y2", 0.0),
            charge=_int_attr(atom_el, "formalCharge", 0),
            isotope=_optional_int_attr(atom_el, "isotopeNumber"),
            radical_electrons=_int_attr(atom_el, _cq("radicalElectrons"), 0),
            explicit_h=_optional_int_attr(atom_el, _cq("explicitH")),
            is_explicit=_bool_attr(atom_el, _cq("isExplicit"), element != "C"),
            no_implicit=_bool_attr(atom_el, _cq("noImplicit"), False),
            r_group_substituents=_split_csv(atom_el.attrib.get(_cq("rGroupSubstituents"))),
        )
        atom_id_map[cml_id] = atom.id

    for bond_el in _elements(root, "bond"):
        refs = str(bond_el.attrib.get("atomRefs2", "")).split()
        if len(refs) != 2 or refs[0] not in atom_id_map or refs[1] not in atom_id_map:
            continue
        order_text = str(bond_el.attrib.get("order", "1"))
        is_aromatic = order_text.upper() == "A"
        order = 1 if is_aromatic else _safe_int(order_text, 1)
        style = _parse_enum(BondStyle, bond_el.attrib.get(_cq("style")), BondStyle.PLAIN)
        stereo = _parse_enum(BondStereo, bond_el.attrib.get(_cq("stereo")), BondStereo.NONE)
        donor_ref = bond_el.attrib.get(_cq("donorAtomId"))
        donor_atom_id = atom_id_map.get(donor_ref) if donor_ref else None
        graph.add_bond(
            atom_id_map[refs[0]],
            atom_id_map[refs[1]],
            order=order,
            style=style,
            stereo=stereo,
            is_aromatic=is_aromatic,
            display_order=_optional_int_attr(bond_el, _cq("displayOrder")),
            donor_atom_id=donor_atom_id,
        )
    return graph


def _q(name: str) -> str:
    return f"{{{_CML_NS}}}{name}"


def _cq(name: str) -> str:
    return f"{{{_CHEMUSON_NS}}}{name}"


def _elements(root: ET.Element, local_name: str) -> list[ET.Element]:
    return [element for element in root.iter() if element.tag.split("}")[-1] == local_name]


def _num(value: object) -> str:
    return f"{float(value):.6g}"


def _safe_int(value: object, default: int = 0) -> int:
    try:
        return int(value)
    except Exception:
        return int(default)


def _int_attr(element: ET.Element, name: str, default: int = 0) -> int:
    return _safe_int(element.attrib.get(name), default)


def _optional_int_attr(element: ET.Element, name: str) -> int | None:
    if name not in element.attrib:
        return None
    return _safe_int(element.attrib.get(name), 0)


def _float_attr(element: ET.Element, name: str, default: float = 0.0) -> float:
    try:
        return float(element.attrib.get(name, default))
    except Exception:
        return float(default)


def _bool_attr(element: ET.Element, name: str, default: bool = False) -> bool:
    if name not in element.attrib:
        return bool(default)
    return str(element.attrib.get(name, "")).strip().lower() in {"1", "true", "yes"}


def _split_csv(value: str | None) -> tuple[str, ...]:
    if not value:
        return ()
    return tuple(item.strip() for item in value.split(",") if item.strip())


def _parse_enum(enum_cls, value: str | None, default):
    try:
        return enum_cls(value)
    except Exception:
        return default
