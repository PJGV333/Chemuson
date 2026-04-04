from __future__ import annotations

import math
from typing import Any, Callable, Iterable
from uuid import uuid4

from chemuson.core.model import BondStereo, BondStyle, MolGraph
from chemuson.gui.templates import (
    build_chair_template,
    build_haworth_template,
    build_linear_chain_template,
)


def ring_graph(
    elements: Iterable[str],
    aromatic: bool = True,
    bond_length: float = 40.0,
) -> MolGraph:
    """Construye un anillo regular simple para plantillas predefinidas."""
    symbols = list(elements)
    n = len(symbols)
    graph = MolGraph()
    if n < 3:
        return graph
    radius = bond_length / (2.0 * math.sin(math.pi / n))
    start = math.pi / 6.0
    atom_ids: list[int] = []
    for i, symbol in enumerate(symbols):
        theta = start + (2.0 * math.pi * i / n)
        x = radius * math.cos(theta)
        y = radius * math.sin(theta)
        atom = graph.add_atom(symbol, x, y, is_explicit=(symbol != "C"))
        atom_ids.append(atom.id)
    for i in range(n):
        a1 = atom_ids[i]
        a2 = atom_ids[(i + 1) % n]
        graph.add_bond(
            a1,
            a2,
            order=1,
            style=BondStyle.PLAIN,
            stereo=BondStereo.NONE,
            is_aromatic=aromatic,
        )
    return graph


def builtin_template_specs() -> list[tuple[str, str, MolGraph]]:
    """Devuelve las plantillas base que siempre deben existir."""
    return [
        ("Benceno", "Aromáticos", ring_graph(["C", "C", "C", "C", "C", "C"], aromatic=True)),
        ("Piridina", "Heterociclos", ring_graph(["N", "C", "C", "C", "C", "C"], aromatic=True)),
        ("Cadena lineal", "Bioquímicos", build_linear_chain_template(40.0)),
        ("Haworth β", "Bioquímicos", build_haworth_template(40.0, anomeric_up=True, bold_front=True)),
        ("Silla β", "Bioquímicos", build_chair_template(40.0, anomeric_up=True, bold_front=True)),
    ]


def default_data(
    *,
    library_version: int,
    default_category_user: str,
    molgraph_to_molfile_fn: Callable[[MolGraph], str],
    safe_smiles_fn: Callable[[MolGraph], str],
) -> dict[str, Any]:
    """Genera la biblioteca inicial con categorías y plantillas base."""
    categories = ["Aromáticos", "Heterociclos", "Bioquímicos", default_category_user]
    templates: list[dict[str, Any]] = []
    for name, category, graph in builtin_template_specs():
        try:
            molblock = molgraph_to_molfile_fn(graph)
        except Exception:
            continue
        templates.append(
            {
                "id": f"tpl_{uuid4().hex}",
                "name": name,
                "category": category,
                "molblock": molblock,
                "smiles": safe_smiles_fn(graph),
            }
        )
    return {
        "version": library_version,
        "categories": categories,
        "templates": templates,
    }


def sync_builtin_templates(
    data: dict[str, Any],
    *,
    clean_name_fn: Callable[[str, str], str],
    molgraph_to_molfile_fn: Callable[[MolGraph], str],
    safe_smiles_fn: Callable[[MolGraph], str],
) -> None:
    """Sincroniza plantillas predefinidas por nombre/categoría."""
    builtin_by_key: dict[tuple[str, str], dict[str, str]] = {}
    for name, category, graph in builtin_template_specs():
        key = (clean_name_fn(name, name), clean_name_fn(category, category))
        try:
            molblock = molgraph_to_molfile_fn(graph)
        except Exception:
            continue
        builtin_by_key[key] = {
            "name": key[0],
            "category": key[1],
            "molblock": molblock,
            "smiles": safe_smiles_fn(graph),
        }
    if not builtin_by_key:
        return

    existing = list(data.get("templates", []))
    keep_templates: list[dict[str, Any]] = []
    matched_ids: dict[tuple[str, str], str] = {}
    for template in existing:
        key = (
            clean_name_fn(template.get("name", ""), ""),
            clean_name_fn(template.get("category", ""), ""),
        )
        if key in builtin_by_key:
            matched_ids.setdefault(key, str(template.get("id", "")).strip() or f"tpl_{uuid4().hex}")
            continue
        keep_templates.append(template)

    used_ids = {
        str(template.get("id", "")).strip()
        for template in keep_templates
        if str(template.get("id", "")).strip()
    }
    builtin_templates: list[dict[str, Any]] = []
    for key, payload in builtin_by_key.items():
        template_id = matched_ids.get(key, f"tpl_{uuid4().hex}")
        while template_id in used_ids:
            template_id = f"tpl_{uuid4().hex}"
        used_ids.add(template_id)
        builtin_templates.append(
            {
                "id": template_id,
                "name": payload["name"],
                "category": payload["category"],
                "molblock": payload["molblock"],
                "smiles": payload["smiles"],
            }
        )
        if payload["category"] not in data["categories"]:
            data["categories"].append(payload["category"])

    data["templates"] = keep_templates + builtin_templates
