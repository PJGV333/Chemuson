"""Persistencia y gestión de biblioteca de plantillas del usuario."""

from __future__ import annotations

import json
import math
import os
from copy import deepcopy
from pathlib import Path
from typing import Any, Dict, Iterable, Optional
from uuid import uuid4

from chemio.rdkit_io import (
    molfile_to_molgraph,
    molgraph_to_molfile,
    molgraph_to_smiles,
    smiles_to_molgraph,
)
from core.model import BondStereo, BondStyle, MolGraph
from gui.templates import build_chair_template, build_haworth_template, build_linear_chain_template


LIBRARY_VERSION = 1
DEFAULT_CATEGORY_USER = "Usuario"


def _default_library_path() -> Path:
    """Devuelve la ruta por defecto de la biblioteca de plantillas."""
    base = os.environ.get("XDG_CONFIG_HOME")
    if not base:
        base = os.path.join(os.path.expanduser("~"), ".config")
    return Path(base) / "Chemuson" / "template_library.json"


def _clean_name(value: str, fallback: str) -> str:
    """Normaliza nombres para categorías y plantillas."""
    cleaned = " ".join(str(value).strip().split())
    return cleaned or fallback


def _safe_smiles(graph: MolGraph) -> str:
    """Exporta SMILES en modo tolerante."""
    try:
        return molgraph_to_smiles(graph)
    except Exception:
        return ""


def _ring_graph(elements: Iterable[str], aromatic: bool = True, bond_length: float = 40.0) -> MolGraph:
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


def _default_data() -> dict[str, Any]:
    """Genera la biblioteca inicial con categorías y plantillas base."""
    categories = ["Aromáticos", "Heterociclos", "Bioquímicos", DEFAULT_CATEGORY_USER]
    templates: list[dict[str, Any]] = []

    builtins = [
        ("Benceno", "Aromáticos", _ring_graph(["C", "C", "C", "C", "C", "C"], aromatic=True)),
        ("Piridina", "Heterociclos", _ring_graph(["N", "C", "C", "C", "C", "C"], aromatic=True)),
        ("Cadena lineal", "Bioquímicos", build_linear_chain_template(40.0)),
        ("Haworth β", "Bioquímicos", build_haworth_template(40.0, anomeric_up=True, bold_front=True)),
        ("Silla β", "Bioquímicos", build_chair_template(40.0, anomeric_up=True, bold_front=True)),
    ]

    for name, category, graph in builtins:
        try:
            molblock = molgraph_to_molfile(graph)
        except Exception:
            continue
        templates.append(
            {
                "id": f"tpl_{uuid4().hex}",
                "name": name,
                "category": category,
                "molblock": molblock,
                "smiles": _safe_smiles(graph),
            }
        )

    return {
        "version": LIBRARY_VERSION,
        "categories": categories,
        "templates": templates,
    }


class TemplateLibrary:
    """Maneja categorías, plantillas y persistencia JSON del usuario."""

    def __init__(self, path: Optional[str | Path] = None) -> None:
        self._path = Path(path) if path is not None else _default_library_path()
        self._data: dict[str, Any] = {}
        self.load()

    @property
    def path(self) -> Path:
        """Ruta física de almacenamiento."""
        return self._path

    def load(self) -> None:
        """Carga biblioteca desde disco o inicializa valores por defecto."""
        if self._path.exists():
            with self._path.open("r", encoding="utf-8") as fh:
                raw = json.load(fh)
            self._data = self._normalize(raw)
            self.save()
            return
        self._data = _default_data()
        self.save()

    def save(self) -> None:
        """Guarda biblioteca a disco."""
        self._path.parent.mkdir(parents=True, exist_ok=True)
        with self._path.open("w", encoding="utf-8") as fh:
            json.dump(self._data, fh, ensure_ascii=False, indent=2)

    def _normalize(self, raw: dict[str, Any]) -> dict[str, Any]:
        """Normaliza un diccionario de biblioteca potencialmente incompleto."""
        normalized = {
            "version": int(raw.get("version", LIBRARY_VERSION)),
            "categories": [],
            "templates": [],
        }
        seen_categories: set[str] = set()
        for category in raw.get("categories", []):
            name = _clean_name(category, "")
            if not name or name in seen_categories:
                continue
            seen_categories.add(name)
            normalized["categories"].append(name)
        for template in raw.get("templates", []):
            if not isinstance(template, dict):
                continue
            name = _clean_name(template.get("name", ""), "")
            category = _clean_name(template.get("category", DEFAULT_CATEGORY_USER), DEFAULT_CATEGORY_USER)
            molblock = str(template.get("molblock", "")).strip()
            smiles = str(template.get("smiles", "")).strip()
            if not name:
                continue
            if not molblock and not smiles:
                continue
            template_id = str(template.get("id", "")).strip() or f"tpl_{uuid4().hex}"
            normalized["templates"].append(
                {
                    "id": template_id,
                    "name": name,
                    "category": category,
                    "molblock": molblock,
                    "smiles": smiles,
                }
            )
            if category not in seen_categories:
                seen_categories.add(category)
                normalized["categories"].append(category)
        if DEFAULT_CATEGORY_USER not in seen_categories:
            normalized["categories"].append(DEFAULT_CATEGORY_USER)
        return normalized

    def as_dict(self) -> dict[str, Any]:
        """Devuelve una copia serializable del estado actual."""
        return deepcopy(self._data)

    def categories(self) -> list[str]:
        """Lista categorías existentes."""
        return list(self._data.get("categories", []))

    def grouped_templates(self) -> list[dict[str, Any]]:
        """Devuelve plantillas agrupadas por categoría."""
        by_category: dict[str, list[dict[str, Any]]] = {name: [] for name in self.categories()}
        for template in self._data.get("templates", []):
            category = template.get("category", DEFAULT_CATEGORY_USER)
            by_category.setdefault(category, []).append(deepcopy(template))
        grouped: list[dict[str, Any]] = []
        for category in self.categories():
            templates = sorted(
                by_category.get(category, []),
                key=lambda tpl: str(tpl.get("name", "")).lower(),
            )
            grouped.append({"name": category, "templates": templates})
        for category, templates in sorted(by_category.items()):
            if category in {group["name"] for group in grouped}:
                continue
            grouped.append(
                {
                    "name": category,
                    "templates": sorted(
                        templates,
                        key=lambda tpl: str(tpl.get("name", "")).lower(),
                    ),
                }
            )
        return grouped

    def add_category(self, name: str) -> str:
        """Crea una categoría nueva o devuelve la existente."""
        category = _clean_name(name, DEFAULT_CATEGORY_USER)
        if category not in self._data["categories"]:
            self._data["categories"].append(category)
            self.save()
        return category

    def rename_category(self, old_name: str, new_name: str) -> str:
        """Renombra una categoría y mueve sus plantillas."""
        old_category = _clean_name(old_name, "")
        if not old_category or old_category not in self._data["categories"]:
            raise ValueError("Categoría no encontrada")
        new_category = _clean_name(new_name, old_category)
        if new_category == old_category:
            return old_category
        if new_category not in self._data["categories"]:
            self._data["categories"].append(new_category)
        self._data["categories"] = [
            category for category in self._data["categories"] if category != old_category
        ]
        for template in self._data["templates"]:
            if template.get("category") == old_category:
                template["category"] = new_category
        self.save()
        return new_category

    def delete_category(self, category_name: str, fallback_category: str = DEFAULT_CATEGORY_USER) -> None:
        """Elimina una categoría moviendo sus plantillas a una de respaldo."""
        category = _clean_name(category_name, "")
        if not category or category not in self._data["categories"]:
            raise ValueError("Categoría no encontrada")
        fallback = _clean_name(fallback_category, DEFAULT_CATEGORY_USER)
        if fallback == category:
            fallback = DEFAULT_CATEGORY_USER if category != DEFAULT_CATEGORY_USER else "General"
        self.add_category(fallback)
        for template in self._data["templates"]:
            if template.get("category") == category:
                template["category"] = fallback
        self._data["categories"] = [name for name in self._data["categories"] if name != category]
        self.save()

    def list_templates(self) -> list[dict[str, Any]]:
        """Devuelve una copia de todas las plantillas."""
        return [deepcopy(template) for template in self._data.get("templates", [])]

    def get_template(self, template_id: str) -> dict[str, Any]:
        """Recupera una plantilla por ID."""
        target = str(template_id).strip()
        for template in self._data.get("templates", []):
            if str(template.get("id", "")).strip() == target:
                return deepcopy(template)
        raise ValueError("Plantilla no encontrada")

    def add_template(
        self,
        name: str,
        category: str,
        molblock: str,
        smiles: str = "",
    ) -> dict[str, Any]:
        """Agrega una plantilla cruda (molblock + smiles opcional)."""
        template_name = _clean_name(name, "Plantilla")
        template_category = self.add_category(category)
        clean_molblock = str(molblock or "").strip()
        clean_smiles = str(smiles or "").strip()
        if not clean_molblock and not clean_smiles:
            raise ValueError("La plantilla requiere molblock o smiles")
        entry = {
            "id": f"tpl_{uuid4().hex}",
            "name": template_name,
            "category": template_category,
            "molblock": clean_molblock,
            "smiles": clean_smiles,
        }
        self._data["templates"].append(entry)
        self.save()
        return deepcopy(entry)

    def add_template_from_graph(self, name: str, category: str, graph: MolGraph) -> dict[str, Any]:
        """Agrega una plantilla a partir de un `MolGraph`."""
        if not graph.atoms:
            raise ValueError("La plantilla está vacía")
        molblock = molgraph_to_molfile(graph)
        smiles = _safe_smiles(graph)
        return self.add_template(name, category, molblock, smiles)

    def rename_template(self, template_id: str, new_name: str) -> str:
        """Renombra una plantilla existente."""
        clean_name = _clean_name(new_name, "Plantilla")
        target = str(template_id).strip()
        for template in self._data["templates"]:
            if str(template.get("id", "")).strip() == target:
                template["name"] = clean_name
                self.save()
                return clean_name
        raise ValueError("Plantilla no encontrada")

    def delete_template(self, template_id: str) -> None:
        """Elimina una plantilla por ID."""
        target = str(template_id).strip()
        before = len(self._data["templates"])
        self._data["templates"] = [
            template
            for template in self._data["templates"]
            if str(template.get("id", "")).strip() != target
        ]
        if len(self._data["templates"]) == before:
            raise ValueError("Plantilla no encontrada")
        self.save()

    def graph_from_template(self, template_id: str) -> MolGraph:
        """Convierte una plantilla almacenada a `MolGraph`."""
        template = self.get_template(template_id)
        molblock = str(template.get("molblock", "")).strip()
        if molblock:
            return molfile_to_molgraph(molblock)
        smiles = str(template.get("smiles", "")).strip()
        if smiles:
            return smiles_to_molgraph(smiles)
        raise ValueError("Plantilla sin contenido químico")

    def export_to_file(self, output_path: str | Path) -> None:
        """Exporta la biblioteca actual a un archivo JSON."""
        out_path = Path(output_path)
        out_path.parent.mkdir(parents=True, exist_ok=True)
        with out_path.open("w", encoding="utf-8") as fh:
            json.dump(self._data, fh, ensure_ascii=False, indent=2)

    def import_from_file(self, input_path: str | Path, merge: bool = True) -> int:
        """Importa biblioteca JSON.

        Returns:
            Número de plantillas incorporadas.
        """
        in_path = Path(input_path)
        with in_path.open("r", encoding="utf-8") as fh:
            raw = json.load(fh)
        imported = self._normalize(raw)
        if not merge:
            self._data = imported
            self.save()
            return len(imported.get("templates", []))

        existing_signatures = {
            (
                template.get("name", ""),
                template.get("category", ""),
                template.get("molblock", ""),
                template.get("smiles", ""),
            )
            for template in self._data.get("templates", [])
        }
        added = 0
        for category in imported.get("categories", []):
            self.add_category(category)
        existing_ids = {template.get("id", "") for template in self._data.get("templates", [])}
        for template in imported.get("templates", []):
            signature = (
                template.get("name", ""),
                template.get("category", ""),
                template.get("molblock", ""),
                template.get("smiles", ""),
            )
            if signature in existing_signatures:
                continue
            new_entry = deepcopy(template)
            template_id = str(new_entry.get("id", "")).strip() or f"tpl_{uuid4().hex}"
            while template_id in existing_ids:
                template_id = f"tpl_{uuid4().hex}"
            new_entry["id"] = template_id
            new_entry["category"] = self.add_category(
                new_entry.get("category", DEFAULT_CATEGORY_USER)
            )
            self._data["templates"].append(new_entry)
            existing_ids.add(template_id)
            existing_signatures.add(signature)
            added += 1
        self.save()
        return added
