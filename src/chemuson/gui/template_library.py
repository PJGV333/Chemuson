"""Persistencia y gestión de biblioteca de plantillas del usuario."""

from __future__ import annotations

import json
from pathlib import Path
from typing import Any, Optional

from chemuson.chemio.rdkit_io import (
    molfile_to_molgraph,
    molgraph_to_molfile,
    molgraph_to_smiles,
    normalize_molblock_header,
    smiles_to_molgraph,
)
from chemuson.core.model import MolGraph
from chemuson.gui import template_service
from chemuson.gui.template_builtins import default_data as build_default_data
from chemuson.gui.template_builtins import ring_graph, sync_builtin_templates
from chemuson.gui.template_conversion import graph_from_template_payload, safe_smiles
from chemuson.gui.template_store import (
    clean_name,
    default_config_dir,
    default_library_path,
    is_windows_platform,
    normalize_library_data,
    save_library_data,
)


LIBRARY_VERSION = 1
DEFAULT_CATEGORY_USER = "Usuario"


def _is_windows_platform() -> bool:
    """Indica si la plataforma de ejecución es Windows."""
    return is_windows_platform()


def _default_config_dir() -> Path:
    """Devuelve directorio de configuración por plataforma."""
    return default_config_dir(is_windows_fn=_is_windows_platform)


def _default_library_path() -> Path:
    """Devuelve la ruta por defecto de la biblioteca de plantillas."""
    return default_library_path(default_config_dir_fn=_default_config_dir)


def _clean_name(value: str, fallback: str) -> str:
    """Normaliza nombres para categorías y plantillas."""
    return clean_name(value, fallback)


def _safe_smiles(graph: MolGraph) -> str:
    """Exporta SMILES en modo tolerante."""
    return safe_smiles(graph, molgraph_to_smiles_fn=molgraph_to_smiles)


def _ring_graph(
    elements,
    aromatic: bool = True,
    bond_length: float = 40.0,
) -> MolGraph:
    """Construye un anillo regular simple para plantillas predefinidas."""
    return ring_graph(elements, aromatic=aromatic, bond_length=bond_length)


def _default_data() -> dict[str, Any]:
    """Genera la biblioteca inicial con categorías y plantillas base."""
    return build_default_data(
        library_version=LIBRARY_VERSION,
        default_category_user=DEFAULT_CATEGORY_USER,
        molgraph_to_molfile_fn=molgraph_to_molfile,
        safe_smiles_fn=_safe_smiles,
    )


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
            self._sync_builtin_templates()
            self.save()
            return
        self._data = _default_data()
        self.save()

    def save(self) -> None:
        """Guarda biblioteca a disco."""
        save_library_data(self._path, self._data)

    def _normalize(self, raw: dict[str, Any]) -> dict[str, Any]:
        """Normaliza un diccionario de biblioteca potencialmente incompleto."""
        return normalize_library_data(
            raw,
            library_version=LIBRARY_VERSION,
            default_category_user=DEFAULT_CATEGORY_USER,
            clean_name_fn=_clean_name,
            normalize_molblock_header_fn=normalize_molblock_header,
        )

    def _sync_builtin_templates(self) -> None:
        """Reemplaza built-ins desfasadas manteniendo IDs existentes."""
        sync_builtin_templates(
            self._data,
            clean_name_fn=_clean_name,
            molgraph_to_molfile_fn=molgraph_to_molfile,
            safe_smiles_fn=_safe_smiles,
        )

    def as_dict(self) -> dict[str, Any]:
        """Devuelve una copia serializable del estado actual."""
        return template_service.as_dict(self._data)

    def categories(self) -> list[str]:
        """Lista categorías existentes."""
        return template_service.categories(self._data)

    def grouped_templates(self) -> list[dict[str, Any]]:
        """Devuelve plantillas agrupadas por categoría."""
        return template_service.grouped_templates(
            self._data,
            default_category_user=DEFAULT_CATEGORY_USER,
        )

    def add_category(self, name: str) -> str:
        """Crea una categoría nueva o devuelve la existente."""
        category = template_service.add_category(
            self._data,
            name,
            clean_name_fn=_clean_name,
            default_category_user=DEFAULT_CATEGORY_USER,
        )
        self.save()
        return category

    def rename_category(self, old_name: str, new_name: str) -> str:
        """Renombra una categoría y mueve sus plantillas."""
        category = template_service.rename_category(
            self._data,
            old_name,
            new_name,
            clean_name_fn=_clean_name,
        )
        self.save()
        return category

    def delete_category(
        self,
        category_name: str,
        fallback_category: str = DEFAULT_CATEGORY_USER,
    ) -> None:
        """Elimina una categoría moviendo sus plantillas a una de respaldo."""
        template_service.delete_category(
            self._data,
            category_name,
            clean_name_fn=_clean_name,
            default_category_user=DEFAULT_CATEGORY_USER,
            fallback_category=fallback_category,
        )
        self.save()

    def list_templates(self) -> list[dict[str, Any]]:
        """Devuelve una copia de todas las plantillas."""
        return template_service.list_templates(self._data)

    def get_template(self, template_id: str) -> dict[str, Any]:
        """Recupera una plantilla por ID."""
        return template_service.get_template(self._data, template_id)

    def add_template(
        self,
        name: str,
        category: str,
        molblock: str,
        smiles: str = "",
    ) -> dict[str, Any]:
        """Agrega una plantilla cruda (molblock + smiles opcional)."""
        template = template_service.add_template(
            self._data,
            name,
            category,
            molblock,
            smiles=smiles,
            clean_name_fn=_clean_name,
            default_category_user=DEFAULT_CATEGORY_USER,
            normalize_molblock_header_fn=normalize_molblock_header,
        )
        self.save()
        return template

    def add_template_from_graph(self, name: str, category: str, graph: MolGraph) -> dict[str, Any]:
        """Agrega una plantilla a partir de un `MolGraph`."""
        if not graph.atoms:
            raise ValueError("La plantilla está vacía")
        molblock = molgraph_to_molfile(graph)
        return self.add_template(name, category, molblock, _safe_smiles(graph))

    def rename_template(self, template_id: str, new_name: str) -> str:
        """Renombra una plantilla existente."""
        name = template_service.rename_template(
            self._data,
            template_id,
            new_name,
            clean_name_fn=_clean_name,
        )
        self.save()
        return name

    def delete_template(self, template_id: str) -> None:
        """Elimina una plantilla por ID."""
        template_service.delete_template(self._data, template_id)
        self.save()

    def graph_from_template(self, template_id: str) -> MolGraph:
        """Convierte una plantilla almacenada a `MolGraph`."""
        template = self.get_template(template_id)
        return graph_from_template_payload(
            template,
            molfile_to_molgraph_fn=molfile_to_molgraph,
            smiles_to_molgraph_fn=smiles_to_molgraph,
        )

    def export_to_file(self, output_path: str | Path) -> None:
        """Exporta la biblioteca actual a un archivo JSON."""
        save_library_data(output_path, self._data)

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
        added = template_service.merge_imported_library(
            self._data,
            imported,
            default_category_user=DEFAULT_CATEGORY_USER,
            clean_name_fn=_clean_name,
        )
        self.save()
        return added
