"""Persistencia y gestión de biblioteca de plantillas del usuario."""

from __future__ import annotations

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
from chemuson.gui.template_chem_adapter import TemplateChemAdapter
from chemuson.gui.template_builtins import default_data as build_default_data
from chemuson.gui.template_builtins import ring_graph, sync_builtin_templates
from chemuson.gui.template_conversion import safe_smiles
from chemuson.gui.template_repository import TemplateRepository
from chemuson.gui.template_store import (
    clean_name,
    default_config_dir,
    default_library_path,
    is_windows_platform,
    normalize_library_data,
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
        resolved_path = Path(path) if path is not None else _default_library_path()
        self._repository = TemplateRepository(
            resolved_path,
            default_data_fn=_default_data,
            normalize_fn=self._normalize,
            sync_builtins_fn=self._sync_builtin_templates,
            default_category_user=DEFAULT_CATEGORY_USER,
            clean_name_fn=_clean_name,
            normalize_molblock_header_fn=normalize_molblock_header,
        )
        self._chem_adapter = TemplateChemAdapter(
            molgraph_to_molfile_fn=molgraph_to_molfile,
            safe_smiles_fn=_safe_smiles,
            molfile_to_molgraph_fn=lambda molblock: molfile_to_molgraph(molblock),
            smiles_to_molgraph_fn=lambda smiles: smiles_to_molgraph(smiles),
        )

    @property
    def path(self) -> Path:
        """Ruta física de almacenamiento."""
        return self._repository.path

    def load(self) -> None:
        """Carga biblioteca desde disco o inicializa valores por defecto."""
        self._repository.load()

    def save(self) -> None:
        """Guarda biblioteca a disco."""
        self._repository.save()

    def _normalize(self, raw: dict[str, Any]) -> dict[str, Any]:
        """Normaliza un diccionario de biblioteca potencialmente incompleto."""
        return normalize_library_data(
            raw,
            library_version=LIBRARY_VERSION,
            default_category_user=DEFAULT_CATEGORY_USER,
            clean_name_fn=_clean_name,
            normalize_molblock_header_fn=normalize_molblock_header,
        )

    def _sync_builtin_templates(self, data: Optional[dict[str, Any]] = None) -> None:
        """Reemplaza built-ins desfasadas manteniendo IDs existentes."""
        target_data = data if data is not None else self._repository.data
        sync_builtin_templates(
            target_data,
            clean_name_fn=_clean_name,
            molgraph_to_molfile_fn=molgraph_to_molfile,
            safe_smiles_fn=_safe_smiles,
        )

    def as_dict(self) -> dict[str, Any]:
        """Devuelve una copia serializable del estado actual."""
        return self._repository.as_dict()

    def categories(self) -> list[str]:
        """Lista categorías existentes."""
        return self._repository.categories()

    def grouped_templates(self) -> list[dict[str, Any]]:
        """Devuelve plantillas agrupadas por categoría."""
        return self._repository.grouped_templates()

    def add_category(self, name: str) -> str:
        """Crea una categoría nueva o devuelve la existente."""
        return self._repository.add_category(name, clean_name_fn=_clean_name)

    def rename_category(self, old_name: str, new_name: str) -> str:
        """Renombra una categoría y mueve sus plantillas."""
        return self._repository.rename_category(
            old_name,
            new_name,
            clean_name_fn=_clean_name,
        )

    def delete_category(
        self,
        category_name: str,
        fallback_category: str = DEFAULT_CATEGORY_USER,
    ) -> None:
        """Elimina una categoría moviendo sus plantillas a una de respaldo."""
        self._repository.delete_category(
            category_name,
            clean_name_fn=_clean_name,
            fallback_category=fallback_category,
        )

    def list_templates(self) -> list[dict[str, Any]]:
        """Devuelve una copia de todas las plantillas."""
        return self._repository.list_templates()

    def get_template(self, template_id: str) -> dict[str, Any]:
        """Recupera una plantilla por ID."""
        return self._repository.get_template(template_id)

    def add_template(
        self,
        name: str,
        category: str,
        molblock: str,
        smiles: str = "",
    ) -> dict[str, Any]:
        """Agrega una plantilla cruda (molblock + smiles opcional)."""
        return self._repository.add_template(
            name,
            category,
            molblock,
            smiles=smiles,
            clean_name_fn=_clean_name,
            normalize_molblock_header_fn=normalize_molblock_header,
        )

    def add_template_from_graph(self, name: str, category: str, graph: MolGraph) -> dict[str, Any]:
        """Agrega una plantilla a partir de un `MolGraph`."""
        return self._chem_adapter.add_template_from_graph(self._repository, name, category, graph)

    def rename_template(self, template_id: str, new_name: str) -> str:
        """Renombra una plantilla existente."""
        return self._repository.rename_template(
            template_id,
            new_name,
            clean_name_fn=_clean_name,
        )

    def delete_template(self, template_id: str) -> None:
        """Elimina una plantilla por ID."""
        self._repository.delete_template(template_id)

    def graph_from_template(self, template_id: str) -> MolGraph:
        """Convierte una plantilla almacenada a `MolGraph`."""
        return self._chem_adapter.graph_from_template(self.get_template(template_id))

    def export_to_file(self, output_path: str | Path) -> None:
        """Exporta la biblioteca actual a un archivo JSON."""
        self._repository.export_to_file(output_path)

    def import_from_file(self, input_path: str | Path, merge: bool = True) -> int:
        """Importa biblioteca JSON.

        Returns:
            Número de plantillas incorporadas.
        """
        return self._repository.import_from_file(
            input_path,
            merge=merge,
            clean_name_fn=_clean_name,
        )
