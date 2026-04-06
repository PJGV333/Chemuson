"""Persistencia y CRUD de la biblioteca de plantillas."""

from __future__ import annotations

import json
from pathlib import Path
from typing import Any, Callable, Optional

from chemuson.gui import template_service
from chemuson.gui.template_store import save_library_data


class TemplateRepository:
    """Mantiene el estado JSON persistido de la biblioteca de plantillas."""

    def __init__(
        self,
        path: str | Path,
        *,
        default_data_fn: Callable[[], dict[str, Any]],
        normalize_fn: Callable[[dict[str, Any]], dict[str, Any]],
        sync_builtins_fn: Callable[[dict[str, Any]], None],
        default_category_user: str,
        clean_name_fn: Callable[[str, str], str],
        normalize_molblock_header_fn: Callable[[str], str],
    ) -> None:
        self._path = Path(path)
        self._default_data_fn = default_data_fn
        self._normalize_fn = normalize_fn
        self._sync_builtins_fn = sync_builtins_fn
        self._default_category_user = default_category_user
        self._clean_name_fn = clean_name_fn
        self._normalize_molblock_header_fn = normalize_molblock_header_fn
        self._data: dict[str, Any] = {}
        self.load()

    @property
    def path(self) -> Path:
        return self._path

    @property
    def data(self) -> dict[str, Any]:
        return self._data

    def load(self) -> None:
        if self._path.exists():
            with self._path.open("r", encoding="utf-8") as fh:
                raw = json.load(fh)
            self._data = self._normalize_fn(raw)
            self._sync_builtins_fn(self._data)
            self.save()
            return
        self._data = self._default_data_fn()
        self.save()

    def save(self) -> None:
        save_library_data(self._path, self._data)

    def as_dict(self) -> dict[str, Any]:
        return template_service.as_dict(self._data)

    def categories(self) -> list[str]:
        return template_service.categories(self._data)

    def grouped_templates(self) -> list[dict[str, Any]]:
        return template_service.grouped_templates(
            self._data,
            default_category_user=self._default_category_user,
        )

    def add_category(self, name: str, *, clean_name_fn: Callable[[str, str], str]) -> str:
        category = template_service.add_category(
            self._data,
            name,
            clean_name_fn=clean_name_fn,
            default_category_user=self._default_category_user,
        )
        self.save()
        return category

    def rename_category(
        self,
        old_name: str,
        new_name: str,
        *,
        clean_name_fn: Callable[[str, str], str],
    ) -> str:
        category = template_service.rename_category(
            self._data,
            old_name,
            new_name,
            clean_name_fn=clean_name_fn,
        )
        self.save()
        return category

    def delete_category(
        self,
        category_name: str,
        *,
        clean_name_fn: Callable[[str, str], str],
        fallback_category: str,
    ) -> None:
        template_service.delete_category(
            self._data,
            category_name,
            clean_name_fn=clean_name_fn,
            default_category_user=self._default_category_user,
            fallback_category=fallback_category,
        )
        self.save()

    def list_templates(self) -> list[dict[str, Any]]:
        return template_service.list_templates(self._data)

    def get_template(self, template_id: str) -> dict[str, Any]:
        return template_service.get_template(self._data, template_id)

    def add_template(
        self,
        name: str,
        category: str,
        molblock: str,
        *,
        smiles: str,
        clean_name_fn: Optional[Callable[[str, str], str]] = None,
        normalize_molblock_header_fn: Optional[Callable[[str], str]] = None,
    ) -> dict[str, Any]:
        template = template_service.add_template(
            self._data,
            name,
            category,
            molblock,
            smiles=smiles,
            clean_name_fn=clean_name_fn or self._clean_name_fn,
            default_category_user=self._default_category_user,
            normalize_molblock_header_fn=(
                normalize_molblock_header_fn or self._normalize_molblock_header_fn
            ),
        )
        self.save()
        return template

    def rename_template(
        self,
        template_id: str,
        new_name: str,
        *,
        clean_name_fn: Callable[[str, str], str],
    ) -> str:
        name = template_service.rename_template(
            self._data,
            template_id,
            new_name,
            clean_name_fn=clean_name_fn,
        )
        self.save()
        return name

    def delete_template(self, template_id: str) -> None:
        template_service.delete_template(self._data, template_id)
        self.save()

    def export_to_file(self, output_path: str | Path) -> None:
        save_library_data(output_path, self._data)

    def import_from_file(
        self,
        input_path: str | Path,
        *,
        merge: bool,
        clean_name_fn: Callable[[str, str], str],
    ) -> int:
        in_path = Path(input_path)
        with in_path.open("r", encoding="utf-8") as fh:
            raw = json.load(fh)
        imported = self._normalize_fn(raw)
        if not merge:
            self._data = imported
            self.save()
            return len(imported.get("templates", []))
        added = template_service.merge_imported_library(
            self._data,
            imported,
            default_category_user=self._default_category_user,
            clean_name_fn=clean_name_fn,
        )
        self.save()
        return added
