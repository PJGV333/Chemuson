from __future__ import annotations

import json
import os
from pathlib import Path
from typing import Any, Callable
from uuid import uuid4


def is_windows_platform() -> bool:
    """Indica si la plataforma de ejecución es Windows."""
    return os.name == "nt"


def default_config_dir(*, is_windows_fn: Callable[[], bool] = is_windows_platform) -> Path:
    """Devuelve directorio de configuración por plataforma."""
    override = os.environ.get("CHEMUSON_CONFIG_HOME")
    if override:
        return Path(override)

    if is_windows_fn():
        base = os.environ.get("APPDATA") or os.environ.get("LOCALAPPDATA")
        if base:
            return Path(base)
        home = Path(os.path.expanduser("~"))
        return home / "AppData" / "Roaming"

    base = os.environ.get("XDG_CONFIG_HOME")
    if base:
        return Path(base)
    return Path(os.path.expanduser("~")) / ".config"


def default_library_path(
    *,
    default_config_dir_fn: Callable[[], Path] = default_config_dir,
) -> Path:
    """Devuelve la ruta por defecto de la biblioteca de plantillas."""
    return default_config_dir_fn() / "Chemuson" / "template_library.json"


def clean_name(value: str, fallback: str) -> str:
    """Normaliza nombres para categorías y plantillas."""
    cleaned = " ".join(str(value).strip().split())
    return cleaned or fallback


def save_library_data(path: str | Path, data: dict[str, Any]) -> None:
    """Guarda biblioteca a disco."""
    target = Path(path)
    target.parent.mkdir(parents=True, exist_ok=True)
    with target.open("w", encoding="utf-8") as fh:
        json.dump(data, fh, ensure_ascii=False, indent=2)


def normalize_library_data(
    raw: dict[str, Any],
    *,
    library_version: int,
    default_category_user: str,
    clean_name_fn: Callable[[str, str], str] = clean_name,
    normalize_molblock_header_fn: Callable[[str], str],
) -> dict[str, Any]:
    """Normaliza un diccionario de biblioteca potencialmente incompleto."""
    normalized = {
        "version": int(raw.get("version", library_version)),
        "categories": [],
        "templates": [],
    }
    seen_categories: set[str] = set()
    for category in raw.get("categories", []):
        name = clean_name_fn(category, "")
        if not name or name in seen_categories:
            continue
        seen_categories.add(name)
        normalized["categories"].append(name)
    for template in raw.get("templates", []):
        if not isinstance(template, dict):
            continue
        name = clean_name_fn(template.get("name", ""), "")
        category = clean_name_fn(template.get("category", default_category_user), default_category_user)
        molblock = normalize_molblock_header_fn(str(template.get("molblock", ""))).rstrip()
        smiles = str(template.get("smiles", "")).strip()
        if not name or (not molblock and not smiles):
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
    if default_category_user not in seen_categories:
        normalized["categories"].append(default_category_user)
    return normalized
