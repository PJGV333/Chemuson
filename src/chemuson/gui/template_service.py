from __future__ import annotations

from copy import deepcopy
from typing import Any, Callable
from uuid import uuid4


def as_dict(data: dict[str, Any]) -> dict[str, Any]:
    """Devuelve una copia serializable del estado actual."""
    return deepcopy(data)


def categories(data: dict[str, Any]) -> list[str]:
    """Lista categorías existentes."""
    return list(data.get("categories", []))


def grouped_templates(
    data: dict[str, Any],
    *,
    default_category_user: str,
) -> list[dict[str, Any]]:
    """Devuelve plantillas agrupadas por categoría."""
    by_category: dict[str, list[dict[str, Any]]] = {name: [] for name in categories(data)}
    for template in data.get("templates", []):
        category = template.get("category", default_category_user)
        by_category.setdefault(category, []).append(deepcopy(template))
    grouped: list[dict[str, Any]] = []
    known_categories = categories(data)
    for category in known_categories:
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


def add_category(
    data: dict[str, Any],
    name: str,
    *,
    clean_name_fn: Callable[[str, str], str],
    default_category_user: str,
) -> str:
    """Crea una categoría nueva o devuelve la existente."""
    category = clean_name_fn(name, default_category_user)
    if category not in data["categories"]:
        data["categories"].append(category)
    return category


def rename_category(
    data: dict[str, Any],
    old_name: str,
    new_name: str,
    *,
    clean_name_fn: Callable[[str, str], str],
) -> str:
    """Renombra una categoría y mueve sus plantillas."""
    old_category = clean_name_fn(old_name, "")
    if not old_category or old_category not in data["categories"]:
        raise ValueError("Categoría no encontrada")
    new_category = clean_name_fn(new_name, old_category)
    if new_category == old_category:
        return old_category
    if new_category not in data["categories"]:
        data["categories"].append(new_category)
    data["categories"] = [category for category in data["categories"] if category != old_category]
    for template in data["templates"]:
        if template.get("category") == old_category:
            template["category"] = new_category
    return new_category


def delete_category(
    data: dict[str, Any],
    category_name: str,
    *,
    clean_name_fn: Callable[[str, str], str],
    default_category_user: str,
    fallback_category: str,
) -> None:
    """Elimina una categoría moviendo sus plantillas a una de respaldo."""
    category = clean_name_fn(category_name, "")
    if not category or category not in data["categories"]:
        raise ValueError("Categoría no encontrada")
    fallback = clean_name_fn(fallback_category, default_category_user)
    if fallback == category:
        fallback = default_category_user if category != default_category_user else "General"
    add_category(
        data,
        fallback,
        clean_name_fn=clean_name_fn,
        default_category_user=default_category_user,
    )
    for template in data["templates"]:
        if template.get("category") == category:
            template["category"] = fallback
    data["categories"] = [name for name in data["categories"] if name != category]


def list_templates(data: dict[str, Any]) -> list[dict[str, Any]]:
    """Devuelve una copia de todas las plantillas."""
    return [deepcopy(template) for template in data.get("templates", [])]


def get_template(data: dict[str, Any], template_id: str) -> dict[str, Any]:
    """Recupera una plantilla por ID."""
    target = str(template_id).strip()
    for template in data.get("templates", []):
        if str(template.get("id", "")).strip() == target:
            return deepcopy(template)
    raise ValueError("Plantilla no encontrada")


def add_template(
    data: dict[str, Any],
    name: str,
    category: str,
    molblock: str,
    *,
    smiles: str,
    clean_name_fn: Callable[[str, str], str],
    default_category_user: str,
    normalize_molblock_header_fn: Callable[[str], str],
) -> dict[str, Any]:
    """Agrega una plantilla cruda (molblock + smiles opcional)."""
    template_name = clean_name_fn(name, "Plantilla")
    template_category = add_category(
        data,
        category,
        clean_name_fn=clean_name_fn,
        default_category_user=default_category_user,
    )
    clean_molblock = normalize_molblock_header_fn(str(molblock or "")).rstrip()
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
    data["templates"].append(entry)
    return deepcopy(entry)


def rename_template(
    data: dict[str, Any],
    template_id: str,
    new_name: str,
    *,
    clean_name_fn: Callable[[str, str], str],
) -> str:
    """Renombra una plantilla existente."""
    clean_template_name = clean_name_fn(new_name, "Plantilla")
    target = str(template_id).strip()
    for template in data["templates"]:
        if str(template.get("id", "")).strip() == target:
            template["name"] = clean_template_name
            return clean_template_name
    raise ValueError("Plantilla no encontrada")


def delete_template(data: dict[str, Any], template_id: str) -> None:
    """Elimina una plantilla por ID."""
    target = str(template_id).strip()
    before = len(data["templates"])
    data["templates"] = [
        template
        for template in data["templates"]
        if str(template.get("id", "")).strip() != target
    ]
    if len(data["templates"]) == before:
        raise ValueError("Plantilla no encontrada")


def merge_imported_library(
    data: dict[str, Any],
    imported: dict[str, Any],
    *,
    default_category_user: str,
    clean_name_fn: Callable[[str, str], str],
) -> int:
    """Fusiona una biblioteca ya normalizada y devuelve cuántas plantillas añadió."""
    existing_signatures = {
        (
            template.get("name", ""),
            template.get("category", ""),
            template.get("molblock", ""),
            template.get("smiles", ""),
        )
        for template in data.get("templates", [])
    }
    added = 0
    for category in imported.get("categories", []):
        add_category(
            data,
            category,
            clean_name_fn=clean_name_fn,
            default_category_user=default_category_user,
        )
    existing_ids = {template.get("id", "") for template in data.get("templates", [])}
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
        new_entry["category"] = add_category(
            data,
            new_entry.get("category", default_category_user),
            clean_name_fn=clean_name_fn,
            default_category_user=default_category_user,
        )
        data["templates"].append(new_entry)
        existing_ids.add(template_id)
        existing_signatures.add(signature)
        added += 1
    return added
