"""Stubs para detección de módulos especiales (fase 4, MVP pendiente)."""

from __future__ import annotations

from typing import Any

from .molview import MolView


def detect_special_templates(view: MolView) -> tuple[str, dict[int, int], dict[str, Any]] | None:
    """Hook de plantillas especiales (carbohidratos/esteroides).

    TODO(PR28/Phase4):
    - Integrar plantillas en `chemname/templates/special/`.
    - Soportar match exacto sin RDKit.
    - Soportar subestructura con RDKit cuando esté disponible.
    """
    _ = view
    return None
