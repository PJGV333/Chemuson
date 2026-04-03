"""Compositor liviano de acciones de proyecto."""

from .numbering_actions import create_numbering_actions
from .structure_actions import create_structure_actions
from .view_actions import create_view_actions


def create_project_actions(window) -> None:
    """Inicializa acciones de vista/estructura/proyecto en ``window``."""
    create_view_actions(window)
    create_numbering_actions(window)
    create_structure_actions(window)
