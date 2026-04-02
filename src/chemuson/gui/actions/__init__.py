"""Factory helpers para crear grupos de acciones de GUI."""

from chemuson.gui.actions.edit_actions import create_edit_actions
from chemuson.gui.actions.file_actions import create_file_actions
from chemuson.gui.actions.project_actions import create_project_actions
from chemuson.gui.actions.update_actions import create_update_actions

__all__ = [
    "create_file_actions",
    "create_edit_actions",
    "create_project_actions",
    "create_update_actions",
]
