from __future__ import annotations

from .canvas_context_menu import CanvasContextMenuMixin
from .canvas_keyboard import CanvasKeyboardMixin
from .canvas_selection_input import CanvasSelectionInputMixin
from .canvas_tools_annotations import CanvasToolsAnnotationsMixin
from .canvas_tools_bonding import CanvasToolsBondingMixin
from .canvas_tools_rings_chains import CanvasToolsRingsChainsMixin


class CanvasInputMixin(
    CanvasSelectionInputMixin,
    CanvasKeyboardMixin,
    CanvasContextMenuMixin,
    CanvasToolsAnnotationsMixin,
    CanvasToolsBondingMixin,
    CanvasToolsRingsChainsMixin,
):
    """Compatibilidad histórica para la API de input del canvas."""
