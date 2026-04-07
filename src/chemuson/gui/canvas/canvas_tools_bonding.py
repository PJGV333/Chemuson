from __future__ import annotations

from .canvas_bond_state import CanvasBondStateMixin
from .canvas_bond_model_ops import CanvasBondModelOpsMixin
from .canvas_bond_drag import CanvasBondDragMixin
from .canvas_bond_geometry import CanvasBondGeometryMixin
from .canvas_bond_hit_testing import CanvasBondHitTestingMixin


class CanvasToolsBondingMixin(
    CanvasBondStateMixin,
    CanvasBondModelOpsMixin,
    CanvasBondDragMixin,
    CanvasBondHitTestingMixin,
    CanvasBondGeometryMixin,
):
    """Fachada compatible para las herramientas de bonding del canvas."""
