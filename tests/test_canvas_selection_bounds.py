"""Pruebas funcionales del módulo selection_bounds."""

from __future__ import annotations

from PyQt6.QtCore import QRectF, Qt

import pytest
from unittest.mock import MagicMock

from chemuson.gui.canvas.selection_bounds import (
    resolve_selected_atom_ids,
    selection_bounds,
)
from chemuson.gui.canvas.canvas_selection import CanvasSelectionMixin
from chemuson.gui.composite_diagram_item import CompositeDiagramItem
from chemuson.gui.items import (
    ArrowItem,
    BracketItem,
    EnergyDiagramItem,
    ImageAnnotationItem,
    OrbitalAnnotationItem,
    TextAnnotationItem,
)



# ---------------------------------------------------------------------------
# Fixtures helpers
# ---------------------------------------------------------------------------


class _ModelStub:
    """Stub para model.bonds con soporte de 'in' y get_bond()."""

    def __init__(self) -> None:
        self.bonds: dict[int, _BondStub] = {}

    def get_bond(self, bond_id: int) -> _BondStub:
        return self.bonds[bond_id]


class _BondStub:
    def __init__(self, a1_id: int, a2_id: int) -> None:
        self.a1_id = a1_id
        self.a2_id = a2_id


class _SceneStub:
    """Stub de escena (un solo objeto compartido por 'is')."""

    pass


class _ItemStub:
    """Stub para atom/bond/graphic items."""

    def __init__(
        self,
        scene: object,
        pen_style: Qt.PenStyle = Qt.PenStyle.SolidLine,
        brush_style: Qt.BrushStyle = Qt.BrushStyle.SolidPattern,
        rect: QRectF | None = None,
        label_visible: bool = True,
        charge_label_visible: bool = True,
        label_rect: QRectF | None = None,
        charge_label_rect: QRectF | None = None,
    ) -> None:
        self._scene_obj = scene
        self._pen_style = pen_style
        self._brush_style = brush_style
        self._rect = rect or QRectF(0, 0, 10, 10)
        self._label_visible = label_visible
        self._charge_label_visible = charge_label_visible
        self._label_rect = label_rect or QRectF(0, 0, 5, 5)
        self._charge_label_rect = charge_label_rect or QRectF(0, 0, 5, 5)

    def scene(self) -> object:
        return self._scene_obj

    def pen(self) -> _PenStub:
        return _PenStub(self._pen_style)

    def brush(self) -> _BrushStub:
        return _BrushStub(self._brush_style)

    def sceneBoundingRect(self) -> QRectF:
        return self._rect

    @property
    def label(self) -> _LabelStub:
        return _LabelStub(self._label_visible, self._label_rect)

    @property
    def charge_label(self) -> _LabelStub:
        return _LabelStub(self._charge_label_visible, self._charge_label_rect)


class _PenStub:
    def __init__(self, style: Qt.PenStyle) -> None:
        self._style = style

    def style(self) -> Qt.PenStyle:
        return self._style


class _BrushStub:
    def __init__(self, style: Qt.BrushStyle) -> None:
        self._style = style

    def style(self) -> Qt.BrushStyle:
        return self._style


class _LabelStub:
    def __init__(self, visible: bool, rect: QRectF) -> None:
        self._visible = visible
        self._rect = rect

    def isVisible(self) -> bool:
        return self._visible

    def sceneBoundingRect(self) -> QRectF:
        return self._rect


class _OverlayLineStub:
    def __init__(self, scene: object | None, visible: bool, rect: QRectF) -> None:
        self._scene_obj = scene
        self._visible = visible
        self._rect = rect

    def scene(self) -> object | None:
        return self._scene_obj

    def isVisible(self) -> bool:
        return self._visible

    def sceneBoundingRect(self) -> QRectF:
        return self._rect


# ---------------------------------------------------------------------------
# Tests de resolve_selected_atom_ids
# ---------------------------------------------------------------------------


def test_resolve_selected_atom_ids_with_atoms_only() -> None:
    model = _ModelStub()
    result = resolve_selected_atom_ids(
        selected_atom_ids=[1, 2, 3],
        selected_bond_ids=[],
        model=model,
    )
    assert result == {1, 2, 3}


def test_resolve_selected_atom_ids_via_bonds() -> None:
    model = _ModelStub()
    model.bonds[10] = _BondStub(4, 5)
    result = resolve_selected_atom_ids(
        selected_atom_ids=[],
        selected_bond_ids=[10],
        model=model,
    )
    assert result == {4, 5}


def test_resolve_selected_atom_ids_with_missing_bond() -> None:
    model = _ModelStub()
    result = resolve_selected_atom_ids(
        selected_atom_ids=[1],
        selected_bond_ids=[999],
        model=model,
    )
    assert result == {1}


def test_resolve_selected_atom_ids_no_mutation() -> None:
    model = _ModelStub()
    model.bonds[1] = _BondStub(10, 20)
    atoms_in = [1, 2]
    bonds_in = [1]
    atoms_copy = list(atoms_in)
    bonds_copy = list(bonds_in)
    resolve_selected_atom_ids(
        selected_atom_ids=atoms_in,
        selected_bond_ids=bonds_in,
        model=model,
    )
    assert atoms_in == atoms_copy
    assert bonds_in == bonds_copy


def test_resolve_selected_atom_ids_empty() -> None:
    model = _ModelStub()
    result = resolve_selected_atom_ids(
        selected_atom_ids=[],
        selected_bond_ids=[],
        model=model,
    )
    assert result == set()


# ---------------------------------------------------------------------------
# Tests de selection_bounds
# ---------------------------------------------------------------------------


def test_selection_bounds_empty_returns_none() -> None:
    scene = _SceneStub()
    result = selection_bounds(
        scene=scene,
        atom_items={},
        bond_items={},
        implicit_h_overlays={},
        atom_ids=[],
        bond_ids=[],
        graphic_items=[],
    )
    assert result is None


def test_selection_bounds_invalid_rects_ignored() -> None:
    """Un QRectF nulo no debe romper la unión."""
    scene = _SceneStub()

    class _NullRectItem:
        def __init__(self, scene_obj: object) -> None:
            self._scene_obj = scene_obj
        def scene(self) -> object:
            return self._scene_obj
        def sceneBoundingRect(self) -> QRectF:
            # QRectF() es nulo (x=0, y=0, w=-1, h=-1 → invalid)
            return QRectF()
        def pen(self):
            class _Pen:
                def style(self):
                    return Qt.PenStyle.SolidLine
            return _Pen()
        def brush(self):
            class _Brush:
                def style(self):
                    return Qt.BrushStyle.SolidPattern
            return _Brush()
        @property
        def label(self):
            class _Label:
                def isVisible(self): return False
                def sceneBoundingRect(self): return QRectF()
            return _Label()
        @property
        def charge_label(self):
            class _Label:
                def isVisible(self): return False
                def sceneBoundingRect(self): return QRectF()
            return _Label()

    item = _NullRectItem(scene)
    result = selection_bounds(
        scene=scene,
        atom_items={1: item},
        bond_items={},
        implicit_h_overlays={},
        atom_ids=[1],
        bond_ids=[],
        graphic_items=[],
    )
    assert result is None


def test_selection_bounds_union_of_multiple_rects() -> None:
    scene = _SceneStub()
    item1 = _ItemStub(scene=scene, rect=QRectF(0, 0, 10, 10))
    item2 = _ItemStub(scene=scene, rect=QRectF(20, 0, 10, 10))
    result = selection_bounds(
        scene=scene,
        atom_items={1: item1, 2: item2},
        bond_items={},
        implicit_h_overlays={},
        atom_ids=[1, 2],
        bond_ids=[],
        graphic_items=[],
    )
    assert result is not None
    assert result.x() == pytest.approx(0.0)
    assert result.y() == pytest.approx(0.0)
    assert result.width() == pytest.approx(30.0)
    assert result.height() == pytest.approx(10.0)


def test_selection_bounds_atom_visible_body() -> None:
    scene = _SceneStub()
    # Ocultar label y charge_label para que solo el body cuente
    item = _ItemStub(
        scene=scene,
        rect=QRectF(5, 5, 15, 15),
        label_visible=False,
        charge_label_visible=False,
    )
    result = selection_bounds(
        scene=scene,
        atom_items={1: item},
        bond_items={},
        implicit_h_overlays={},
        atom_ids=[1],
        bond_ids=[],
        graphic_items=[],
    )
    assert result is not None
    assert result.x() == pytest.approx(5.0)
    assert result.y() == pytest.approx(5.0)


def test_selection_bounds_atom_hidden_no_pen_no_brush() -> None:
    scene = _SceneStub()
    item = _ItemStub(
        scene=scene,
        pen_style=Qt.PenStyle.NoPen,
        brush_style=Qt.BrushStyle.NoBrush,
        rect=QRectF(5, 5, 15, 15),
        label_visible=False,
        charge_label_visible=False,
    )
    result = selection_bounds(
        scene=scene,
        atom_items={1: item},
        bond_items={},
        implicit_h_overlays={},
        atom_ids=[1],
        bond_ids=[],
        graphic_items=[],
    )
    assert result is None


def test_selection_bounds_label_visible_included() -> None:
    scene = _SceneStub()
    item = _ItemStub(
        scene=scene,
        rect=QRectF(0, 0, 10, 10),
        charge_label_visible=False,
        label_visible=True,
        label_rect=QRectF(10, 0, 5, 5),
    )
    result = selection_bounds(
        scene=scene,
        atom_items={1: item},
        bond_items={},
        implicit_h_overlays={},
        atom_ids=[1],
        bond_ids=[],
        graphic_items=[],
    )
    assert result is not None
    assert result.width() >= 15.0


def test_selection_bounds_charge_label_visible_included() -> None:
    scene = _SceneStub()
    item = _ItemStub(
        scene=scene,
        rect=QRectF(0, 0, 10, 10),
        label_visible=False,
        charge_label_visible=True,
        charge_label_rect=QRectF(0, 10, 5, 5),
    )
    result = selection_bounds(
        scene=scene,
        atom_items={1: item},
        bond_items={},
        implicit_h_overlays={},
        atom_ids=[1],
        bond_ids=[],
        graphic_items=[],
    )
    assert result is not None
    assert result.height() >= 15.0


def test_selection_bounds_implicit_h_visible() -> None:
    scene = _SceneStub()
    item = _ItemStub(
        scene=scene,
        rect=QRectF(0, 0, 10, 10),
        label_visible=False,
        charge_label_visible=False,
    )
    overlay_line = _OverlayLineStub(scene=scene, visible=True, rect=QRectF(20, 20, 5, 5))
    overlay_text = _OverlayLineStub(scene=scene, visible=True, rect=QRectF(25, 25, 5, 5))
    result = selection_bounds(
        scene=scene,
        atom_items={1: item},
        bond_items={},
        implicit_h_overlays={1: [(overlay_line, overlay_text)]},
        atom_ids=[1],
        bond_ids=[],
        graphic_items=[],
    )
    assert result is not None
    assert result.x() == pytest.approx(0.0)
    assert result.y() == pytest.approx(0.0)
    assert result.width() >= 30.0
    assert result.height() >= 30.0


def test_selection_bounds_overlays_hidden_or_removed() -> None:
    scene = _SceneStub()
    item = _ItemStub(
        scene=scene,
        rect=QRectF(0, 0, 10, 10),
        label_visible=False,
        charge_label_visible=False,
    )
    removed_line = _OverlayLineStub(scene=None, visible=True, rect=QRectF(100, 100, 5, 5))
    hidden_line = _OverlayLineStub(scene=scene, visible=False, rect=QRectF(200, 200, 5, 5))
    result = selection_bounds(
        scene=scene,
        atom_items={1: item},
        bond_items={},
        implicit_h_overlays={1: [(removed_line, hidden_line)]},
        atom_ids=[1],
        bond_ids=[],
        graphic_items=[],
    )
    assert result is not None
    assert result.x() == pytest.approx(0.0)
    assert result.width() == pytest.approx(10.0)


def test_selection_bounds_bond_item() -> None:
    scene = _SceneStub()
    bond_item = _ItemStub(scene=scene, rect=QRectF(30, 0, 20, 5))
    result = selection_bounds(
        scene=scene,
        atom_items={},
        bond_items={1: bond_item},
        implicit_h_overlays={},
        atom_ids=[],
        bond_ids=[1],
        graphic_items=[],
    )
    assert result is not None
    assert result.x() == pytest.approx(30.0)


def test_selection_bounds_graphic_item() -> None:
    scene = _SceneStub()
    gitem = _ItemStub(scene=scene, rect=QRectF(50, 50, 10, 10))
    result = selection_bounds(
        scene=scene,
        atom_items={},
        bond_items={},
        implicit_h_overlays={},
        atom_ids=[],
        bond_ids=[],
        graphic_items=[gitem],
    )
    assert result is not None
    assert result.x() == pytest.approx(50.0)


def test_selection_bounds_graphic_item_removed_from_scene() -> None:
    scene = _SceneStub()
    gitem = _ItemStub(scene=None, rect=QRectF(50, 50, 10, 10))
    result = selection_bounds(
        scene=scene,
        atom_items={},
        bond_items={},
        implicit_h_overlays={},
        atom_ids=[],
        bond_ids=[],
        graphic_items=[gitem],
    )
    assert result is None


def test_selection_bounds_combination() -> None:
    scene = _SceneStub()
    atom_item = _ItemStub(
        scene=scene,
        rect=QRectF(0, 0, 10, 10),
        label_visible=False,
        charge_label_visible=False,
    )
    bond_item = _ItemStub(scene=scene, rect=QRectF(20, 0, 10, 10))
    gitem = _ItemStub(scene=scene, rect=QRectF(40, 0, 10, 10))
    result = selection_bounds(
        scene=scene,
        atom_items={1: atom_item},
        bond_items={1: bond_item},
        implicit_h_overlays={},
        atom_ids=[1],
        bond_ids=[1],
        graphic_items=[gitem],
    )
    assert result is not None
    assert result.x() == pytest.approx(0.0)
    assert result.width() == pytest.approx(50.0)


def test_selection_bounds_null_rect_ignored() -> None:
    """Un QRectF nulo no debe romper la unión."""
    scene = _SceneStub()

    class _NullRectItem:
        def scene(self):
            return scene
        def sceneBoundingRect(self):
            return QRectF()

    item = _NullRectItem()
    result = selection_bounds(
        scene=scene,
        atom_items={},
        bond_items={},
        implicit_h_overlays={},
        atom_ids=[],
        bond_ids=[],
        graphic_items=[item],
    )
    assert result is None


# ---------------------------------------------------------------------------
# Regresiones: resolve_selected_atom_ids con solo enlaces
# ---------------------------------------------------------------------------


def test_resolve_selected_atom_ids_bond_only_includes_end_atoms() -> None:
    """Seleccionar solo un enlace debe incluir sus dos átomos extremos."""
    model = _ModelStub()
    model.bonds[10] = _BondStub(100, 200)
    result = resolve_selected_atom_ids(
        selected_atom_ids=[],
        selected_bond_ids=[10],
        model=model,
    )
    assert result == {100, 200}


def test_resolve_selected_atom_ids_mixed_bond_and_atoms() -> None:
    """Selección mixta de átomos y enlaces debe unir ambos."""
    model = _ModelStub()
    model.bonds[10] = _BondStub(100, 200)
    model.bonds[20] = _BondStub(300, 400)
    result = resolve_selected_atom_ids(
        selected_atom_ids=[1, 2],
        selected_bond_ids=[10, 20],
        model=model,
    )
    assert result == {1, 2, 100, 200, 300, 400}


def test_resolve_selected_atom_ids_bond_nonexistent_ignored() -> None:
    """Enlaces inexistentes en la selección no deben afectar el resultado."""
    model = _ModelStub()
    model.bonds[10] = _BondStub(100, 200)
    result = resolve_selected_atom_ids(
        selected_atom_ids=[],
        selected_bond_ids=[10, 999, 888],
        model=model,
    )
    assert result == {100, 200}


# ---------------------------------------------------------------------------
# Regresiones: selección_bounds numérica
# ---------------------------------------------------------------------------


def test_selection_bounds_numerical_union_precision() -> None:
    """La unión de rects debe producir valores numéricos exactos."""
    scene = _SceneStub()
    item1 = _ItemStub(scene=scene, rect=QRectF(0.0, 0.0, 10.5, 20.5))
    item2 = _ItemStub(scene=scene, rect=QRectF(10.5, 0.0, 15.3, 10.0))
    result = selection_bounds(
        scene=scene,
        atom_items={1: item1, 2: item2},
        bond_items={},
        implicit_h_overlays={},
        atom_ids=[1, 2],
        bond_ids=[],
        graphic_items=[],
    )
    assert result is not None
    assert result.x() == pytest.approx(0.0, abs=1e-9)
    assert result.y() == pytest.approx(0.0, abs=1e-9)
    assert result.width() == pytest.approx(25.8, abs=1e-9)
    assert result.height() == pytest.approx(20.5, abs=1e-9)


def test_selection_bounds_bond_item_only() -> None:
    """Solo bond items debe devolver el bbox del enlace."""
    scene = _SceneStub()
    bond_item = _ItemStub(scene=scene, rect=QRectF(50.0, 60.0, 25.0, 15.0))
    result = selection_bounds(
        scene=scene,
        atom_items={},
        bond_items={1: bond_item},
        implicit_h_overlays={},
        atom_ids=[],
        bond_ids=[1],
        graphic_items=[],
    )
    assert result is not None
    assert result.x() == pytest.approx(50.0)
    assert result.y() == pytest.approx(60.0)
    assert result.width() == pytest.approx(25.0)
    assert result.height() == pytest.approx(15.0)


def test_selection_bounds_graphic_items_only() -> None:
    """Solo graphic items debe devolver el bbox unido."""
    scene = _SceneStub()
    g1 = _ItemStub(scene=scene, rect=QRectF(0, 0, 10, 10))
    g2 = _ItemStub(scene=scene, rect=QRectF(10, 10, 20, 20))
    result = selection_bounds(
        scene=scene,
        atom_items={},
        bond_items={},
        implicit_h_overlays={},
        atom_ids=[],
        bond_ids=[],
        graphic_items=[g1, g2],
    )
    assert result is not None
    assert result.x() == pytest.approx(0.0)
    assert result.y() == pytest.approx(0.0)
    assert result.width() == pytest.approx(30.0)
    assert result.height() == pytest.approx(30.0)


# ---------------------------------------------------------------------------
# Regresiones: wrappers reales de CanvasSelectionMixin
# ---------------------------------------------------------------------------


def _make_mock_graphic_item(item_class: type, scene_obj: object, rect: QRectF) -> MagicMock:
    """Crear un MagicMock que pase isinstance(item_class) con un rect dado."""
    item = MagicMock()
    item.__class__ = item_class
    item.scene.return_value = scene_obj
    item.sceneBoundingRect.return_value = rect
    return item


def test_wrapper_selected_atom_ids_bond_only_includes_end_atoms() -> None:
    """Seleccionar solo un enlace en el wrapper real debe incluir sus dos átomos extremos."""
    mixin = CanvasSelectionMixin()
    mixin.model = _ModelStub()
    mixin.model.bonds[10] = _BondStub(100, 200)
    mixin.state = MagicMock()
    mixin.state.selected_atoms = []
    mixin.state.selected_bonds = [10]
    result = mixin._selected_atom_ids_for_transform()
    assert result == {100, 200}


def test_wrapper_selected_items_bbox_includes_selected_annotations() -> None:
    """_selected_items_bbox del wrapper real incluye TextAnnotationItem seleccionado."""
    mixin = CanvasSelectionMixin()
    mixin.model = _ModelStub()
    mixin.scene = MagicMock()
    mixin.atom_items = {}
    mixin.bond_items = {}
    mixin._implicit_h_overlays = {}
    mixin.state = MagicMock()
    mixin.state.selected_atoms = []
    mixin.state.selected_bonds = []

    # Crear una anotación que pase isinstance(TextAnnotationItem)
    # y que devuelva el mismo objeto mixin.scene desde .scene() para el check "is"
    annotation = _make_mock_graphic_item(TextAnnotationItem, mixin.scene, QRectF(100, 100, 30, 15))
    mixin.scene.selectedItems.return_value = [annotation]

    bbox = mixin._selected_items_bbox()
    assert bbox is not None
    assert bbox.x() == pytest.approx(100.0)
    assert bbox.width() == pytest.approx(30.0)

    # También debe incluir un átomo seleccionado además de la anotación
    atom_item = _ItemStub(scene=mixin.scene, rect=QRectF(0, 0, 10, 10))
    mixin.atom_items = {1: atom_item}
    mixin.state.selected_atoms = [1]
    mixin.state.selected_bonds = []
    mixin.scene.selectedItems.return_value = [annotation]

    bbox = mixin._selected_items_bbox()
    assert bbox is not None
    assert bbox.x() == pytest.approx(0.0)
    assert bbox.width() == pytest.approx(130.0)


def test_wrapper_targets_bbox_preserves_all_groups() -> None:
    """_targets_bbox del wrapper real conserva todos los grupos explícitos."""
    mixin = CanvasSelectionMixin()
    mixin.scene = MagicMock()
    mixin.atom_items = {}
    mixin.bond_items = {}
    mixin._implicit_h_overlays = {}

    text_item = _make_mock_graphic_item(TextAnnotationItem, mixin.scene, QRectF(0, 0, 10, 10))
    arrow_item = _make_mock_graphic_item(ArrowItem, mixin.scene, QRectF(10, 0, 20, 5))
    bracket_item = _make_mock_graphic_item(BracketItem, mixin.scene, QRectF(30, 0, 15, 10))
    energy_item = _make_mock_graphic_item(EnergyDiagramItem, mixin.scene, QRectF(45, 0, 25, 15))
    semantic_item = _make_mock_graphic_item(CompositeDiagramItem, mixin.scene, QRectF(70, 0, 10, 10))
    orbital_item = _make_mock_graphic_item(OrbitalAnnotationItem, mixin.scene, QRectF(80, 0, 15, 10))
    image_item = _make_mock_graphic_item(ImageAnnotationItem, mixin.scene, QRectF(95, 0, 20, 20))

    result = mixin._targets_bbox(
        atom_ids=[],
        bond_ids=[],
        text_items=[text_item],
        arrow_items=[arrow_item],
        bracket_items=[bracket_item],
        energy_diagram_items=[energy_item],
        semantic_diagram_items=[semantic_item],
        orbital_items=[orbital_item],
        image_items=[image_item],
    )
    assert result is not None
    assert result.x() == pytest.approx(0.0)
    assert result.width() == pytest.approx(115.0)


def test_wrapper_results_match_selection_bounds_direct() -> None:
    """Los resultados del wrapper real deben coincidir con llamada directa a selection_bounds."""
    scene = MagicMock()
    atom_item = _ItemStub(scene=scene, rect=QRectF(5, 5, 20, 20))
    bond_item = _ItemStub(scene=scene, rect=QRectF(25, 5, 10, 10))

    direct_result = selection_bounds(
        scene=scene,
        atom_items={1: atom_item},
        bond_items={1: bond_item},
        implicit_h_overlays={},
        atom_ids=[1],
        bond_ids=[1],
        graphic_items=[],
    )

    mixin = CanvasSelectionMixin()
    mixin.model = _ModelStub()
    mixin.scene = scene
    mixin.atom_items = {1: atom_item}
    mixin.bond_items = {1: bond_item}
    mixin._implicit_h_overlays = {}
    mixin.state = MagicMock()
    mixin.state.selected_atoms = [1]
    mixin.state.selected_bonds = [1]
    wrapper_result = mixin._selected_items_bbox()

    assert direct_result is not None
    assert wrapper_result is not None
    assert direct_result.x() == pytest.approx(wrapper_result.x())
    assert direct_result.y() == pytest.approx(wrapper_result.y())
    assert direct_result.width() == pytest.approx(wrapper_result.width())
    assert direct_result.height() == pytest.approx(wrapper_result.height())
