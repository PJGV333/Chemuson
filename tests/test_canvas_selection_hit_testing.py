"""Functional regressions for consultative canvas selection hit testing."""

from __future__ import annotations

from dataclasses import dataclass

from PyQt6.QtCore import QPointF, QRectF

from chemuson.gui.canvas.selection_hit_testing import (
    get_item_at,
    resolve_click_item,
    selected_annotation_item_at,
    semantic_diagram_parent,
)


class Composite:
    pass


class Atom:
    def __init__(self, atom_id: int) -> None:
        self.atom_id = atom_id

    def parentItem(self):
        return None


class Bond:
    def parentItem(self):
        return None


class Text:
    def __init__(self, parent=None) -> None:
        self._parent = parent

    def parentItem(self):
        return self._parent


class Node:
    def __init__(self, parent=None) -> None:
        self._parent = parent

    def parentItem(self):
        return self._parent


@dataclass
class Scene:
    scene_items: list[object]

    def items(self, _point: QPointF) -> list[object]:
        return list(self.scene_items)


class Annotation:
    def __init__(
        self,
        scene,
        rect: QRectF,
        z: float,
        *,
        visible: bool = True,
        contains: bool = True,
    ) -> None:
        self._scene = scene
        self._rect = rect
        self._z = z
        self._visible = visible
        self._contains = contains

    def scene(self):
        return self._scene

    def isVisible(self):
        return self._visible

    def sceneBoundingRect(self):
        return self._rect

    def mapFromScene(self, _point: QPointF):
        return QPointF(0.0, 0.0)

    def contains(self, _point: QPointF):
        return self._contains

    def zValue(self):
        return self._z


def _get_item_at(scene, point):
    return get_item_at(
        scene=scene,
        scene_pos=point,
        composite_type=Composite,
        atom_type=Atom,
        text_type=Text,
        selectable_types=(Atom, Bond),
        is_disposable_orphan_atom=lambda atom_id: atom_id == 99,
    )


def test_semantic_parent_promotes_descendant_and_root() -> None:
    root = Composite()
    child = Node(Node(root))

    assert semantic_diagram_parent(child, composite_type=Composite) is root
    assert semantic_diagram_parent(root, composite_type=Composite) is root
    assert semantic_diagram_parent(Node(), composite_type=Composite) is None


def test_scene_query_maps_atom_text_child_to_atom() -> None:
    atom = Atom(1)
    assert _get_item_at(Scene([Text(atom)]), QPointF()) is atom


def test_scene_query_skips_disposable_orphan_atom() -> None:
    assert _get_item_at(Scene([Atom(99)]), QPointF()) is None


def test_scene_query_promotes_composite_before_selectable_item() -> None:
    root = Composite()
    child = Node(root)
    assert _get_item_at(Scene([child]), QPointF()) is root


def test_scene_query_returns_selectable_item_and_none_for_unknown() -> None:
    bond = Bond()
    assert _get_item_at(Scene([bond]), QPointF()) is bond
    assert _get_item_at(Scene([Node()]), QPointF()) is None


def test_resolve_click_prioritizes_annotation_then_atom_then_scene_then_bond() -> None:
    annotation = object()
    scene_item = object()
    atom_item = object()
    bond_item = object()

    assert resolve_click_item(
        scene_item=annotation,
        scene_pos=QPointF(),
        annotation_types=(object,),
        pick_hover_target=lambda _point: (1, 2),
        atom_items={1: atom_item},
        bond_items={2: bond_item},
    ) is annotation
    assert resolve_click_item(
        scene_item=scene_item,
        scene_pos=QPointF(),
        annotation_types=(),
        pick_hover_target=lambda _point: (1, 2),
        atom_items={1: atom_item},
        bond_items={2: bond_item},
    ) is atom_item
    assert resolve_click_item(
        scene_item=scene_item,
        scene_pos=QPointF(),
        annotation_types=(),
        pick_hover_target=lambda _point: (None, 2),
        atom_items={},
        bond_items={2: bond_item},
    ) is scene_item
    assert resolve_click_item(
        scene_item=None,
        scene_pos=QPointF(),
        annotation_types=(),
        pick_hover_target=lambda _point: (None, 2),
        atom_items={},
        bond_items={2: bond_item},
    ) is bond_item


def test_selected_annotation_requires_scene_visibility_bounds_and_contains() -> None:
    active_scene = object()
    point = QPointF(5.0, 5.0)
    winner = Annotation(active_scene, QRectF(0, 0, 10, 10), 4)
    lower = Annotation(active_scene, QRectF(0, 0, 10, 10), 2)
    removed = Annotation(object(), QRectF(0, 0, 10, 10), 9)
    hidden = Annotation(active_scene, QRectF(0, 0, 10, 10), 9, visible=False)
    miss = Annotation(active_scene, QRectF(20, 20, 10, 10), 9)
    contains_miss = Annotation(active_scene, QRectF(0, 0, 10, 10), 9, contains=False)

    result = selected_annotation_item_at(
        scene=active_scene,
        scene_pos=point,
        items=[lower, winner, removed, hidden, miss, contains_miss],
    )

    assert result is winner


def test_selected_annotation_returns_none_without_valid_target() -> None:
    assert selected_annotation_item_at(scene=object(), scene_pos=QPointF(), items=[]) is None
