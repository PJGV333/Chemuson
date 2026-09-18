"""Consultative hit-testing policy for canvas selection.

The module deliberately receives item classes and callbacks from its owner. It
knows how to query scene/items, but it does not own event dispatch, molecular
picking, commands, or scene mutation.
"""

from __future__ import annotations

from collections.abc import Callable, Iterable, Mapping


def semantic_diagram_parent(item: object | None, *, composite_type: type) -> object | None:
    """Return the nearest semantic diagram root in an item's parent chain."""
    current = item
    while current is not None:
        if isinstance(current, composite_type):
            return current
        current = current.parentItem()
    return None


def get_item_at(
    *,
    scene: object,
    scene_pos: object,
    composite_type: type,
    atom_type: type,
    text_type: type,
    selectable_types: tuple[type, ...],
    is_disposable_orphan_atom: Callable[[int], bool],
) -> object | None:
    """Resolve a selectable item under a scene point in scene order."""
    for item in scene.items(scene_pos):
        semantic_parent = semantic_diagram_parent(item, composite_type=composite_type)
        if semantic_parent is not None:
            return semantic_parent
        if isinstance(item, atom_type) and is_disposable_orphan_atom(item.atom_id):
            continue
        if isinstance(item, selectable_types):
            return item
        if isinstance(item, text_type):
            parent = item.parentItem()
            if isinstance(parent, atom_type):
                if is_disposable_orphan_atom(parent.atom_id):
                    continue
                return parent
            semantic_parent = semantic_diagram_parent(parent, composite_type=composite_type)
            if semantic_parent is not None:
                return semantic_parent
    return None


def resolve_click_item(
    *,
    scene_item: object | None,
    scene_pos: object,
    annotation_types: tuple[type, ...],
    pick_hover_target: Callable[[object], tuple[int | None, int | None]],
    atom_items: Mapping[int, object],
    bond_items: Mapping[int, object],
) -> object | None:
    """Resolve a click using the canvas's existing priority policy."""
    if isinstance(scene_item, annotation_types):
        return scene_item

    atom_id, bond_id = pick_hover_target(scene_pos)
    if atom_id is not None:
        atom_item = atom_items.get(atom_id)
        if atom_item is not None:
            return atom_item

    if scene_item is not None:
        return scene_item

    if bond_id is not None:
        return bond_items.get(bond_id)
    return None


def selected_annotation_item_at(
    *,
    scene: object,
    scene_pos: object,
    items: Iterable[object],
) -> object | None:
    """Return the topmost selected annotation containing a scene point."""
    best_item = None
    best_z = float("-inf")
    for item in items:
        if item.scene() is not scene:
            continue
        try:
            if not item.isVisible():
                continue
            if not item.sceneBoundingRect().contains(scene_pos):
                continue
            local_pos = item.mapFromScene(scene_pos)
            if not item.contains(local_pos):
                continue
            z_value = float(item.zValue())
        except RuntimeError:
            continue
        if best_item is None or z_value >= best_z:
            best_item = item
            best_z = z_value
    return best_item
