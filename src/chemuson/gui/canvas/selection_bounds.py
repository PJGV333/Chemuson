"""Lógica de consulta de límites visuales de selección del canvas.

Este módulo proporciona funciones puras de consulta para resolver el conjunto
de IDs de átomos seleccionados y calcular los límites visuales de una
selección en el canvas. No importa QtGui, QtWidgets ni ningún módulo de
renderizado, comandos, diálogos o controladores.
"""

from __future__ import annotations

from collections.abc import Iterable, Mapping

from PyQt6.QtCore import QRectF, Qt


def resolve_selected_atom_ids(
    *,
    selected_atom_ids: Iterable[int],
    selected_bond_ids: Iterable[int],
    model: object,
) -> set[int]:
    """Resolver el conjunto de IDs de átomos para transformaciones.

    Copia los IDs de átomos directamente seleccionados y añade los extremos
    de cada enlace seleccionado que todavía exista en el modelo.

    Args:
        selected_atom_ids: IDs de átomos directamente seleccionados.
        selected_bond_ids: IDs de enlaces seleccionados.
        model: Objeto con atributo ``bonds`` mapeable de
            ``{bond_id: bond}`` donde cada ``bond`` tiene ``a1_id`` y
            ``a2_id``.

    Returns:
        Conjunto estable de IDs de átomos.
    """
    result: set[int] = set(selected_atom_ids)
    for bond_id in selected_bond_ids:
        if bond_id in model.bonds:
            bond = model.get_bond(bond_id)
            result.add(bond.a1_id)
            result.add(bond.a2_id)
    return result


def selection_bounds(
    *,
    scene: object,
    atom_items: Mapping[int, object],
    bond_items: Mapping[int, object],
    implicit_h_overlays: Mapping[int, Iterable[tuple[object, object]]],
    atom_ids: Iterable[int] = (),
    bond_ids: Iterable[int] = (),
    graphic_items: Iterable[object] = (),
) -> QRectF | None:
    """Calcular el bounding box visual de una selección.

    Reúne límites de átomos (cuerpo, etiqueta, carga, hidrógenos implícitos),
    enlaces y graphic items, uniendo únicamente los QRectF válidos y no nulos.

    Args:
        scene: Objeto escena Qt para verificar pertenencia.
        atom_items: Mapeo de atom_id al item gráfico correspondiente.
        bond_items: Mapeo de bond_id al item de enlace.
        implicit_h_overlays: Mapeo de atom_id a pares (line_item, text_item).
        atom_ids: IDs de átomos a incluir.
        bond_ids: IDs de enlaces a incluir.
        graphic_items: Items gráficos adicionales a incluir.

    Returns:
        QRectF con el bounding box unido, o None si no hay límites válidos.
    """
    rect: QRectF | None = None

    def _extend(candidate: QRectF) -> None:
        """Amiar ``rect`` con ``candidate`` si es válido y no nulo."""
        nonlocal rect
        if not candidate.isValid() or candidate.isNull():
            return
        if rect is None:
            rect = candidate
        else:
            rect = rect.united(candidate)

    def _extend_atom_bounds(atom_id: int) -> None:
        """Recopilar límites de un átomo y sus subelementos visibles."""
        item = atom_items.get(atom_id)
        if item is None:
            return
        if item.scene() is not scene:
            return

        # Cuerpo del átomo: solo si tiene pen o brush visibles
        pen_style = item.pen().style()
        brush_style = item.brush().style()

        if pen_style != Qt.PenStyle.NoPen or brush_style != Qt.BrushStyle.NoBrush:
            _extend(item.sceneBoundingRect())

        # Etiqueta visible
        label = item.label
        if label is not None and label.isVisible():
            _extend(label.sceneBoundingRect())

        # Etiqueta de carga visible
        charge_label = item.charge_label
        if charge_label is not None and charge_label.isVisible():
            _extend(charge_label.sceneBoundingRect())

        # Hidrógenos implícitos
        overlays = implicit_h_overlays.get(atom_id)
        if overlays:
            for line_item, text_item in overlays:
                if line_item.scene() is scene and line_item.isVisible():
                    _extend(line_item.sceneBoundingRect())
                if text_item.scene() is scene and text_item.isVisible():
                    _extend(text_item.sceneBoundingRect())

    # Átomos
    for atom_id in atom_ids:
        _extend_atom_bounds(atom_id)

    # Enlaces
    for bond_id in bond_ids:
        bond_item = bond_items.get(bond_id)
        if bond_item is not None:
            _extend(bond_item.sceneBoundingRect())

    # Graphic items adicionales
    for item in graphic_items:
        if item.scene() is scene:
            _extend(item.sceneBoundingRect())

    return rect
