"""Adapter that calibrates and exports real ChemusonCanvas graphics items."""
from __future__ import annotations

from copy import deepcopy
from dataclasses import dataclass

from PyQt6.QtCore import QRectF, Qt
from PyQt6.QtGui import QBrush, QColor, QImage, QPainter, QPen, QTransform
from PyQt6.QtWidgets import QApplication, QGraphicsItem, QGraphicsRectItem

from chemuson.core.model import MolGraph
from chemuson.gui.canvas import ChemusonCanvas
from chemuson.gui.style import CHEMDOODLE_LIKE

from .parameters import DepictionStyleParameters
from .scenes import clone_graph, graph_snapshot

PAPER_COLOR = QColor("#FFFFFF")
THEME_BACKGROUNDS = {"light": QColor("#F1F5F9"), "dark": QColor("#111827")}


def graph_scene_rect(graph: MolGraph, padding: float = 56.0) -> QRectF:
    if not graph.atoms:
        return QRectF(-padding, -padding, padding * 2, padding * 2)
    xs = [float(atom.x) for atom in graph.atoms.values()]
    ys = [float(atom.y) for atom in graph.atoms.values()]
    return QRectF(
        min(xs) - padding,
        min(ys) - padding,
        max(1.0, max(xs) - min(xs) + padding * 2),
        max(1.0, max(ys) - min(ys) + padding * 2),
    )


@dataclass
class ProductionCanvasAdapter:
    """An isolated production canvas and its graph snapshot."""

    canvas: ChemusonCanvas
    source_snapshot: tuple
    paper_item: QGraphicsRectItem
    scene_rect: QRectF
    theme: str
    diagnostic_items: list[QGraphicsItem] | None = None

    @classmethod
    def create(
        cls,
        graph: MolGraph,
        *,
        parameters: DepictionStyleParameters | None = None,
        theme: str = "light",
        parent=None,
    ) -> "ProductionCanvasAdapter":
        if theme not in THEME_BACKGROUNDS:
            raise ValueError(f"unsupported theme: {theme}")
        app = QApplication.instance()
        if app is None:
            raise RuntimeError("QApplication must exist before creating a depiction canvas")
        copied_graph = clone_graph(graph)
        before = graph_snapshot(copied_graph)
        canvas = ChemusonCanvas(parent)
        canvas.setInteractive(False)
        canvas.setFocusPolicy(Qt.FocusPolicy.NoFocus)
        canvas.setHorizontalScrollBarPolicy(Qt.ScrollBarPolicy.ScrollBarAsNeeded)
        canvas.setVerticalScrollBarPolicy(Qt.ScrollBarPolicy.ScrollBarAsNeeded)
        canvas.setAlignment(Qt.AlignmentFlag.AlignCenter)
        canvas.model = copied_graph
        canvas.drawing_style = (
            parameters or DepictionStyleParameters.current()
        ).to_drawing_style(CHEMDOODLE_LIKE)
        if getattr(canvas, "paper", None) is not None:
            canvas.paper.setVisible(False)
        canvas.show_grid = False
        canvas.show_rulers = False
        for item_name in ("_grid_minor_item", "_grid_major_item"):
            item = getattr(canvas, item_name, None)
            if item is not None:
                item.setVisible(False)
        for item in getattr(canvas, "_numbering_overlay_items", ()):
            item.setVisible(False)

        # Production validation currently fills an `implicit_h` cache while
        # rebuilding items. Restore every model field immediately: this lab is
        # presentation-only and must not retain that canvas-side mutation.
        atom_state = {atom_id: deepcopy(vars(atom)) for atom_id, atom in copied_graph.atoms.items()}
        bond_state = {bond_id: deepcopy(vars(bond)) for bond_id, bond in copied_graph.bonds.items()}
        canvas._rebuild_items_from_model()
        for atom_id, state in atom_state.items():
            atom = copied_graph.atoms[atom_id]
            vars(atom).clear()
            vars(atom).update(state)
        for bond_id, state in bond_state.items():
            bond = copied_graph.bonds[bond_id]
            vars(bond).clear()
            vars(bond).update(state)
        rect = graph_scene_rect(copied_graph)
        canvas.scene.setSceneRect(rect)
        canvas.scene.setBackgroundBrush(QBrush(THEME_BACKGROUNDS[theme]))
        paper = QGraphicsRectItem(rect)
        paper.setBrush(QBrush(PAPER_COLOR))
        paper.setPen(QPen(Qt.PenStyle.NoPen))
        paper.setZValue(-1000.0)
        paper.setAcceptedMouseButtons(Qt.MouseButton.NoButton)
        paper.setFlag(QGraphicsItem.GraphicsItemFlag.ItemIsSelectable, False)
        canvas.scene.addItem(paper)
        canvas.centerOn(rect.center())
        adapter = cls(canvas, before, paper, rect, theme)
        adapter.assert_unchanged()
        return adapter

    def assert_unchanged(self) -> None:
        if graph_snapshot(self.canvas.model) != self.source_snapshot:
            raise AssertionError("visual renderer modified its isolated chemical graph")

    def apply_parameters(self, parameters: DepictionStyleParameters) -> None:
        """Refresh actual BondItem/AtomItem styles without rebuilding the scene."""
        style = parameters.to_drawing_style(CHEMDOODLE_LIKE)
        canvas = self.canvas
        canvas.drawing_style = style
        for item in canvas.atom_items.values():
            item.set_style(style)

        ring_pairs = canvas._compute_ring_bond_pairs()
        for bond_id, item in canvas.bond_items.items():
            bond = canvas.model.get_bond(bond_id)
            atom1 = canvas.model.get_atom(bond.a1_id)
            atom2 = canvas.model.get_atom(bond.a2_id)
            item.set_style(style, atom1, atom2)
            canvas._configure_bond_rendering(bond, item)
            canvas._set_bond_item_join_context(bond, item, ring_pairs)
            item.set_endpoint_trim(
                canvas._bond_endpoint_trim(bond, bond.a1_id),
                canvas._bond_endpoint_trim(bond, bond.a2_id),
            )
            item.set_endpoint_extend(
                canvas._bond_endpoint_extend(bond, bond.a1_id),
                canvas._bond_endpoint_extend(bond, bond.a2_id),
            )
            item.update_positions(atom1, atom2)
        canvas._update_bond_label_shrinks(set(canvas.model.atoms))
        self.assert_unchanged()

    def set_theme(self, theme: str) -> None:
        if theme not in THEME_BACKGROUNDS:
            raise ValueError(f"unsupported theme: {theme}")
        self.theme = theme
        self.canvas.scene.setBackgroundBrush(QBrush(THEME_BACKGROUNDS[theme]))

    def set_view(self, zoom: float, angle: float) -> None:
        transform = QTransform()
        transform.rotate(float(angle))
        transform.scale(float(zoom), float(zoom))
        self.canvas.setTransform(transform)
        self.canvas.centerOn(self.scene_rect.center())

    def render_image(
        self,
        *,
        width: int = 1200,
        height: int = 820,
        include_diagnostics: bool = False,
    ) -> QImage:
        """Render the existing scene at reproducible dimensions."""
        hidden_items: list[QGraphicsItem] = []
        if not include_diagnostics:
            hidden_items = list(self.diagnostic_items or ())
            for item in hidden_items:
                item.setVisible(False)
        image = QImage(width, height, QImage.Format.Format_ARGB32_Premultiplied)
        image.fill(THEME_BACKGROUNDS[self.theme])
        painter = QPainter(image)
        painter.setRenderHint(QPainter.RenderHint.Antialiasing, True)
        painter.setRenderHint(QPainter.RenderHint.TextAntialiasing, True)
        self.canvas.scene.render(
            painter,
            QRectF(0.0, 0.0, float(width), float(height)),
            self.scene_rect,
            Qt.AspectRatioMode.KeepAspectRatio,
        )
        painter.end()
        for item in hidden_items:
            item.setVisible(True)
        self.assert_unchanged()
        return image

    def close(self) -> None:
        self.canvas.close()
        self.canvas.deleteLater()


def create_production_canvas(
    graph: MolGraph,
    parameters: DepictionStyleParameters | None = None,
    *,
    theme: str = "light",
    parent=None,
) -> ProductionCanvasAdapter:
    """Convenience entry point used by UI, batch export, and tests."""
    return ProductionCanvasAdapter.create(
        graph,
        parameters=parameters,
        theme=theme,
        parent=parent,
    )
