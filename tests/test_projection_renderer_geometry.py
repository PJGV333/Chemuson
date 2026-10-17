from __future__ import annotations

from PyQt6.QtWidgets import QApplication

from chemuson.core.model import BondStereo
from chemuson.gui.projections import (
    ChairParameters,
    ChairSubstituent,
    FischerParameters,
    build_chair_geometry,
    build_fischer_geometry,
)
from tools.visual.depiction_lab_app.renderer import ProductionCanvasAdapter


def _assert_visible(adapter: ProductionCanvasAdapter) -> None:
    canvas = adapter.canvas
    target = canvas.mapFromScene(adapter.scene_rect.center())
    center = canvas.viewport().rect().center()
    assert abs(target.x() - center.x()) <= 10
    assert abs(target.y() - center.y()) <= 10
    image = canvas.viewport().grab().toImage()
    assert image.width() >= 380 and image.height() >= 300
    foreground = 0
    for y in range(image.height()):
        for x in range(image.width()):
            color = image.pixelColor(x, y)
            foreground += max(color.red(), color.green(), color.blue()) < 170
    assert foreground > 20
    assert {bond.stereo for bond in canvas.model.bonds.values()} == {BondStereo.NONE}
    bond_shapes = {
        bond_id: item.mapToScene(item.shape())
        for bond_id, item in canvas.bond_items.items()
    }
    label_shapes = {
        atom_id: item.label.mapToScene(item.label.shape())
        for atom_id, item in canvas.atom_items.items()
        if item.label.isVisible() and item.label.toPlainText().strip()
    }
    assert all(
        not bond_shape.intersects(label_shape)
        for bond_shape in bond_shapes.values()
        for label_shape in label_shapes.values()
    )
    bonds = list(canvas.model.bonds.values())
    for index, bond in enumerate(bonds):
        for other in bonds[index + 1 :]:
            if {bond.a1_id, bond.a2_id}.isdisjoint({other.a1_id, other.a2_id}):
                assert not bond_shapes[bond.id].intersects(bond_shapes[other.id])


def test_fischer_convention_renders_through_production_canvas_at_multiple_scales() -> None:
    app = QApplication.instance() or QApplication([])
    geometry = build_fischer_geometry(
        FischerParameters(
            substituent_sides=("right", "left", "right", "left"),
            center_configurations=(None, None, None, None),
        )
    )
    adapter = ProductionCanvasAdapter.create(geometry.to_projection_instance().render_graph)
    adapter.canvas.resize(700, 520)
    adapter.canvas.show()
    app.processEvents()

    try:
        assert geometry.semantics.horizontal_toward_observer
        assert geometry.semantics.vertical_away_from_observer
        for zoom in (0.75, 1.0, 1.2):
            adapter.set_view(zoom, 0.0)
            app.processEvents()
            _assert_visible(adapter)
    finally:
        adapter.close()
        app.processEvents()


def test_chair_flip_semantics_render_with_explicit_axial_equatorial_relations() -> None:
    app = QApplication.instance() or QApplication([])
    parameters = ChairParameters(
        substituents=(ChairSubstituent(0, "OH", "up", "axial"),)
    )
    geometries = (build_chair_geometry(parameters), build_chair_geometry(parameters.flipped()))
    adapters = [
        ProductionCanvasAdapter.create(geometry.to_projection_instance().render_graph)
        for geometry in geometries
    ]
    try:
        for adapter in adapters:
            adapter.canvas.resize(680, 500)
            adapter.canvas.show()
        app.processEvents()
        assert geometries[0].semantics.conformation == "chair_a"
        assert geometries[1].semantics.conformation == "chair_b"
        assert [item.orientation for item in geometries[0].semantics.substituents] == [
            item.orientation for item in geometries[1].semantics.substituents
        ]
        assert [item.placement for item in geometries[1].semantics.substituents] == [
            "equatorial"
        ]
        for zoom, angle in ((0.75, 0.0), (1.0, 0.0), (1.25, 25.0)):
            for adapter in adapters:
                adapter.set_view(zoom, angle)
            app.processEvents()
            for adapter in adapters:
                _assert_visible(adapter)
    finally:
        for adapter in adapters:
            adapter.close()
        app.processEvents()
