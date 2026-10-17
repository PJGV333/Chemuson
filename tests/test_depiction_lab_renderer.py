from __future__ import annotations

import hashlib

import pytest


from chemuson.core.model import BondStyle
from PyQt6.QtGui import QImage

from tools.visual.depiction_lab_app.diagnostics import GeometryDiagnostics
from tools.visual.depiction_lab_app.parameters import DepictionStyleParameters
from tools.visual.depiction_lab_app.renderer import ProductionCanvasAdapter
from tools.visual.depiction_lab_app.scenes import CASES, case_graph, graph_snapshot


def _signature(item):
    path = item.path()
    points = tuple(
        (round(path.elementAt(index).x, 6), round(path.elementAt(index).y, 6))
        for index in range(path.elementCount())
    )
    pen = item.pen()
    return points, round(pen.widthF(), 6), pen.capStyle(), pen.joinStyle()


def _image_digest(image) -> str:
    image = image.convertToFormat(QImage.Format.Format_ARGB32)
    bits = image.bits()
    bits.setsize(image.sizeInBytes())
    return hashlib.sha256(bytes(bits)).hexdigest()


def _find_item(adapter, *, style=None, order=None):
    for bond_id, bond in adapter.canvas.model.bonds.items():
        item = adapter.canvas.bond_items[bond_id]
        if style is not None and bond.style != style:
            continue
        if order is not None and bond.order != order:
            continue
        return item
    raise AssertionError(f"fixture has no bond matching style={style!r}, order={order!r}")


@pytest.mark.parametrize("case_id", [case.case_id for case in CASES])
def test_gallery_case_builders_are_deterministic(case_id) -> None:
    first = case_graph(case_id)
    second = case_graph(case_id)
    assert graph_snapshot(first) == graph_snapshot(second)


def test_gallery_has_twelve_stable_cases_and_marks_geometric_fixture() -> None:
    ids = [case.case_id for case in CASES]
    assert len(ids) >= 12
    assert len(set(ids)) == len(ids)
    assert ids == [f"A{index}" for index in range(1, 5)] + [f"B{index}" for index in range(1, 5)] + [f"C{index}" for index in range(1, 5)]
    eligible = [case.case_id for case in CASES if case.editorial_calibration_eligible]
    assert eligible == ["A1", "A2", "A3", "A4", "B1", "B2", "B3"]
    b4 = next(case for case in CASES if case.case_id == "B4")
    assert b4.kind == "geometric-fixture"
    assert not b4.editorial_calibration_eligible
    assert "geometric fixture" in b4.editorial_exclusion_reason
    assert "non-editorial" in b4.editorial_exclusion_reason
    assert "excluded from preset selection / style calibration" in b4.editorial_exclusion_reason


@pytest.mark.parametrize(
    "case_id,bond_style,bond_order,key,value",
    (
        ("B3", BondStyle.PLAIN, 2, "bond_stroke_ratio", 0.09),
        ("B1", BondStyle.PLAIN, 2, "double_bond_offset_ratio", 0.17),
        ("B1", BondStyle.PLAIN, 2, "double_bond_inner_trim_ratio", 0.31),
        ("B4", BondStyle.BOLD, 2, "bold_double_bond_offset_ratio", 0.32),
        ("B4", BondStyle.PLAIN, 3, "triple_bond_offset_ratio", 0.33),
        ("A1", BondStyle.WEDGE, 1, "solid_wedge_width_ratio", 0.34),
        ("B4", BondStyle.WEDGE, 1, "solid_wedge_tip_width_ratio", 0.38),
        ("A3", BondStyle.HASHED, 1, "hashed_wedge_count", 9),
        ("A3", BondStyle.HASHED, 1, "hashed_wedge_min_width_ratio", 0.09),
        ("A3", BondStyle.HASHED, 1, "hashed_wedge_max_width_ratio", 0.34),
        ("A3", BondStyle.HASHED, 1, "hashed_wedge_stroke_ratio", 0.08),
        ("B1", BondStyle.HASHED, 1, "hashed_wedge_terminal_shift_ratio", 0.45),
        ("A1", BondStyle.WEDGE, 1, "label_clearance_ratio", 0.07),
        ("B3", BondStyle.PLAIN, 2, "cap_style", "flat"),
        ("B3", BondStyle.PLAIN, 2, "join_style", "bevel"),
    ),
)
def test_each_enabled_parameter_changes_production_item(case_id, bond_style, bond_order, key, value) -> None:
    graph = case_graph(case_id)
    original_snapshot = graph_snapshot(graph)
    adapter = ProductionCanvasAdapter.create(graph)
    try:
        item = _find_item(adapter, style=bond_style, order=bond_order)
        before = _signature(item)
        values = DepictionStyleParameters.current().to_mapping()
        values[key] = value
        working = DepictionStyleParameters.from_mapping(values)
        adapter.apply_parameters(working)
        after = _signature(item)

        assert after != before, f"{key} did not alter the effective BondItem output"
        assert graph_snapshot(adapter.canvas.model) == original_snapshot
        assert graph_snapshot(graph) == original_snapshot
        assert adapter.canvas.drawing_style == working.to_drawing_style()
    finally:
        adapter.close()
        from PyQt6.QtWidgets import QApplication
        QApplication.processEvents()


def test_working_update_leaves_current_panel_and_chemical_data_unchanged() -> None:
    graph = case_graph("A1")
    snapshot = graph_snapshot(graph)
    current = ProductionCanvasAdapter.create(graph)
    working = ProductionCanvasAdapter.create(graph)
    try:
        current_item = _find_item(current, style=BondStyle.WEDGE)
        working_item = _find_item(working, style=BondStyle.WEDGE)
        stable_signature = _signature(current_item)
        stable_image = _image_digest(current.render_image(width=480, height=360))
        original_working_signature = _signature(working_item)

        values = DepictionStyleParameters.current().to_mapping()
        values["solid_wedge_width_ratio"] = 0.36
        working.apply_parameters(DepictionStyleParameters.from_mapping(values))

        assert _signature(current_item) == stable_signature
        assert _signature(working_item) != original_working_signature
        assert _image_digest(current.render_image(width=480, height=360)) == stable_image
        assert graph_snapshot(current.canvas.model) == snapshot
        assert graph_snapshot(working.canvas.model) == snapshot
    finally:
        current.close()
        working.close()
        from PyQt6.QtWidgets import QApplication
        QApplication.processEvents()


def test_all_fixture_rendering_preserves_atom_bond_and_stereo_snapshots() -> None:
    for case in CASES:
        graph = case_graph(case.case_id)
        expected = graph_snapshot(graph)
        adapter = ProductionCanvasAdapter.create(graph)
        try:
            adapter.apply_parameters(DepictionStyleParameters.current())
            adapter.set_view(1.5, 45.0)
            assert graph_snapshot(graph) == expected
            assert graph_snapshot(adapter.canvas.model) == expected
            assert adapter.canvas.scene.sceneRect() == adapter.scene_rect
        finally:
            adapter.close()
            from PyQt6.QtWidgets import QApplication
            QApplication.processEvents()


def test_diagnostic_overlay_is_scene_only_and_can_be_omitted_from_export() -> None:
    graph = case_graph("A1")
    expected = graph_snapshot(graph)
    adapter = ProductionCanvasAdapter.create(graph)
    try:
        diagnostics = GeometryDiagnostics(adapter)
        baseline = _image_digest(adapter.render_image(width=600, height=420))
        diagnostics.show()
        assert adapter.diagnostic_items
        with_overlay = _image_digest(
            adapter.render_image(width=600, height=420, include_diagnostics=True)
        )
        without_overlay = _image_digest(
            adapter.render_image(width=600, height=420, include_diagnostics=False)
        )
        assert with_overlay != baseline
        assert without_overlay == baseline
        assert graph_snapshot(adapter.canvas.model) == expected
        diagnostics.clear()
        assert not adapter.diagnostic_items
    finally:
        adapter.close()
        from PyQt6.QtWidgets import QApplication
        QApplication.processEvents()
