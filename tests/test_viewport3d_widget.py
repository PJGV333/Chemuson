from __future__ import annotations

import pytest
from PyQt6.QtWidgets import QApplication

from chemuson.geometry3d import SceneAtom3D, SceneMolecule3D
from chemuson.gui.main_window import ChemusonWindow
from chemuson.gui.viewport3d import MolecularViewportDock, MolecularViewportWidget


@pytest.fixture(scope="module", autouse=True)
def _qapp():
    return QApplication.instance() or QApplication([])


def test_viewport_widget_keeps_selection_local() -> None:
    widget = MolecularViewportWidget()
    scene = SceneMolecule3D(
        atoms=(SceneAtom3D(1, "C", 0.0, 0.0, 0.0, 1.7, "#404040"),),
        bonds=(),
    )
    widget.set_scene(scene)
    assert widget.scene is scene
    assert widget.selected_atom_id is None
    assert widget.selected_bond_id is None
    assert not widget.gl_ready


def test_viewport_dock_is_hidden_until_explicitly_opened() -> None:
    dock = MolecularViewportDock()
    assert not dock.isVisible()
    assert dock.viewport.scene.atoms == ()
    dock.close()


def test_main_window_exposes_one_viewport_action_without_side_panel_duplication() -> None:
    window = ChemusonWindow()
    try:
        assert not window.viewport3d_dock.isVisible()
        assert window.action_toggle_viewport3d.text() == "Visor molecular 3D"
        assert window.action_toggle_viewport3d in window.view_menu.actions()
        assert window.side_panel.PAGE_KEYS == (
            "inspector",
            "validation",
            "properties",
            "templates",
            "appearance",
            "spectroscopy",
            "compchem",
        )
        window.action_toggle_viewport3d.trigger()
        assert not window.viewport3d_dock.isHidden()
        window.action_toggle_viewport3d.trigger()
        assert window.viewport3d_dock.isHidden()
    finally:
        window.close()
