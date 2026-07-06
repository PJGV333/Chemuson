from __future__ import annotations

import pytest
from PyQt6.QtCore import Qt
from PyQt6.QtWidgets import QApplication

from chemuson.gui.docks import PlantillasDock


@pytest.fixture(scope="module", autouse=True)
def _qapp():
    return QApplication.instance() or QApplication([])


def test_plantillas_dock_renders_empty_placeholder() -> None:
    dock = PlantillasDock()

    dock.set_templates([])

    assert dock.tree.topLevelItemCount() == 1
    placeholder = dock.tree.topLevelItem(0)
    assert placeholder.text(0) == "(Sin plantillas)"
    assert not (placeholder.flags() & Qt.ItemFlag.ItemIsSelectable)


def test_plantillas_dock_sets_category_and_template_payloads() -> None:
    dock = PlantillasDock()

    dock.set_templates(
        [
            {
                "name": "Personal",
                "templates": [{"id": "tpl_1", "name": "Uno"}],
            }
        ]
    )

    category = dock.tree.topLevelItem(0)
    template = category.child(0)
    assert category.data(0, Qt.ItemDataRole.UserRole) == {
        "kind": "category",
        "name": "Personal",
    }
    assert template.data(0, Qt.ItemDataRole.UserRole) == {
        "kind": "template",
        "id": "tpl_1",
        "name": "Uno",
        "category": "Personal",
    }


def test_plantillas_dock_emits_template_selection() -> None:
    dock = PlantillasDock()
    emitted: list[dict] = []
    dock.template_selected.connect(emitted.append)
    dock.set_templates(
        [
            {
                "name": "Personal",
                "templates": [{"id": "tpl_1", "name": "Uno"}],
            }
        ]
    )

    dock._emit_template(dock.tree.topLevelItem(0).child(0))

    assert emitted == [
        {
            "kind": "template",
            "id": "tpl_1",
            "name": "Uno",
            "category": "Personal",
        }
    ]
