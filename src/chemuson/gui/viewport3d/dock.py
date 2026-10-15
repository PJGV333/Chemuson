"""Dock discreto para el visor molecular 3D."""

from __future__ import annotations

from PyQt6.QtCore import Qt
from PyQt6.QtWidgets import QDockWidget, QLabel, QVBoxLayout, QWidget

from chemuson.geometry3d import SceneMolecule3D

from .widget import MolecularViewportWidget


class MolecularViewportDock(QDockWidget):
    """Contenedor integrable que no expone el modelo químico al renderer."""

    def __init__(self, parent=None) -> None:
        super().__init__("Visor molecular 3D", parent)
        self.setObjectName("molecularViewportDock")
        self.setAllowedAreas(Qt.DockWidgetArea.LeftDockWidgetArea | Qt.DockWidgetArea.RightDockWidgetArea)
        self.setMinimumWidth(360)
        container = QWidget(self)
        layout = QVBoxLayout(container)
        layout.setContentsMargins(6, 6, 6, 6)
        layout.setSpacing(4)
        self.viewport = MolecularViewportWidget(container)
        self.status_label = QLabel("Sin conformero 3D")
        self.status_label.setObjectName("molecularViewportStatus")
        self.status_label.setWordWrap(True)
        layout.addWidget(self.viewport, 1)
        layout.addWidget(self.status_label)
        self.setWidget(container)
        self.viewport.statusChanged.connect(self.set_status)

    def set_scene(self, scene: SceneMolecule3D | None) -> None:
        self.viewport.set_scene(scene)
        if scene is None or not scene.atoms:
            self.set_status("Sin conformero 3D")
        else:
            self.set_status(f"{len(scene.atoms)} átomos · {len(scene.bonds)} enlaces")

    def clear_scene(self) -> None:
        self.set_scene(None)

    def set_status(self, message: str) -> None:
        self.status_label.setText(str(message))
