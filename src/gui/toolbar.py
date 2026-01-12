"""
Chemuson Toolbar
Vertical toolbar with palette-style tool buttons.
"""
from __future__ import annotations

from PyQt6.QtWidgets import QToolBar, QToolButton, QMenu, QInputDialog
from PyQt6.QtGui import QAction, QActionGroup
from PyQt6.QtCore import pyqtSignal, Qt, QSize

from core.model import BondStyle, BondStereo
from gui.icons import (
    draw_atom_icon,
    draw_bond_icon,
    draw_generic_icon,
    draw_ring_icon,
)


from gui.styles import TOOL_PALETTE_STYLESHEET

class ChemusonToolbar(QToolBar):
    """
    Vertical toolbar for selecting drawing tools.
    Organized into tools + palettes (bonds, rings, labels).
    """

    tool_changed = pyqtSignal(str)
    bond_palette_changed = pyqtSignal(object)
    ring_palette_changed = pyqtSignal(object)
    element_palette_changed = pyqtSignal(str)
    periodic_table_requested = pyqtSignal()

    def __init__(self, parent=None) -> None:
        super().__init__("Herramientas de Dibujo", parent)
        self.setOrientation(Qt.Orientation.Vertical)
        self.setMovable(False)
        self.setFloatable(False)
        self.setIconSize(QSize(32, 32))

        self.setStyleSheet(TOOL_PALETTE_STYLESHEET)

        self.action_group = QActionGroup(self)
        self.action_group.setExclusive(True)

        self._add_tool_action(draw_generic_icon("pointer"), "Seleccionar", "tool_select")
        self._add_tool_action(draw_generic_icon("eraser"), "Borrar", "tool_erase")
        self.addSeparator()

        self.bond_button, self.bond_action = self._add_palette_button(
            draw_bond_icon("single"), "Enlaces", "tool_bond"
        )
        self._build_bond_palette(self.bond_button.menu())

        self.ring_button, self.ring_action = self._add_palette_button(
            draw_ring_icon(), "Anillos", "tool_ring"
        )
        self._build_ring_palette(self.ring_button.menu())

        self.label_button, self.label_action = self._add_palette_button(
            draw_atom_icon("C"), "Labels", "tool_atom"
        )
        self._build_label_palette(self.label_button.menu())

        default_bond_spec = {
            "order": 1,
            "style": BondStyle.PLAIN,
            "stereo": BondStereo.NONE,
            "mode": "increment",
        }
        default_ring_spec = {"size": 6, "aromatic": True}
        default_element = "C"

        self._current_bond_spec = dict(default_bond_spec)
        self._current_ring_spec = dict(default_ring_spec)
        self._current_element = default_element

        self._select_bond_palette(
            self.bond_button,
            draw_bond_icon("single"),
            "Enlace sencillo (incremental)",
            default_bond_spec,
        )
        self._select_ring_palette(
            self.ring_button,
            draw_ring_icon(),
            "Benceno",
            default_ring_spec,
        )
        self._select_element_palette(default_element)

    def _add_tool_action(self, icon, tooltip: str, internal_id: str) -> QAction:
        action = QAction(icon, "", self)
        action.setObjectName(internal_id)
        action.setToolTip(tooltip)
        action.setCheckable(True)
        self.action_group.addAction(action)
        self.addAction(action)
        action.triggered.connect(lambda checked, id=internal_id: self.tool_changed.emit(id))
        return action

    def _add_palette_button(self, icon, tooltip: str, internal_id: str):
        action = QAction(icon, "", self)
        action.setObjectName(internal_id)
        action.setToolTip(tooltip)
        action.setCheckable(True)
        self.action_group.addAction(action)
        action.triggered.connect(lambda checked, id=internal_id: self.tool_changed.emit(id))

        button = QToolButton(self)
        button.setDefaultAction(action)
        button.setPopupMode(QToolButton.ToolButtonPopupMode.MenuButtonPopup)
        button.setMenu(QMenu(button))
        self.addWidget(button)
        return button, action

    def _build_bond_palette(self, menu: QMenu) -> None:
        self._add_bond_action(
            menu,
            draw_bond_icon("single"),
            "Enlace sencillo (incremental)",
            {"order": 1, "style": BondStyle.PLAIN, "stereo": BondStereo.NONE, "mode": "increment"},
        )
        self._add_bond_action(
            menu,
            draw_bond_icon("single"),
            "Enlace sencillo",
            {"order": 1, "style": BondStyle.PLAIN, "stereo": BondStereo.NONE, "mode": "set"},
        )
        self._add_bond_action(
            menu,
            draw_bond_icon("double"),
            "Enlace doble",
            {"order": 2, "style": BondStyle.PLAIN, "stereo": BondStereo.NONE, "mode": "set"},
        )
        self._add_bond_action(
            menu,
            draw_bond_icon("triple"),
            "Enlace triple",
            {"order": 3, "style": BondStyle.PLAIN, "stereo": BondStereo.NONE, "mode": "set"},
        )
        menu.addSeparator()
        self._add_bond_action(
            menu,
            draw_bond_icon("wedge"),
            "Wedge",
            {"order": 1, "style": BondStyle.WEDGE, "stereo": BondStereo.UP, "mode": "set"},
        )
        self._add_bond_action(
            menu,
            draw_bond_icon("hashed"),
            "Wedge hashed",
            {"order": 1, "style": BondStyle.HASHED, "stereo": BondStereo.DOWN, "mode": "set"},
        )
        self._add_bond_action(
            menu,
            draw_bond_icon("wavy"),
            "Wavy",
            {"order": 1, "style": BondStyle.WAVY, "stereo": BondStereo.EITHER, "mode": "set"},
        )

    def _add_bond_action(self, menu: QMenu, icon, text: str, spec: dict) -> None:
        action = QAction(icon, text, self)
        action.triggered.connect(
            lambda checked=False, a=action, s=spec: self._select_bond_palette(
                self.bond_button, a.icon(), a.text(), s
            )
        )
        menu.addAction(action)

    def _build_ring_palette(self, menu: QMenu) -> None:
        for size in range(3, 9):
            self._add_ring_action(
                menu,
                draw_ring_icon(),
                f"Anillo {size}",
                {"size": size, "aromatic": False},
            )
        menu.addSeparator()
        self._add_ring_action(
            menu,
            draw_ring_icon(),
            "Benceno",
            {"size": 6, "aromatic": True},
        )
        menu.addSeparator()
        custom = QAction("Tamaño personalizado...", self)
        custom.triggered.connect(self._select_custom_ring_size)
        menu.addAction(custom)

    def _add_ring_action(self, menu: QMenu, icon, text: str, spec: dict) -> None:
        action = QAction(icon, text, self)
        action.triggered.connect(
            lambda checked=False, a=action, s=spec: self._select_ring_palette(
                self.ring_button, a.icon(), a.text(), s
            )
        )
        menu.addAction(action)

    def _build_label_palette(self, menu: QMenu) -> None:
        elements = ["C", "N", "O", "S", "P", "F", "Cl", "Br", "I", "H"]
        for element in elements:
            icon = draw_atom_icon(element)
            action = QAction(icon, element, self)
            action.triggered.connect(
                lambda checked=False, el=element: self._select_element_palette(el)
            )
            menu.addAction(action)
        menu.addSeparator()
        periodic_action = QAction("Tabla periódica...", self)
        periodic_action.triggered.connect(self.periodic_table_requested.emit)
        menu.addAction(periodic_action)

    def _select_bond_palette(self, button: QToolButton, icon, text: str, spec: dict) -> None:
        self.bond_action.setIcon(icon)
        self.bond_action.setToolTip(text)
        button.setToolTip(text)
        self._current_bond_spec = dict(spec)
        self.bond_palette_changed.emit(spec)
        self.tool_changed.emit("tool_bond")
        self.bond_action.setChecked(True)

    def _select_ring_palette(self, button: QToolButton, icon, text: str, spec: dict) -> None:
        self.ring_action.setIcon(icon)
        self.ring_action.setToolTip(text)
        button.setToolTip(text)
        self._current_ring_spec = dict(spec)
        self.ring_palette_changed.emit(spec)
        self.tool_changed.emit("tool_ring")
        self.ring_action.setChecked(True)

    def _select_element_palette(self, element: str) -> None:
        icon = draw_atom_icon(element)
        self.label_action.setIcon(icon)
        self.label_action.setToolTip(f"Elemento {element}")
        self.label_button.setToolTip(f"Elemento {element}")
        self._current_element = element
        self.element_palette_changed.emit(element)
        self.tool_changed.emit("tool_atom")
        self.label_action.setChecked(True)

    def _select_custom_ring_size(self) -> None:
        size, ok = QInputDialog.getInt(self, "Tamaño de anillo", "Número de miembros:", 6, 3, 12)
        if not ok:
            return
        self._select_ring_palette(
            self.ring_button,
            draw_ring_icon(),
            f"Anillo {size}",
            {"size": size, "aromatic": False},
        )

    def select_element(self, element: str) -> None:
        self._select_element_palette(element)

    def current_bond_spec(self) -> dict:
        return dict(self._current_bond_spec)

    def current_ring_spec(self) -> dict:
        return dict(self._current_ring_spec)

    def current_element(self) -> str:
        return self._current_element
