"""
Chemuson Toolbar
Vertical toolbar with palette-style tool buttons.
"""
from __future__ import annotations

from PyQt6.QtWidgets import (
    QToolBar,
    QToolButton,
    QMenu,
    QInputDialog,
    QWidget,
    QGridLayout,
    QWidgetAction,
)
from PyQt6.QtGui import QAction, QActionGroup
from PyQt6.QtCore import pyqtSignal, Qt, QSize

from core.model import BondStyle, BondStereo
from gui.icons import (
    draw_arrow_icon,
    draw_bond_icon,
    draw_generic_icon,
    draw_glyph_icon,
    draw_ring_icon,
    draw_charge_icon,
    draw_electron_icon,
    draw_radical_charge_icon,
    draw_wavy_anchor_icon,
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
        self.setIconSize(QSize(28, 28))

        self.setStyleSheet(TOOL_PALETTE_STYLESHEET)

        self.action_group = QActionGroup(self)
        self.action_group.setExclusive(True)

        self._selection_meta = {
            "tool_select": (draw_generic_icon("pointer"), "Seleccionar"),
            "tool_select_lasso": (draw_generic_icon("lasso"), "Seleccion libre"),
        }
        self._bracket_meta = {
            "tool_brackets_round": (draw_glyph_icon("()"), "Parentesis ()"),
            "tool_brackets_square": (draw_glyph_icon("[]"), "Corchetes []"),
            "tool_brackets_curly": (draw_glyph_icon("{}"), "Llaves {}"),
        }
        
        self._arrow_meta = {
            "tool_arrow_forward": (draw_arrow_icon("forward"), "Flecha directa"),
            "tool_arrow_forward_open": (draw_arrow_icon("forward_open"), "Flecha directa abierta"),
            "tool_arrow_forward_dashed": (
                draw_arrow_icon("forward_dashed"),
                "Flecha directa discontinua",
            ),
            "tool_arrow_retro": (draw_arrow_icon("retro"), "Flecha retro"),
            "tool_arrow_retro_open": (draw_arrow_icon("retro_open"), "Flecha retro abierta"),
            "tool_arrow_retro_dashed": (
                draw_arrow_icon("retro_dashed"),
                "Flecha retro discontinua",
            ),
            "tool_arrow_both": (draw_arrow_icon("both"), "Flecha doble"),
            "tool_arrow_both_open": (draw_arrow_icon("both_open"), "Flecha doble abierta"),
            "tool_arrow_both_dashed": (
                draw_arrow_icon("both_dashed"),
                "Flecha doble discontinua",
            ),
            "tool_arrow_equilibrium": (draw_arrow_icon("equilibrium"), "Equilibrio"),
            "tool_arrow_equilibrium_dashed": (
                draw_arrow_icon("equilibrium_dashed"),
                "Equilibrio discontinuo",
            ),
            "tool_arrow_retrosynthetic": (
                draw_arrow_icon("retrosynthetic"),
                "Flecha retrosintesis",
            ),
            "tool_arrow_curved": (draw_arrow_icon("curved"), "Flecha curva"),
            "tool_arrow_curved_fishhook": (
                draw_arrow_icon("curved_fishhook"),
                "Flecha curva (1 e-)",
            ),
        }

        self._current_select_tool_id = "tool_select"
        self._current_bracket_tool_id = "tool_brackets_square"
        self._current_arrow_tool_id = "tool_arrow_forward"

        select_icon, select_tip = self._selection_meta[self._current_select_tool_id]
        self.select_button, self.select_action = self._add_palette_button(
            select_icon,
            select_tip,
            "tool_select",
            trigger_callback=self._emit_current_selection_tool,
        )
        self._build_select_palette(self.select_button.menu())
        self._add_tool_action(draw_generic_icon("eraser"), "Borrar", "tool_erase")
        self.addSeparator()

        self.bond_button, self.bond_action = self._add_palette_button(
            draw_bond_icon("single"), "Enlaces", "tool_bond"
        )
        self._build_bond_palette(self.bond_button.menu())
        self._add_tool_action(draw_generic_icon("chain"), "Cadena", "tool_chain")
        self.addSeparator()

        self.ring_button, self.ring_action = self._add_palette_button(
            draw_ring_icon(6, aromatic=True), "Anillos", "tool_ring"
        )
        self._build_ring_palette(self.ring_button.menu())
        self.addSeparator()

        self.label_button, self.label_action = self._add_palette_button(
            draw_glyph_icon("C"), "Labels", "tool_atom"
        )
        self._build_label_palette(self.label_button.menu())
        self.addSeparator()
        
        # --- Text Tool ---
        self.text_button, self.text_action = self._add_palette_button(
            draw_glyph_icon("T"), "Texto", "tool_text"
        )
        
        # --- Brackets Tool ---
        bracket_icon, bracket_tip = self._bracket_meta[self._current_bracket_tool_id]
        self.bracket_button, self.bracket_action = self._add_palette_button(
            bracket_icon,
            bracket_tip,
            "tool_brackets",
            trigger_callback=self._emit_current_bracket_tool
        )
        self._build_bracket_palette(self.bracket_button.menu())
        
        self.addSeparator()

        # --- Arrows (Annotation) Tool ---
        arrow_icon, arrow_tip = self._arrow_meta[self._current_arrow_tool_id]
        self.annotation_button, self.annotation_action = self._add_palette_button(
            arrow_icon,
            arrow_tip,
            "tool_annotation",
            trigger_callback=self._emit_current_arrow_tool,
        )
        self._build_arrow_palette(self.annotation_button.menu())

        default_bond_spec = {
            "order": 1,
            "style": BondStyle.PLAIN,
            "stereo": BondStereo.NONE,
            "mode": "increment",
            "aromatic": False,
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
            draw_ring_icon(6, aromatic=True),
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

    def _add_palette_button(
        self,
        icon,
        tooltip: str,
        internal_id: str,
        trigger_callback=None,
    ):
        action = QAction(icon, "", self)
        action.setObjectName(internal_id)
        action.setToolTip(tooltip)
        action.setCheckable(True)
        self.action_group.addAction(action)
        if trigger_callback is None:
            action.triggered.connect(lambda checked, id=internal_id: self.tool_changed.emit(id))
        else:
            action.triggered.connect(trigger_callback)

        button = QToolButton(self)
        button.setDefaultAction(action)
        button.setPopupMode(QToolButton.ToolButtonPopupMode.MenuButtonPopup)
        button.setMenu(QMenu(button))
        self.addWidget(button)
        return button, action

    def _emit_current_selection_tool(self, checked: bool = False) -> None:
        self.tool_changed.emit(self._current_select_tool_id)

    def _emit_current_bracket_tool(self, checked: bool = False) -> None:
        self.tool_changed.emit(self._current_bracket_tool_id)

    def _emit_current_arrow_tool(self, checked: bool = False) -> None:
        self.tool_changed.emit(self._current_arrow_tool_id)

    def _make_palette_entry(self, icon, tooltip: str, callback, enabled: bool = True) -> dict:
        return {
            "icon": icon,
            "tooltip": tooltip,
            "callback": callback,
            "enabled": enabled,
        }

    def _trigger_palette_action(self, callback, menu: QMenu) -> None:
        callback()
        menu.close()

    def set_text_menu(self, actions: list[QAction], color_actions: list[QAction]) -> None:
        if not hasattr(self, "text_button"):
            return
        menu = self.text_button.menu()
        menu.clear()
        for action in actions:
            if action is None:
                menu.addSeparator()
            else:
                menu.addAction(action)
        if color_actions:
            color_menu = menu.addMenu("Color de etiquetas")
            for action in color_actions:
                color_menu.addAction(action)

    def _populate_grid_menu(self, menu: QMenu, entries: list[dict], columns: int) -> None:
        container = QWidget(menu)
        container.setObjectName("palette_grid")
        layout = QGridLayout(container)
        layout.setContentsMargins(6, 6, 6, 6)
        layout.setSpacing(4)
        icon_size = QSize(22, 22)

        for index, entry in enumerate(entries):
            button = QToolButton(container)
            button.setIcon(entry["icon"])
            button.setIconSize(icon_size)
            button.setToolTip(entry["tooltip"])
            button.setEnabled(entry.get("enabled", True))
            button.setAutoRaise(True)
            callback = entry.get("callback")
            if callback is not None:
                button.clicked.connect(
                    lambda checked=False, cb=callback, m=menu: self._trigger_palette_action(cb, m)
                )
            layout.addWidget(button, index // columns, index % columns)

        action = QWidgetAction(menu)
        action.setDefaultWidget(container)
        menu.addAction(action)

    def _build_select_palette(self, menu: QMenu) -> None:
        entries = []
        for tool_id, (icon, tooltip) in self._selection_meta.items():
            entries.append(
                self._make_palette_entry(
                    icon,
                    tooltip,
                    lambda tid=tool_id, ic=icon, tip=tooltip: self._select_selection_palette(
                        tid, ic, tip
                    ),
                )
            )
        self._populate_grid_menu(menu, entries, columns=2)

    def _build_bond_palette(self, menu: QMenu) -> None:
        icon_single = draw_bond_icon("single")
        icon_double = draw_bond_icon("double")
        icon_triple = draw_bond_icon("triple")
        icon_aromatic = draw_bond_icon("aromatic")
        icon_interaction = draw_bond_icon("interaction")
        icon_wedge = draw_bond_icon("wedge")
        icon_hashed = draw_bond_icon("hashed")
        icon_wavy = draw_bond_icon("wavy")

        entries = [
            self._make_palette_entry(
                icon_single,
                "Enlace sencillo (incremental)",
                lambda ic=icon_single: self._select_bond_palette(
                    self.bond_button,
                    ic,
                    "Enlace sencillo (incremental)",
                    {
                        "order": 1,
                        "style": BondStyle.PLAIN,
                        "stereo": BondStereo.NONE,
                        "mode": "increment",
                    },
                ),
            ),
            self._make_palette_entry(
                icon_double,
                "Enlace doble",
                lambda ic=icon_double: self._select_bond_palette(
                    self.bond_button,
                    ic,
                    "Enlace doble",
                    {"order": 2, "style": BondStyle.PLAIN, "stereo": BondStereo.NONE, "mode": "set"},
                ),
            ),
            self._make_palette_entry(
                icon_triple,
                "Enlace triple",
                lambda ic=icon_triple: self._select_bond_palette(
                    self.bond_button,
                    ic,
                    "Enlace triple",
                    {"order": 3, "style": BondStyle.PLAIN, "stereo": BondStereo.NONE, "mode": "set"},
                ),
            ),
            self._make_palette_entry(
                icon_aromatic,
                "Enlace aromatico",
                lambda ic=icon_aromatic: self._select_bond_palette(
                    self.bond_button,
                    ic,
                    "Enlace aromatico",
                    {
                        "order": 1,
                        "style": BondStyle.PLAIN,
                        "stereo": BondStereo.NONE,
                        "mode": "set",
                        "aromatic": True,
                    },
                ),
            ),
            self._make_palette_entry(
                icon_wedge,
                "Wedge",
                lambda ic=icon_wedge: self._select_bond_palette(
                    self.bond_button,
                    ic,
                    "Wedge",
                    {"order": 1, "style": BondStyle.WEDGE, "stereo": BondStereo.UP, "mode": "set"},
                ),
            ),
            self._make_palette_entry(
                icon_hashed,
                "Wedge hashed",
                lambda ic=icon_hashed: self._select_bond_palette(
                    self.bond_button,
                    ic,
                    "Wedge hashed",
                    {"order": 1, "style": BondStyle.HASHED, "stereo": BondStereo.DOWN, "mode": "set"},
                ),
            ),
            self._make_palette_entry(
                icon_wavy,
                "Wavy",
                lambda ic=icon_wavy: self._select_bond_palette(
                    self.bond_button,
                    ic,
                    "Wavy",
                    {"order": 1, "style": BondStyle.WAVY, "stereo": BondStereo.EITHER, "mode": "set"},
                ),
            ),
            self._make_palette_entry(
                icon_interaction,
                "Enlace intermolecular",
                lambda ic=icon_interaction: self._select_bond_palette(
                    self.bond_button,
                    ic,
                    "Enlace intermolecular",
                    {
                        "order": 1,
                        "style": BondStyle.INTERACTION,
                        "stereo": BondStereo.NONE,
                        "mode": "set",
                    },
                ),
            ),
        ]
        self._populate_grid_menu(menu, entries, columns=3)

    def _build_ring_palette(self, menu: QMenu) -> None:
        entries = []
        icon_benzene = draw_ring_icon(6, aromatic=True)
        entries.append(
            self._make_palette_entry(
                icon_benzene,
                "Benceno",
                lambda ic=icon_benzene: self._select_ring_palette(
                    self.ring_button, ic, "Benceno", {"size": 6, "aromatic": True}
                ),
            )
        )
        for size in range(3, 13):
            icon = draw_ring_icon(size, aromatic=False)
            entries.append(
                self._make_palette_entry(
                    icon,
                    f"Anillo {size}",
                    lambda s=size, ic=icon: self._select_ring_palette(
                        self.ring_button, ic, f"Anillo {s}", {"size": s, "aromatic": False}
                    ),
                )
            )
        self._populate_grid_menu(menu, entries, columns=4)
        menu.addSeparator()
        custom = QAction("Tamaño personalizado...", self)
        custom.triggered.connect(self._select_custom_ring_size)
        menu.addAction(custom)

    def _build_label_palette(self, menu: QMenu) -> None:
        elements = ["C", "N", "O", "S", "P", "F", "Cl", "Br", "I", "H"]
        entries = []
        for element in elements:
            icon = draw_glyph_icon(element)
            entries.append(
                self._make_palette_entry(
                    icon,
                    element,
                    lambda el=element: self._select_element_palette(el),
                )
            )
        self._populate_grid_menu(menu, entries, columns=5)
        menu.addSeparator()
        periodic_action = QAction("Tabla periódica...", self)
        periodic_action.triggered.connect(self.periodic_table_requested.emit)
        menu.addAction(periodic_action)

    def _build_bracket_palette(self, menu: QMenu) -> None:
        entries = []
        for tool_id in sorted(self._bracket_meta.keys()):
            icon, tooltip = self._bracket_meta[tool_id]
            entries.append(
                self._make_palette_entry(
                    icon,
                    tooltip,
                    lambda tid=tool_id: self._select_bracket_tool(tid),
                )
            )
        self._populate_grid_menu(menu, entries, columns=3)

    def _build_arrow_palette(self, menu: QMenu) -> None:
        entries = []
        # Defined order for arrow tools
        tool_order = [
            "tool_arrow_forward",
            "tool_arrow_forward_open",
            "tool_arrow_forward_dashed",
            "tool_arrow_retro",
            "tool_arrow_retro_open",
            "tool_arrow_retro_dashed",
            "tool_arrow_both",
            "tool_arrow_both_open",
            "tool_arrow_both_dashed",
            "tool_arrow_equilibrium",
            "tool_arrow_equilibrium_dashed",
            "tool_arrow_retrosynthetic",
            "tool_arrow_curved",
            "tool_arrow_curved_fishhook",
        ]
        
        for tool_id in tool_order:
            icon, tooltip = self._arrow_meta[tool_id]
            entries.append(
                self._make_palette_entry(
                    icon,
                    tooltip,
                    lambda tid=tool_id: self._select_arrow_tool(tid),
                )
            )
        self._populate_grid_menu(menu, entries, columns=4)

    def _select_selection_palette(self, tool_id: str, icon, tooltip: str) -> None:
        self._current_select_tool_id = tool_id
        self.select_action.setIcon(icon)
        self.select_action.setToolTip(tooltip)
        self.select_button.setToolTip(tooltip)
        self.tool_changed.emit(tool_id)
        self.select_action.setChecked(True)

    def _select_bracket_tool(self, tool_id: str) -> None:
        icon, tooltip = self._bracket_meta[tool_id]
        self._current_bracket_tool_id = tool_id
        self.bracket_action.setIcon(icon)
        self.bracket_action.setToolTip(tooltip)
        self.bracket_button.setToolTip(tooltip)
        self.bracket_action.setChecked(True)
        self.tool_changed.emit(tool_id)

    def _select_arrow_tool(self, tool_id: str) -> None:
        icon, tooltip = self._arrow_meta[tool_id]
        self._current_arrow_tool_id = tool_id
        self.annotation_action.setIcon(icon)
        self.annotation_action.setToolTip(tooltip)
        self.annotation_button.setToolTip(tooltip)
        self.annotation_action.setChecked(True)
        self.tool_changed.emit(tool_id)

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
        icon = draw_glyph_icon(element)
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
            draw_ring_icon(size, aromatic=False),
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


class SymbolPaletteToolbar(QToolBar):
    """
    Right-side toolbar for chemical symbol tools (charges, radicals, electrons, etc.).
    """

    tool_changed = pyqtSignal(str)

    def __init__(self, action_group: QActionGroup, parent=None) -> None:
        super().__init__("Simbolismos químicos", parent)
        self.setOrientation(Qt.Orientation.Vertical)
        self.setMovable(False)
        self.setFloatable(False)
        self.setIconSize(QSize(28, 28))
        self.setStyleSheet(TOOL_PALETTE_STYLESHEET)

        self._action_group = action_group
        self._current_symbol_tool_id = "tool_charge_plus"

        icon, tip = self._symbol_meta()[self._current_symbol_tool_id]
        self.symbol_button, self.symbol_action = self._add_palette_button(
            icon,
            tip,
            "tool_symbol_palette",
            trigger_callback=self._emit_current_symbol_tool,
        )
        self._build_symbol_palette(self.symbol_button.menu())

    def _add_palette_button(
        self,
        icon,
        tooltip: str,
        internal_id: str,
        trigger_callback=None,
    ):
        action = QAction(icon, "", self)
        action.setObjectName(internal_id)
        action.setToolTip(tooltip)
        action.setCheckable(True)
        if self._action_group is not None:
            self._action_group.addAction(action)
        if trigger_callback is None:
            action.triggered.connect(lambda checked, id=internal_id: self.tool_changed.emit(id))
        else:
            action.triggered.connect(trigger_callback)

        button = QToolButton(self)
        button.setDefaultAction(action)
        button.setPopupMode(QToolButton.ToolButtonPopupMode.MenuButtonPopup)
        button.setMenu(QMenu(button))
        self.addWidget(button)
        return button, action

    def _symbol_meta(self) -> dict[str, tuple]:
        return {
            "tool_charge_plus": (draw_charge_icon("+"), "Carga positiva"),
            "tool_charge_minus": (draw_charge_icon("-"), "Carga negativa"),
            "tool_charge": (draw_glyph_icon("±"), "Carga alterna (+/-)"),
            "tool_symbol_plus": (draw_glyph_icon("+"), "Signo más"),
            "tool_symbol_minus": (draw_glyph_icon("-"), "Signo menos"),
            "tool_symbol_radical": (draw_electron_icon(1), "Electrón desapareado"),
            "tool_symbol_lone_pair": (draw_electron_icon(2, spread=5.0), "Par solitario"),
            "tool_symbol_wavy_anchor": (draw_wavy_anchor_icon(), "Ancla ondulada"),
            "tool_symbol_radical_cation": (draw_radical_charge_icon("+"), "Radical catión"),
            "tool_symbol_radical_anion": (draw_radical_charge_icon("-"), "Radical anión"),
            "tool_symbol_partial_plus": (draw_glyph_icon("δ+"), "Carga parcial (δ+)"),
            "tool_symbol_partial_minus": (draw_glyph_icon("δ-"), "Carga parcial (δ-)"),
        }

    def _make_palette_entry(self, icon, tooltip: str, callback, enabled: bool = True) -> dict:
        return {
            "icon": icon,
            "tooltip": tooltip,
            "callback": callback,
            "enabled": enabled,
        }

    def _populate_grid_menu(self, menu: QMenu, entries: list[dict], columns: int) -> None:
        container = QWidget(menu)
        container.setObjectName("palette_grid")
        layout = QGridLayout(container)
        layout.setContentsMargins(6, 6, 6, 6)
        layout.setSpacing(4)
        icon_size = QSize(22, 22)

        for index, entry in enumerate(entries):
            button = QToolButton(container)
            button.setIcon(entry["icon"])
            button.setIconSize(icon_size)
            button.setToolTip(entry["tooltip"])
            button.setEnabled(entry.get("enabled", True))
            button.setAutoRaise(True)
            callback = entry.get("callback")
            if callback is not None:
                button.clicked.connect(lambda checked=False, cb=callback, m=menu: self._trigger_palette_action(cb, m))
            layout.addWidget(button, index // columns, index % columns)

        action = QWidgetAction(menu)
        action.setDefaultWidget(container)
        menu.addAction(action)

    def _trigger_palette_action(self, callback, menu: QMenu) -> None:
        callback()
        menu.close()

    def _build_symbol_palette(self, menu: QMenu) -> None:
        entries = []
        for tool_id, (icon, tooltip) in self._symbol_meta().items():
            entries.append(
                self._make_palette_entry(
                    icon,
                    tooltip,
                    lambda tid=tool_id: self._select_symbol_tool(tid),
                )
            )
        self._populate_grid_menu(menu, entries, columns=4)

    def _emit_current_symbol_tool(self, checked: bool = False) -> None:
        self.tool_changed.emit(self._current_symbol_tool_id)

    def _select_symbol_tool(self, tool_id: str) -> None:
        icon, tooltip = self._symbol_meta()[tool_id]
        self._current_symbol_tool_id = tool_id
        self.symbol_action.setIcon(icon)
        self.symbol_action.setToolTip(tooltip)
        self.symbol_button.setToolTip(tooltip)
        self.symbol_action.setChecked(True)
        self.tool_changed.emit(tool_id)
