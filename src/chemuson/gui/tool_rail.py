"""
Rail de herramientas unificado + flyouts (Fase 4 del plan de modernización).

El ``ToolRail`` es una **superficie visual** (QWidget vertical ~58 px) que
reemplaza visualmente las dos barras históricas de toolbars
(``ChemusonToolbar`` izquierda, ``SymbolPaletteToolbar`` derecha). **No
posee lógica de herramienta propia**: cada interacción delega 1:1 a los
``QAction``, ``QActionGroup``, ``tool_id``, señales y callbacks de los
toolbars originales:

- Botón de acción: dispara el ``QAction``/handler existente
  (``chain_action.trigger()`` etc.).
- Botón de paleta: clic izquierdo activa la herramienta *actual* (mismo
  ``QAction`` que la barra histórica); clic derecho abre el :class:`Flyout`
  construido **leyendo el ``QMenu`` original** del toolbar (mismas celdas;
  el clic en una celda ejecuta ``QToolButton.click()`` del botón del menú
  original → el callback del toolbar original se ejecuta y emite las
  señales originales).

El estado activo (highlight) del rail es **derivado** de las señales
``tool_changed`` de los toolbars (no es la fuente de verdad; el estado
vive en el canvas + los toolbars).

Los atajos de letra simple (``V, A, L, B, R, C, T, N, G, E, O``) se
implementan con :class:`ToolShortcutDispatcher` (event filter contextual),
**sin ``QShortcut``** (evita conflictos de contexto y doble conexión): solo
con modificadores nulos, sin diálogo modal activo y sin el foco en un
widget de entrada de texto.
"""
from __future__ import annotations

from dataclasses import dataclass
from typing import Callable

from PyQt6.QtCore import QObject, QEvent, QSize, Qt
from PyQt6.QtGui import QAction, QCursor, QMouseEvent
from PyQt6.QtWidgets import (
    QAbstractSpinBox,
    QApplication,
    QComboBox,
    QFrame,
    QGridLayout,
    QLabel,
    QLineEdit,
    QMenu,
    QPlainTextEdit,
    QScrollArea,
    QSizePolicy,
    QToolButton,
    QTextEdit,
    QVBoxLayout,
    QWidget,
    QWidgetAction,
)

from chemuson.gui.flyout import Flyout, FlyoutFooter, FlyoutFooterButton, FlyoutItem
from chemuson.gui.toolbar import ChemusonToolbar, SymbolPaletteToolbar

__all__ = [
    "ToolRail",
    "ToolRailButton",
    "ToolShortcutDispatcher",
]

#: Ancho fijo del rail (px) — ``METRICS["railW"]`` (ver ``theme/tokens.py``).
RAIL_WIDTH = 58
#: Tamaño del botón del rail (px) — métrica del spike aprobado (42 px).
_RAIL_BTN_SIZE = 42


@dataclass
class _CellEntry:
    """Celda leída de un ``QMenu`` original (delegación 1:1)."""

    icon: object
    label: str
    enabled: bool
    trigger: Callable[[], None]


@dataclass
class _MenuEntry:
    """Entrada de pie leída de un ``QMenu`` original.

    Args:
        label: Texto del botón de pie.
        trigger: Callback (``QAction.trigger`` o popup del submenú original).
        close_after: Cerrar el flyout tras el callback.
        tooltip: Tooltip del botón de pie.
    """

    label: str
    trigger: Callable[[], None]
    close_after: bool = True
    tooltip: str = ""


@dataclass
class _RailSpec:
    """Especificación de un botón del rail."""

    key: str
    label: str
    tooltip: str
    kbd: str | None
    trigger: Callable[[], None]
    group: str
    menu: QMenu | None = None
    columns: int = 4
    plain_as: str = "cells"  # "cells" | "footer"
    icon_action: QAction | None = None
    flyout_title: str = ""
    footer_labels: dict[str, str] | None = None


class ToolRailButton(QToolButton):
    """Botón del rail: icono 21 px + kbd-hint opcional + estado activo."""

    def __init__(self, parent: QWidget | None = None) -> None:
        super().__init__(parent)
        self.setObjectName("railBtn")
        self.setAutoRaise(True)
        self.setIconSize(QSize(21, 21))
        self.setToolButtonStyle(Qt.ToolButtonStyle.ToolButtonIconOnly)
        self.setFixedSize(_RAIL_BTN_SIZE, _RAIL_BTN_SIZE)
        self.setSizePolicy(QSizePolicy.Policy.Fixed, QSizePolicy.Policy.Fixed)
        self._kbd: QLabel | None = None
        self._open_flyout: Callable[[], None] | None = None

    def set_kbd(self, text: str | None) -> None:
        """Fija la pista de tecla (esquina inferior derecha)."""
        if self._kbd is not None:
            self._kbd.setParent(None)
            self._kbd.deleteLater()
            self._kbd = None
        if text:
            label = QLabel(text, self)
            label.setObjectName("railKbd")
            font = label.font()
            font.setPixelSize(9)
            font.setBold(True)
            label.setFont(font)
            self._kbd = label
            self._place_kbd()

    def set_open_flyout(self, callback: Callable[[], None]) -> None:
        """Fija el callback de abrir flyout (clic derecho)."""
        self._open_flyout = callback

    def mousePressEvent(self, e: QMouseEvent) -> None:  # noqa: N802
        if e.button() == Qt.MouseButton.RightButton and self._open_flyout is not None:
            self._open_flyout()
            e.accept()
            return
        super().mousePressEvent(e)

    def resizeEvent(self, e) -> None:  # noqa: N802
        super().resizeEvent(e)
        self._place_kbd()

    def _place_kbd(self) -> None:
        if self._kbd is None:
            return
        self._kbd.adjustSize()
        self._kbd.move(self.width() - self._kbd.width() - 4, self.height() - self._kbd.height() - 3)

    def set_active(self, on: bool) -> None:
        """Marca/desmarca el estado activo (property QSS ``active``)."""
        self.setProperty("active", on)
        style = self.style()
        style.unpolish(self)
        style.polish(self)


class ToolRail(QWidget):
    """Rail vertical de herramientas (superficie visual, delegación 1:1).

    Args:
        toolbar: ``ChemusonToolbar`` histórico (fuente de verdad de
            selección/enlace/anillo/átomo/acciones simples).
        symbols_toolbar: ``SymbolPaletteToolbar`` histórico (texto, corchetes,
            flechas, placas, símbolos, energía, orbitales).
        window: Ventana (para las acciones de limpiado/validación/numeración).
        parent: Widget padre opcional.
    """

    def __init__(
        self,
        toolbar: ChemusonToolbar,
        symbols_toolbar: SymbolPaletteToolbar,
        window: QWidget | None = None,
        parent: QWidget | None = None,
    ) -> None:
        super().__init__(parent)
        self.setObjectName("toolRail")
        self.setFixedWidth(RAIL_WIDTH)
        self._toolbar = toolbar
        self._symbols = symbols_toolbar
        self._window = window
        self._buttons: dict[str, ToolRailButton] = {}
        self._flyouts: dict[str, Flyout] = {}
        self._specs: list[_RailSpec] = []
        self._current_group: str | None = None
        self._last_theme: str | None = None

        # El contenido vive en un ``QScrollArea`` compacto (sin marco y sin
        # scrollbar visible): a 1440×900 todos los botones caben sin
        # scroll; en ventanas pequeñas (980×600) el rail se desplaza con la
        # rueda manteniendo botones de 42 px e iconos de 21 px (sin
        # micro-iconos).
        outer = QVBoxLayout(self)
        outer.setContentsMargins(0, 0, 0, 0)
        outer.setSpacing(0)
        scroll = QScrollArea(self)
        scroll.setObjectName("railScroll")
        scroll.setWidgetResizable(True)
        scroll.setFrameShape(QFrame.Shape.NoFrame)
        scroll.setHorizontalScrollBarPolicy(Qt.ScrollBarPolicy.ScrollBarAlwaysOff)
        scroll.setVerticalScrollBarPolicy(Qt.ScrollBarPolicy.ScrollBarAlwaysOff)
        outer.addWidget(scroll)
        self._scroll = scroll

        inner = QWidget(scroll)
        scroll.setWidget(inner)
        layout = QVBoxLayout(inner)
        layout.setContentsMargins(0, 10, 0, 10)
        layout.setSpacing(3)
        self._content_layout = layout

        t, s = toolbar, symbols_toolbar
        w = window

        def _trigger(action: QAction | None) -> Callable[[], None]:
            if action is None:
                return lambda: None
            return action.trigger

        specs: list[_RailSpec] = [
            _RailSpec(
                key="select",
                label="Seleccionar",
                tooltip="Seleccionar (V)",
                kbd="V",
                trigger=_trigger(t.select_action),
                group="select",
                menu=t.select_button.menu(),
                columns=2,
                icon_action=t.select_action,
                flyout_title="Selección",
            ),
            _RailSpec(
                key="bond",
                label="Enlaces",
                tooltip="Enlaces (B)",
                kbd="B",
                trigger=_trigger(t.bond_action),
                group="bond",
                menu=t.bond_button.menu(),
                columns=3,
                icon_action=t.bond_action,
                flyout_title="Enlaces",
            ),
            _RailSpec(
                key="chain",
                label="Cadena",
                tooltip="Cadena lineal (L)",
                kbd="L",
                trigger=_trigger(t.chain_action),
                group="chain",
                icon_action=t.chain_action,
            ),
            _RailSpec(
                key="ring",
                label="Anillos",
                tooltip="Anillos (R)",
                kbd="R",
                trigger=_trigger(t.ring_action),
                group="ring",
                menu=t.ring_button.menu(),
                columns=4,
                plain_as="footer",
                icon_action=t.ring_action,
                flyout_title="Anillos",
                footer_labels={"Tamaño personalizado...": "Tamaño personalizado…"},
            ),
            _RailSpec(
                key="atom",
                label="Átomos",
                tooltip="Átomo / tabla periódica (C)",
                kbd="C",
                trigger=_trigger(t.label_action),
                group="atom",
                menu=t.label_button.menu(),
                columns=5,
                plain_as="footer",
                icon_action=t.label_action,
                flyout_title="Átomos",
                footer_labels={"Tabla periódica...": "Tabla periódica…"},
            ),
            _RailSpec(
                key="coord",
                label="Coordinación",
                tooltip="Centro de coordinación (esfera)",
                kbd=None,
                trigger=_trigger(t.coord_action),
                group="coord",
                icon_action=t.coord_action,
            ),
            _RailSpec(
                key="text",
                label="Texto",
                tooltip="Texto (T)",
                kbd="T",
                trigger=_trigger(s.text_action),
                group="text",
                menu=s.text_button.menu(),
                columns=4,
                icon_action=s.text_action,
                flyout_title="Texto",
            ),
            _RailSpec(
                key="arrows",
                label="Flechas",
                tooltip="Flechas de anotación (N)",
                kbd="N",
                trigger=_trigger(s.annotation_action),
                group="arrows",
                menu=s.annotation_button.menu(),
                columns=4,
                icon_action=s.annotation_action,
                flyout_title="Flechas",
            ),
            _RailSpec(
                key="brackets",
                label="Corchetes",
                tooltip="Corchetes / paréntesis (G)",
                kbd="G",
                trigger=_trigger(s.bracket_action),
                group="brackets",
                menu=s.bracket_button.menu(),
                columns=2,
                icon_action=s.bracket_action,
                flyout_title="Corchetes",
            ),
            _RailSpec(
                key="symbols",
                label="Símbolos",
                tooltip="Símbolos, cargas y pares de electrones",
                kbd=None,
                trigger=_trigger(s.symbol_action),
                group="symbols",
                menu=s.symbol_button.menu(),
                columns=4,
                icon_action=s.symbol_action,
                flyout_title="Símbolos",
            ),
            _RailSpec(
                key="plates",
                label="Placas",
                tooltip="Placas (TLC / gel de electroforesis)",
                kbd=None,
                trigger=_trigger(s.plate_action),
                group="plates",
                menu=s.plate_button.menu(),
                columns=2,
                icon_action=s.plate_action,
                flyout_title="Placas",
            ),
            _RailSpec(
                key="energy",
                label="Energía",
                tooltip="Diagramas de energía (E)",
                kbd="E",
                trigger=_trigger(s.energy_diagram_action),
                group="energy",
                menu=s.energy_diagram_button.menu(),
                columns=2,
                icon_action=s.energy_diagram_action,
                flyout_title="Diagramas de energía",
            ),
            _RailSpec(
                key="orbitals",
                label="Orbitales",
                tooltip="Orbitales (O)",
                kbd="O",
                trigger=_trigger(s.orbital_action),
                group="orbitals",
                menu=s.orbital_button.menu(),
                columns=4,
                icon_action=s.orbital_action,
                flyout_title="Orbitales",
            ),
        ]
        # Clean2D / Validar / Numerar no son botones permanentes del rail
        # (convergencia visual con el spike): sus ``QAction`` siguen
        # disponibles en menús, atajos (Ctrl+K) y flyouts de contexto.
        self._specs = specs

        separators_after = {"symbols", "plates"}
        for spec in self._specs:
            if spec.key in separators_after:
                layout.addSpacing(2)
                sep = QFrame(inner)
                sep.setObjectName("railSep")
                sep.setFixedHeight(1)
                sep.setFixedWidth(26)
                layout.addWidget(sep, 0, Qt.AlignmentFlag.AlignHCenter)
                layout.addSpacing(2)
            self._add_button(spec)
        layout.addStretch(1)
        self._build_all_flyouts()

    # ------------------------------------------------------------------
    # Construcción
    # ------------------------------------------------------------------
    def _add_button(self, spec: _RailSpec) -> ToolRailButton:
        button = ToolRailButton(self)
        # Sin badges kbd visibles: el atajo se mantiene en el tooltip
        # (p. ej. "Enlaces (B)") y sigue funcional vía
        # ``ToolShortcutDispatcher``.
        button.setToolTip(spec.tooltip)
        if spec.icon_action is not None:
            button.setIcon(spec.icon_action.icon())
        button.clicked.connect(lambda checked=False, s=spec: self._on_button_clicked(s))
        if spec.menu is not None:
            button.set_open_flyout(lambda k=spec.key: self.open_flyout(k))
        self._buttons[spec.key] = button
        # Centrado horizontal (métrica del spike: botones de 42 px centrados
        # en el rail de 58 px).
        self._layout().addWidget(button, 0, Qt.AlignmentFlag.AlignHCenter)
        return button

    def _on_button_clicked(self, spec: _RailSpec) -> None:
        """Clic izquierdo: primer clic activa la herramienta actual de la
        categoría; un segundo clic sobre la misma categoría *activa* abre
        su flyout. El clic derecho abre el flyout directamente
        (``ToolRailButton.mousePressEvent``).
        """
        if spec.group == self._current_group and spec.menu is not None:
            self.open_flyout(spec.key)
            return
        spec.trigger()

    def _layout(self) -> QVBoxLayout:
        """Layout de contenido del rail (dentro del scroll compacto)."""
        return self._content_layout

    @staticmethod
    def _introspect_menu(
        menu: QMenu, plain_as: str
    ) -> tuple[list[_CellEntry], list[_MenuEntry]]:
        """Lee un ``QMenu`` original (celdas de grid + entradas de pie).

        Las celdas de grid reutilizan el ``QToolButton`` del menú original
        (``click()`` → el callback del toolbar original). Las ``QAction``
        simples y submenús se convierten en botones de pie según
        ``plain_as`` (los submenús se re-abren vía ``popup()``: mismo menú,
        mismas ``QAction``).
        """
        cells: list[_CellEntry] = []
        menus: list[_MenuEntry] = []
        for action in menu.actions():
            if action.isSeparator():
                continue
            if action.menu() is not None:
                submenu = action.menu()
                menus.append(
                    _MenuEntry(
                        label="…",
                        trigger=lambda m=submenu: m.popup(QCursor.pos()),
                        close_after=False,
                        tooltip=action.text() or submenu.title() or "Menú",
                    )
                )
                continue
            widget_action = _as_widget_action(action)
            if widget_action is not None:
                container = widget_action.defaultWidget()
                lay = container.layout() if container is not None else None
                if isinstance(lay, QGridLayout):
                    seen: set = set()
                    for row in range(lay.rowCount()):
                        for col in range(lay.columnCount()):
                            item = lay.itemAtPosition(row, col)
                            if item is None:
                                continue
                            w = item.widget()
                            if (
                                isinstance(w, QToolButton)
                                and w not in seen
                                # Slots de alineación del grid original (sin
                                # etiqueta/icono): no representan una función.
                                and not (not w.isEnabled() and not (w.toolTip() or w.text()).strip())
                            ):
                                seen.add(w)
                                cells.append(
                                    _CellEntry(
                                        icon=w.icon(),
                                        label=w.toolTip() or w.text() or "",
                                        enabled=w.isEnabled(),
                                        trigger=w.click,
                                    )
                                )
            elif plain_as == "footer":
                menus.append(
                    _MenuEntry(
                        label=action.text() or "…",
                        trigger=action.trigger,
                        close_after=True,
                        tooltip=action.text() or "",
                    )
                )
            else:
                cells.append(
                    _CellEntry(
                        icon=action.icon(),
                        label=action.text() or "",
                        enabled=action.isEnabled(),
                        trigger=action.trigger,
                    )
                )
        return cells, menus

    def _build_all_flyouts(self) -> None:
        for spec in self._specs:
            if spec.menu is None:
                continue
            flyout = self._flyouts.get(spec.key)
            if flyout is None:
                flyout = Flyout()
                self._flyouts[spec.key] = flyout
            self._populate_flyout(spec, flyout)

    def _populate_flyout(self, spec: _RailSpec, flyout: Flyout) -> None:
        cells, menu_entries = self._introspect_menu(spec.menu, spec.plain_as)
        items = [
            FlyoutItem(
                item_id=f"{spec.group}-{index}",
                icon=cell.icon,
                label=cell.label,
                enabled=cell.enabled,
                tooltip=cell.label,
                trigger=cell.trigger,
            )
            for index, cell in enumerate(cells)
        ]
        footer = None
        if menu_entries:
            buttons: list[FlyoutFooterButton] = []
            for entry in menu_entries:
                text = entry.label
                if entry.tooltip and entry.label == "…":
                    # Submenú: etiqueta corta + tooltip original.
                    short = _short_menu_label(entry.tooltip)
                    text = short
                if spec.footer_labels and entry.label in spec.footer_labels:
                    text = spec.footer_labels[entry.label]
                buttons.append(
                    FlyoutFooterButton(
                        text=text,
                        callback=entry.trigger,
                        close_after=entry.close_after,
                    )
                )
            footer = FlyoutFooter(buttons=buttons)
        flyout.populate(
            title=spec.flyout_title or spec.label,
            items=items,
            columns=spec.columns,
            on_select=None,
            close_on_select=True,
            footer=footer,
        )

    # ------------------------------------------------------------------
    # API pública
    # ------------------------------------------------------------------
    def open_flyout(self, key: str) -> Flyout | None:
        """Abre el flyout de ``key`` junto al botón del rail."""
        spec = next((s for s in self._specs if s.key == key), None)
        if spec is None or spec.menu is None:
            return None
        button = self._buttons.get(key)
        flyout = self._flyouts.get(key)
        if flyout is None or button is None or self.window() is None:
            return None
        if flyout.isVisible():
            flyout.close_with(None)
            return flyout
        theme = self._last_theme or "light"
        flyout.set_theme(theme)
        flyout.show_near(button, self.window())
        return flyout

    def set_active_tool(self, tool_id: str) -> None:
        """Sincroniza el highlight activo con el ``tool_id`` emitido."""
        group = _group_for_tool(tool_id)
        for button in self._buttons.values():
            button.set_active(False)
        self._current_group = group
        if group is None:
            return

        if group == "select":
            # Icono pointer/lasso: re-leído del toolbar (fuente de verdad).
            select_btn = self._buttons.get("select")
            if select_btn is not None:
                select_btn.setIcon(self._toolbar.select_action.icon())
        button = self._buttons.get(group)
        if button is not None:
            if group in {"bond", "ring", "atom", "brackets", "arrows",
                        "plates", "symbols", "energy", "orbitals", "text"}:
                action = self._action_for_group(group)
                if action is not None:
                    button.setIcon(action.icon())
            button.set_active(True)

    def clear_active(self) -> None:
        """Quita el highlight (cambio de pestaña; canvas → ``tool_none``)."""
        for button in self._buttons.values():
            button.set_active(False)
        self._current_group = None

    def refresh_icons(self, theme_name: str | None = None) -> None:
        """Regenera iconos de botones y flyouts (tras cambio de tema).

        Re-lee los iconos de los ``QAction`` de los toolbars (ya
        regenerados por sus ``refresh_icons()``) y reconstruye los flyouts
        (los menús originales son reconstruidos por los toolbars).
        """
        if theme_name:
            self._last_theme = theme_name
        for spec in self._specs:
            button = self._buttons.get(spec.key)
            if button is None or spec.icon_action is None:
                continue
            button.setIcon(spec.icon_action.icon())
        self._build_all_flyouts()

    def shortcut_map(self) -> dict[int, Callable[[], None]]:
        """Mapeo de teclas de atajo → callback (mismo que el clic en el rail)."""
        t, s = self._toolbar, self._symbols

        def _lasso() -> None:
            icon, tip = t._selection_meta.get(
                "tool_select_lasso", (None, "Seleccion libre")
            )
            t._select_selection_palette("tool_select_lasso", icon, tip)

        return {
            int(Qt.Key.Key_V): t.select_action.trigger,
            int(Qt.Key.Key_A): _lasso,
            int(Qt.Key.Key_L): t.chain_action.trigger,
            int(Qt.Key.Key_B): t.bond_action.trigger,
            int(Qt.Key.Key_R): t.ring_action.trigger,
            int(Qt.Key.Key_C): t.label_action.trigger,
            int(Qt.Key.Key_T): s.text_action.trigger,
            int(Qt.Key.Key_N): s.annotation_action.trigger,
            int(Qt.Key.Key_G): s.bracket_action.trigger,
            int(Qt.Key.Key_E): s.energy_diagram_action.trigger,
            int(Qt.Key.Key_O): s.orbital_action.trigger,
        }

    def button_count(self) -> int:
        """Número de botones del rail (auditoría de paridad)."""
        return len(self._buttons)

    # ------------------------------------------------------------------
    # Internos
    # ------------------------------------------------------------------
    _last_theme: str | None = None

    def _action_for_group(self, group: str) -> QAction | None:
        t, s = self._toolbar, self._symbols
        return {
            "select": t.select_action,
            "bond": t.bond_action,
            "ring": t.ring_action,
            "atom": t.label_action,
            "text": s.text_action,
            "arrows": s.annotation_action,
            "brackets": s.bracket_action,
            "plates": s.plate_action,
            "symbols": s.symbol_action,
            "energy": s.energy_diagram_action,
            "orbitals": s.orbital_action,
        }.get(group)


def _short_menu_label(title: str) -> str:
    """Etiqueta corta para botones de pie de submenús."""
    mapping = {
        "Electronic Diagrams": "Diagrams ▾",
        "Electronic Diagram Presets": "Presets ▾",
        "Color de etiquetas": "Colores ▾",
    }
    return mapping.get(title, title if len(title) <= 14 else f"{title[:13]}…")


def _group_for_tool(tool_id: str) -> str | None:
    """Resuelve el grupo del rail para un ``tool_id`` (señal o normalizado)."""
    tool_id = str(tool_id or "tool_none")
    if tool_id in {"tool_select", "tool_select_lasso", "tool_rotate_3d_precise"}:
        return "select"
    if tool_id == "tool_bond" or tool_id.startswith("bond_"):
        return "bond"
    if tool_id == "tool_ring":
        return "ring"
    if tool_id == "tool_atom" or tool_id.startswith("atom_"):
        return "atom"
    if tool_id == "tool_chain":
        return "chain"
    if tool_id == "tool_coordination_center" or tool_id.startswith("coord_"):
        return "coord"
    if tool_id == "tool_text":
        return "text"
    if tool_id == "tool_arrow" or tool_id.startswith("tool_arrow_"):
        return "arrows"
    if tool_id == "tool_brackets" or tool_id.startswith("tool_brackets_"):
        return "brackets"
    if tool_id in {"tool_tlc", "tool_electrophoresis"}:
        return "plates"
    if tool_id.startswith("tool_charge") or tool_id.startswith("tool_symbol_"):
        return "symbols"
    if (
        tool_id == "tool_energy_diagram"
        or tool_id.startswith("tool_energy_diagram_")
    ):
        return "energy"
    if tool_id == "tool_orbital" or tool_id.startswith("tool_orbital_"):
        return "orbitals"
    return None


def _as_widget_action(action: QAction):
    """Devuelve ``action`` como ``QWidgetAction`` o ``None``."""
    if isinstance(action, QWidgetAction):
        return action
    return None


class ToolShortcutDispatcher(QObject):
    """Event filter contextual para atajos de letra simple (sin ``QShortcut``).

    Filtros (KeyPress sobre la ventana):

    1. Tecla en el mapeo y ``modifiers() == NoModifier``.
    2. Sin diálogo modal activo.
    3. El widget con foco no es de entrada de texto
       (``QLineEdit``/``QTextEdit``/``QPlainTextEdit``/``QComboBox``/
       ``QAbstractSpinBox``).
    """

    def __init__(
        self,
        window: QWidget,
        mappings: dict[int, Callable[[], None]],
        parent: QObject | None = None,
        suppress_predicate: Callable[[], bool] | None = None,
    ) -> None:
        """Crea el dispatcher.

        Args:
            window: Ventana sobre la que se filtran los eventos.
            mappings: Tecla (int) → callback de herramienta.
            parent: Parent QObject opcional.
            suppress_predicate: Devuelve ``True`` cuando el atajo NO debe
                consumirse (p. ej. edición de texto en el canvas). El
                dispatcher simplemente deja de consumir el evento: no
                termina la edición ni cambia de herramienta; la tecla sigue
                su curso normal hacia el receptor con foco.
        """
        super().__init__(parent)
        self._window = window
        self._mappings = mappings
        self._suppress_predicate = suppress_predicate
        app = QApplication.instance()
        if app is not None:
            app.installEventFilter(self)

    def _top_level_of(self, obj: QObject) -> QObject | None:
        node: QObject | None = obj
        top = node
        while top is not None and top.parent() is not None:
            top = top.parent()
        return top

    def eventFilter(self, obj: QObject, event: QEvent) -> bool:
        if event.type() != QEvent.Type.KeyPress:
            return False
        key_event = event  # QKeyEvent
        callback = self._mappings.get(int(key_event.key()))
        if callback is None:
            return False
        if (
            key_event.modifiers()
            != Qt.KeyboardModifier.NoModifier
        ):
            return False
        app = QApplication.instance()
        if app is not None and app.activeModalWidget() is not None:
            return False
        focus = app.focusWidget() if app is not None else None
        if focus is not None and isinstance(
            focus,
            (
                QLineEdit,
                QTextEdit,
                QPlainTextEdit,
                QComboBox,
                QAbstractSpinBox,
            ),
        ):
            return False
        if self._top_level_of(obj) is not self._window:
            return False
        if self._suppress_predicate is not None and self._suppress_predicate():
            return False
        callback()
        event.accept()
        return True
