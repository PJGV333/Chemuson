"""Paleta de comandos (Ctrl+K) de la Fase 6 del plan de modernización.

Responsabilidad: **presentar, filtrar y ejecutar** `QAction` existentes de la
ventana. La `QAction` es la fuente de verdad: la paleta conserva su identidad,
estado ``enabled``, ``checkable``/``checked``, shortcut mostrado, icono,
``triggered`` y handlers, y ejecuta la acción por delegación a ``trigger`` (o
``activate``) sin generar una segunda ``QAction`` ni duplicar handlers.

La paleta NO implementa lógica química, NO conoce internals de ``gui.canvas``
ni llama a Clean2D directamente: depende únicamente de ``QAction``/metadatos,
``theme`` (tokens/QSS/IconProvider) y el shell (la ventana que la monta).

Lenguaje visual del spike aprobado
(``docs/ui-modernization/pyqt6-spike``): overlay hijo de la ventana con tarjeta
centrada (~560 px, limitado al ancho de la ventana), input superior y lista de
resultados con icono/título/sección/atajo. QSS por tokens en ``theme/qss.py``
(``#commandPalette``). Ver OpenSpec
``2026-10-01-modernize-ui-command-palette``.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Sequence

from PyQt6.QtCore import Qt, QEvent, pyqtSignal
from PyQt6.QtGui import QAction
from PyQt6.QtWidgets import (
    QFrame,
    QHBoxLayout,
    QLabel,
    QLineEdit,
    QScrollArea,
    QVBoxLayout,
    QWidget,
)

from chemuson.gui.theme import METRICS
from chemuson.gui.theme.icon_provider import IconProvider

__all__ = ["CommandEntry", "CommandRegistry", "CommandPalette"]

_PALLETTE_ICON_SIZE = 18
_PALLETTE_SEARCH_ICON = 16
_MAX_RESULTS = 12
_ROW_MIN_HEIGHT = 34


@dataclass
class CommandEntry:
    """Una entrada presentacional que apunta a una ``QAction`` existente.

    Solo añade metadata de presentación (``section``, ``keywords``, ``icon``);
    no copia estado funcional. ``action`` es la fuente de verdad (identidad,
    enabled, checkable/checked, shortcut, icono, ``triggered``).
    """

    action: QAction
    section: str = "General"
    keywords: tuple[str, ...] = ()
    icon: str = ""

    @property
    def title(self) -> str:
        """Título mostrado: el texto de la ``QAction`` (fuente de verdad)."""
        return self.action.text()

    def _fields(self, query: str) -> tuple[bool, bool, bool, bool]:
        """Flags (prefijo/substring) sobre título, keywords y sección.

        Devuelve ``(title_prefix, other_prefix, title_substr, other_substr)``.
        ``query`` debe venir ya normalizado (lower, sin espacios).
        """
        title = self.title.lower().strip()
        kws = [k.lower().strip() for k in self.keywords if k]
        section = self.section.lower().strip()
        title_prefix = title.startswith(query)
        other_prefix = any(k.startswith(query) for k in kws) or section.startswith(query)
        title_substr = query in title
        # substring en keywords o sección
        other_substr = any(query in k for k in kws) or (query in section)
        return title_prefix, other_prefix, title_substr, other_substr

    def rank(self, query: str) -> int:
        """Tier de ranking (menor = mejor): 0 título-prefijo, 1 otra-prefijo,
        2 título-substring, 3 otra-substring; 999 si no coincide."""
        tp, op, ts, osu = self._fields(query)
        if tp:
            return 0
        if op:
            return 1
        if ts:
            return 2
        if osu:
            return 3
        return 999


class CommandRegistry:
    """Registro de comandos con deduplicación por identidad de ``QAction``.

    ``register`` indexa por ``id(action)``: registrar la misma ``QAction`` más
    de una vez es no-op (se conserva la primera inscripción), de modo que no
    existe una segunda ``QAction``/entrada por comando.
    """

    def __init__(self) -> None:
        self._entries: list[CommandEntry] = []
        self._by_action: dict[int, CommandEntry] = {}

    def register(
        self,
        action: QAction,
        section: str = "General",
        keywords: Sequence[str] = (),
        icon: str = "",
    ) -> CommandEntry | None:
        """Registra una ``QAction``; devuelve la entrada o ``None`` si ya existía."""
        if action is None:
            return None
        key = id(action)
        if key in self._by_action:
            return None  # dedup por identidad: no duplicar la misma QAction
        entry = CommandEntry(
            action=action, section=section, keywords=tuple(keywords), icon=icon
        )
        self._entries.append(entry)
        self._by_action[key] = entry
        return entry

    def entries(self) -> list[CommandEntry]:
        return list(self._entries)

    def unique_action_ids(self) -> list[int]:
        """Identidades únicas de ``QAction`` representadas (sin repetidos)."""
        return [id(e.action) for e in self._entries]

    def count(self) -> int:
        return len(self._entries)

    def find_by_action(self, action: QAction) -> CommandEntry | None:
        return self._by_action.get(id(action))


class CommandPalette(QFrame):
    """Overlay de comandos: input + lista filtrada, centrado sobre la ventana."""

    #: Emite la ``QAction`` ejecutada (para tests y para que la ventana sepa).
    commandExecuted = pyqtSignal(object)

    def __init__(
        self,
        registry: CommandRegistry,
        *,
        parent: QWidget | None = None,
    ) -> None:
        super().__init__(parent)
        self.setObjectName("commandPalette")
        self.setAttribute(Qt.WidgetAttribute.WA_StyledBackground, True)
        self._icons = IconProvider()
        self._theme_name = "light"
        self.registry = registry
        self._filtered: list[CommandEntry] = []
        self._idx = 0
        self._rows: list[QFrame] = []
        self._last_action: QAction | None = None
        self._open = False

        layout = QVBoxLayout(self)
        layout.setContentsMargins(0, 0, 0, 0)
        layout.setSpacing(0)

        self.card = QFrame(self)
        self.card.setObjectName("paletteCard")
        self.card_layout = QVBoxLayout(self.card)
        self.card_layout.setContentsMargins(0, 0, 0, 0)
        self.card_layout.setSpacing(0)
        layout.addWidget(self.card)

        # --- fila de input -------------------------------------------------
        self.input_row = QFrame(self.card)
        self.input_row.setObjectName("paletteInputRow")
        ir = QHBoxLayout(self.input_row)
        ir.setContentsMargins(14, 11, 14, 11)
        ir.setSpacing(10)
        self._search_ic = QLabel(self.input_row)
        self.input = QLineEdit(self.input_row)
        self.input.setObjectName("paletteInput")
        self.input.setFrame(False)
        self.input.setPlaceholderText("Buscar o ejecutar un comando…")
        self.input.setTextMargins(0, 0, 0, 0)
        self.input.setClearButtonEnabled(False)
        self._esc_kbd = QLabel("Esc", self.input_row)
        self._esc_kbd.setObjectName("paletteKbd")
        ir.addWidget(self._search_ic)
        ir.addWidget(self.input, 1)
        ir.addWidget(self._esc_kbd)
        self.card_layout.addWidget(self.input_row)

        # --- lista ---------------------------------------------------------
        self.scroll_area = QScrollArea(self.card)
        self.scroll_area.setObjectName("paletteScroll")
        self.scroll_area.setWidgetResizable(True)
        self.scroll_area.setFrameShape(QFrame.Shape.NoFrame)
        self.scroll_area.setHorizontalScrollBarPolicy(Qt.ScrollBarPolicy.ScrollBarAlwaysOff)
        self._list_host = QWidget()
        self._list_host.setObjectName("paletteList")
        self._list_lay = QVBoxLayout(self._list_host)
        self._list_lay.setContentsMargins(6, 4, 6, 6)
        self._list_lay.setSpacing(2)
        self._list_lay.addStretch(1)
        self.scroll_area.setWidget(self._list_host)
        self.card_layout.addWidget(self.scroll_area, 1)

        self.empty = QLabel("Sin resultados", self.card)
        self.empty.setObjectName("paletteEmpty")
        self.empty.setAlignment(Qt.AlignmentFlag.AlignCenter)
        self.empty.setContentsMargins(20, 18, 20, 18)
        self.card_layout.addWidget(self.empty)

        self.input.textChanged.connect(self._on_query)
        # El foco está en el ``QLineEdit`` al abrir la paleta; intercepta
        # ↑/↓/Enter/Esc desde el input para navegar/ejecutar/cerrar (de lo
        # contrario esas teclas no llegarían a ``keyPressEvent`` de la paleta).
        self.input.installEventFilter(self)
        self._refresh_theme_icons()
        self.hide()

    # ------------------------------------------------------------------
    # Tema (tokens; sin colores hardcodeados fuera del sistema)
    # ------------------------------------------------------------------
    def _refresh_theme_icons(self) -> None:
        self._search_ic.setPixmap(
            self._icons.pixmap("search", self._color("text3"), _PALLETTE_SEARCH_ICON)
        )

    def _color(self, role: str) -> str:
        return self._icons.theme_color(self._theme_name, role)

    def refresh_theme(self, theme_name: str) -> None:
        """Re-estiliza la paleta para el tema resuelto (``"light"``/``"dark"``)."""
        resolved = "dark" if theme_name == "dark" else "light"
        self._theme_name = resolved
        self._refresh_theme_icons()
        for row in self._rows:
            row.setProperty("cls", "paletteItem")
            row.setProperty("selected", row.property("_is_selected"))
            row.style().unpolish(row)
            row.style().polish(row)

    # ------------------------------------------------------------------
    # Apertura / cierre
    # ------------------------------------------------------------------
    def open(self) -> None:
        """Muestra la paleta centrada sobre la ventana con foco en el input."""
        self._open = True
        self.setGeometry(self.parent().rect())
        self.input.clear()
        self._idx = 0
        self._apply_filter("")
        self._rebuild()
        self.show()
        self.raise_()
        self._position_card()
        self.input.setFocus()

    def close_overlay(self) -> None:
        self._open = False
        self.hide()

    def is_open(self) -> bool:
        return self._open

    def _position_card(self) -> None:
        """Centra la tarjeta: ancho min(paletteW, 92% de la ventana), dentro
        del rect de la ventana (no se sale en ventanas estrechas)."""
        pr = self.parent().rect()
        pw, ph = pr.width(), pr.height()
        target_w = METRICS.get("paletteW", 560)
        cw = min(target_w, max(320, int(pw * 0.92)))
        cw = min(cw, max(320, pw - 24))
        # altura acotada al rect de la ventana (11% de margen superior)
        hint_h = self.card.sizeHint().height()
        max_h = max(140, int(ph * 0.78))
        ch = max(140, min(hint_h, max_h))
        x = max(12, (pw - cw) // 2)
        y = max(12, int(ph * 0.12))
        # clamp inferior para no salirse por el borde
        if y + ch > ph - 12:
            y = max(12, ph - ch - 12)
        self.card.setGeometry(x, y, cw, ch)
        self.scroll_area.setVisible(bool(self._filtered))
        self.empty.setVisible(not self._filtered)

    # ------------------------------------------------------------------
    # Filtro + ranking
    # ------------------------------------------------------------------
    def _on_query(self, text: str) -> None:
        self._apply_filter(text)
        self._idx = 0
        self._rebuild()

    def _apply_filter(self, text: str) -> None:
        query = " ".join(text.lower().split())
        all_entries = self.registry.entries()
        if not query:
            # query vacío: último comando usado (sesión) primero, luego catálogo
            ordered = list(all_entries)
            if self._last_action is not None:
                last_entry = self.registry.find_by_action(self._last_action)
                if last_entry is not None:
                    ordered = [e for e in ordered if e is not last_entry]
                    ordered.insert(0, last_entry)
            self._filtered = ordered
            return
        scored = [
            (e.rank(query), i, e)
            for i, e in enumerate(all_entries)
            if e.rank(query) != 999
        ]

        def _sort_key(item):
            tier, i, e = item
            return (tier, e.section, e.title, i)

        scored.sort(key=_sort_key)
        self._filtered = [e for (_t, _i, e) in scored]

    def filter_count(self) -> int:
        return len(self._filtered)

    def _rebuild(self) -> None:
        # limpiar filas anteriores
        for row in self._rows:
            row.setParent(None)
            row.deleteLater()
        self._rows = []
        # limpiar labels de sección residuales. ``setParent(None)`` es síncrono
        # (el header sale del árbol de ``_list_host`` al instante); solo con
        # ``deleteLater()`` (asíncrono) cabría que headers viejos se pintaran
        # junto a la lista nueva hasta que el event loop procese la eliminación.
        for lab in self._list_host.findChildren(QLabel):
            lab.setParent(None)
            lab.deleteLater()

        t = self._theme_name
        last_section: str | None = None
        for i, entry in enumerate(self._filtered):
            if entry.section != last_section:
                sec = QLabel(entry.section.upper(), self._list_host)
                sec.setObjectName("paletteSection")
                sec.setContentsMargins(10, 8, 10, 3)
                self._list_lay.insertWidget(self._list_lay.count() - 1, sec)
                last_section = entry.section
            row = self._build_row(entry, i)
            self._list_lay.insertWidget(self._list_lay.count() - 1, row)
            self._rows.append(row)

        self.empty.setText(f'Sin resultados para "{self.input.text()}"')
        self.scroll_area.setVisible(bool(self._filtered))
        self.empty.setVisible(not self._filtered)
        self._scroll_selected_into_view()
        if self._open:
            self._position_card()

    def _build_row(self, entry: CommandEntry, index: int) -> QFrame:
        row = QFrame(self._list_host)
        row.setProperty("cls", "paletteItem")
        row.setProperty("selected", index == self._idx)
        row.setProperty("_is_selected", index == self._idx)
        row.setMinimumHeight(_ROW_MIN_HEIGHT)
        row.setEnabled(entry.action.isEnabled())
        row.setCursor(Qt.CursorShape.PointingHandCursor)
        rl = QHBoxLayout(row)
        rl.setContentsMargins(10, 7, 10, 7)
        rl.setSpacing(10)

        icon_label = QLabel(row)
        icon_label.setFixedSize(_PALLETTE_ICON_SIZE, _PALLETTE_ICON_SIZE)
        icon_name = entry.icon or "flask"
        tint = self._color("text2")
        icon_label.setPixmap(self._icons.pixmap(icon_name, tint, _PALLETTE_ICON_SIZE))
        rl.addWidget(icon_label)

        title = QLabel(entry.title, row)
        title.setObjectName("paletteTitle")
        title.setToolTip(self._action_tooltip(entry))
        rl.addWidget(title, 1)

        section_label = QLabel(entry.section, row)
        section_label.setObjectName("paletteSectionCell")
        rl.addWidget(section_label)

        shortcut = entry.action.shortcut().toString()
        if shortcut:
            kbd = QLabel(shortcut, row)
            kbd.setObjectName("paletteKbd")
            rl.addWidget(kbd)

        row.mousePressEvent = lambda _e, i=index: self._execute_index(i)
        return row

    @staticmethod
    def _action_tooltip(entry: CommandEntry) -> str:
        sc = entry.action.shortcut().toString()
        base = entry.title
        if entry.action.isEnabled() is False:
            return f"{base} (no disponible)"
        if sc:
            return f"{base}  ({sc})"
        return base

    # ------------------------------------------------------------------
    # Navegación
    # ------------------------------------------------------------------
    def _move(self, delta: int) -> None:
        if not self._filtered:
            return
        n = len(self._filtered)
        self._idx = (self._idx + delta) % n
        for i, row in enumerate(self._rows):
            on = i == self._idx
            row.setProperty("selected", on)
            row.setProperty("_is_selected", on)
            row.style().unpolish(row)
            row.style().polish(row)
        self._scroll_selected_into_view()

    def _scroll_selected_into_view(self) -> None:
        if 0 <= self._idx < len(self._rows):
            row = self._rows[self._idx]
            vh = self.scroll_area.viewport().height()
            sb = self.scroll_area.verticalScrollBar()
            sb.setValue(max(0, row.y() - vh // 2))

    # ------------------------------------------------------------------
    # Ejecución
    # ------------------------------------------------------------------
    def _execute_index(self, index: int) -> None:
        if 0 <= index < len(self._filtered):
            self._idx = index
            self._execute_selected()

    def _execute_selected(self) -> None:
        if not (0 <= self._idx < len(self._filtered)):
            return
        entry = self._filtered[self._idx]
        if not entry.action.isEnabled():
            # QAction disabled: no se ejecuta (se queda la paleta abierta).
            return
        self._last_action = entry.action
        self.close_overlay()
        entry.action.trigger()
        self.commandExecuted.emit(entry.action)

    # ------------------------------------------------------------------
    # Teclado
    # ------------------------------------------------------------------
    def _handle_navigation_key(self, key: int) -> bool:
        """Procesa una tecla de navegación (↑/↓/Enter/Esc). True si la consumió."""
        if key == Qt.Key.Key_Escape:
            self.close_overlay()
            return True
        if key == Qt.Key.Key_Down:
            self._move(1)
            return True
        if key == Qt.Key.Key_Up:
            self._move(-1)
            return True
        if key in (Qt.Key.Key_Return, Qt.Key.Key_Enter):
            self._execute_selected()
            return True
        return False

    def eventFilter(self, obj, event) -> bool:  # noqa: N802
        """Intercepta ↑/↓/Enter/Esc desde el ``QLineEdit`` de la paleta."""
        if self._open and obj is self.input and event.type() == QEvent.Type.KeyPress:
            if self._handle_navigation_key(event.key()):
                event.accept()
                return True
        return super().eventFilter(obj, event)

    def keyPressEvent(self, event) -> None:  # noqa: N802
        if not self._open:
            super().keyPressEvent(event)
            return
        if self._handle_navigation_key(event.key()):
            event.accept()
            return
        super().keyPressEvent(event)
