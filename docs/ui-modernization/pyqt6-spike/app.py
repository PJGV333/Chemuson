"""Spike visual PyQt6 de la propuesta de modernización de UI de Chemuson.

Maqueta independiente: NO importa nada de `src/chemuson/`, no se conecta con la
lógica real y no modifica producción. Reproduce la apariencia de
`docs/ui-modernization/mockup-ui.html` para evaluar si PyQt6/QtWidgets sostiene
el lenguaje visual propuesto en PLAN.md.

Ejecución:
    python docs/ui-modernization/pyqt6-spike/app.py
    python docs/ui-modernization/pyqt6-spike/app.py --theme dark --size 1440x900
    python docs/ui-modernization/pyqt6-spike/app.py --smoke /tmp/spike_shots
"""
from __future__ import annotations

import sys
from pathlib import Path

from PyQt6.QtCore import QEvent, QSize, Qt, QTimer
from PyQt6.QtGui import QKeySequence, QShortcut
from PyQt6.QtWidgets import (
    QApplication, QFrame, QGridLayout, QHBoxLayout, QLabel, QLineEdit,
    QPushButton, QScrollArea, QStackedWidget, QTabBar, QToolButton,
    QVBoxLayout, QWidget,
)

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))

from canvas_demo import CanvasDemo, ZOOM_STEPS  # noqa: E402
from icons import IconProvider  # noqa: E402
from palette import CommandPalette, PaletteAction  # noqa: E402
from theme import METRICS, Theme  # noqa: E402
from widgets import (  # noqa: E402
    ChipToggle, Flyout, Kbd, RailButton, SearchPill, SideTabRow,
    make_bar, make_issue, make_kv, make_pill, make_row, make_tpl_card, shadow,
)

# ---------------------------------------------------------------------------
# datos de demostración (espejo del mockup)
# ---------------------------------------------------------------------------
FLYOUTS = {
    "bond": dict(title="Enlace · estilo (11)", cols=4, icon="bond",
                 items=[("single", "Simple"), ("double", "Doble"), ("triple", "Triple"),
                        ("partial", "Parcial"), ("coordinate", "Coordinación"),
                        ("wavy", "Ondulado"), ("hashed", "Truncado"), ("bold", "Grueso"),
                        ("wedge", "Cuña"), ("wedgeHashed", "Cuña sombreada"),
                        ("dotted", "Punteado")]),
    "ring": dict(title="Anillos", cols=5, icon="ring",
                 foot=("customRing", "Anillo personalizado…"),
                 items=[("benzene", "Benceno"), ("pyridine", "Piridina"), ("r3", "3"),
                        ("r4", "4"), ("r5", "5"), ("r6", "6"), ("r7", "7"), ("r8", "8"),
                        ("r9", "9"), ("r10", "10"), ("r11", "11"), ("r12", "12")]),
    "atom": dict(title="Átomos (10)", cols=5, icon="atom",
                 foot=("periodic", "Tabla periódica…"),
                 items=[("C", "C"), ("N", "N"), ("O", "O"), ("S", "S"), ("P", "P"),
                        ("F", "F"), ("Cl", "Cl"), ("Br", "Br"), ("I", "I"), ("H", "H")]),
    "arrow": dict(title="Flechas (15)", cols=4, icon="arrow",
                 items=[("forward", "Hacia adelante"), ("half", "Media"),
                        ("curved", "Curvada"), ("curvedHalf", "Curva media"),
                        ("fishhook", "Fishhook"), ("resonance", "Resonancia"),
                        ("equilibrium", "Equilibrio"), ("pushing", "Pushing"),
                        ("pushingHalf", "Pushing media"), ("retro", "Retro"),
                        ("doubleHead", "Doble cabeza"), ("left", "Izquierda"),
                        ("right", "Derecha"), ("up", "Arriba"), ("down", "Abajo")]),
    "bracket": dict(title="Corchetes (10)", cols=5, icon="bracket",
                 items=[("b1", "+"), ("b2", "2+"), ("b3", "3+"), ("b-1", "−"),
                        ("b-2", "2−"), ("b-3", "3−"), ("br1", "•"), ("br2", "2•"),
                        ("bn", "n"), ("b0", "∅")]),
    "symbol": dict(title="Símbolos", cols=5, icon="symbol",
                 items=[("s+", "+"), ("s-", "−"), ("sr", "•"), ("sp1", ":"),
                        ("sp2", "··"), ("sd+", "δ+"), ("sd-", "δ−"), ("sH", "H")]),
    "plate": dict(title="Placas", cols=2, icon="tlc", per_item=True,
                 items=[("tlc", "Cromatografía (TLC)"), ("gel", "Gel de electroforesis")]),
}

# (grupo, nombre, clave icono SVG, atajo, tipo: tool|flyout|action)
RAIL = [
    ("select", "Seleccionar", "pointer", "V", "tool"),
    ("lasso", "Selección área", "lasso", "A", "tool"),
    ("SEP", None, None, None, None),
    ("bond", "Enlace", "bond", "B", "flyout"),
    ("chain", "Cadena", "chain", "L", "tool"),
    ("ring", "Anillo", "ring", "R", "flyout"),
    ("atom", "Átomo", "atom", "C", "flyout"),
    ("coord", "Centro de coordinación", "coord", "", "tool"),
    ("3d", "Rotación 3D", "cube", "", "action"),
    ("SEP", None, None, None, None),
    ("text", "Texto", "text", "T", "tool"),
    ("arrow", "Flecha", "arrow", "N", "flyout"),
    ("bracket", "Corchetes", "bracket", "G", "flyout"),
    ("symbol", "Símbolos", "symbol", "", "flyout"),
    ("nums", "Numeración", "num", "", "action"),
    ("SEP", None, None, None, None),
    ("energy", "Diagrama de energía", "energy", "E", "action"),
    ("orbital", "Orbitales", "orbital", "O", "action"),
    ("plate", "Placas", "tlc", "", "flyout"),
    ("SEP", None, None, None, None),
    ("clean2d", "Limpiar 2D", "clean", "", "action"),
    ("validate", "Validar", "shield", "", "action"),
]


def _fly_items(spec: dict) -> list[tuple[str, str, str]]:
    """(id, label, visual): visual = clave SVG o 'glyph:<texto>'."""
    out = []
    for iid, label in spec["items"]:
        if spec.get("per_item"):
            visual = iid if iid in ("tlc", "gel") else spec["icon"]
        elif spec["icon"] in ("atom", "bracket", "symbol"):
            visual = f"glyph:{label}||13"
        else:
            visual = spec["icon"]
        out.append((iid, label, visual))
    return out


class Toast:
    """Mensaje flotante breve (equivalente al .toast del mockup)."""

    def __init__(self, host: QWidget, theme_getter):
        self.label = QLabel(host)
        self.label.setObjectName("toast")
        self.label.hide()
        self.label.setStyleSheet(
            "QLabel#toast { background: %s; color: %s; font-size: 12px; font-weight: 600; "
            "border-radius: 9px; padding: 8px 14px; }"
            % (theme_getter()["text1"], theme_getter()["bg"]))
        self._tg = theme_getter
        self._timer = QTimer(host)
        self._timer.setSingleShot(True)
        self._timer.timeout.connect(self.label.hide)

    def show(self, msg: str) -> None:
        self.label.setText(msg)
        w, h = self.label.sizeHint().width(), self.label.sizeHint().height()
        host = self.label.parent()
        self.label.move((host.width() - w) // 2, host.height() - h - 46)
        self.label.show()
        self.label.raise_()
        self._timer.start(1800)

    def refresh_theme(self) -> None:
        self.label.setStyleSheet(
            "QLabel#toast { background: %s; color: %s; font-size: 12px; font-weight: 600; "
            "border-radius: 9px; padding: 8px 14px; }"
            % (self._tg()["text1"], self._tg()["bg"]))


class CanvasHost(QWidget):
    """Contenedor del canvas: apila vista, chips, pill de zoom y toast."""

    def __init__(self, canvas: CanvasDemo, icons: IconProvider, theme_getter,
                 toast: Toast, on_toggle_grid, on_toggle_nums,
                 on_zoom, on_zoom_fit):
        super().__init__()
        self._zoom_btns: list = []
        self.icons = icons
        self._tg = theme_getter
        self.grid = QGridLayout(self)
        self.grid.setContentsMargins(0, 0, 0, 0)
        self.grid.addWidget(canvas.view, 0, 0)

        chips = QWidget()
        ch = QHBoxLayout(chips)
        ch.setContentsMargins(0, 0, 14, 0)
        ch.setSpacing(6)
        ch.addStretch(1)
        self.chip_grid = ChipToggle("Rejilla", icons, theme_getter, "grid", checked=True)
        self.chip_nums = ChipToggle("Números", icons, theme_getter, "num", checked=False)
        self.chip_grid.toggled.connect(on_toggle_grid)
        self.chip_nums.toggled.connect(on_toggle_nums)
        ch.addWidget(self.chip_grid); ch.addWidget(self.chip_nums)
        self.grid.addWidget(chips, 0, 0, Qt.AlignmentFlag.AlignTop | Qt.AlignmentFlag.AlignRight)

        self.pill = QFrame()
        self.pill.setObjectName("zoomPill")
        pl = QHBoxLayout(self.pill)
        pl.setContentsMargins(0, 0, 0, 0)
        pl.setSpacing(0)
        self.btn_out = self._zoom_btn("minus", on_zoom, -1, "Alejar")
        self.val = QLabel("100 %")
        self.val.setObjectName("zoomVal")
        self.val.setAlignment(Qt.AlignmentFlag.AlignCenter)
        self.val.setFixedWidth(44)
        self.btn_in = self._zoom_btn("plus", on_zoom, 1, "Acercar")
        self.btn_fit = self._zoom_btn("fit", on_zoom_fit, 0, "Ajustar a la hoja")
        pl.addWidget(self.btn_out); pl.addWidget(self.val)
        pl.addWidget(self.btn_in); pl.addWidget(self.btn_fit)
        self.pill_wrapper = QWidget()
        pw = QHBoxLayout(self.pill_wrapper)
        pw.setContentsMargins(0, 0, 14, 14)
        pw.addWidget(self.pill)
        self.grid.addWidget(self.pill_wrapper, 0, 0,
                            Qt.AlignmentFlag.AlignBottom | Qt.AlignmentFlag.AlignRight)

        wrap = QWidget()
        wl = QVBoxLayout(wrap)
        wl.setContentsMargins(0, 0, 0, 0)
        wl.addStretch(1)
        wl.addWidget(toast.label, 0, Qt.AlignmentFlag.AlignHCenter)
        wl.addSpacing(46)
        self.grid.addWidget(wrap, 0, 0, Qt.AlignmentFlag.AlignBottom)

    def _zoom_btn(self, icon: str, fn, delta: int, tip: str) -> QToolButton:
        from theme import Theme
        b = QToolButton(self.pill)
        b.setProperty("cls", "zoomBtn")
        b.setIconSize(QSize(13, 13))
        b.setFixedSize(26, 30)
        b.setToolTip(tip)
        b.clicked.connect(lambda _=False: fn(delta))
        self._zoom_btns.append((b, icon))
        return b

    def refresh_theme(self) -> None:
        self.chip_grid.refresh_theme()
        self.chip_nums.refresh_theme()
        for b, icon in self._zoom_btns:
            b.setIcon(self.icons.icon(icon, self._tg()["text2"], 13))


class SpikeWindow(QWidget):
    def __init__(self, app: QApplication):
        super().__init__()
        self.setObjectName("rootWindow")
        self.setWindowTitle("Chemuson — propuesta de UI moderna (spike PyQt6)")
        self.dpr = app.primaryScreen().devicePixelRatio()
        self.icons = IconProvider(self.dpr)
        self._theme_name = "light"
        self.theme = Theme(self._theme_name)
        self._tg = lambda: self.theme

        # estado de herramientas (contrato conceptual con M09: tool_id)
        self.tool_group = "select"
        self.tool_name = "Seleccionar"
        self.tool_sub: str | None = None
        self.tool_key = "V"
        self._fly_active: dict[str, str | None] = {g: None for g in FLYOUTS}
        self._fly_open_group: str | None = None
        self._abar_icons: dict = {}

        self._build()
        self.set_theme("light")
        self._min_size = QSize(900, 560)
        self.setMinimumSize(self._min_size)

    # ------------------------------------------------------------------
    def _build(self) -> None:
        root = QVBoxLayout(self)
        root.setContentsMargins(0, 0, 0, 0)
        root.setSpacing(0)
        root.addWidget(self._build_appbar())

        main = QHBoxLayout()
        main.setContentsMargins(0, 0, 0, 0)
        main.setSpacing(0)
        main.addWidget(self._build_rail())
        self.toast = Toast(self, self._tg)
        self.canvas = CanvasDemo(self.theme, self.icons)
        self.host = CanvasHost(
            self.canvas, self.icons, self._tg, self.toast,
            on_toggle_grid=lambda on: self.canvas.set_grid(on),
            on_toggle_nums=lambda on: self.canvas.set_numbers(on),
            on_zoom=self._zoom_delta, on_zoom_fit=self._zoom_fit)
        main.addWidget(self.host, 1)
        main.addWidget(self._build_side())
        main_wrap = QWidget()
        main_wrap.setLayout(main)
        root.addWidget(main_wrap, 1)

        root.addWidget(self._build_statusbar())

        # flyout + paleta + atajos
        self.flyout = Flyout(self.icons, self._tg, parent=self)
        self.flyout.on_select = self._flyout_select
        self.palette = CommandPalette(self.icons, self._tg, self._make_actions(), parent=self)
        QShortcut(QKeySequence("Ctrl+K"), self).activated.connect(self.palette.show_overlay)
        app = QApplication.instance()
        if app is not None:
            app.installEventFilter(self)

        self.canvas.view.cursor_moved.connect(self._on_cursor)

    # ------------------------------------------------------------------
    # app bar
    # ------------------------------------------------------------------
    def _build_appbar(self) -> QWidget:
        bar = QFrame(); bar.setObjectName("appbar")
        bar.setFixedHeight(METRICS["appbarH"])
        lay = QHBoxLayout(bar); lay.setContentsMargins(12, 0, 12, 0); lay.setSpacing(12)

        brand_ic = QLabel()
        brand_ic.setPixmap(self.icons.pixmap("flask", self.theme["accent"], 24))
        brand_name = QLabel("Chemuson"); brand_name.setObjectName("brandName")
        ver = QLabel("0.4"); ver.setObjectName("verPill")
        self._brand_ic = brand_ic
        lay.addWidget(brand_ic); lay.addWidget(brand_name); lay.addWidget(ver)

        self.tabs = QTabBar(); self.tabs.setObjectName("docTabs")
        self.tabs.setExpanding(False)
        self.tabs.setDrawBase(False)
        self.tabs.setElideMode(Qt.TextElideMode.ElideRight)
        self._tab_counter = 2
        self.tabs.currentChanged.connect(self._on_tab_changed)
        self._add_tab("cafeina.cmsn", dirty=True)
        self._add_tab("Sin título 2", dirty=False)
        lay.addWidget(self.tabs)

        self.tab_new = QToolButton(); self.tab_new.setObjectName("tabNew")
        self.tab_new.setIconSize(QSize(14, 14))
        self.tab_new.setFixedSize(28, 28)
        self.tab_new.setToolTip("Nuevo documento (Ctrl+N)")
        self.tab_new.clicked.connect(lambda: self._add_tab(None))
        lay.addWidget(self.tab_new)

        self.search = SearchPill(self.icons, self._tg, on_click=lambda: self.palette.show_overlay())
        self.search.setFixedWidth(250)
        lay.addStretch(1)
        lay.addWidget(self.search)

        self.btn_undo = self._abar_btn("undo", "Deshacer (Ctrl+Z)", enabled=False)
        self.btn_redo = self._abar_btn("redo", "Rehacer (Ctrl+Y)")
        sep = QFrame(); sep.setObjectName("abarSep"); sep.setFixedHeight(22)
        self.btn_theme = self._abar_btn("moon", "Cambiar tema claro/oscuro")
        self.btn_theme.clicked.connect(self._toggle_theme)
        self.btn_settings = self._abar_btn("sliders", "Apariencia / ajustes")
        self.btn_settings.clicked.connect(lambda: self._open_side(4))
        lay.addWidget(self.btn_undo); lay.addWidget(self.btn_redo)
        lay.addWidget(sep)
        lay.addWidget(self.btn_theme); lay.addWidget(self.btn_settings)
        return bar

    def _abar_btn(self, icon: str, tip: str, enabled: bool = True) -> QToolButton:
        b = QToolButton(self)
        b.setProperty("cls", "abar")
        b.setIconSize(QSize(18, 18))
        b.setFixedSize(32, 32)
        b.setToolTip(tip)
        b.setEnabled(enabled)
        self._abar_icons[b] = icon
        return b

    # ------------------------------------------------------------------
    # pestañas de documento
    # ------------------------------------------------------------------
    def _add_tab(self, name: str | None, dirty: bool = True) -> int:
        if name is None:
            self._tab_counter += 1
            name = f"Sin título {self._tab_counter}"
        idx = self.tabs.addTab(self.icons.icon("mol", self.theme["text2"], 14), name)
        w = QWidget()
        wl = QHBoxLayout(w); wl.setContentsMargins(2, 0, 2, 0); wl.setSpacing(5)
        dot = QLabel(); dot.setFixedSize(7, 7)
        dot.setStyleSheet("background: %s; border-radius: 4px;" % self.theme["accent"])
        if not dirty:
            dot.hide()
        close = QToolButton(w)
        close.setProperty("cls", "tabClose")
        close.setIconSize(QSize(11, 11))
        close.setFixedSize(18, 18)
        close.setIcon(self.icons.icon("x", self.theme["text3"], 11))
        close.setToolTip("Cerrar pestaña")
        close.clicked.connect(lambda _=False, i=idx: self._close_tab(i))
        wl.addWidget(dot); wl.addWidget(close)
        self.tabs.setTabButton(idx, QTabBar.ButtonPosition.RightSide, w)
        self.tabs.setCurrentIndex(idx)
        return idx

    def _close_tab(self, idx: int) -> None:
        self.tabs.removeTab(idx)
        if self.tabs.count() == 0:
            self._add_tab(None, dirty=False)

    def _on_tab_changed(self, _idx: int) -> None:
        pass

    # ------------------------------------------------------------------
    # rail
    # ------------------------------------------------------------------
    def _build_rail(self) -> QWidget:
        rail = QFrame(); rail.setObjectName("rail")
        rail.setFixedWidth(METRICS["railW"])
        lay = QVBoxLayout(rail)
        lay.setContentsMargins(0, 10, 0, 10)
        lay.setSpacing(3)
        lay.addStretch(0)
        self._rail_buttons: dict[str, RailButton] = {}
        for group, name, icon, key, kind in RAIL:
            if kind is None:
                sep = QFrame(); sep.setObjectName("railSep")
                sep.setFixedSize(26, 1)
                lay.addWidget(sep, 0, Qt.AlignmentFlag.AlignHCenter)
                continue
            rb = RailButton(self.icons, self._tg, icon, key or "", name or "")
            rb._kind = kind  # type: ignore[attr-defined]
            rb.clicked.connect(lambda _=False, g=group: self._rail_clicked(g))
            lay.addWidget(rb, 0, Qt.AlignmentFlag.AlignHCenter)
            self._rail_buttons[group] = rb
        lay.addStretch(0)
        if self.tool_group in self._rail_buttons:
            self._rail_buttons[self.tool_group].set_active(True)
        return rail

    def _rail_clicked(self, group: str) -> None:
        kind = next(k for g, _n, _i, _k, k in RAIL if g == group)
        if kind == "flyout":
            if self.flyout.isVisible() and self._fly_open_group == group:
                self.flyout.hide()
                self._fly_open_group = None
                return
            spec = FLYOUTS[group]
            self._fly_open_group = group
            self.flyout.show_near(
                self._rail_buttons[group], group, spec["title"], spec["cols"],
                _fly_items(spec), foot=spec.get("foot"),
                active_id=self._fly_active.get(group))
            return
        if kind == "tool":
            self._activate_tool(group)
            return
        self._do_action(group)

    def _activate_tool(self, group: str, sub: str | None = None, sub_label: str | None = None) -> None:
        self.tool_group = group
        # sub_label es la etiqueta visible (mockup: 'Enlace · Doble'); sub es el id
        self.tool_sub = sub_label or sub
        self._fly_active[group] = sub
        name = next((n for g, n, _i, _k, _x in RAIL if g == group), group)
        self.tool_name = name
        self.tool_key = next((k for g, _n, _i, k, _x in RAIL if g == group), "") or ""
        for g, b in self._rail_buttons.items():
            b.set_active(g == group)
        self._update_status()
        if sub_label and group in FLYOUTS:
            self.toast.show(f"{name} · {sub_label}")
        self.flyout.set_active(sub if group == getattr(self, "_fly_open_group", None) else None)

    def _flyout_select(self, iid: str, label: str) -> None:
        group = self.flyout._group
        if iid == "periodic":
            self.toast.show("Tabla periódica (demo)")
            self.flyout.hide(); self._fly_open_group = None
            return
        if iid == "customRing":
            self.toast.show("Anillo personalizado (demo)")
            self.flyout.hide(); self._fly_open_group = None
            return
        self._activate_tool(group, iid, label)
        self.flyout.hide()
        self._fly_open_group = None

    def _do_action(self, action: str) -> None:
        if action == "3d":
            self.toast.show("Modo rotación 3D (demo)")
        elif action == "nums":
            self.host.chip_nums.toggle()
        elif action == "energy":
            self.toast.show("Diagramas de energía (demo)")
        elif action == "orbital":
            self.toast.show("Orbitales (demo)")
        elif action == "clean2d":
            self.toast.show("Clean2D: limpieza 2D (demo)")
        elif action == "validate":
            self._open_side(1)
            self.toast.show("Validación ejecutada (demo)")

    # ------------------------------------------------------------------
    # canvas
    # ------------------------------------------------------------------
    def _zoom_delta(self, d: int) -> None:
        if d > 0:
            self.canvas.zoom_in()
        else:
            self.canvas.zoom_out()
        self.host.val.setText(self.canvas.zoom_label())

    def _zoom_fit(self) -> None:
        self.canvas.zoom_fit()
        self.host.val.setText(self.canvas.zoom_label())

    def _on_cursor(self, p) -> None:
        self.lbl_cursor.setText(f"{round(p.x())} , {round(p.y())}")

    # ------------------------------------------------------------------
    # panel derecho
    # ------------------------------------------------------------------
    def _build_side(self) -> QWidget:
        wrap = QFrame(); wrap.setObjectName("sideWrap")
        wrap.setFixedWidth(METRICS["sideW"])
        lay = QVBoxLayout(wrap); lay.setContentsMargins(0, 0, 0, 0); lay.setSpacing(0)

        self.side_tabs = SideTabRow(self._tg)
        for label in ("Inspector", "Validación", "Propiedades", "Plantillas", "Apariencia"):
            self.side_tabs.add_tab(label)
        self.side_tabs.on_change = self._on_side_tab
        lay.addWidget(self.side_tabs)

        self.side_stack = QStackedWidget()
        self.side_stack.setObjectName("sideBody")
        lay.addWidget(self.side_stack, 1)

        self.side_stack.addWidget(self._page_inspector())
        self.side_stack.addWidget(self._page_validation())
        self.side_stack.addWidget(self._page_properties())
        self.side_stack.addWidget(self._page_templates())
        self.side_stack.addWidget(self._page_appearance())
        return wrap

    def _on_side_tab(self, i: int) -> None:
        self.side_stack.setCurrentIndex(i)

    def _open_side(self, i: int) -> None:
        self.side_tabs.set_active(i)

    def _scroll_page(self, inner: QWidget) -> QWidget:
        sc = QScrollArea(); sc.setObjectName("sideScroll")
        sc.setWidgetResizable(True)
        sc.setFrameShape(QFrame.Shape.NoFrame)
        host = QWidget(); host.setAttribute(Qt.WidgetAttribute.WA_TranslucentBackground, True)
        v = QVBoxLayout(host); v.setContentsMargins(14, 14, 14, 14); v.setSpacing(12)
        v.addWidget(inner)
        v.addStretch(1)
        sc.setWidget(host)
        return sc

    def _sec_title(self, text: str, pill: tuple[str, str] | None = None) -> QWidget:
        w = QWidget(); w.setAttribute(Qt.WidgetAttribute.WA_TranslucentBackground, True)
        h = QHBoxLayout(w); h.setContentsMargins(2, 4, 2, 8); h.setSpacing(8)
        t = QLabel(text.upper()); t.setObjectName("secTitle")
        h.addWidget(t)
        if pill:
            h.addWidget(make_pill(pill[0], pill[1]))
        h.addStretch(1)
        return w

    def _page_inspector(self) -> QWidget:
        inner = QWidget(); inner.setAttribute(Qt.WidgetAttribute.WA_TranslucentBackground, True)
        v = QVBoxLayout(inner); v.setContentsMargins(0, 0, 0, 0); v.setSpacing(14)

        v.addWidget(self._sec_title("Selección", ("1 grupo", "accent")))
        kv = QGridLayout(); kv.setSpacing(7)
        cards = [
            ("Átomo central", "N", "nitrógeno · Z 7"),
            ("Carga formal", "0", ""),
            ("Hibridación", "sp²", ""),
            ("Oxidación", "+3", ""),
        ]
        for i, (k, vv, s) in enumerate(cards):
            kv.addWidget(make_kv(k, vv, s), i // 2, i % 2)
        wrap_kv = QWidget(); wrap_kv.setLayout(kv)
        wrap_kv.setAttribute(Qt.WidgetAttribute.WA_TranslucentBackground, True)
        v.addWidget(wrap_kv)

        v.addWidget(self._sec_title("Geometría"))
        for k, val in (("Enlace C–N", "1.47 Å"), ("Enlace N=O", "1.20 Å"),
                       ("Ángulo C–N–O", "118.4°")):
            v.addWidget(make_row(k, val))
        v.addWidget(make_row("Densidad de enlaces", make_bar(0.62)))

        v.addWidget(self._sec_title("Acciones"))
        btns = QWidget(); bh = QHBoxLayout(btns); bh.setContentsMargins(0, 0, 0, 0); bh.setSpacing(8)
        b1 = self._btn("Limpiar 2D", primary=False)
        b1.clicked.connect(lambda: self.toast.show("Clean2D (demo)"))
        b2 = self._btn("Validar", primary=True)
        b2.clicked.connect(lambda: (self._open_side(1), self.toast.show("Validación (demo)")))
        bh.addWidget(b1); bh.addWidget(b2)
        v.addWidget(btns)
        v.addWidget(self._hint("Equivalente al dock «Inspector» actual, reutilizado sin reescribir su contenido (PLAN.md, fase 5)."))
        return self._scroll_page(inner)

    def _page_validation(self) -> QWidget:
        inner = QWidget(); inner.setAttribute(Qt.WidgetAttribute.WA_TranslucentBackground, True)
        v = QVBoxLayout(inner); v.setContentsMargins(0, 0, 0, 0); v.setSpacing(14)
        v.addWidget(self._sec_title("Resultados", ("2 avisos", "warn")))
        self._issues: list = []
        for sev, ttl, sub, chip, selected in (
            ("err", "Carga formal no neutra", "N con 4 enlaces en el grupo NO₂ (formal +1)", "Error", True),
            ("warn", "Longitud de enlace fuera de rango", "N=O 1.20 Å (esperado 1.15–1.19 Å)", "Aviso", False),
            ("ok", "Coherente con el SMILES de referencia", "Cc1ccc(cc1)[N+](=O)[O-]", "OK", False),
        ):
            issue = make_issue(sev, ttl, sub, chip, self._tg, selected=selected)
            issue.clicked.connect(lambda _=False, it=issue: self._select_issue(it))
            self._issues.append(issue)
            v.addWidget(issue)
        btns = QWidget(); bh = QHBoxLayout(btns); bh.setContentsMargins(0, 0, 0, 0); bh.setSpacing(8)
        bh.addWidget(self._btn("Revalidar todo", primary=True))
        bh.addWidget(self._btn("Filtrar"))
        v.addWidget(btns)
        v.addWidget(self._hint("En la app, pulsar un aviso selecciona el átomo afectado en el lienzo."))
        return self._scroll_page(inner)

    def _select_issue(self, it) -> None:
        for x in self._issues:
            x.setProperty("selected", x is it)
            st = x.style(); st.unpolish(x); st.polish(x)

    def _page_properties(self) -> QWidget:
        inner = QWidget(); inner.setAttribute(Qt.WidgetAttribute.WA_TranslucentBackground, True)
        v = QVBoxLayout(inner); v.setContentsMargins(0, 0, 0, 0); v.setSpacing(14)
        v.addWidget(self._sec_title("Molécula"))
        for k, val, mono in (
            ("Fórmula", "C7H7NO2", True), ("Masa molar", "137.14 g/mol", False),
            ("Carga neta", "0", False), ("Átomos / enlaces", "11 / 10", False),
            ("SMILES", "Cc1ccc(cc1)[N+](=O)[O-]", True),
        ):
            v.addWidget(make_row(k, val, mono=mono))
        v.addWidget(self._sec_title("Propiedades calculadas"))
        for k, val in (("LogP", "1.98"), ("Densidad", "1.12 g/cm³"),
                       ("Punto de fusión", "51 °C"), ("Momento dipolar", "4.10 D")):
            v.addWidget(make_row(k, val))
        v.addWidget(self._btn("Recalcular", primary=False))
        return self._scroll_page(inner)

    def _page_templates(self) -> QWidget:
        inner = QWidget(); inner.setAttribute(Qt.WidgetAttribute.WA_TranslucentBackground, True)
        v = QVBoxLayout(inner); v.setContentsMargins(0, 0, 0, 0); v.setSpacing(12)
        self.tpl_search = QLineEdit(); self.tpl_search.setObjectName("tplSearch")
        self.tpl_search.setPlaceholderText("Buscar plantilla…")
        self.tpl_search.textChanged.connect(self._filter_templates)
        v.addWidget(self.tpl_search)
        grid = QGridLayout(); grid.setSpacing(8)
        self._tpl_cards: list[tuple[QWidget, str]] = []
        for i, (icon, nm, meta) in enumerate(
            (("ring", "Benceno", "Anillo · C6"), ("ring", "Piridina", "Anillo · C5N"),
             ("tpl", "Glicina", "Aminoácido"), ("mol", "Nitro (–NO₂)", "Sustituyente"))):
            card = make_tpl_card(self.icons, self._tg, icon, nm, meta,
                                 on_click=lambda n=nm: self.toast.show(f"Plantilla «{n}» insertada (demo)"))
            grid.addWidget(card, i // 2, i % 2)
            self._tpl_cards.append((card, nm))
        wrap = QWidget(); wrap.setLayout(grid)
        wrap.setAttribute(Qt.WidgetAttribute.WA_TranslucentBackground, True)
        v.addWidget(wrap)
        v.addWidget(self._hint("Equivalente al dock «Plantillas» (TemplateBrowserService), sin reescribir su contenido."))
        return self._scroll_page(inner)

    def _filter_templates(self, q: str) -> None:
        qq = q.strip().lower()
        for card, nm in self._tpl_cards:
            card.setVisible(not qq or qq in nm.lower())

    def _page_appearance(self) -> QWidget:
        inner = QWidget(); inner.setAttribute(Qt.WidgetAttribute.WA_TranslucentBackground, True)
        v = QVBoxLayout(inner); v.setContentsMargins(0, 0, 0, 0); v.setSpacing(14)
        v.addWidget(self._sec_title("Tema"))
        for label, fn in (("Claro", lambda: self.set_theme("light")),
                          ("Oscuro", lambda: self.set_theme("dark"))):
            b = self._btn(label, primary=False)
            b.clicked.connect(fn)
            v.addWidget(b)
        v.addWidget(self._sec_title("Lienzo"))
        v.addWidget(self._btn("Alternar rejilla", primary=False))
        v.addWidget(self._hint("Hoy el tema es manual (claro/oscuro); la propuesta añade «seguir sistema» persistido vía platform.settings (M21)."))
        return self._scroll_page(inner)

    def _btn(self, text: str, primary: bool = False) -> QPushButton:
        b = QPushButton(text)
        b.setProperty("cls", "btnPrimary" if primary else "btn")
        return b

    def _hint(self, text: str) -> QLabel:
        l = QLabel(text); l.setObjectName("hint")
        l.setWordWrap(True)
        return l

    # ------------------------------------------------------------------
    # status bar
    # ------------------------------------------------------------------
    def _build_statusbar(self) -> QWidget:
        bar = QFrame(); bar.setObjectName("statusbar")
        bar.setFixedHeight(METRICS["statusH"])
        lay = QHBoxLayout(bar); lay.setContentsMargins(14, 0, 14, 0); lay.setSpacing(14)
        dot = QLabel(); dot.setObjectName("toolDot"); dot.setFixedSize(8, 8)
        dot.setStyleSheet("background: %s; border-radius: 4px;" % self.theme["accent"])
        self._status_dot = dot
        self.lbl_tool = QLabel("Seleccionar (V)"); self.lbl_tool.setObjectName("toolName")
        self.lbl_cursor = QLabel("— , —"); self.lbl_cursor.setObjectName("cursorPos")
        lay.addWidget(dot); lay.addWidget(self.lbl_tool); lay.addWidget(self.lbl_cursor)
        lay.addStretch(1)
        self.lbl_formula = QLabel("C₇H₇NO₂"); self.lbl_formula.setObjectName("stFormula")
        self.lbl_iupac = QLabel("1-metil-4-nitrobenzeno"); self.lbl_iupac.setObjectName("stIupac")
        self.lbl_charge = QLabel("Carga 0"); self.lbl_charge.setObjectName("stCharge")
        aw = QWidget(); ahl = QHBoxLayout(aw)
        ahl.setContentsMargins(0, 0, 0, 0); ahl.setSpacing(5)
        ai = QLabel(); ai.setPixmap(self.icons.pixmap("check", self.theme["ok"], 12))
        al = QLabel("autosave 14:32"); al.setObjectName("stAutosave")
        ahl.addWidget(ai); ahl.addWidget(al)
        self._autosave_ic = ai
        lay.addWidget(self.lbl_formula); lay.addWidget(self.lbl_iupac)
        lay.addWidget(self.lbl_charge); lay.addWidget(aw)
        return bar

    def _update_status(self) -> None:
        parts = [self.tool_name]
        if self.tool_sub:
            parts.append(self.tool_sub)
        txt = " · ".join(parts)
        if self.tool_key:
            txt += f" ({self.tool_key})"
        self.lbl_tool.setText(txt)

    # ------------------------------------------------------------------
    # tema
    # ------------------------------------------------------------------
    def set_theme(self, name: str) -> None:
        self._theme_name = name
        self.theme = Theme(name)
        app = QApplication.instance()
        if app is not None:
            self.theme.apply(app)
        self.canvas.apply_theme(self.theme)
        self._refresh_themed_widgets()
        self._update_status()

    def _toggle_theme(self) -> None:
        self.set_theme("dark" if self._theme_name == "light" else "light")

    def _refresh_themed_widgets(self) -> None:
        t = self.theme
        self._brand_ic.setPixmap(self.icons.pixmap("flask", t["accent"], 24))
        self.search.refresh_theme()
        # iconos de la app bar
        icon = "sun" if self._theme_name == "dark" else "moon"
        self.btn_theme.setIcon(self.icons.icon(icon, t["icon"], 18))
        for b, key in list(getattr(self, "_abar_icons", {}).items()):
            if b is self.btn_theme:
                continue
            color = t["iconHover"] if b.underMouse() else t["icon"]
            b.setIcon(self.icons.icon(key, color, 18))
        self.tab_new.setIcon(self.icons.icon("plus", t["text2"], 14))
        # pestañas
        for i in range(self.tabs.count()):
            self.tabs.setTabIcon(i, self.icons.icon("mol", t["text2"], 14))
        # rail + chips + canvas + panel
        for b in self._rail_buttons.values():
            b.refresh_theme = b._refresh_icon  # type: ignore[attr-defined]
            b._refresh_icon(); b.update()
        self.host.refresh_theme()
        self.side_tabs.refresh_theme()
        for it in getattr(self, "_issues", []):
            it.refresh_theme()
        for card, _nm in getattr(self, "_tpl_cards", []):
            card.refresh_theme()
        # status
        self._status_dot.setStyleSheet("background: %s; border-radius: 4px;" % t["accent"])
        self._autosave_ic.setPixmap(self.icons.pixmap("check", t["ok"], 12))
        # paleta + flyout + toast
        self.palette.refresh_theme()
        self.flyout.refresh_theme()
        self.toast.refresh_theme()

    # ------------------------------------------------------------------
    # acciones de la paleta
    # ------------------------------------------------------------------
    def _make_actions(self) -> list:
        A = PaletteAction
        return [
            A("Archivo", "Nuevo documento", "doc-new", "Ctrl N", lambda: self._add_tab(None, dirty=False)),
            A("Archivo", "Abrir…", "doc-open", "Ctrl O", lambda: self.toast.show("Abrir… (demo)")),
            A("Archivo", "Guardar", "doc-save", "Ctrl S", lambda: self.toast.show("Guardado (demo)")),
            A("Archivo", "Exportar como PNG…", "export", "", lambda: self.toast.show("Exportar PNG (demo)")),
            A("Edición", "Deshacer", "undo", "Ctrl Z", lambda: self.toast.show("Deshacer (demo)")),
            A("Edición", "Rehacer", "redo", "Ctrl Y", lambda: self.toast.show("Rehacer (demo)")),
            A("Edición", "Copiar", "copy", "Ctrl C", lambda: self.toast.show("Copiado (demo)")),
            A("Edición", "Pegar", "paste", "Ctrl V", lambda: self.toast.show("Pegado (demo)")),
            A("Herramientas", "Seleccionar", "pointer", "V", lambda: self._activate_tool("select")),
            A("Herramientas", "Enlace", "bond", "B", lambda: self._activate_tool("bond")),
            A("Herramientas", "Anillo", "ring", "R", lambda: self._activate_tool("ring")),
            A("Herramientas", "Átomo", "atom", "C", lambda: self._activate_tool("atom")),
            A("Herramientas", "Texto", "text", "T", lambda: self._activate_tool("text")),
            A("Herramientas", "Flechas", "arrow", "N", lambda: self._activate_tool("arrow")),
            A("Herramientas", "Corchetes", "bracket", "G", lambda: self._activate_tool("bracket")),
            A("Vista", "Cambiar tema claro/oscuro", "moon", "", self._toggle_theme),
            A("Vista", "Rejilla del lienzo", "grid", "", lambda: self.host.chip_grid.toggle()),
            A("Vista", "Numeración de átomos", "num", "", lambda: self.host.chip_nums.toggle()),
            A("Vista", "Ajustar zoom a la hoja", "fit", "", self._zoom_fit),
            A("Análisis", "Validar estructura", "shield", "",
              lambda: (self._open_side(1), self.toast.show("Validación (demo)"))),
            A("Análisis", "Limpiar 2D (Clean2D)", "clean", "", lambda: self.toast.show("Clean2D (demo)")),
            A("Análisis", "Propiedades químicas", "mol", "", lambda: self._open_side(2)),
            A("Análisis", "Espectroscopía", "spec", "", lambda: self.toast.show("Espectroscopía (demo)")),
            A("Análisis", "CompChem 3D", "cube", "", lambda: self.toast.show("CompChem 3D (demo)")),
            A("Plantillas", "Benceno", "ring", "", lambda: self.toast.show("Plantilla «Benceno» insertada (demo)")),
            A("Plantillas", "Piridina", "ring", "", lambda: self.toast.show("Plantilla «Piridina» insertada (demo)")),
            A("Plantillas", "Glicina", "tpl", "", lambda: self.toast.show("Plantilla «Glicina» insertada (demo)")),
        ]

    # ------------------------------------------------------------------
    # event filter: clic fuera cierra el flyout
    # ------------------------------------------------------------------
    def eventFilter(self, obj, event) -> bool:  # noqa: N802
        if (event.type() == QEvent.Type.MouseButtonPress and self.flyout.isVisible()
                and obj is not self.flyout and obj not in self.flyout.children()):
            if not self._is_in_rail_btn(event):
                self.flyout.hide()
                self._fly_open_group = None
        return False

    def _is_in_rail_btn(self, event) -> bool:
        w = event.widget() if hasattr(event, "widget") else getattr(event, "object", None)
        if w is None:
            return False
        for b in self._rail_buttons.values():
            if w is b or b in w.ancestors():
                return True
        return False

    # ------------------------------------------------------------------
    def refresh_theme_all(self) -> None:  # API para el smoke
        self._refresh_themed_widgets()


# ---------------------------------------------------------------------------
# smoke
# ---------------------------------------------------------------------------
def run_smoke(win: SpikeWindow, app: QApplication, out: Path) -> int:
    """Secuencia de interacción + capturas; devuelve nº de fallos."""
    out.mkdir(parents=True, exist_ok=True)
    fails: list[str] = []

    def check(name: str, cond: bool) -> None:
        print(("  ✓ " if cond else "  ✗ FALLO: ") + name)
        if not cond:
            fails.append(name)

    def shot(name: str) -> None:
        app.processEvents()
        pm = win.grab()
        pm.save(str(out / name))

    win.resize(1440, 900)
    win.show(); app.processEvents()
    shot("01_light_1440.png")

    # estado inicial
    check("tema claro", win._theme_name == "light")
    check("2 pestañas de documento", win.tabs.count() == 2)
    check("selector activo", win._rail_buttons["select"].property("active"))
    check("status tool", win.lbl_tool.text() == "Seleccionar (V)")

    # tema oscuro
    win._toggle_theme(); app.processEvents()
    check("tema oscuro", win._theme_name == "dark")
    shot("02_dark_1440.png")
    win._toggle_theme(); app.processEvents()

    # rail + flyout
    win._rail_buttons["bond"].click(); app.processEvents()
    check("flyout enlace visible", win.flyout.isVisible())
    check("flyout 11 ítems", len(win.flyout._items) == 11)
    shot("03_flyout_bond.png")
    win.flyout._items[1].click()  # Doble
    app.processEvents()
    check("flyout cerrado tras elegir", not win.flyout.isVisible())
    check("tool = enlace·doble", win.tool_group == "bond" and win.tool_sub == "Doble")
    check("status", win.lbl_tool.text() == "Enlace · Doble (B)")

    win._rail_buttons["atom"].click(); app.processEvents()
    check("flyout átomos 10 ítems", len(win.flyout._items) == 10)
    shot("04_flyout_atom.png")
    app.processEvents()
    win.flyout.keyPressEvent(_key(Qt.Key.Key_Escape))

    # tabs del panel derecho
    win._open_side(1); app.processEvents()
    check("tab Validación activo", win.side_stack.currentIndex() == 1)
    shot("05_validation.png")
    win._open_side(3); app.processEvents()
    check("tab Plantillas activo", win.side_stack.currentIndex() == 3)
    win.tpl_search.setText("piri")
    visible = sum(1 for c, _ in win._tpl_cards if c.isVisible())
    check("filtro plantillas", visible == 1)

    # zoom / rejilla / números
    win._zoom_delta(1); win._zoom_delta(1); app.processEvents()
    check("zoom 135 %", win.host.val.text() == "135 %")
    win._zoom_fit(); app.processEvents()
    check("fit 100 %", win.host.val.text() == "100 %")
    win.host.chip_grid.toggle(); app.processEvents()
    check("rejilla off", not win.canvas.grid_on)
    win.host.chip_nums.toggle(); app.processEvents()
    check("números on", win.canvas.num_items[0].isVisible())
    shot("06_zoom_grid.png")

    # pestañas de documento
    n = win.tabs.count()
    win._add_tab(None, dirty=False); app.processEvents()
    check("nueva pestaña", win.tabs.count() == n + 1)
    win._close_tab(0); app.processEvents()
    check("cerrar pestaña", win.tabs.count() == n)

    # paleta Ctrl+K
    win.palette.show_overlay(); app.processEvents()
    check("paleta visible", win.palette.is_visible())
    check("27 acciones", len(win.palette.actions) == 27)
    shot("07_palette.png")
    win.palette.input.setText("benc")
    app.processEvents()
    check("filtro benc", len(win.palette._filtered) == 1)
    win.palette.input.setText("")
    app.processEvents()
    win.palette.keyPressEvent(_key(Qt.Key.Key_Escape))
    check("Esc cierra", not win.palette.is_visible())

    # ventana pequeña
    win.resize(980, 600); app.processEvents()
    shot("08_light_980.png")
    win.set_theme("dark"); app.processEvents()
    shot("09_dark_980.png")

    print(f"\nSMOKE: {'OK' if not fails else 'CON FALLOS'} — {len(fails)} fallos")
    return len(fails)


def _key(key: int):
    from PyQt6.QtCore import QEvent, Qt as _Qt
    from PyQt6.QtGui import QKeyEvent
    return QKeyEvent(QEvent.Type.KeyPress, key, _Qt.KeyboardModifier.NoModifier)


def main(argv: list[str]) -> int:
    import argparse
    ap = argparse.ArgumentParser(description="Spike PyQt6 de la propuesta de UI")
    ap.add_argument("--theme", choices=["light", "dark"], default="light")
    ap.add_argument("--size", default="1440x900")
    ap.add_argument("--smoke", metavar="DIR", help="ejecuta la secuencia de verificación y guarda capturas en DIR")
    args = ap.parse_args(argv)

    app = QApplication([])
    win = SpikeWindow(app)
    w, h = (int(x) for x in args.size.lower().split("x"))
    win.resize(w, h)

    if args.smoke:
        win.set_theme(args.theme)
        rc = run_smoke(win, app, Path(args.smoke))
        return rc

    win.set_theme(args.theme)
    win.show()
    return app.exec()


if __name__ == "__main__":
    raise SystemExit(main(sys.argv[1:]))
