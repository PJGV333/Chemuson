"""Design tokens y generador de QSS para el spike de modernización de UI.

Traducción a Qt de los design tokens del mockup (docs/ui-modernization/mockup-ui.html,
:root / html[data-theme="dark"]) y de PLAN.md §2.2. NO es una copia literal del CSS:
es la misma tabla de tokens expresada como QSS generado a partir de constantes.

Uso:
    from theme import Theme
    theme = Theme("light")
    theme.apply(app)          # setStyleSheet + fuente + paleta base
"""
from __future__ import annotations

from typing import Callable

from PyQt6.QtCore import Qt
from PyQt6.QtGui import QColor, QFont, QPalette

# callable que devuelve el Theme *actual* (para widgets pintados a mano)
ThemeGetter = Callable[[], "Theme"]


# ---------------------------------------------------------------------------
# Tokens (misma tabla que el mockup; claridad es el objetivo, no el tamaño)
# ---------------------------------------------------------------------------
TOKENS = {
    "light": {
        "bg":          "#F1F5F9",
        "surface":     "#FFFFFF",
        "surface2":    "#F8FAFC",
        "surface3":    "#EEF2F7",
        "border":      "#E2E8F0",
        "borderStrong":"#CBD5E1",
        "text1":       "#0F172A",
        "text2":       "#475569",
        "text3":       "#94A3B8",
        "accent":      "#0E7490",
        "accentStrong":"#155E75",
        "accentHover": "#0891B2",
        "accentSoft":  "#ECFEFF",
        "accentBorder":"#A5F3FC",
        "onAccent":    "#FFFFFF",
        "danger":      "#DC2626",
        "dangerSoft":  "#FEF2F2",
        "warn":        "#D97706",
        "warnSoft":    "#FFFBEB",
        "ok":          "#059669",
        "okSoft":      "#ECFDF5",
        "sheet":       "#FFFFFF",
        "sheetGrid":   "#E9EEF4",
        "canvasBg":    "#E9EEF4",
        # sombras (para QGraphicsDropShadowEffect; QSS no tiene box-shadow)
        "shadow1":     QColor(15, 23, 42, 25),
        "shadow2":     QColor(15, 23, 42, 38),
        "sheetShadow": QColor(15, 23, 42, 31),
        # color de tinte de iconos por estado
        "icon":        "#475569",
        "iconHover":   "#0F172A",
        "iconActive":  "#155E75",
    },
    "dark": {
        "bg":          "#0B1120",
        "surface":     "#0F172A",
        "surface2":    "#16213A",
        "surface3":    "#1E293B",
        "border":      "#263349",
        "borderStrong":"#334155",
        "text1":       "#F1F5F9",
        "text2":       "#C3CEDF",
        "text3":       "#64748B",
        "accent":      "#22D3EE",
        "accentStrong":"#67E8F9",
        "accentHover": "#06B6D4",
        "accentSoft":  "rgba(34,211,238,31)",
        "accentBorder":"rgba(34,211,238,115)",
        "onAccent":    "#083344",
        "danger":      "#F87171",
        "dangerSoft":  "rgba(248,113,113,33)",
        "warn":        "#FBBF24",
        "warnSoft":    "rgba(251,191,36,33)",
        "ok":          "#34D399",
        "okSoft":      "rgba(52,211,153,33)",
        "sheet":       "#FFFFFF",
        "sheetGrid":   "#E9EEF4",
        "canvasBg":    "#0A0F1A",
        "shadow1":     QColor(0, 0, 0, 115),
        "shadow2":     QColor(0, 0, 0, 140),
        "sheetShadow": QColor(0, 0, 0, 140),
        "icon":        "#C3CEDF",
        "iconHover":   "#F1F5F9",
        "iconActive":  "#67E8F9",
    },
}

# Métricas del mockup (px) — equivalen a la grilla de 8 px de PLAN.md §2.2
METRICS = {
    "appbarH": 54,
    "statusH": 34,
    "railW": 58,
    "railBtn": 42,
    "railIcon": 21,
    "sideW": 324,
    "radiusSurface": 10,
    "radiusBtn": 9,
    "radiusRail": 11,
    "radiusChip": 8,
    "fontBase": 13,
    "fontSmall": 12,
    "fontTiny": 11,
}


def pick_font_family(candidates) -> str:
    from PyQt6.QtGui import QFontDatabase
    avail = set(QFontDatabase.families())
    for fam in candidates:
        if fam in avail:
            return fam
    return QFontDatabase().families()[0]


_FONT_CANDIDATES = ("Inter", "Segoe UI", "SF Pro Text", "Roboto", "DejaVu Sans")


QSS_TEMPLATE = """
/* ===== base ===== */
QWidget { color: @@text1@@; }
QLabel { background: transparent; }

#rootWindow { background: @@surface@@; }

/* ===== barra de aplicación ===== */
#appbar { background: @@surface@@; border-bottom: 1px solid @@border@@; }
#brandName { font-size: 14px; font-weight: 700; }
#brandIcon { }
#verPill {
  color: @@text3@@; font-size: 10px; font-weight: 600;
  border: 1px solid @@border@@; border-radius: 5px;
  padding: 1px 5px;
}

/* pestañas de documento */
QTabBar#docTabs { background: transparent; qproperty-drawBase: 0; }
QTabBar#docTabs::tab {
  background: transparent;
  border: 1px solid transparent;
  border-radius: 8px;
  padding: 5px 8px 5px 10px;
  margin-right: 2px;
  color: @@text2@@;
  font-size: 12px;
  min-width: 86px;
  max-width: 190px;
}
QTabBar#docTabs::tab:hover { background: @@surface2@@; color: @@text1@@; }
QTabBar#docTabs::tab:selected {
  background: @@surface3@@;
  color: @@text1@@;
  font-weight: 600;
  border-bottom: 2px solid @@accent@@;
}
QTabBar#docTabs::tab:selected:hover { background: @@surface3@@; }
#tabNew {
  background: transparent; border: 1px dashed @@borderStrong@@;
  border-radius: 8px; color: @@text2@@;
}
#tabNew:hover { color: @@accentStrong@@; border: 1px dashed @@accentBorder@@; background: @@accentSoft@@; }
QToolButton[cls="tabClose"] { background: transparent; border: none; border-radius: 5px; }
QToolButton[cls="tabClose"]:hover { background: @@border@@; }

/* píldora de búsqueda */
#searchPill {
  background: @@surface2@@; border: 1px solid @@border@@;
  border-radius: 9px; padding: 6px 10px;
}
#searchPill:hover { border: 1px solid @@borderStrong@@; }
#searchPillTxt { color: @@text3@@; font-size: 12px; }
#searchPill:hover #searchPillTxt { color: @@text2@@; }
#kbdK {
  color: @@text2@@; font-size: 10px; font-weight: 700;
  background: @@surface@@; border: 1px solid @@borderStrong@@;
  border-bottom: 2px solid @@borderStrong@@; border-radius: 5px;
  padding: 1px 5px;
}

/* botones de la app bar (32 px) */
QToolButton[cls="abar"] {
  background: transparent; border: none; border-radius: 9px;
  color: @@text2@@;
}
QToolButton[cls="abar"]:hover { background: @@surface2@@; color: @@text1@@; }
QToolButton[cls="abar"]:pressed { background: @@surface3@@; }
QToolButton[cls="abar"]:disabled { opacity: 0.38; }
#abarSep { background: @@border@@; max-width: 1px; min-width: 1px; }

/* ===== rail ===== */
#rail { background: @@surface@@; border-right: 1px solid @@border@@; }
#railSep { background: @@border@@; }
QToolButton[cls="rail"] {
  background: transparent; border: 1px solid transparent;
  border-radius: @@radiusRail@@px; color: @@text2@@;
}
QToolButton[cls="rail"]:hover { background: @@surface2@@; color: @@text1@@; }
QToolButton[cls="rail"]:pressed { background: @@surface3@@; }
QToolButton[cls="rail"][active="true"] {
  background: @@accentSoft@@;
  border: 1.5px solid @@accentBorder@@;
  color: @@accentStrong@@;
}

/* ===== flyout ===== */
QFrame#flyout, QFrame#paletteCard {
  background: @@surface@@; border: 1px solid @@border@@;
}
QFrame#flyout { border-radius: 13px; }
#flyoutTitle { color: @@text3@@; font-size: 10px; font-weight: 700; }
QToolButton[cls="flyItem"] {
  background: transparent; border: 1px solid transparent;
  border-radius: 9px; padding: 7px 2px 6px;
  color: @@text2@@; font-size: 10px;
}
QToolButton[cls="flyItem"]:hover { background: @@surface2@@; color: @@text1@@; }
QToolButton[cls="flyItem"][active="true"] {
  background: @@accentSoft@@; border: 1px solid @@accentBorder@@;
  color: @@accentStrong@@; font-weight: 600;
}
#flyoutFootTxt { color: @@text3@@; font-size: 10px; }
QToolButton[cls="flyFoot"] {
  background: transparent; border: none;
  color: @@accentStrong@@; font-size: 10px; font-weight: 600;
}
QToolButton[cls="flyFoot"]:hover { text-decoration: underline; }
#flyoutSep { background: @@border@@; }

/* ===== canvas ===== */
QGraphicsView#canvas { border: none; background: transparent; }
QGraphicsView#canvas { margin: 0; }

/* chips Rejilla / Números */
QToolButton[cls="chip"] {
  background: @@surface@@; border: 1px solid @@border@@;
  border-radius: 8px; padding: 5px 10px;
  color: @@text2@@; font-size: 11px; font-weight: 600;
}
QToolButton[cls="chip"]:hover { border: 1px solid @@borderStrong@@; color: @@text1@@; }
QToolButton[cls="chip"]:checked {
  background: @@accentSoft@@; border: 1px solid @@accentBorder@@;
  color: @@accentStrong@@;
}

/* pill de zoom */
QFrame#zoomPill {
  background: @@surface@@; border: 1px solid @@border@@;
  border-radius: 9px;
}
QToolButton[cls="zoomBtn"] {
  background: transparent; border: none; color: @@text2@@;
}
QToolButton[cls="zoomBtn"]:hover { color: @@text1@@; }
#zoomVal { color: @@text2@@; font-size: 11px; font-weight: 600; min-width: 40px; }

/* ===== panel derecho ===== */
#sideWrap { background: @@surface@@; border-left: 1px solid @@border@@; }
QToolButton[cls="sideTab"] {
  background: transparent; border: none;
  padding: 9px 10px 11px; color: @@text2@@; font-size: 12px; font-weight: 500;
}
QToolButton[cls="sideTab"]:hover { color: @@text1@@; }
QToolButton[cls="sideTab"][active="true"] { color: @@text1@@; font-weight: 600; }
#sideTabSep { background: @@border@@; }
#sideBody { background: @@surface@@; }
QScrollArea#sideScroll { border: none; background: transparent; }
QScrollArea#sideScroll > QWidget > QWidget { background: transparent; }

/* contenido de paneles */
#secTitle { color: @@text3@@; font-size: 10px; font-weight: 700; }
QFrame[cls="kv"] {
  background: @@surface2@@; border: 1px solid @@border@@;
  border-radius: 9px;
}
#kvK { color: @@text3@@; font-size: 10px; font-weight: 600; }
#kvV { color: @@text1@@; font-size: 13px; font-weight: 700; }
#kvSmall { color: @@text3@@; font-size: 10px; font-weight: 600; }
QFrame[cls="row"] { border-bottom: 1px dashed @@border@@; }
#rowK { color: @@text2@@; font-size: 12px; }
#rowV { color: @@text1@@; font-size: 12px; font-weight: 600; }
#rowV.mono { font-family: monospace; font-size: 11px; font-weight: 500; }

QPushButton[cls="btn"] {
  background: @@surface@@; color: @@text1@@;
  border: 1px solid @@borderStrong@@; border-radius: 8px;
  padding: 7px 12px; font-size: 12px; font-weight: 600;
}
QPushButton[cls="btn"]:hover { background: @@surface2@@; }
QPushButton[cls="btn"]:pressed { background: @@surface3@@; }
QPushButton[cls="btn"]:disabled { opacity: 0.45; }
QPushButton[cls="btnPrimary"] {
  background: @@accent@@; color: @@onAccent@@;
  border: 1px solid @@accent@@; border-radius: 8px;
  padding: 7px 12px; font-size: 12px; font-weight: 600;
}
QPushButton[cls="btnPrimary"]:hover { background: @@accentHover@@; border: 1px solid @@accentHover@@; }
QPushButton[cls="btnPrimary"]:pressed { background: @@accentHover@@; }

QFrame[cls="pill"] {
  border-radius: 999px; padding: 2px 9px;
  font-size: 10px; font-weight: 600;
  background: @@surface2@@; color: @@text2@@;
}
QFrame[cls="pill"][sev="warn"] { background: @@warnSoft@@; color: @@warn@@; }
QFrame[cls="pill"][sev="err"] { background: @@dangerSoft@@; color: @@danger@@; }
QFrame[cls="pill"][sev="ok"] { background: @@okSoft@@; color: @@ok@@; }
QFrame[cls="pill"][sev="accent"] { background: @@accentSoft@@; color: @@accentStrong@@; border: 1px solid @@accentBorder@@; }

QFrame[cls="issue"] {
  background: @@surface@@; border: 1px solid @@border@@;
  border-radius: 10px;
}
QFrame[cls="issue"]:hover { background: @@surface2@@; }
QFrame[cls="issue"][selected="true"] {
  background: @@accentSoft@@; border: 1px solid @@accentBorder@@;
}
#issueTtl { font-size: 12px; font-weight: 600; }
#issueSub { color: @@text3@@; font-size: 10px; }

QProgressBar {
  background: @@surface3@@; border: none; border-radius: 3px;
  min-height: 5px; max-height: 5px; text-align: center;
}
QProgressBar::chunk { background: @@accent@@; border-radius: 3px; }

QLineEdit#tplSearch {
  background: @@surface2@@; border: 1px solid @@border@@;
  border-radius: 9px; padding: 7px 10px;
  color: @@text1@@; font-size: 12px;
}
QLineEdit#tplSearch:focus { border: 1px solid @@accentBorder@@; }

QFrame[cls="tpl"] {
  background: @@surface@@; border: 1px solid @@border@@;
  border-radius: 11px;
}
QFrame[cls="tpl"]:hover { background: @@accentSoft@@; border: 1px solid @@accentBorder@@; }
#tplNm { font-size: 12px; font-weight: 600; }
#tplMeta { color: @@text3@@; font-size: 10px; }

#hint { color: @@text3@@; font-size: 10px; }

/* ===== barra de estado ===== */
#statusbar { background: @@surface@@; border-top: 1px solid @@border@@; }
#toolDot { background: @@accent@@; border-radius: 4px; }
#toolName { font-size: 12px; font-weight: 600; }
#cursorPos { color: @@text3@@; font-size: 11px; }
#stFormula { font-size: 12px; font-weight: 700; }
#stIupac { color: @@text3@@; font-size: 11px; }
#stCharge { color: @@text2@@; font-size: 11px; }
#stAutosave { color: @@ok@@; font-size: 11px; font-weight: 600; }

/* ===== paleta de comandos ===== */
#paletteOverlay { background: rgba(15, 23, 42, 107); }
QFrame#paletteCard { border-radius: 14px; }
#palInputRow { border-bottom: 1px solid @@border@@; }
QLineEdit#palInput {
  border: none; background: transparent;
  color: @@text1@@; font-size: 14px;
}
QScrollArea#palScroll { border: none; background: transparent; }
QScrollArea#palScroll > QWidget > QWidget { background: transparent; }
#palSec { color: @@text3@@; font-size: 10px; font-weight: 700; }
QFrame[cls="palItem"] { border-radius: 9px; }
QFrame[cls="palItem"]:hover { background: @@surface2@@; }
QFrame[cls="palItem"][selected="true"] { background: @@accentSoft@@; }
#palTitle { font-size: 13px; }
QFrame[cls="palItem"][selected="true"] #palTitle { color: @@accentStrong@@; font-weight: 600; }
#palKbd { color: @@text3@@; font-size: 10px; font-weight: 600; }
#palEmpty { color: @@text3@@; font-size: 12px; }

/* ===== scrollbars finos ===== */
QScrollArea, QAbstractScrollArea { background: transparent; }
QScrollBar:vertical { background: transparent; width: 8px; margin: 0; }
QScrollBar::handle:vertical {
  background: @@borderStrong@@; border-radius: 4px; min-height: 24px;
}
QScrollBar::handle:vertical:hover { background: @@text3@@; }
QScrollBar::add-line:vertical, QScrollBar::sub-line:vertical { height: 0; }
QScrollBar::add-page:vertical, QScrollBar::sub-page:vertical { background: transparent; }
QScrollBar:horizontal { background: transparent; height: 8px; margin: 0; }
QScrollBar::handle:horizontal {
  background: @@borderStrong@@; border-radius: 4px; min-width: 24px;
}
QScrollBar::add-line:horizontal, QScrollBar::sub-line:horizontal { width: 0; }
QScrollBar::add-page:horizontal, QScrollBar::sub-page:horizontal { background: transparent; }

QToolTip {
  background: @@surface@@; color: @@text1@@;
  border: 1px solid @@border@@; border-radius: 6px; padding: 4px 7px;
  font-size: 11px;
}
"""


class Theme:
    """Un tema resuelto: tokens + QSS generado."""

    def __init__(self, name: str = "light"):
        if name not in TOKENS:
            raise ValueError(f"tema desconocido: {name!r}")
        self.name = name
        self.tokens = TOKENS[name]

    def __getitem__(self, key):
        return self.tokens[key]

    def color(self, key: str) -> QColor:
        return QColor(self.tokens[key])

    def build_qss(self) -> str:
        qss = QSS_TEMPLATE
        for k, v in self.tokens.items():
            qss = qss.replace(f"@@{k}@@", str(v))
        qss = qss.replace("@@radiusRail@@", str(METRICS["radiusRail"]))
        # propiedades numéricas de métricas que usan placeholder
        return qss

    def font(self, size_px: int | None = None, bold: bool = False) -> QFont:
        f = QFont()
        f.setFamily(pick_font_family(_FONT_CANDIDATES))
        f.setPixelSize(METRICS["fontBase"] if size_px is None else size_px)
        f.setBold(bold)
        return f

    def apply(self, app) -> None:
        """Aplica QSS + fuente base + paleta (base/texto) a toda la app."""
        app.setFont(self.font())
        app.setStyleSheet(self.build_qss())
        pal = QPalette()
        t = self.tokens
        pal.setColor(QPalette.ColorRole.Window, QColor(t["surface"]))
        pal.setColor(QPalette.ColorRole.WindowText, QColor(t["text1"]))
        pal.setColor(QPalette.ColorRole.Base, QColor(t["surface2"]))
        pal.setColor(QPalette.ColorRole.AlternateBase, QColor(t["surface2"]))
        pal.setColor(QPalette.ColorRole.Text, QColor(t["text1"]))
        pal.setColor(QPalette.ColorRole.Button, QColor(t["surface"]))
        pal.setColor(QPalette.ColorRole.ButtonText, QColor(t["text1"]))
        pal.setColor(QPalette.ColorRole.PlaceholderText, QColor(t["text3"]))
        pal.setColor(QPalette.ColorRole.Highlight, QColor(t["accent"]))
        pal.setColor(QPalette.ColorRole.HighlightedText, QColor(t["onAccent"]))
        app.setPalette(pal)
