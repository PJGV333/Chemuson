"""Generadores de hojas de estilo Qt a partir de los design tokens.

Fase 1 del plan de modernización de la UI. Sustituye la generación QSS de
``chemuson.gui.styles`` (que pasa a ser fachada de compatibilidad): las
mismas dos hojas públicas —``get_main_stylesheet`` y
``get_tool_palette_stylesheet``— se construyen exclusivamente desde la
tabla de tokens del tema.

Solo se usan propiedades QSS nativas de Qt (sin ``box-shadow``,
``transitions`` ni ``opacity`` QSS): los estados ``disabled`` se pintan con
colores explícitos de tokens.
"""

from __future__ import annotations

from chemuson.gui.theme.tokens import METRICS, get_tokens

__all__ = [
    "get_main_stylesheet",
    "get_tool_palette_stylesheet",
]

# Números de métrica interpolados en el QSS (px).
_RADIUS_SURFACE = METRICS["radiusSurface"]
_RADIUS_BTN = METRICS["radiusBtn"]
_RADIUS_CHIP = METRICS["radiusChip"]
_FONT_BASE = METRICS["fontBase"]
_FONT_SMALL = METRICS["fontSmall"]
_FONT_TINY = METRICS["fontTiny"]


def get_main_stylesheet(theme_name: str) -> str:
    """Genera la hoja de estilo principal según el tema (100 % tokens)."""
    t = get_tokens(theme_name)

    def c(key: str) -> str:
        return str(t[key])

    return f"""
/* =========================================================== */
/* Main Window (tokens: tema {theme_name})                     */
/* =========================================================== */
QMainWindow {{
    background-color: {c('bg')};
}}

QWidget {{
    color: {c('text1')};
}}

QLabel {{
    color: {c('text1')};
    background: transparent;
}}

/* =========================================================== */
/* Menu Bar (dirección del spike: superficie + borde inferior) */
/* =========================================================== */
QMenuBar {{
    background-color: {c('surface')};
    color: {c('text1')};
    border: none;
    border-bottom: 1px solid {c('border')};
    padding: 4px 8px;
    spacing: 4px;
    font-size: {_FONT_BASE}px;
    min-height: 34px;
}}

QMenuBar::item {{
    background: transparent;
    color: {c('text2')};
    padding: 5px 12px;
    border-radius: {_RADIUS_CHIP}px;
    margin: 3px 2px;
}}

QMenuBar::item:hover {{
    background-color: {c('surface2')};
    color: {c('text1')};
}}

QMenuBar::item:selected {{
    background-color: {c('accentSoft')};
    color: {c('accentStrong')};
}}

QMenuBar::item:pressed {{
    background-color: {c('surface3')};
}}

/* =========================================================== */
/* Dropdown Menus                                              */
/* =========================================================== */
QMenu {{
    background-color: {c('surface')};
    color: {c('text1')};
    border: 1px solid {c('border')};
    border-radius: {_RADIUS_SURFACE}px;
    padding: 8px 4px;
    margin: 4px;
}}

QMenu::item {{
    padding: 7px 28px 7px 12px;
    border-radius: {_RADIUS_CHIP}px;
    margin: 2px 4px;
    color: {c('text1')};
}}

QMenu::item:hover,
QMenu::item:selected {{
    background-color: {c('accentSoft')};
    color: {c('accentStrong')};
}}

QMenu::item:disabled {{
    color: {c('text3')};
}}

QMenu::separator {{
    height: 1px;
    background-color: {c('border')};
    margin: 6px 12px;
}}

QMenu::icon {{
    margin-left: 8px;
}}

QMenu::indicator {{
    width: 18px;
    height: 18px;
    margin-left: 8px;
}}

/* =========================================================== */
/* Toolbars (superficie + borde inferior, como el spike)       */
/* =========================================================== */
QToolBar {{
    background-color: {c('surface')};
    color: {c('text1')};
    border: none;
    border-bottom: 1px solid {c('border')};
    spacing: 6px;
    padding: 8px 10px;
}}

QToolBar::separator {{
    width: 1px;
    background-color: {c('border')};
    margin: 8px 10px;
}}

QToolBar::handle {{
    background-color: transparent;
}}

/* =========================================================== */
/* Tool Buttons                                                */
/* =========================================================== */
QToolButton {{
    background-color: transparent;
    border: 1px solid transparent;
    border-radius: {_RADIUS_CHIP}px;
    padding: 8px;
    color: {c('text2')};
}}

QToolButton:hover {{
    background-color: {c('surface2')};
    color: {c('text1')};
}}

QToolButton:pressed {{
    background-color: {c('surface3')};
}}

QToolButton:checked {{
    background-color: {c('accentSoft')};
    border: 1px solid {c('accentBorder')};
    color: {c('accentStrong')};
}}

QToolButton:disabled {{
    color: {c('text3')};
    background-color: transparent;
    border: 1px solid transparent;
}}

QToolButton::menu-indicator {{
    image: none;
    subcontrol-position: right bottom;
    subcontrol-origin: padding;
    width: 8px;
    height: 8px;
}}

/* =========================================================== */
/* Dock Widgets                                                */
/* =========================================================== */
QDockWidget {{
    titlebar-close-icon: none;
    titlebar-normal-icon: none;
    font-weight: 600;
    color: {c('text1')};
}}

QDockWidget::title {{
    background-color: {c('surface2')};
    padding: 12px 14px;
    border-bottom: 1px solid {c('border')};
    text-align: left;
    font-size: {_FONT_BASE}px;
}}

QDockWidget::close-button,
QDockWidget::float-button {{
    border: none;
    background: transparent;
    padding: 4px;
    border-radius: 4px;
}}

QDockWidget::close-button:hover,
QDockWidget::float-button:hover {{
    background-color: {c('border')};
    border-radius: 6px;
}}

/* =========================================================== */
/* Status Bar (dirección del spike: superficie + borde)        */
/* =========================================================== */
QStatusBar {{
    background-color: {c('surface')};
    color: {c('text2')};
    border-top: 1px solid {c('border')};
    padding: 6px 14px;
    font-size: {_FONT_SMALL}px;
    min-height: 28px;
}}

QStatusBar::item {{
    border: none;
}}

QStatusBar QLabel {{
    color: {c('text2')};
    padding: 0 4px;
}}

/* =========================================================== */
/* Scrollbars (finos, como el spike)                           */
/* =========================================================== */
QScrollBar:vertical {{
    background-color: transparent;
    width: 8px;
    margin: 0;
}}

QScrollBar::handle:vertical {{
    background-color: {c('borderStrong')};
    border-radius: 4px;
    min-height: 24px;
}}

QScrollBar::handle:vertical:hover {{
    background-color: {c('text3')};
}}

QScrollBar::add-line:vertical,
QScrollBar::sub-line:vertical {{
    height: 0px;
    background: transparent;
}}

QScrollBar::add-page:vertical,
QScrollBar::sub-page:vertical {{
    background: transparent;
}}

QScrollBar:horizontal {{
    background-color: transparent;
    height: 8px;
    margin: 0;
}}

QScrollBar::handle:horizontal {{
    background-color: {c('borderStrong')};
    border-radius: 4px;
    min-width: 24px;
}}

QScrollBar::handle:horizontal:hover {{
    background-color: {c('text3')};
}}

QScrollBar::add-line:horizontal,
QScrollBar::sub-line:horizontal {{
    width: 0px;
    background: transparent;
}}

QScrollBar::add-page:horizontal,
QScrollBar::sub-page:horizontal {{
    background: transparent;
}}

/* =========================================================== */
/* Tables (Inspector Dock, etc.)                               */
/* =========================================================== */
QTableView,
QTableWidget {{
    background-color: {c('surface')};
    alternate-background-color: {c('surface2')};
    border: 1px solid {c('border')};
    border-radius: {_RADIUS_CHIP}px;
    gridline-color: {c('border')};
}}

QTableView::item,
QTableWidget::item {{
    padding: 6px 10px;
    color: {c('text1')};
}}

QTableView::item:selected,
QTableWidget::item:selected {{
    background-color: {c('accentSoft')};
    color: {c('text1')};
}}

QTableView::item:hover,
QTableWidget::item:hover {{
    background-color: {c('surface2')};
}}

QHeaderView::section {{
    background-color: {c('surface2')};
    color: {c('text2')};
    padding: 8px 10px;
    border: none;
    border-bottom: 2px solid {c('borderStrong')};
    font-weight: 600;
    font-size: {_FONT_SMALL}px;
}}

/* =========================================================== */
/* Tree Widgets (Plantillas, etc.)                             */
/* =========================================================== */
QTreeView,
QTreeWidget {{
    background-color: {c('surface')};
    border: 1px solid {c('border')};
    border-radius: {_RADIUS_CHIP}px;
}}

QTreeView::item,
QTreeWidget::item {{
    padding: 8px 6px;
    border-radius: 6px;
    margin: 2px 4px;
    color: {c('text1')};
}}

QTreeView::item:hover,
QTreeWidget::item:hover {{
    background-color: {c('surface2')};
}}

QTreeView::item:selected,
QTreeWidget::item:selected {{
    background-color: {c('accentSoft')};
    color: {c('text1')};
}}

QTreeWidget::branch:has-children:!has-siblings:closed,
QTreeWidget::branch:closed:has-children:has-siblings {{
    border-image: none;
}}

QTreeWidget::branch:open:has-children:!has-siblings,
QTreeWidget::branch:open:has-children:has-siblings {{
    border-image: none;
}}

/* =========================================================== */
/* Item views genéricos (listas de diálogos)                   */
/* =========================================================== */
QAbstractItemView {{
    background-color: {c('surface')};
    alternate-background-color: {c('surface2')};
    color: {c('text1')};
    border: 1px solid {c('border')};
    border-radius: {_RADIUS_CHIP}px;
    selection-background-color: {c('accentSoft')};
    selection-color: {c('text1')};
}}

QAbstractItemView::item:hover {{
    background-color: {c('surface2')};
}}

QAbstractItemView::item:selected {{
    background-color: {c('accentSoft')};
    color: {c('text1')};
}}

QAbstractScrollArea {{
    background-color: {c('surface')};
}}

/* =========================================================== */
/* Dialogs                                                     */
/* =========================================================== */
QDialog {{
    background-color: {c('surface')};
}}

QDialog QLabel {{
    color: {c('text1')};
}}

/* =========================================================== */
/* Push Buttons                                                */
/* =========================================================== */
QPushButton {{
    background-color: {c('accent')};
    color: {c('onAccent')};
    border: none;
    border-radius: {_RADIUS_BTN}px;
    padding: 10px 22px;
    font-weight: 600;
    min-width: 80px;
}}

QPushButton:hover {{
    background-color: {c('accentHover')};
}}

QPushButton:pressed {{
    background-color: {c('accentStrong')};
}}

QPushButton:disabled {{
    background-color: {c('borderStrong')};
    color: {c('text3')};
}}

/* Botón secundario (QSS por propiedad, como la versión anterior) */
QPushButton[flat="true"] {{
    background-color: transparent;
    color: {c('accentStrong')};
    border: 1px solid {c('borderStrong')};
}}

QPushButton[flat="true"]:hover {{
    background-color: {c('accentSoft')};
    border: 1px solid {c('accentBorder')};
}}

QPushButton[flat="true"]:disabled {{
    color: {c('text3')};
    border: 1px solid {c('border')};
}}

/* =========================================================== */
/* Line Edits                                                  */
/* =========================================================== */
QLineEdit {{
    background-color: {c('surface')};
    border: 1px solid {c('borderStrong')};
    border-radius: {_RADIUS_CHIP}px;
    padding: 8px 12px;
    color: {c('text1')};
    selection-background-color: {c('accentSoft')};
    selection-color: {c('text1')};
}}

QLineEdit:focus {{
    border: 1px solid {c('accent')};
    background-color: {c('surface')};
}}

QLineEdit:disabled {{
    background-color: {c('surface2')};
    color: {c('text3')};
    border-color: {c('border')};
}}

/* =========================================================== */
/* Combo Boxes                                                 */
/* =========================================================== */
QComboBox {{
    background-color: {c('surface')};
    border: 1px solid {c('borderStrong')};
    border-radius: {_RADIUS_CHIP}px;
    padding: 8px 12px;
    color: {c('text1')};
    min-width: 100px;
}}

QComboBox:hover {{
    border-color: {c('accent')};
}}

QComboBox:focus {{
    border: 1px solid {c('accent')};
}}

QComboBox::drop-down {{
    border: none;
    width: 24px;
}}

QComboBox QAbstractItemView {{
    background-color: {c('surface')};
    color: {c('text1')};
    border: 1px solid {c('border')};
    border-radius: {_RADIUS_CHIP}px;
    selection-background-color: {c('accentSoft')};
    selection-color: {c('accentStrong')};
    padding: 4px;
}}

QComboBox QAbstractItemView::item {{
    color: {c('text1')};
    background-color: {c('surface')};
}}

QComboBox QAbstractItemView::item:hover,
QComboBox QAbstractItemView::item:selected {{
    color: {c('accentStrong')};
    background-color: {c('accentSoft')};
}}

/* =========================================================== */
/* Spin Boxes                                                  */
/* =========================================================== */
QSpinBox, QDoubleSpinBox {{
    background-color: {c('surface')};
    border: 1px solid {c('borderStrong')};
    border-radius: {_RADIUS_CHIP}px;
    padding: 8px 12px;
    color: {c('text1')};
}}

QSpinBox:focus, QDoubleSpinBox:focus {{
    border: 1px solid {c('accent')};
}}

QSpinBox:disabled, QDoubleSpinBox:disabled {{
    background-color: {c('surface2')};
    color: {c('text3')};
    border-color: {c('border')};
}}

/* =========================================================== */
/* Checkboxes / Radios                                         */
/* =========================================================== */
QCheckBox {{
    color: {c('text1')};
    spacing: 10px;
}}

QCheckBox:disabled {{
    color: {c('text3')};
}}

QCheckBox::indicator {{
    width: 18px;
    height: 18px;
    border: 2px solid {c('borderStrong')};
    border-radius: 5px;
    background-color: {c('surface')};
}}

QCheckBox::indicator:hover {{
    border-color: {c('accent')};
}}

QCheckBox::indicator:checked {{
    background-color: {c('accent')};
    border-color: {c('accent')};
}}

QCheckBox::indicator:disabled {{
    background-color: {c('surface2')};
    border-color: {c('border')};
}}

QRadioButton {{
    color: {c('text1')};
    spacing: 10px;
}}

QRadioButton:disabled {{
    color: {c('text3')};
}}

QRadioButton::indicator {{
    width: 18px;
    height: 18px;
    border: 2px solid {c('borderStrong')};
    border-radius: 9px;
    background-color: {c('surface')};
}}

QRadioButton::indicator:hover {{
    border-color: {c('accent')};
}}

QRadioButton::indicator:checked {{
    background-color: {c('accent')};
    border-color: {c('accent')};
    border-width: 5px;
}}

QRadioButton::indicator:disabled {{
    background-color: {c('surface2')};
    border-color: {c('border')};
}}

/* =========================================================== */
/* Tab Widgets (lenguaje de pestañas del spike)                */
/* =========================================================== */
QTabWidget::pane {{
    border: 1px solid {c('border')};
    border-radius: {_RADIUS_CHIP}px;
    background-color: {c('surface')};
}}

QTabBar::tab {{
    background-color: transparent;
    color: {c('text2')};
    padding: 9px 18px;
    border: 1px solid transparent;
    border-bottom: 1px solid {c('border')};
    border-top-left-radius: {_RADIUS_CHIP}px;
    border-top-right-radius: {_RADIUS_CHIP}px;
    margin-right: 2px;
    font-size: {_FONT_SMALL}px;
}}

QTabBar::tab:hover {{
    background-color: {c('surface2')};
    color: {c('text1')};
}}

QTabBar::tab:selected {{
    background-color: {c('surface')};
    color: {c('text1')};
    font-weight: 600;
    border: 1px solid {c('border')};
    border-bottom: 2px solid {c('accent')};
}}

QTabBar::tab:selected:hover {{
    background-color: {c('surface')};
}}

/* =========================================================== */
/* Group Boxes                                                 */
/* =========================================================== */
QGroupBox {{
    border: 1px solid {c('border')};
    border-radius: {_RADIUS_CHIP}px;
    margin-top: 12px;
    padding-top: 16px;
    font-weight: 600;
    color: {c('text1')};
    background-color: {c('surface')};
}}

QGroupBox::title {{
    subcontrol-origin: margin;
    subcontrol-position: top left;
    padding: 0 8px;
    color: {c('text2')};
}}

/* =========================================================== */
/* Progress Bars                                               */
/* =========================================================== */
QProgressBar {{
    background-color: {c('surface3')};
    border: none;
    border-radius: 4px;
    text-align: center;
    color: {c('text2')};
    font-size: {_FONT_TINY}px;
}}

QProgressBar::chunk {{
    background-color: {c('accent')};
    border-radius: 4px;
}}

/* =========================================================== */
/* ToolTip (lenguaje del spike: superficie clara)              */
/* =========================================================== */
QToolTip {{
    background-color: {c('surface')};
    color: {c('text1')};
    border: 1px solid {c('border')};
    border-radius: 6px;
    padding: 6px 10px;
    font-size: {_FONT_SMALL}px;
}}

/* =========================================================== */
/* Palette Grid Buttons (paletas de herramientas)              */
/* =========================================================== */
#palette_grid QToolButton {{
    background-color: {c('surface')};
    border: 1px solid {c('border')};
    border-radius: 6px;
    padding: 4px;
    min-width: 28px;
    min-height: 28px;
}}

#palette_grid QToolButton:hover {{
    background-color: {c('surface2')};
    border: 1px solid {c('accentBorder')};
}}

#palette_grid QToolButton:checked {{
    background-color: {c('accentSoft')};
    border: 1px solid {c('accentBorder')};
    color: {c('accentStrong')};
}}

#palette_grid QToolButton:disabled {{
    color: {c('text3')};
    background-color: {c('surface2')};
    border: 1px solid {c('border')};
}}
"""


def get_tool_palette_stylesheet(theme_name: str) -> str:
    """Genera la hoja de estilo de las paletas verticales (100 % tokens)."""
    t = get_tokens(theme_name)

    def c(key: str) -> str:
        return str(t[key])

    return f"""
/* =========================================================== */
/* Tool palettes verticales (tokens: tema {theme_name})        */
/* =========================================================== */
QToolBar {{
    background-color: {c('surface')};
    border: none;
    border-right: 1px solid {c('border')};
    spacing: 6px;
    padding: 10px 8px;
}}

QToolBar::separator {{
    height: 1px;
    background-color: {c('border')};
    margin: 10px 8px;
}}

QToolButton {{
    background-color: {c('surface')};
    border: 1px solid {c('border')};
    border-radius: {_RADIUS_BTN}px;
    padding: 8px;
    min-width: 32px;
    min-height: 32px;
    color: {c('text2')};
}}

QToolButton:hover {{
    background-color: {c('surface2')};
    border: 1px solid {c('accentBorder')};
    color: {c('text1')};
}}

QToolButton:pressed {{
    background-color: {c('surface3')};
}}

QToolButton:checked {{
    background-color: {c('accentSoft')};
    border: 2px solid {c('accentBorder')};
    color: {c('accentStrong')};
}}

QToolButton:disabled {{
    color: {c('text3')};
    background-color: {c('surface2')};
    border: 1px solid {c('border')};
}}

#palette_grid QToolButton {{
    background-color: {c('surface')};
    border: 1px solid {c('border')};
    border-radius: 6px;
    padding: 6px;
    min-width: 30px;
    min-height: 30px;
}}

#palette_grid QToolButton:hover {{
    background-color: {c('surface2')};
    border: 1px solid {c('accentBorder')};
}}

#palette_grid QToolButton:disabled {{
    color: {c('text3')};
    background-color: {c('surface2')};
    border: 1px solid {c('border')};
}}
"""
