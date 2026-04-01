"""Sistema de estilos de Chemuson basado en tokens por tema."""

from __future__ import annotations

from dataclasses import dataclass


ThemeName = str


@dataclass(frozen=True)
class ThemeTokens:
    """Colección de tokens de color para construir paletas y QSS."""

    window: str
    surface: str
    surface_alt: str
    panel: str
    border: str
    border_strong: str
    text: str
    text_muted: str
    text_inverse: str
    accent: str
    accent_hover: str
    accent_pressed: str
    selection: str
    selection_text: str
    tooltip_bg: str
    tooltip_text: str
    menubar_bg: str
    menubar_text: str
    status_bg: str
    status_text: str
    scrollbar_bg: str
    scrollbar_handle: str
    scrollbar_handle_hover: str


LIGHT_TOKENS = ThemeTokens(
    window="#F4F7FB",
    surface="#FFFFFF",
    surface_alt="#EEF3F8",
    panel="#ECF2F8",
    border="#D4DEE8",
    border_strong="#A8B6C8",
    text="#0F172A",
    text_muted="#55657A",
    text_inverse="#F8FAFC",
    accent="#0B84B8",
    accent_hover="#0E9AD5",
    accent_pressed="#086489",
    selection="#D5EEFF",
    selection_text="#0F172A",
    tooltip_bg="#1E293B",
    tooltip_text="#F8FAFC",
    menubar_bg="#1E293B",
    menubar_text="#F8FAFC",
    status_bg="#1E293B",
    status_text="#F8FAFC",
    scrollbar_bg="#E8EEF5",
    scrollbar_handle="#B3C0D1",
    scrollbar_handle_hover="#8D9AAF",
)


DARK_TOKENS = ThemeTokens(
    window="#111926",
    surface="#172233",
    surface_alt="#1D2A3D",
    panel="#1B2638",
    border="#2F415B",
    border_strong="#49607E",
    text="#E4ECF6",
    text_muted="#A9B8CC",
    text_inverse="#EAF2FF",
    accent="#46A6E3",
    accent_hover="#67B8EC",
    accent_pressed="#328EC8",
    selection="#2A4666",
    selection_text="#EAF2FF",
    tooltip_bg="#0D1522",
    tooltip_text="#EAF2FF",
    menubar_bg="#0F1728",
    menubar_text="#EAF2FF",
    status_bg="#0F1728",
    status_text="#EAF2FF",
    scrollbar_bg="#1B2638",
    scrollbar_handle="#41546C",
    scrollbar_handle_hover="#576F8C",
)


THEME_TOKENS: dict[ThemeName, ThemeTokens] = {
    "light": LIGHT_TOKENS,
    "dark": DARK_TOKENS,
}


def theme_tokens(theme: ThemeName) -> ThemeTokens:
    """Devuelve tokens válidos para el tema solicitado."""
    return THEME_TOKENS.get(str(theme).strip().lower(), LIGHT_TOKENS)


def build_main_stylesheet(tokens: ThemeTokens) -> str:
    """Construye QSS principal a partir de tokens."""
    return f"""
QMainWindow, QWidget {{ background-color: {tokens.window}; color: {tokens.text}; }}
QDialog {{ background-color: {tokens.surface}; color: {tokens.text}; }}
QMenuBar {{ background-color: {tokens.menubar_bg}; color: {tokens.menubar_text}; padding: 4px 8px; spacing: 4px; border: none; }}
QMenuBar::item {{ background: transparent; padding: 6px 12px; border-radius: 6px; }}
QMenuBar::item:selected {{ background-color: {tokens.accent}; color: {tokens.text_inverse}; }}
QMenu {{ background-color: {tokens.surface}; color: {tokens.text}; border: 1px solid {tokens.border}; border-radius: 8px; padding: 6px; }}
QMenu::item {{ padding: 7px 24px 7px 12px; border-radius: 5px; }}
QMenu::item:selected {{ background: {tokens.accent}; color: {tokens.text_inverse}; }}
QMenu::separator {{ height: 1px; background: {tokens.border}; margin: 5px 10px; }}
QToolBar {{ background-color: {tokens.panel}; border: none; border-bottom: 1px solid {tokens.border}; spacing: 6px; padding: 6px; }}
QToolBar::separator {{ width: 1px; background-color: {tokens.border}; margin: 6px; }}
QToolButton {{ background: transparent; border: 1px solid transparent; border-radius: 7px; padding: 7px; color: {tokens.text}; }}
QToolButton:hover {{ background-color: {tokens.surface_alt}; border-color: {tokens.border_strong}; }}
QToolButton:checked {{ background-color: {tokens.selection}; border-color: {tokens.accent}; }}
QToolButton:pressed {{ background-color: {tokens.surface_alt}; border-color: {tokens.border_strong}; }}
QStatusBar {{ background: {tokens.status_bg}; color: {tokens.status_text}; border: none; padding: 4px 10px; }}
QStatusBar QLabel {{ color: {tokens.status_text}; }}
QTabWidget::pane {{ border: 1px solid {tokens.border}; background: {tokens.surface}; border-radius: 8px; }}
QTabBar::tab {{ background: {tokens.panel}; color: {tokens.text_muted}; padding: 8px 14px; border: 1px solid {tokens.border}; border-bottom: none; border-top-left-radius: 8px; border-top-right-radius: 8px; }}
QTabBar::tab:selected {{ background: {tokens.surface}; color: {tokens.text}; border-bottom: 2px solid {tokens.accent}; }}
QDockWidget {{ color: {tokens.text}; }}
QDockWidget::title {{ background: {tokens.panel}; color: {tokens.text}; padding: 8px 10px; border-bottom: 1px solid {tokens.border}; }}
QPushButton {{ background: {tokens.accent}; color: {tokens.text_inverse}; border: none; border-radius: 7px; padding: 8px 16px; font-weight: 600; }}
QPushButton:hover {{ background: {tokens.accent_hover}; }}
QPushButton:pressed {{ background: {tokens.accent_pressed}; }}
QLineEdit, QComboBox, QSpinBox, QDoubleSpinBox, QTextEdit, QPlainTextEdit {{ background: {tokens.surface}; color: {tokens.text}; border: 1px solid {tokens.border_strong}; border-radius: 7px; padding: 6px 9px; selection-background-color: {tokens.selection}; selection-color: {tokens.selection_text}; }}
QLineEdit:focus, QComboBox:focus, QSpinBox:focus, QDoubleSpinBox:focus, QTextEdit:focus, QPlainTextEdit:focus {{ border-color: {tokens.accent}; }}
QComboBox QAbstractItemView {{ background: {tokens.surface}; color: {tokens.text}; border: 1px solid {tokens.border}; selection-background-color: {tokens.selection}; selection-color: {tokens.selection_text}; }}
QCheckBox, QRadioButton {{ color: {tokens.text}; spacing: 8px; }}
QCheckBox::indicator, QRadioButton::indicator {{ width: 17px; height: 17px; border: 1px solid {tokens.border_strong}; background: {tokens.surface}; }}
QRadioButton::indicator {{ border-radius: 8px; }}
QCheckBox::indicator {{ border-radius: 4px; }}
QCheckBox::indicator:checked, QRadioButton::indicator:checked {{ background: {tokens.accent}; border-color: {tokens.accent}; }}
QTableWidget, QTreeWidget {{ background: {tokens.surface}; color: {tokens.text}; border: 1px solid {tokens.border}; gridline-color: {tokens.border}; alternate-background-color: {tokens.surface_alt}; }}
QTableWidget::item:selected, QTreeWidget::item:selected {{ background: {tokens.selection}; color: {tokens.selection_text}; }}
QHeaderView::section {{ background: {tokens.panel}; color: {tokens.text_muted}; border: none; border-bottom: 1px solid {tokens.border}; padding: 6px 8px; font-weight: 600; }}
QScrollBar:vertical {{ background: {tokens.scrollbar_bg}; width: 12px; margin: 2px; border-radius: 6px; }}
QScrollBar:horizontal {{ background: {tokens.scrollbar_bg}; height: 12px; margin: 2px; border-radius: 6px; }}
QScrollBar::handle:vertical, QScrollBar::handle:horizontal {{ background: {tokens.scrollbar_handle}; border-radius: 5px; min-height: 24px; min-width: 24px; margin: 2px; }}
QScrollBar::handle:vertical:hover, QScrollBar::handle:horizontal:hover {{ background: {tokens.scrollbar_handle_hover}; }}
QScrollBar::add-line, QScrollBar::sub-line, QScrollBar::add-page, QScrollBar::sub-page {{ background: transparent; border: none; }}
QToolTip {{ background: {tokens.tooltip_bg}; color: {tokens.tooltip_text}; border: 1px solid {tokens.border}; border-radius: 6px; padding: 5px 8px; }}
#palette_grid QToolButton {{ background: {tokens.surface}; border: 1px solid {tokens.border}; border-radius: 6px; padding: 4px; min-width: 28px; min-height: 28px; }}
#palette_grid QToolButton:hover {{ background: {tokens.surface_alt}; border-color: {tokens.accent}; }}
"""


def build_tool_palette_stylesheet(tokens: ThemeTokens) -> str:
    """Construye QSS específico para paletas laterales."""
    return f"""
QToolBar {{ background: {tokens.panel}; border: none; border-right: 1px solid {tokens.border}; spacing: 6px; padding: 8px 6px; }}
QToolButton {{ background: {tokens.surface}; border: 1px solid {tokens.border}; border-radius: 8px; padding: 8px; min-width: 32px; min-height: 32px; }}
QToolButton:hover {{ background: {tokens.surface_alt}; border-color: {tokens.accent}; }}
QToolButton:checked {{ background: {tokens.selection}; border: 2px solid {tokens.accent}; }}
#palette_grid QToolButton {{ min-width: 30px; min-height: 30px; padding: 5px; }}
"""


# Compatibilidad con importaciones antiguas.
MAIN_STYLESHEET = build_main_stylesheet(LIGHT_TOKENS)
TOOL_PALETTE_STYLESHEET = build_tool_palette_stylesheet(LIGHT_TOKENS)
