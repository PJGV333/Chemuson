"""
Estilos visuales de la aplicación Chemuson.

Define paletas de color y hojas de estilo Qt para la ventana principal,
menús, barras de herramientas y docks.
"""

# ============================================================================
# Paletas de Colores (Modo Claro y Oscuro)
# ============================================================================

# Colores Modo Claro
LIGHT_COLORS = {
    "primary_dark": "#1E293B",
    "primary_medium": "#334155",
    "accent_primary": "#0891B2",
    "accent_hover": "#06B6D4",
    "accent_pressed": "#0E7490",
    "bg_main": "#F8FAFC",
    "bg_elevated": "#FFFFFF",
    "bg_toolbar": "#F1F5F9",
    "bg_dock": "#F8FAFC",
    "border_light": "#E2E8F0",
    "border_medium": "#CBD5E1",
    "border_dark": "#94A3B8",
    "text_primary": "#0F172A",
    "text_secondary": "#475569",
    "text_muted": "#64748B",
    "text_inverse": "#F8FAFC",
    "palette_bg": "#F1F5F9",
    "palette_border": "#E2E8F0",
    "palette_button_bg": "#FFFFFF",
    "palette_button_border": "#CBD5E1",
    "palette_button_hover": "#ECFEFF",
    "palette_selected_bg": "#CFFAFE",
    "palette_selected_border": "#0891B2",
}

# Colores Modo Oscuro
DARK_COLORS = {
    "primary_dark": "#0F172A",
    "primary_medium": "#1E293B",
    "accent_primary": "#22D3EE",
    "accent_hover": "#67E8F9",
    "accent_pressed": "#06B6D4",
    "bg_main": "#0F172A",
    "bg_elevated": "#1E293B",
    "bg_toolbar": "#1E293B",
    "bg_dock": "#1E293B",
    "border_light": "#334155",
    "border_medium": "#475569",
    "border_dark": "#64748B",
    "text_primary": "#F1F5F9",
    "text_secondary": "#CBD5E1",
    "text_muted": "#94A3B8",
    "text_inverse": "#F8FAFC",
    "palette_bg": "#1E293B",
    "palette_border": "#334155",
    "palette_button_bg": "#334155",
    "palette_button_border": "#475569",
    "palette_button_hover": "#0F172A",
    "palette_selected_bg": "#164E63",
    "palette_selected_border": "#22D3EE",
}

# ============================================================================
# Generación de Hojas de Estilo Dinámicas
# ============================================================================

DEFAULT_THEME = "light"


def _get_theme_colors(theme_name: str = DEFAULT_THEME) -> dict[str, str]:
    """Devuelve la paleta asociada al tema solicitado."""
    return DARK_COLORS if theme_name == "dark" else LIGHT_COLORS


def get_main_stylesheet(theme_name: str = DEFAULT_THEME) -> str:
    """Genera la hoja de estilo principal según el tema."""
    colors = _get_theme_colors(theme_name)
    
    return f"""
/* =========================================== */
/* Main Window                                 */
/* =========================================== */
QMainWindow {{
    background-color: {colors['bg_main']};
}}

/* =========================================== */
/* Menu Bar - Modern Dark Header               */
/* =========================================== */
QMenuBar {{
    background-color: {colors['primary_dark']};
    color: {colors['text_inverse']};
    border: none;
    padding: 4px 10px;
    spacing: 4px;
    font-size: 13px;
    font-weight: 500;
    min-height: 32px;
}}

QMenuBar::item {{
    background: transparent;
    color: {colors['text_inverse']};
    padding: 6px 14px;
    border-radius: 6px;
    margin: 2px 2px;
}}

QMenuBar::item:selected {{
    background-color: {colors['accent_primary']};
    color: {colors['text_inverse']};
}}

QMenuBar::item:pressed {{
    background-color: {colors['accent_pressed']};
}}

/* =========================================== */
/* Dropdown Menus                              */
/* =========================================== */
QMenu {{
    background-color: {colors['bg_elevated']};
    color: {colors['text_primary']};
    border: 1px solid {colors['border_medium']};
    border-radius: 10px;
    padding: 8px 4px;
    margin: 4px;
}}

QMenu::item {{
    padding: 8px 28px 8px 14px;
    border-radius: 6px;
    margin: 2px 4px;
    color: {colors['text_primary']};
}}

QMenu::item:selected {{
    background-color: {colors['accent_primary']};
    color: {colors['text_inverse']};
}}

QMenu::item:disabled {{
    color: {colors['text_muted']};
}}

QMenu::separator {{
    height: 1px;
    background-color: {colors['border_light']};
    margin: 6px 14px;
}}

QMenu::icon {{
    margin-left: 8px;
}}

QMenu::indicator {{
    width: 18px;
    height: 18px;
    margin-left: 8px;
}}

/* =========================================== */
/* Toolbars                                    */
/* =========================================== */
QToolBar {{
    background-color: {colors['bg_toolbar']};
    color: {colors['text_primary']};
    border: none;
    border-bottom: 1px solid {colors['border_light']};
    spacing: 6px;
    padding: 8px 10px;
}}

QToolBar::separator {{
    width: 1px;
    background-color: {colors['border_medium']};
    margin: 8px 10px;
}}

/* =========================================== */
/* Tool Buttons                                */
/* =========================================== */
QToolButton {{
    background-color: transparent;
    border: 1px solid transparent;
    border-radius: 8px;
    padding: 8px;
    color: {colors['text_primary']};
}}

QToolButton:hover {{
    background-color: {colors['palette_button_hover']};
    border: 1px solid {colors['border_medium']};
}}

QToolButton:pressed {{
    background-color: {colors['border_light']};
    border: 1px solid {colors['border_dark']};
}}

QToolButton:checked {{
    background-color: {colors['palette_selected_bg']};
    border: 1px solid {colors['palette_selected_border']};
}}

QToolButton::menu-indicator {{
    image: none;
    subcontrol-position: right bottom;
    subcontrol-origin: padding;
    width: 8px;
    height: 8px;
}}

/* =========================================== */
/* Dock Widgets                                */
/* =========================================== */
QDockWidget {{
    titlebar-close-icon: none;
    titlebar-normal-icon: none;
    font-weight: 600;
    color: {colors['text_primary']};
}}

QDockWidget::title {{
    background: qlineargradient(x1:0, y1:0, x2:0, y2=1,
        stop:0 {colors['bg_toolbar']}, stop:1 {colors['bg_elevated']});
    padding: 12px 14px;
    border-bottom: 1px solid {colors['border_light']};
    text-align: left;
    font-size: 13px;
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
    background-color: {colors['palette_button_hover']};
    border-radius: 6px;
}}

/* =========================================== */
/* Status Bar                                  */
/* =========================================== */
QStatusBar {{
    background-color: {colors['primary_dark']};
    color: {colors['text_inverse']};
    border-top: none;
    padding: 6px 14px;
    font-size: 12px;
    min-height: 28px;
}}

QStatusBar::item {{
    border: none;
}}

QStatusBar QLabel {{
    color: {colors['text_inverse']};
    padding: 0 4px;
}}

/* =========================================== */
/* Scrollbars                                  */
/* =========================================== */
QScrollBar:vertical {{
    background-color: {colors['bg_toolbar']};
    width: 10px;
    border-radius: 5px;
    margin: 2px;
}}

QScrollBar::handle:vertical {{
    background-color: {colors['border_medium']};
    min-height: 30px;
    border-radius: 4px;
    margin: 2px;
}}

QScrollBar::handle:vertical:hover {{
    background-color: {colors['text_muted']};
}}

QScrollBar::handle:vertical:pressed {{
    background-color: {colors['text_secondary']};
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
    background-color: {colors['bg_toolbar']};
    height: 10px;
    border-radius: 5px;
    margin: 2px;
}}

QScrollBar::handle:horizontal {{
    background-color: {colors['border_medium']};
    min-width: 30px;
    border-radius: 4px;
    margin: 2px;
}}

QScrollBar::handle:horizontal:hover {{
    background-color: {colors['text_muted']};
}}

QScrollBar::handle:horizontal:pressed {{
    background-color: {colors['text_secondary']};
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

/* =========================================== */
/* Tables (Inspector Dock)                     */
/* =========================================== */
QTableWidget {{
    background-color: {colors['bg_elevated']};
    alternate-background-color: {colors['bg_dock']};
    border: 1px solid {colors['border_light']};
    border-radius: 8px;
    gridline-color: {colors['border_light']};
}}

QTableWidget::item {{
    padding: 6px 10px;
    color: {colors['text_primary']};
}}

QTableWidget::item:selected {{
    background-color: {colors['palette_selected_bg']};
    color: {colors['text_primary']};
}}

QTableWidget::item:hover {{
    background-color: {colors['palette_button_hover']};
}}

QHeaderView::section {{
    background-color: {colors['bg_toolbar']};
    color: {colors['text_secondary']};
    padding: 8px 10px;
    border: none;
    border-bottom: 2px solid {colors['border_medium']};
    font-weight: 600;
    font-size: 12px;
}}

/* =========================================== */
/* Tree Widgets (Templates Dock)               */
/* =========================================== */
QTreeWidget {{
    background-color: {colors['bg_elevated']};
    border: 1px solid {colors['border_light']};
    border-radius: 8px;
}}

QTreeWidget::item {{
    padding: 8px 6px;
    border-radius: 6px;
    margin: 2px 4px;
}}

QTreeWidget::item:hover {{
    background-color: {colors['palette_button_hover']};
}}

QTreeWidget::item:selected {{
    background-color: {colors['palette_selected_bg']};
    color: {colors['text_primary']};
    border: 1px solid {colors['palette_selected_border']};
}}

QTreeWidget::branch:has-children:!has-siblings:closed,
QTreeWidget::branch:closed:has-children:has-siblings {{
    border-image: none;
}}

QTreeWidget::branch:open:has-children:!has-siblings,
QTreeWidget::branch:open:has-children:has-siblings {{
    border-image: none;
}}

/* =========================================== */
/* Labels                                      */
/* =========================================== */
QLabel {{
    color: {colors['text_primary']};
}}

/* =========================================== */
/* Dialogs                                     */
/* =========================================== */
QDialog {{
    background-color: {colors['bg_elevated']};
}}

QDialog QLabel {{
    color: {colors['text_primary']};
}}

/* =========================================== */
/* Push Buttons                                */
/* =========================================== */
QPushButton {{
    background-color: {colors['accent_primary']};
    color: {colors['text_inverse']};
    border: none;
    border-radius: 8px;
    padding: 10px 22px;
    font-weight: 600;
    min-width: 80px;
}}

QPushButton:hover {{
    background-color: {colors['accent_hover']};
}}

QPushButton:pressed {{
    background-color: {colors['accent_pressed']};
}}

QPushButton:disabled {{
    background-color: {colors['border_medium']};
    color: {colors['bg_toolbar']};
}}

/* Secondary button style */
QPushButton[flat="true"] {{
    background-color: transparent;
    color: {colors['accent_primary']};
    border: 1px solid {colors['accent_primary']};
}}

QPushButton[flat="true"]:hover {{
    background-color: {colors['palette_button_hover']};
}}

/* =========================================== */
/* Line Edits                                  */
/* =========================================== */
QLineEdit {{
    background-color: {colors['bg_elevated']};
    border: 1px solid {colors['border_medium']};
    border-radius: 8px;
    padding: 8px 12px;
    color: {colors['text_primary']};
    selection-background-color: {colors['palette_selected_bg']};
}}

QLineEdit:focus {{
    border: 1px solid {colors['accent_primary']};
    background-color: {colors['bg_elevated']};
}}

QLineEdit:disabled {{
    background-color: {colors['bg_dock']};
    color: {colors['text_muted']};
    border-color: {colors['border_light']};
}}

/* =========================================== */
/* Combo Boxes                                 */
/* =========================================== */
QComboBox {{
    background-color: {colors['bg_elevated']};
    border: 1px solid {colors['border_medium']};
    border-radius: 8px;
    padding: 8px 12px;
    color: {colors['text_primary']};
    min-width: 100px;
}}

QComboBox:hover {{
    border-color: {colors['accent_primary']};
}}

QComboBox:focus {{
    border: 1px solid {colors['accent_primary']};
}}

QComboBox::drop-down {{
    border: none;
    width: 24px;
}}

QComboBox QAbstractItemView {{
    background-color: {colors['bg_elevated']};
    color: {colors['text_primary']};
    border: 1px solid {colors['border_medium']};
    border-radius: 8px;
    selection-background-color: {colors['accent_primary']};
    selection-color: {colors['text_inverse']};
    padding: 4px;
}}

QComboBox QAbstractItemView::item {{
    color: {colors['text_primary']};
    background-color: {colors['bg_elevated']};
}}

QComboBox QAbstractItemView::item:selected {{
    color: {colors['text_inverse']};
    background-color: {colors['accent_primary']};
}}

/* =========================================== */
/* Spin Boxes                                  */
/* =========================================== */
QSpinBox, QDoubleSpinBox {{
    background-color: {colors['bg_elevated']};
    border: 1px solid {colors['border_medium']};
    border-radius: 8px;
    padding: 8px 12px;
    color: {colors['text_primary']};
}}

QSpinBox:focus, QDoubleSpinBox:focus {{
    border: 1px solid {colors['accent_primary']};
}}

QSpinBox:disabled, QDoubleSpinBox:disabled {{
    background-color: {colors['bg_dock']};
    color: {colors['text_muted']};
    border-color: {colors['border_light']};
}}

/* =========================================== */
/* Checkboxes                                  */
/* =========================================== */
QCheckBox {{
    color: {colors['text_primary']};
    spacing: 10px;
}}

QCheckBox::indicator {{
    width: 20px;
    height: 20px;
    border: 2px solid {colors['border_dark']};
    border-radius: 5px;
    background-color: {colors['bg_elevated']};
}}

QCheckBox::indicator:hover {{
    border-color: {colors['accent_primary']};
}}

QCheckBox::indicator:checked {{
    background-color: {colors['accent_primary']};
    border-color: {colors['accent_primary']};
}}

QCheckBox::indicator:disabled {{
    background-color: {colors['bg_dock']};
    border-color: {colors['border_light']};
}}

/* =========================================== */
/* Radio Buttons                               */
/* =========================================== */
QRadioButton {{
    color: {colors['text_primary']};
    spacing: 10px;
}}

QRadioButton::indicator {{
    width: 20px;
    height: 20px;
    border: 2px solid {colors['border_dark']};
    border-radius: 10px;
    background-color: {colors['bg_elevated']};
}}

QRadioButton::indicator:hover {{
    border-color: {colors['accent_primary']};
}}

QRadioButton::indicator:checked {{
    background-color: {colors['accent_primary']};
    border-color: {colors['accent_primary']};
    border-width: 5px;
}}

/* =========================================== */
/* Tab Widgets                                 */
/* =========================================== */
QTabWidget::pane {{
    border: 1px solid {colors['border_light']};
    border-radius: 8px;
    background-color: {colors['bg_elevated']};
}}

QTabBar::tab {{
    background-color: {colors['bg_toolbar']};
    color: {colors['text_secondary']};
    padding: 10px 18px;
    border: 1px solid {colors['border_light']};
    border-bottom: none;
    border-top-left-radius: 8px;
    border-top-right-radius: 8px;
    margin-right: 2px;
}}

QTabBar::tab:selected {{
    background-color: {colors['bg_elevated']};
    color: {colors['text_primary']};
    border-bottom: 2px solid {colors['accent_primary']};
}}

QTabBar::tab:hover:!selected {{
    background-color: {colors['palette_button_hover']};
    color: {colors['text_primary']};
}}

/* =========================================== */
/* Group Boxes                                 */
/* =========================================== */
QGroupBox {{
    border: 1px solid {colors['border_light']};
    border-radius: 8px;
    margin-top: 12px;
    padding-top: 16px;
    font-weight: 600;
    color: {colors['text_primary']};
}}

QGroupBox::title {{
    subcontrol-origin: margin;
    subcontrol-position: top left;
    padding: 0 8px;
    color: {colors['text_secondary']};
}}

/* =========================================== */
/* ToolTip                                     */
/* =========================================== */
QToolTip {{
    background-color: {colors['primary_dark']};
    color: {colors['text_inverse']};
    border: 1px solid {colors['primary_medium']};
    border-radius: 6px;
    padding: 6px 10px;
    font-size: 12px;
}}

/* =========================================== */
/* Palette Grid Buttons                        */
/* =========================================== */
#palette_grid QToolButton {{
    background-color: {colors['palette_button_bg']};
    border: 1px solid {colors['palette_button_border']};
    border-radius: 6px;
    padding: 4px;
    min-width: 28px;
    min-height: 28px;
}}

#palette_grid QToolButton:hover {{
    background-color: {colors['palette_button_hover']};
    border: 1px solid {colors['palette_selected_border']};
}}

#palette_grid QToolButton:disabled {{
    color: {colors['text_muted']};
    background-color: {colors['bg_dock']};
    border: 1px solid {colors['border_light']};
}}
"""

def get_tool_palette_stylesheet(theme_name: str = DEFAULT_THEME) -> str:
    """Genera la hoja de estilo de la barra de herramientas."""
    colors = _get_theme_colors(theme_name)
    
    return f"""
QToolBar {{
    background-color: {colors['palette_bg']};
    border: none;
    border-right: 1px solid {colors['border_light']};
    spacing: 6px;
    padding: 10px 8px;
}}

QToolBar::separator {{
    height: 1px;
    background-color: {colors['border_medium']};
    margin: 10px 8px;
}}

QToolButton {{
    background-color: {colors['palette_button_bg']};
    border: 1px solid {colors['palette_button_border']};
    border-radius: 8px;
    padding: 8px;
    min-width: 32px;
    min-height: 32px;
}}

QToolButton:hover {{
    background-color: {colors['palette_button_hover']};
    border: 1px solid {colors['accent_primary']};
}}

QToolButton:pressed {{
    background-color: {colors['border_light']};
}}

QToolButton:checked {{
    background-color: {colors['palette_selected_bg']};
    border: 2px solid {colors['accent_primary']};
}}

QToolButton[popupMode="1"] {{
    padding-right: 14px;
}}

QToolButton::menu-indicator {{
    subcontrol-position: right bottom;
    subcontrol-origin: padding;
    right: 3px;
    bottom: 3px;
}}

#palette_grid QToolButton {{
    background-color: {colors['palette_button_bg']};
    border: 1px solid {colors['palette_button_border']};
    border-radius: 6px;
    padding: 6px;
    min-width: 30px;
    min-height: 30px;
}}

#palette_grid QToolButton:hover {{
    background-color: {colors['palette_button_hover']};
    border: 1px solid {colors['accent_primary']};
}}

#palette_grid QToolButton:disabled {{
    color: {colors['text_muted']};
    background-color: {colors['bg_dock']};
    border: 1px solid {colors['border_light']};
}}
"""


# Compatibilidad retroactiva para módulos que aún importan constantes.
MAIN_STYLESHEET = get_main_stylesheet(DEFAULT_THEME)
TOOL_PALETTE_STYLESHEET = get_tool_palette_stylesheet(DEFAULT_THEME)
