"""
Estilos visuales de la aplicación Chemuson.

Define paletas de color y hojas de estilo Qt para la ventana principal,
menús, barras de herramientas y docks.
"""

# ============================================================================
# Paleta moderna de colores
# ============================================================================

# Colores principales
PRIMARY_DARK = "#2D3E50"        # Deep blue-gray for headers
PRIMARY_MEDIUM = "#34495E"      # Slightly lighter variant
ACCENT_PRIMARY = "#00BCD4"      # Cyan/teal accent
ACCENT_HOVER = "#26C6DA"        # Lighter cyan for hover
ACCENT_PRESSED = "#00ACC1"      # Darker cyan for pressed

# Colores de fondo
BG_MAIN = "#ECEFF1"             # Light blue-gray background
BG_ELEVATED = "#FFFFFF"         # White for elevated surfaces
BG_TOOLBAR = "#F5F7FA"          # Very light gray for toolbars
BG_DOCK = "#FAFBFC"             # Almost white for dock panels

# Colores de bordes y separadores
BORDER_LIGHT = "#E0E4E8"        # Subtle border
BORDER_MEDIUM = "#CFD8DC"       # Medium border
BORDER_DARK = "#B0BEC5"         # Darker border

# Colores de texto
TEXT_PRIMARY = "#212121"        # Near-black for primary text
TEXT_SECONDARY = "#546E7A"      # Blue-gray for secondary text
TEXT_MUTED = "#78909C"          # Muted text
TEXT_INVERSE = "#FFFFFF"        # White text on dark backgrounds

# Paleta específica de herramientas
PALETTE_BG = "#F5F7FA"
PALETTE_BORDER = "#E0E4E8"
PALETTE_BUTTON_BG = "#FFFFFF"
PALETTE_BUTTON_BORDER = "#CFD8DC"
PALETTE_BUTTON_HOVER = "#E3F2FD"
PALETTE_SELECTED_BG = "#B2EBF2"
PALETTE_SELECTED_BORDER = "#00BCD4"

# Sombra para efecto de profundidad (Qt no soporta box-shadow estándar)
SHADOW_COLOR = "rgba(0, 0, 0, 0.1)"

# ============================================================================
# Hoja de estilos principal
# ============================================================================

MAIN_STYLESHEET = f"""
/* =========================================== */
/* Main Window                                 */
/* =========================================== */
QMainWindow {{
    background-color: {BG_MAIN};
}}

/* =========================================== */
/* Menu Bar - Modern Dark Header               */
/* =========================================== */
QMenuBar {{
    background-color: {PRIMARY_DARK};
    color: {TEXT_INVERSE};
    border: none;
    padding: 4px 8px;
    spacing: 4px;
    font-size: 13px;
    font-weight: 500;
    min-height: 28px;
}}

QMenuBar::item {{
    background: transparent;
    padding: 6px 12px;
    border-radius: 4px;
    margin: 1px 2px;
}}

QMenuBar::item:selected {{
    background-color: {ACCENT_PRIMARY};
    color: {TEXT_INVERSE};
}}

QMenuBar::item:pressed {{
    background-color: {ACCENT_PRESSED};
}}

/* =========================================== */
/* Dropdown Menus                              */
/* =========================================== */
QMenu {{
    background-color: {BG_ELEVATED};
    border: 1px solid {BORDER_MEDIUM};
    border-radius: 8px;
    padding: 6px 4px;
    margin: 2px;
}}

QMenu::item {{
    padding: 8px 24px 8px 12px;
    border-radius: 4px;
    margin: 2px 4px;
    color: {TEXT_PRIMARY};
}}

QMenu::item:selected {{
    background-color: {ACCENT_PRIMARY};
    color: {TEXT_INVERSE};
}}

QMenu::item:disabled {{
    color: {TEXT_MUTED};
}}

QMenu::separator {{
    height: 1px;
    background-color: {BORDER_LIGHT};
    margin: 6px 12px;
}}

QMenu::icon {{
    margin-left: 8px;
}}

QMenu::indicator {{
    width: 16px;
    height: 16px;
    margin-left: 8px;
}}

/* =========================================== */
/* Toolbars                                    */
/* =========================================== */
QToolBar {{
    background-color: {BG_TOOLBAR};
    border: none;
    border-bottom: 1px solid {BORDER_LIGHT};
    spacing: 6px;
    padding: 6px 8px;
}}

QToolBar::separator {{
    width: 1px;
    background-color: {BORDER_MEDIUM};
    margin: 6px 8px;
}}

/* =========================================== */
/* Tool Buttons                                */
/* =========================================== */
QToolButton {{
    background-color: transparent;
    border: 1px solid transparent;
    border-radius: 6px;
    padding: 6px;
    color: {TEXT_PRIMARY};
}}

QToolButton:hover {{
    background-color: {PALETTE_BUTTON_HOVER};
    border: 1px solid {BORDER_MEDIUM};
}}

QToolButton:pressed {{
    background-color: {BORDER_LIGHT};
    border: 1px solid {BORDER_DARK};
}}

QToolButton:checked {{
    background-color: {PALETTE_SELECTED_BG};
    border: 1px solid {ACCENT_PRIMARY};
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
    color: {TEXT_PRIMARY};
}}

QDockWidget::title {{
    background: qlineargradient(x1:0, y1:0, x2:0, y2:1,
        stop:0 {BG_TOOLBAR}, stop:1 {BG_DOCK});
    padding: 10px 12px;
    border-bottom: 1px solid {BORDER_LIGHT};
    text-align: left;
    font-size: 13px;
}}

QDockWidget::close-button,
QDockWidget::float-button {{
    border: none;
    background: transparent;
    padding: 2px;
}}

QDockWidget::close-button:hover,
QDockWidget::float-button:hover {{
    background-color: {PALETTE_BUTTON_HOVER};
    border-radius: 4px;
}}

/* =========================================== */
/* Status Bar                                  */
/* =========================================== */
QStatusBar {{
    background-color: {PRIMARY_DARK};
    color: {TEXT_INVERSE};
    border-top: none;
    padding: 4px 12px;
    font-size: 12px;
    min-height: 24px;
}}

QStatusBar::item {{
    border: none;
}}

QStatusBar QLabel {{
    color: {TEXT_INVERSE};
    padding: 0 4px;
}}

/* =========================================== */
/* Scrollbars                                  */
/* =========================================== */
QScrollBar:vertical {{
    background-color: {BG_TOOLBAR};
    width: 12px;
    border-radius: 6px;
    margin: 2px;
}}

QScrollBar::handle:vertical {{
    background-color: {BORDER_MEDIUM};
    min-height: 30px;
    border-radius: 5px;
    margin: 2px;
}}

QScrollBar::handle:vertical:hover {{
    background-color: {TEXT_MUTED};
}}

QScrollBar::handle:vertical:pressed {{
    background-color: {TEXT_SECONDARY};
}}

QScrollBar::add-line:vertical,
QScrollBar::sub-line:vertical {{
    height: 14px;
    background-color: {BG_TOOLBAR};
    border: 1px solid {BORDER_LIGHT};
    border-radius: 4px;
    subcontrol-origin: margin;
}}

QScrollBar::sub-line:vertical {{
    subcontrol-position: top;
}}

QScrollBar::add-line:vertical {{
    subcontrol-position: bottom;
}}

QScrollBar::add-page:vertical,
QScrollBar::sub-page:vertical {{
    background: transparent;
}}

QScrollBar:horizontal {{
    background-color: {BG_TOOLBAR};
    height: 12px;
    border-radius: 6px;
    margin: 2px;
}}

QScrollBar::handle:horizontal {{
    background-color: {BORDER_MEDIUM};
    min-width: 30px;
    border-radius: 5px;
    margin: 2px;
}}

QScrollBar::handle:horizontal:hover {{
    background-color: {TEXT_MUTED};
}}

QScrollBar::add-line:horizontal,
QScrollBar::sub-line:horizontal {{
    width: 14px;
    background-color: {BG_TOOLBAR};
    border: 1px solid {BORDER_LIGHT};
    border-radius: 4px;
    subcontrol-origin: margin;
}}

QScrollBar::sub-line:horizontal {{
    subcontrol-position: left;
}}

QScrollBar::add-line:horizontal {{
    subcontrol-position: right;
}}

/* =========================================== */
/* Tables (Inspector Dock)                     */
/* =========================================== */
QTableWidget {{
    background-color: {BG_ELEVATED};
    alternate-background-color: {BG_DOCK};
    border: 1px solid {BORDER_LIGHT};
    border-radius: 4px;
    gridline-color: {BORDER_LIGHT};
}}

QTableWidget::item {{
    padding: 6px 8px;
    color: {TEXT_PRIMARY};
}}

QTableWidget::item:selected {{
    background-color: {PALETTE_SELECTED_BG};
    color: {TEXT_PRIMARY};
}}

QHeaderView::section {{
    background-color: {BG_TOOLBAR};
    color: {TEXT_SECONDARY};
    padding: 8px;
    border: none;
    border-bottom: 1px solid {BORDER_MEDIUM};
    font-weight: 600;
}}

/* =========================================== */
/* Tree Widgets (Templates Dock)               */
/* =========================================== */
QTreeWidget {{
    background-color: {BG_ELEVATED};
    border: 1px solid {BORDER_LIGHT};
    border-radius: 4px;
}}

QTreeWidget::item {{
    padding: 6px 4px;
    border-radius: 4px;
    margin: 1px 4px;
}}

QTreeWidget::item:hover {{
    background-color: {PALETTE_BUTTON_HOVER};
}}

QTreeWidget::item:selected {{
    background-color: {PALETTE_SELECTED_BG};
    color: {TEXT_PRIMARY};
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
    color: {TEXT_PRIMARY};
}}

/* =========================================== */
/* Dialogs                                     */
/* =========================================== */
QDialog {{
    background-color: {BG_ELEVATED};
}}

QDialog QLabel {{
    color: {TEXT_PRIMARY};
}}

/* =========================================== */
/* Push Buttons                                */
/* =========================================== */
QPushButton {{
    background-color: {ACCENT_PRIMARY};
    color: {TEXT_INVERSE};
    border: none;
    border-radius: 6px;
    padding: 8px 20px;
    font-weight: 600;
    min-width: 80px;
}}

QPushButton:hover {{
    background-color: {ACCENT_HOVER};
}}

QPushButton:pressed {{
    background-color: {ACCENT_PRESSED};
}}

QPushButton:disabled {{
    background-color: {BORDER_MEDIUM};
    color: {TEXT_MUTED};
}}

/* Secondary button style */
QPushButton[flat="true"] {{
    background-color: transparent;
    color: {ACCENT_PRIMARY};
    border: 1px solid {ACCENT_PRIMARY};
}}

QPushButton[flat="true"]:hover {{
    background-color: {PALETTE_BUTTON_HOVER};
}}

/* =========================================== */
/* Line Edits                                  */
/* =========================================== */
QLineEdit {{
    background-color: {BG_ELEVATED};
    border: 1px solid {BORDER_MEDIUM};
    border-radius: 6px;
    padding: 8px 12px;
    color: {TEXT_PRIMARY};
    selection-background-color: {PALETTE_SELECTED_BG};
}}

QLineEdit:focus {{
    border: 2px solid {ACCENT_PRIMARY};
    padding: 7px 11px;
}}

/* =========================================== */
/* Combo Boxes                                 */
/* =========================================== */
QComboBox {{
    background-color: {BG_ELEVATED};
    border: 1px solid {BORDER_MEDIUM};
    border-radius: 6px;
    padding: 6px 12px;
    color: {TEXT_PRIMARY};
    min-width: 100px;
}}

QComboBox:hover {{
    border-color: {ACCENT_PRIMARY};
}}

QComboBox::drop-down {{
    border: none;
    width: 24px;
}}

QComboBox QAbstractItemView {{
    background-color: {BG_ELEVATED};
    border: 1px solid {BORDER_MEDIUM};
    border-radius: 4px;
    selection-background-color: {ACCENT_PRIMARY};
    selection-color: {TEXT_INVERSE};
}}

/* =========================================== */
/* Spin Boxes                                  */
/* =========================================== */
QSpinBox, QDoubleSpinBox {{
    background-color: {BG_ELEVATED};
    border: 1px solid {BORDER_MEDIUM};
    border-radius: 6px;
    padding: 6px 12px;
    color: {TEXT_PRIMARY};
}}

QSpinBox:focus, QDoubleSpinBox:focus {{
    border: 2px solid {ACCENT_PRIMARY};
}}

/* =========================================== */
/* Checkboxes                                  */
/* =========================================== */
QCheckBox {{
    color: {TEXT_PRIMARY};
    spacing: 8px;
}}

QCheckBox::indicator {{
    width: 18px;
    height: 18px;
    border: 2px solid {BORDER_DARK};
    border-radius: 4px;
    background-color: {BG_ELEVATED};
}}

QCheckBox::indicator:hover {{
    border-color: {ACCENT_PRIMARY};
}}

QCheckBox::indicator:checked {{
    background-color: {ACCENT_PRIMARY};
    border-color: {ACCENT_PRIMARY};
}}

/* =========================================== */
/* Palette Grid Buttons                        */
/* =========================================== */
#palette_grid QToolButton {{
    background-color: {PALETTE_BUTTON_BG};
    border: 1px solid {PALETTE_BUTTON_BORDER};
    border-radius: 4px;
    padding: 4px;
    min-width: 28px;
    min-height: 28px;
}}

#palette_grid QToolButton:hover {{
    background-color: {PALETTE_BUTTON_HOVER};
    border: 1px solid {PALETTE_SELECTED_BORDER};
}}

#palette_grid QToolButton:disabled {{
    color: {TEXT_MUTED};
    background-color: {BG_DOCK};
    border: 1px solid {BORDER_LIGHT};
}}
"""

# ============================================================================
# Vertical Tool Palette Stylesheet
# ============================================================================

TOOL_PALETTE_STYLESHEET = f"""
QToolBar {{
    background-color: {PALETTE_BG};
    border: none;
    border-right: 1px solid {BORDER_LIGHT};
    spacing: 4px;
    padding: 8px 6px;
}}

QToolBar::separator {{
    height: 1px;
    background-color: {BORDER_MEDIUM};
    margin: 8px 6px;
}}

QToolButton {{
    background-color: {PALETTE_BUTTON_BG};
    border: 1px solid {PALETTE_BUTTON_BORDER};
    border-radius: 6px;
    padding: 5px;
    min-width: 30px;
    min-height: 30px;
}}

QToolButton:hover {{
    background-color: {PALETTE_BUTTON_HOVER};
    border: 1px solid {ACCENT_PRIMARY};
}}

QToolButton:pressed {{
    background-color: {BORDER_LIGHT};
}}

QToolButton:checked {{
    background-color: {PALETTE_SELECTED_BG};
    border: 2px solid {ACCENT_PRIMARY};
}}

QToolButton[popupMode="1"] {{
    padding-right: 12px;
}}

QToolButton::menu-indicator {{
    subcontrol-position: right bottom;
    subcontrol-origin: padding;
    right: 2px;
    bottom: 2px;
}}

#palette_grid QToolButton {{
    background-color: {PALETTE_BUTTON_BG};
    border: 1px solid {PALETTE_BUTTON_BORDER};
    border-radius: 4px;
    padding: 4px;
    min-width: 28px;
    min-height: 28px;
}}

#palette_grid QToolButton:hover {{
    background-color: {PALETTE_BUTTON_HOVER};
    border: 1px solid {ACCENT_PRIMARY};
}}

#palette_grid QToolButton:disabled {{
    color: {TEXT_MUTED};
    background-color: {BG_DOCK};
    border: 1px solid {BORDER_LIGHT};
}}
"""
