"""
Estilos visuales de la aplicación Chemuson.

Define paletas de color y hojas de estilo Qt para la ventana principal,
menús, barras de herramientas y docks.
"""

# ============================================================================
# Paleta moderna de colores
# ============================================================================

# Colores principales
PRIMARY_DARK = "#1E293B"        # Deep slate for headers
PRIMARY_MEDIUM = "#334155"      # Slightly lighter slate
ACCENT_PRIMARY = "#0891B2"      # Refined cyan accent
ACCENT_HOVER = "#06B6D4"        # Lighter cyan for hover
ACCENT_PRESSED = "#0E7490"      # Darker cyan for pressed

# Colores de fondo
BG_MAIN = "#F8FAFC"             # Clean light background
BG_ELEVATED = "#FFFFFF"         # White for elevated surfaces
BG_TOOLBAR = "#F1F5F9"          # Soft light gray for toolbars
BG_DOCK = "#F8FAFC"             # Almost white for dock panels

# Colores de bordes y separadores
BORDER_LIGHT = "#E2E8F0"        # Subtle border
BORDER_MEDIUM = "#CBD5E1"       # Medium border
BORDER_DARK = "#94A3B8"         # Darker border

# Colores de texto
TEXT_PRIMARY = "#0F172A"        # Deep slate for primary text
TEXT_SECONDARY = "#475569"      # Slate for secondary text
TEXT_MUTED = "#64748B"          # Muted text
TEXT_INVERSE = "#F8FAFC"        # Light text on dark backgrounds

# Paleta específica de herramientas
PALETTE_BG = "#F1F5F9"
PALETTE_BORDER = "#E2E8F0"
PALETTE_BUTTON_BG = "#FFFFFF"
PALETTE_BUTTON_BORDER = "#CBD5E1"
PALETTE_BUTTON_HOVER = "#ECFEFF"
PALETTE_SELECTED_BG = "#CFFAFE"
PALETTE_SELECTED_BORDER = "#0891B2"

# Sombra para efecto de profundidad (Qt no soporta box-shadow estándar)
SHADOW_COLOR = "rgba(0, 0, 0, 0.08)"

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
    padding: 4px 10px;
    spacing: 4px;
    font-size: 13px;
    font-weight: 500;
    min-height: 32px;
}}

QMenuBar::item {{
    background: transparent;
    padding: 6px 14px;
    border-radius: 6px;
    margin: 2px 2px;
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
    border-radius: 10px;
    padding: 8px 4px;
    margin: 4px;
}}

QMenu::item {{
    padding: 8px 28px 8px 14px;
    border-radius: 6px;
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
    background-color: {BG_TOOLBAR};
    border: none;
    border-bottom: 1px solid {BORDER_LIGHT};
    spacing: 6px;
    padding: 8px 10px;
}}

QToolBar::separator {{
    width: 1px;
    background-color: {BORDER_MEDIUM};
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
        stop:0 {BG_TOOLBAR}, stop:1 {BG_ELEVATED});
    padding: 12px 14px;
    border-bottom: 1px solid {BORDER_LIGHT};
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
    background-color: {PALETTE_BUTTON_HOVER};
    border-radius: 6px;
}}

/* =========================================== */
/* Status Bar                                  */
/* =========================================== */
QStatusBar {{
    background-color: {PRIMARY_DARK};
    color: {TEXT_INVERSE};
    border-top: none;
    padding: 6px 14px;
    font-size: 12px;
    min-height: 28px;
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
    width: 10px;
    border-radius: 5px;
    margin: 2px;
}}

QScrollBar::handle:vertical {{
    background-color: {BORDER_MEDIUM};
    min-height: 30px;
    border-radius: 4px;
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
    height: 0px;
    background: transparent;
}}

QScrollBar::add-page:vertical,
QScrollBar::sub-page:vertical {{
    background: transparent;
}}

QScrollBar:horizontal {{
    background-color: {BG_TOOLBAR};
    height: 10px;
    border-radius: 5px;
    margin: 2px;
}}

QScrollBar::handle:horizontal {{
    background-color: {BORDER_MEDIUM};
    min-width: 30px;
    border-radius: 4px;
    margin: 2px;
}}

QScrollBar::handle:horizontal:hover {{
    background-color: {TEXT_MUTED};
}}

QScrollBar::handle:horizontal:pressed {{
    background-color: {TEXT_SECONDARY};
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
    border-radius: 8px;
    gridline-color: {BORDER_LIGHT};
}}

QTableWidget::item {{
    padding: 6px 10px;
    color: {TEXT_PRIMARY};
}}

QTableWidget::item:selected {{
    background-color: {PALETTE_SELECTED_BG};
    color: {TEXT_PRIMARY};
}}

QTableWidget::item:hover {{
    background-color: {PALETTE_BUTTON_HOVER};
}}

QHeaderView::section {{
    background-color: {BG_TOOLBAR};
    color: {TEXT_SECONDARY};
    padding: 8px 10px;
    border: none;
    border-bottom: 2px solid {BORDER_MEDIUM};
    font-weight: 600;
    font-size: 12px;
}}

/* =========================================== */
/* Tree Widgets (Templates Dock)               */
/* =========================================== */
QTreeWidget {{
    background-color: {BG_ELEVATED};
    border: 1px solid {BORDER_LIGHT};
    border-radius: 8px;
}}

QTreeWidget::item {{
    padding: 8px 6px;
    border-radius: 6px;
    margin: 2px 4px;
}}

QTreeWidget::item:hover {{
    background-color: {PALETTE_BUTTON_HOVER};
}}

QTreeWidget::item:selected {{
    background-color: {PALETTE_SELECTED_BG};
    color: {TEXT_PRIMARY};
    border: 1px solid {PALETTE_SELECTED_BORDER};
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
    border-radius: 8px;
    padding: 10px 22px;
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
    color: {BG_TOOLBAR};
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
    border-radius: 8px;
    padding: 8px 12px;
    color: {TEXT_PRIMARY};
    selection-background-color: {PALETTE_SELECTED_BG};
}}

QLineEdit:focus {{
    border: 1px solid {ACCENT_PRIMARY};
    background-color: {BG_ELEVATED};
}}

QLineEdit:disabled {{
    background-color: {BG_DOCK};
    color: {TEXT_MUTED};
    border-color: {BORDER_LIGHT};
}}

/* =========================================== */
/* Combo Boxes                                 */
/* =========================================== */
QComboBox {{
    background-color: {BG_ELEVATED};
    border: 1px solid {BORDER_MEDIUM};
    border-radius: 8px;
    padding: 8px 12px;
    color: {TEXT_PRIMARY};
    min-width: 100px;
}}

QComboBox:hover {{
    border-color: {ACCENT_PRIMARY};
}}

QComboBox:focus {{
    border: 1px solid {ACCENT_PRIMARY};
}}

QComboBox::drop-down {{
    border: none;
    width: 24px;
}}

QComboBox QAbstractItemView {{
    background-color: {BG_ELEVATED};
    border: 1px solid {BORDER_MEDIUM};
    border-radius: 8px;
    selection-background-color: {ACCENT_PRIMARY};
    selection-color: {TEXT_INVERSE};
    padding: 4px;
}}

/* =========================================== */
/* Spin Boxes                                  */
/* =========================================== */
QSpinBox, QDoubleSpinBox {{
    background-color: {BG_ELEVATED};
    border: 1px solid {BORDER_MEDIUM};
    border-radius: 8px;
    padding: 8px 12px;
    color: {TEXT_PRIMARY};
}}

QSpinBox:focus, QDoubleSpinBox:focus {{
    border: 1px solid {ACCENT_PRIMARY};
}}

QSpinBox:disabled, QDoubleSpinBox:disabled {{
    background-color: {BG_DOCK};
    color: {TEXT_MUTED};
    border-color: {BORDER_LIGHT};
}}

/* =========================================== */
/* Checkboxes                                  */
/* =========================================== */
QCheckBox {{
    color: {TEXT_PRIMARY};
    spacing: 10px;
}}

QCheckBox::indicator {{
    width: 20px;
    height: 20px;
    border: 2px solid {BORDER_DARK};
    border-radius: 5px;
    background-color: {BG_ELEVATED};
}}

QCheckBox::indicator:hover {{
    border-color: {ACCENT_PRIMARY};
}}

QCheckBox::indicator:checked {{
    background-color: {ACCENT_PRIMARY};
    border-color: {ACCENT_PRIMARY};
}}

QCheckBox::indicator:disabled {{
    background-color: {BG_DOCK};
    border-color: {BORDER_LIGHT};
}}

/* =========================================== */
/* Radio Buttons                               */
/* =========================================== */
QRadioButton {{
    color: {TEXT_PRIMARY};
    spacing: 10px;
}}

QRadioButton::indicator {{
    width: 20px;
    height: 20px;
    border: 2px solid {BORDER_DARK};
    border-radius: 10px;
    background-color: {BG_ELEVATED};
}}

QRadioButton::indicator:hover {{
    border-color: {ACCENT_PRIMARY};
}}

QRadioButton::indicator:checked {{
    background-color: {ACCENT_PRIMARY};
    border-color: {ACCENT_PRIMARY};
    border-width: 5px;
}}

/* =========================================== */
/* Tab Widgets                                 */
/* =========================================== */
QTabWidget::pane {{
    border: 1px solid {BORDER_LIGHT};
    border-radius: 8px;
    background-color: {BG_ELEVATED};
}}

QTabBar::tab {{
    background-color: {BG_TOOLBAR};
    color: {TEXT_SECONDARY};
    padding: 10px 18px;
    border: 1px solid {BORDER_LIGHT};
    border-bottom: none;
    border-top-left-radius: 8px;
    border-top-right-radius: 8px;
    margin-right: 2px;
}}

QTabBar::tab:selected {{
    background-color: {BG_ELEVATED};
    color: {TEXT_PRIMARY};
    border-bottom: 2px solid {ACCENT_PRIMARY};
}}

QTabBar::tab:hover:!selected {{
    background-color: {PALETTE_BUTTON_HOVER};
    color: {TEXT_PRIMARY};
}}

/* =========================================== */
/* Group Boxes                                 */
/* =========================================== */
QGroupBox {{
    border: 1px solid {BORDER_LIGHT};
    border-radius: 8px;
    margin-top: 12px;
    padding-top: 16px;
    font-weight: 600;
    color: {TEXT_PRIMARY};
}}

QGroupBox::title {{
    subcontrol-origin: margin;
    subcontrol-position: top left;
    padding: 0 8px;
    color: {TEXT_SECONDARY};
}}

/* =========================================== */
/* ToolTip                                     */
/* =========================================== */
QToolTip {{
    background-color: {PRIMARY_DARK};
    color: {TEXT_INVERSE};
    border: 1px solid {PRIMARY_MEDIUM};
    border-radius: 6px;
    padding: 6px 10px;
    font-size: 12px;
}}

/* =========================================== */
/* Palette Grid Buttons                        */
/* =========================================== */
#palette_grid QToolButton {{
    background-color: {PALETTE_BUTTON_BG};
    border: 1px solid {PALETTE_BUTTON_BORDER};
    border-radius: 6px;
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
    spacing: 6px;
    padding: 10px 8px;
}}

QToolBar::separator {{
    height: 1px;
    background-color: {BORDER_MEDIUM};
    margin: 10px 8px;
}}

QToolButton {{
    background-color: {PALETTE_BUTTON_BG};
    border: 1px solid {PALETTE_BUTTON_BORDER};
    border-radius: 8px;
    padding: 8px;
    min-width: 32px;
    min-height: 32px;
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
    padding-right: 14px;
}}

QToolButton::menu-indicator {{
    subcontrol-position: right bottom;
    subcontrol-origin: padding;
    right: 3px;
    bottom: 3px;
}}

#palette_grid QToolButton {{
    background-color: {PALETTE_BUTTON_BG};
    border: 1px solid {PALETTE_BUTTON_BORDER};
    border-radius: 6px;
    padding: 6px;
    min-width: 30px;
    min-height: 30px;
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
