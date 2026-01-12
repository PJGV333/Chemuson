"""
Chemuson Application Styles
Centralized stylesheet definitions for a professional appearance.
"""

# Color palette
ACCENT_COLOR = "#6699CC"
ACCENT_HOVER = "#7FAAD4"
BACKGROUND_LIGHT = "#F7F7F7"
BACKGROUND_DARK = "#E0E0E0"
BORDER_COLOR = "#CCCCCC"
TEXT_COLOR = "#333333"
TEXT_MUTED = "#666666"

# Main application stylesheet
MAIN_STYLESHEET = f"""
QMainWindow {{
    background-color: {BACKGROUND_DARK};
}}

QMenuBar {{
    background-color: {BACKGROUND_LIGHT};
    border-bottom: 1px solid {BORDER_COLOR};
    padding: 2px 4px;
    font-size: 13px;
}}

QMenuBar::item {{
    background: transparent;
    padding: 4px 8px;
    border-radius: 4px;
}}

QMenuBar::item:selected {{
    background-color: {ACCENT_COLOR};
    color: white;
}}

QMenu {{
    background-color: white;
    border: 1px solid {BORDER_COLOR};
    border-radius: 4px;
    padding: 4px;
}}

QMenu::item {{
    padding: 6px 24px 6px 8px;
    border-radius: 3px;
}}

QMenu::item:selected {{
    background-color: {ACCENT_COLOR};
    color: white;
}}

QMenu::separator {{
    height: 1px;
    background-color: {BORDER_COLOR};
    margin: 4px 8px;
}}

QToolBar {{
    background-color: {BACKGROUND_LIGHT};
    border: none;
    spacing: 4px;
    padding: 4px;
}}

QToolBar::separator {{
    width: 1px;
    background-color: {BORDER_COLOR};
    margin: 4px 6px;
}}

QToolButton {{
    border-radius: 4px;
    padding: 4px;
    border: 1px solid transparent;
}}

QToolButton:hover {{
    background-color: #E8E8E8;
    border: 1px solid {BORDER_COLOR};
}}

QToolButton:pressed {{
    background-color: #D0D0D0;
}}

QToolButton:checked {{
    background-color: {ACCENT_COLOR};
    color: white;
    border: 1px solid {ACCENT_COLOR};
}}

QToolButton[popupMode="1"] {{
    padding-right: 16px;
}}

QToolButton::menu-indicator {{
    image: none;
    subcontrol-origin: padding;
    subcontrol-position: right center;
    width: 0px;
}}

QDockWidget {{
    titlebar-close-icon: none;
    titlebar-normal-icon: none;
    font-weight: bold;
}}

QDockWidget::title {{
    background-color: {BACKGROUND_LIGHT};
    padding: 6px;
    border-bottom: 1px solid {BORDER_COLOR};
}}

QStatusBar {{
    background-color: {BACKGROUND_LIGHT};
    border-top: 1px solid {BORDER_COLOR};
    color: {TEXT_MUTED};
    font-size: 12px;
}}

QScrollBar:vertical {{
    background-color: {BACKGROUND_LIGHT};
    width: 12px;
    border-radius: 6px;
}}

QScrollBar::handle:vertical {{
    background-color: {BORDER_COLOR};
    min-height: 20px;
    border-radius: 5px;
    margin: 2px;
}}

QScrollBar::handle:vertical:hover {{
    background-color: {TEXT_MUTED};
}}

QScrollBar::add-line:vertical,
QScrollBar::sub-line:vertical {{
    height: 0px;
}}
"""

# Vertical tool palette specific styles
TOOL_PALETTE_STYLESHEET = f"""
QToolBar {{
    background-color: {BACKGROUND_LIGHT};
    border-right: 1px solid {BORDER_COLOR};
    spacing: 2px;
    padding: 6px;
}}

QToolButton {{
    border-radius: 6px;
    padding: 6px;
    min-width: 32px;
    min-height: 32px;
}}

QToolButton:hover {{
    background-color: #E0E8F0;
    border: 1px solid {ACCENT_COLOR};
}}

QToolButton:checked {{
    background-color: {ACCENT_COLOR};
    color: white;
}}
"""
