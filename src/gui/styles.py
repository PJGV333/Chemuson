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
PALETTE_BG = "#EFEFEF"
PALETTE_BORDER = "#BDBDBD"
PALETTE_BUTTON_BG = "#F8F8F8"
PALETTE_BUTTON_BORDER = "#C6C6C6"
PALETTE_BUTTON_HOVER = "#E3EDF9"
PALETTE_SELECTED_BG = "#C9DBF2"
PALETTE_SELECTED_BORDER = "#4A78B8"

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

#palette_grid QToolButton {{
    background-color: {PALETTE_BUTTON_BG};
    border: 1px solid {PALETTE_BUTTON_BORDER};
    border-radius: 2px;
    padding: 2px;
    min-width: 24px;
    min-height: 24px;
}}

#palette_grid QToolButton:hover {{
    background-color: {PALETTE_BUTTON_HOVER};
    border: 1px solid {PALETTE_SELECTED_BORDER};
}}

#palette_grid QToolButton:disabled {{
    color: {TEXT_MUTED};
    background-color: #F0F0F0;
    border: 1px solid {BORDER_COLOR};
}}
"""

# Vertical tool palette specific styles
TOOL_PALETTE_STYLESHEET = f"""
QToolBar {{
    background-color: {PALETTE_BG};
    border-right: 1px solid {PALETTE_BORDER};
    spacing: 2px;
    padding: 6px 4px;
}}

QToolBar::separator {{
    height: 1px;
    background-color: {PALETTE_BORDER};
    margin: 6px 4px;
}}

QToolButton {{
    background-color: {PALETTE_BUTTON_BG};
    border: 1px solid {PALETTE_BUTTON_BORDER};
    border-radius: 3px;
    padding: 4px;
    min-width: 26px;
    min-height: 26px;
}}

QToolButton:hover {{
    background-color: {PALETTE_BUTTON_HOVER};
    border: 1px solid {PALETTE_SELECTED_BORDER};
}}

QToolButton:checked {{
    background-color: {PALETTE_SELECTED_BG};
    border: 1px solid {PALETTE_SELECTED_BORDER};
}}

QToolButton[popupMode="1"] {{
    padding-right: 8px;
}}

#palette_grid QToolButton {{
    background-color: {PALETTE_BUTTON_BG};
    border: 1px solid {PALETTE_BUTTON_BORDER};
    border-radius: 2px;
    padding: 2px;
    min-width: 24px;
    min-height: 24px;
}}

#palette_grid QToolButton:hover {{
    background-color: {PALETTE_BUTTON_HOVER};
    border: 1px solid {PALETTE_SELECTED_BORDER};
}}

#palette_grid QToolButton:disabled {{
    color: {TEXT_MUTED};
    background-color: #F0F0F0;
    border: 1px solid {BORDER_COLOR};
}}
"""
