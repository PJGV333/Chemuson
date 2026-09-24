"""Tests de la fundación de temas de la UI.

OpenSpec: ``2026-09-24-modernize-ui-theme-foundation`` (Fase 1 del plan de
modernización de la UI).

Cubre:
- resolución de tokens light/dark y validez de colores;
- métricas UI centralizadas;
- resolución de nombres de tema (incl. ``system`` y valores inválidos);
- aplicación de ambos temas sin excepciones (QSS + QPalette coherentes);
- compatibilidad de la fachada ``chemuson.gui.styles``;
- ventana real (``ChemusonWindow``) creada y conmutada entre temas sin
  excepciones;
- ``IconProvider``: caché por clave, HiDPI, glifos dinámicos y fallo
  visible sin excepción.
"""

from __future__ import annotations

from PyQt6.QtGui import QColor
from PyQt6.QtWidgets import QApplication, QMainWindow

import pytest

from chemuson.gui import theme
from chemuson.gui.styles import (
    DARK_COLORS,
    DEFAULT_THEME,
    LIGHT_COLORS,
    MAIN_STYLESHEET,
    TOOL_PALETTE_STYLESHEET,
    get_main_stylesheet,
    get_tool_palette_stylesheet,
)
from chemuson.gui.theme.icon_provider import IconProvider
from chemuson.gui.theme.tokens import (
    DEFAULT_THEME_NAME,
    METRICS,
    THEME_NAMES,
    get_tokens,
    theme_color,
)

# Tokens requeridos por la especificación (mínimo).
REQUIRED_COLOR_TOKENS = (
    "bg",
    "surface",
    "surface2",
    "surface3",
    "border",
    "borderStrong",
    "text1",
    "text2",
    "text3",
    "accent",
    "accentStrong",
    "accentHover",
    "accentSoft",
    "accentBorder",
    "onAccent",
    "danger",
    "dangerSoft",
    "warn",
    "warnSoft",
    "ok",
    "okSoft",
    "sheet",
    "sheetGrid",
    "canvasBg",
)

LEGACY_COLOR_KEYS = (
    "primary_dark",
    "primary_medium",
    "accent_primary",
    "accent_hover",
    "accent_pressed",
    "bg_main",
    "bg_elevated",
    "bg_toolbar",
    "bg_dock",
    "border_light",
    "border_medium",
    "border_dark",
    "text_primary",
    "text_secondary",
    "text_muted",
    "text_inverse",
    "palette_bg",
    "palette_border",
    "palette_button_bg",
    "palette_button_border",
    "palette_button_hover",
    "palette_selected_bg",
    "palette_selected_border",
)


def _app() -> QApplication:
    app = QApplication.instance()
    assert app is not None
    return app


# ---------------------------------------------------------------------------
# Tokens
# ---------------------------------------------------------------------------


class TestTokens:
    def test_theme_names_and_default(self) -> None:
        assert THEME_NAMES == ("light", "dark")
        assert DEFAULT_THEME_NAME == "light"
        assert theme.resolve_theme_name("garbage") == "light"

    def test_tokens_resolve_light_and_dark(self) -> None:
        light = get_tokens("light")
        dark = get_tokens("dark")
        for required in REQUIRED_COLOR_TOKENS:
            assert required in light, f"falta token {required!r} en light"
            assert required in dark, f"falta token {required!r} en dark"
        # Valores de referencia (spike aprobado).
        assert light["bg"] == "#F1F5F9"
        assert dark["bg"] == "#0B1120"
        assert light["accent"] == "#0E7490"
        assert dark["accent"] == "#22D3EE"
        # Los dos temas se distinguen.
        assert light["bg"] != dark["bg"]

    def test_unknown_theme_falls_back_to_default(self) -> None:
        assert get_tokens("nope") is get_tokens("light")

    def test_token_colors_are_valid_qcolors(self) -> None:
        for name in THEME_NAMES:
            tokens = get_tokens(name)
            for key in REQUIRED_COLOR_TOKENS:
                color = QColor(str(tokens[key]))
                assert color.isValid(), f"token {key} ({name}) no es color válido"

    def test_theme_color_helper_returns_qcolor(self) -> None:
        color = theme_color("dark", "accent")
        assert isinstance(color, QColor)
        assert color.isValid()
        assert color.name().lower() == "#22d3ee"

    def test_metrics_present(self) -> None:
        for key in (
            "spacingXs",
            "spacingSm",
            "spacingMd",
            "spacingLg",
            "radiusSurface",
            "radiusBtn",
            "radiusChip",
            "fontBase",
            "fontSmall",
            "fontTiny",
        ):
            assert isinstance(METRICS[key], int), f"falta métrica {key}"
        assert METRICS["spacingSm"] == 8  # grilla de 8 px
        assert METRICS["fontBase"] == 13


# ---------------------------------------------------------------------------
# Resolución de nombres
# ---------------------------------------------------------------------------


class TestResolveThemeName:
    def test_light_and_dark_normalize(self) -> None:
        assert theme.resolve_theme_name("light") == "light"
        assert theme.resolve_theme_name("DARK") == "dark"
        assert theme.resolve_theme_name(" dark ") == "dark"

    def test_unknown_name_never_raises(self) -> None:
        assert theme.resolve_theme_name("neon") == "light"
        assert theme.resolve_theme_name("") == "light"

    def test_system_resolves_to_light_or_dark(self) -> None:
        resolved = theme.resolve_theme_name("system")
        assert resolved in THEME_NAMES
        assert resolved == theme.system_theme_name()


# ---------------------------------------------------------------------------
# Aplicación del tema
# ---------------------------------------------------------------------------


@pytest.fixture
def themed_window():
    """QMainWindow con restauración del estado global de la app al terminar."""
    app = _app()
    window = QMainWindow()
    original_stylesheet = app.styleSheet()
    original_palette = app.palette()
    original_font = app.font()
    yield window
    app.setStyleSheet(original_stylesheet)
    app.setPalette(original_palette)
    app.setFont(original_font)
    window.close()
    window.deleteLater()
    QApplication.processEvents()


class TestApplyTheme:
    def test_apply_light_and_dark_no_exceptions(self, themed_window: QMainWindow) -> None:
        applied = theme.apply_theme(themed_window, "light")
        assert applied == "light"
        assert "#F1F5F9" in themed_window.styleSheet()
        assert (
            themed_window.palette().window().color().name().lower()
            == get_tokens("light")["bg"].lower()
        )

        applied = theme.apply_theme(themed_window, "dark")
        assert applied == "dark"
        assert "#0B1120" in themed_window.styleSheet()
        assert (
            themed_window.palette().window().color().name().lower()
            == get_tokens("dark")["bg"].lower()
        )

    def test_apply_unknown_theme_resolves_light(self, themed_window: QMainWindow) -> None:
        applied = theme.apply_theme(themed_window, "neon")
        assert applied == "light"
        assert "#F1F5F9" in themed_window.styleSheet()

    def test_set_theme_from_system_applies_resolved_theme(
        self, themed_window: QMainWindow
    ) -> None:
        applied = theme.set_theme_from_system(themed_window)
        assert applied in THEME_NAMES
        assert themed_window.styleSheet()

    def test_apply_theme_on_application(self) -> None:
        app = _app()
        original_stylesheet = app.styleSheet()
        original_palette = app.palette()
        original_font = app.font()
        try:
            assert theme.apply_theme(app, "dark") == "dark"
            assert (
                app.palette().window().color().name().lower()
                == get_tokens("dark")["bg"].lower()
            )
        finally:
            app.setStyleSheet(original_stylesheet)
            app.setPalette(original_palette)
            app.setFont(original_font)


class TestStylesheetTokens:
    def test_dark_stylesheet_contains_dark_tokens(self) -> None:
        qss = theme.get_main_stylesheet("dark")
        assert get_tokens("dark")["bg"] in qss
        assert get_tokens("dark")["surface"] in qss
        assert get_tokens("dark")["accent"] in qss
        assert get_tokens("dark")["text1"] in qss

    def test_light_stylesheet_contains_light_tokens(self) -> None:
        qss = theme.get_main_stylesheet("light")
        assert get_tokens("light")["bg"] in qss
        assert get_tokens("light")["surface"] in qss
        assert get_tokens("light")["accent"] in qss

    def test_widget_states_use_tokens(self) -> None:
        qss = theme.get_main_stylesheet("light")
        t = get_tokens("light")
        # hover/pressed/checked/disabled con tokens de superficie y acento.
        assert f"background-color: {t['surface2']}" in qss
        assert f"background-color: {t['surface3']}" in qss
        assert f"background-color: {t['accentSoft']}" in qss
        assert f"border: 1px solid {t['accentBorder']}" in qss
        assert f"background-color: {t['accentHover']}" in qss
        # disabled con tokens atenuados.
        assert f"color: {t['text3']}" in qss

    def test_tool_palette_stylesheet_uses_tokens(self) -> None:
        qss = theme.get_tool_palette_stylesheet("dark")
        t = get_tokens("dark")
        assert f"background-color: {t['surface']}" in qss
        assert f"border-right: 1px solid {t['border']}" in qss


# ---------------------------------------------------------------------------
# Fachada styles.py (compatibilidad)
# ---------------------------------------------------------------------------


class TestStylesFacade:
    def test_generators_are_equivalent_to_theme_system(self) -> None:
        assert get_main_stylesheet("light") == theme.get_main_stylesheet("light")
        assert get_main_stylesheet("dark") == theme.get_main_stylesheet("dark")
        assert (
            get_tool_palette_stylesheet("dark")
            == theme.get_tool_palette_stylesheet("dark")
        )

    def test_module_constants(self) -> None:
        assert DEFAULT_THEME == "light"
        assert isinstance(MAIN_STYLESHEET, str) and len(MAIN_STYLESHEET) > 1000
        assert isinstance(TOOL_PALETTE_STYLESHEET, str) and len(TOOL_PALETTE_STYLESHEET) > 100
        assert MAIN_STYLESHEET == theme.get_main_stylesheet("light")
        assert TOOL_PALETTE_STYLESHEET == theme.get_tool_palette_stylesheet("light")

    def test_legacy_palettes_are_token_aliases(self) -> None:
        for key in LEGACY_COLOR_KEYS:
            assert key in LIGHT_COLORS, f"falta clave legada {key!r}"
            assert key in DARK_COLORS, f"falta clave legada {key!r}"
        # Alias de tokens: valor legado == valor del token al que mapea.
        assert LIGHT_COLORS["accent_primary"] == get_tokens("light")["accent"]
        assert DARK_COLORS["accent_primary"] == get_tokens("dark")["accent"]
        assert LIGHT_COLORS["bg_main"] == get_tokens("light")["bg"]
        assert DARK_COLORS["text_primary"] == get_tokens("dark")["text1"]


# ---------------------------------------------------------------------------
# Ventana real
# ---------------------------------------------------------------------------


class TestRealWindow:
    def test_window_creation_and_theme_toggle_no_exceptions(self) -> None:
        from chemuson.gui.main_window import ChemusonWindow

        window = ChemusonWindow()
        try:
            assert window.current_theme in THEME_NAMES
            # El tema se aplica al ensamblar (sin excepciones).
            assert window.styleSheet()

            first = window.current_theme
            window.toggle_theme()
            assert window.current_theme != first
            window.toggle_theme()
            assert window.current_theme == first
            # QSS de toolbars aplicado.
            assert window.toolbar.styleSheet()
            assert window.symbols_toolbar.styleSheet()
        finally:
            window.close()
            window.deleteLater()
            QApplication.processEvents()

    def test_window_loads_persisted_dark_theme(self, tmp_path, monkeypatch) -> None:
        from chemuson.gui.main_window import ChemusonWindow

        # Aísla QSettings en tmp (no tocar la configuración del usuario).
        monkeypatch.setenv("HOME", str(tmp_path))
        from chemuson.platform.settings import application_settings

        application_settings().setValue("ui/theme", "dark")

        window = ChemusonWindow()
        try:
            assert window.current_theme == "dark"
            assert "#0B1120" in window.styleSheet()
        finally:
            window.close()
            window.deleteLater()
            QApplication.processEvents()


# ---------------------------------------------------------------------------
# IconProvider
# ---------------------------------------------------------------------------


class TestIconProvider:
    def test_cache_returns_same_instance_per_key(self) -> None:
        provider = IconProvider(dpr=1.0)
        first = provider.icon("glyph:C|none|12", "#0F172A", 20)
        second = provider.icon("glyph:C|none|12", "#0F172A", 20)
        assert first is second
        # Diferente color (otro tema/estado) → instancia distinta.
        other = provider.icon("glyph:C|none|12", "#F1F5F9", 20)
        assert other is not first

    def test_glyph_pixmap_renders(self) -> None:
        provider = IconProvider(dpr=1.0)
        pixmap = provider.pixmap("glyph:C|circle|14", "#0F172A", 20)
        assert not pixmap.isNull()
        assert pixmap.width() == 20
        assert pixmap.height() == 20
        # El tinte se aplica: hay píxeles con tinta además del transparente.
        image = pixmap.toImage()
        distinct = {
            image.pixelColor(x, y).rgba()
            for x in range(image.width())
            for y in range(image.height())
        }
        assert len(distinct) > 1

    def test_hidpi_device_pixel_ratio(self) -> None:
        provider = IconProvider(dpr=2.0)
        pixmap = provider.pixmap("glyph:N|none|12", "#111111", 20)
        assert pixmap.devicePixelRatio() == 2.0
        # 20 px lógicos a dpr 2 → 40 px físicos.
        assert pixmap.width() == 40

    def test_missing_icon_fails_visibly_without_exception(self) -> None:
        provider = IconProvider(dpr=1.0)
        icon = provider.icon("definitely-missing-icon-xyz", "#000000", 20)
        assert icon.isNull()
        pixmap = provider.pixmap("definitely-missing-icon-xyz", "#000000", 20)
        assert pixmap.isNull()

    def test_theme_color_helper_maps_roles(self) -> None:
        light_icon = IconProvider.theme_color("light", "icon")
        dark_icon = IconProvider.theme_color("dark", "icon")
        assert light_icon.lower() == str(get_tokens("light")["icon"]).lower()
        assert dark_icon.lower() == str(get_tokens("dark")["icon"]).lower()

    def test_clear_cache_invalidates(self) -> None:
        provider = IconProvider(dpr=1.0)
        first = provider.icon("glyph:O|none|12", "#0F172A", 16)
        provider.clear_cache()
        second = provider.icon("glyph:O|none|12", "#0F172A", 16)
        assert first is not second
