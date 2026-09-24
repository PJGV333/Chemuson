"""Tests del sistema de iconos SVG de la UI.

OpenSpec: ``2026-09-24-modernize-ui-svg-icons`` (Fase 2 del plan de
modernización de la UI).

Cubre:
- inventario 1:1: cada SVG estático válido (parse XML, viewBox 24×24,
  ``currentColor``, sin raster) y paridad con los literals usados por
  ``toolbar.py`` y ``main_window_ui_builder.py``;
- ``IconProvider``: API dinámica (caché por key+params, HiDPI, fallo
  seguro), API estática de la Fase 1 inalterada;
- builders dinámicos puros (anillos, átomos CPK, cargas, electrones,
  radicales, esfera, diagramas de energía);
- fachada ``chemuson.gui.icons``: compatibilidad de nombres/firmas,
  retorno ``QIcon`` no nulo para todo el inventario, contrato de
  ``set_icon_theme`` y colores legacy → tokens;
- cambio de tema light→dark→light sin iconos con color del tema anterior
  y con colores químicos (CPK) independientes del tema;
- ventana real creada y refrescada en ambos temas sin excepciones.
"""

from __future__ import annotations

import inspect
import re
import xml.etree.ElementTree as ET
from pathlib import Path

from PyQt6.QtCore import QSize
from PyQt6.QtGui import QIcon
from PyQt6.QtWidgets import QApplication

import pytest

import chemuson.gui.icons as icons
from chemuson.gui.theme import icon_svg
from chemuson.gui.theme.icon_provider import DEFAULT_ICONS_DIR, IconProvider
from chemuson.gui.theme.tokens import get_tokens

REPO_ROOT = Path(__file__).resolve().parent.parent
TOOLBAR_SOURCE = REPO_ROOT / "src" / "chemuson" / "gui" / "toolbar.py"
UI_BUILDER_SOURCE = REPO_ROOT / "src" / "chemuson" / "gui" / "main_window_ui_builder.py"

STANDARD_SIZES = (16, 20, 24, 28)

#: Nombres públicos que la fachada debe conservar (API histórica).
HISTORICAL_PUBLIC_API = (
    "ICON_SIZE",
    "ATOM_COLORS",
    "set_icon_theme",
    "icon_foreground_color",
    "icon_muted_color",
    "icon_fill_color",
    "icon_paper_color",
    "draw_atom_icon",
    "draw_coordination_sphere_icon",
    "draw_glyph_icon",
    "draw_charge_icon",
    "draw_electron_icon",
    "draw_energy_diagram_icon",
    "draw_energy_levels_icon",
    "draw_molecular_orbital_icon",
    "draw_radical_charge_icon",
    "draw_bond_icon",
    "draw_wavy_anchor_icon",
    "draw_ring_icon",
    "draw_ring_template_icon",
    "draw_arrow_icon",
    "draw_generic_icon",
    "get_pointer_icon",
    "get_eraser_icon",
    "get_single_bond_icon",
    "get_double_bond_icon",
    "get_benzene_icon",
)

#: Símbolos estáticos no cubiertos por los tres mapas (energía, orbital, ancla).
EXTRA_STATIC_NAMES = ("energy-levels", "molecular-orbital", "wavy-anchor")


def _app() -> QApplication:
    app = QApplication.instance()
    assert app is not None, "la fixture de sesión debe crear la QApplication"
    return app


def _require_icon(icon: QIcon, label: str) -> None:
    assert isinstance(icon, QIcon), f"{label}: no devuelve QIcon"
    assert not icon.isNull(), f"{label}: QIcon nulo"


def _pixel_signature(icon: QIcon, size: int = 32) -> frozenset[int]:
    """Firma visual: conjunto de colores RGBA renderizados a ``size`` px."""
    image = icon.pixmap(QSize(size, size)).toImage()
    return frozenset(
        image.pixelColor(x, y).rgba()
        for x in range(image.width())
        for y in range(image.height())
    )


def _has_ink(icon: QIcon, size: int = 32) -> bool:
    """True si el icono renderiza algún píxel no transparente."""
    image = icon.pixmap(QSize(size, size)).toImage()
    for y in range(image.height()):
        for x in range(image.width()):
            if image.pixelColor(x, y).alpha() > 0:
                return True
    return False


def _all_inventory_icons() -> dict[str, QIcon]:
    """Un ``QIcon`` por cada nombre del inventario 1:1."""
    out: dict[str, QIcon] = {}
    for shape, name in icons.GENERIC_ICON_MAP.items():
        out[f"generic:{shape}"] = icons.draw_generic_icon(shape)
    for bond, name in icons.BOND_ICON_MAP.items():
        out[f"bond:{bond}"] = icons.draw_bond_icon(bond)
    for kind, name in icons.ARROW_ICON_MAP.items():
        out[f"arrow:{kind}"] = icons.draw_arrow_icon(kind)
    for extra in EXTRA_STATIC_NAMES:
        out[f"static:{extra}"] = icons._static(extra)
    return out


# ---------------------------------------------------------------------------
# 1. Inventario y validez del set SVG
# ---------------------------------------------------------------------------

class TestSvgInventory:
    def test_all_static_files_parse_with_valid_viewbox(self) -> None:
        files = sorted(DEFAULT_ICONS_DIR.glob("i-*.svg"))
        assert len(files) == 55, f"se esperan 55 SVG, hay {len(files)}"
        for path in files:
            root = ET.fromstring(path.read_text())
            # Con ``xmlns`` el tag lleva namespace: ``{http://...}svg``.
            assert root.tag.endswith("svg"), path.name
            assert root.get("viewBox") == "0 0 24 24", path.name

    def test_static_files_use_currentcolor_and_no_raster(self) -> None:
        for path in DEFAULT_ICONS_DIR.iterdir():
            if path.suffix == ".svg":
                assert "currentColor" in path.read_text(), path.name
        # Sin assets raster en la carpeta (solo SVG + LICENSE).
        rasters = [p.name for p in DEFAULT_ICONS_DIR.iterdir()
                   if p.suffix.lower() in {".png", ".jpg", ".jpeg", ".bmp", ".ico"}]
        assert not rasters, f"raster no permitido: {rasters}"
        assert (DEFAULT_ICONS_DIR / "LICENSE.txt").exists()

    def test_inventory_names_resolve_to_files(self) -> None:
        expected: set[str] = set()
        for mapping in (icons.GENERIC_ICON_MAP, icons.BOND_ICON_MAP, icons.ARROW_ICON_MAP):
            expected.update(mapping.values())
        expected.update(EXTRA_STATIC_NAMES)
        for name in expected:
            assert (DEFAULT_ICONS_DIR / f"i-{name}.svg").exists(), name

    def test_map_values_are_unique(self) -> None:
        """Dos APIs distintas no deben apuntar al mismo SVG por error."""
        for mapping in (icons.GENERIC_ICON_MAP, icons.BOND_ICON_MAP, icons.ARROW_ICON_MAP):
            assert len(set(mapping.values())) == len(mapping)


def _source_literals(path: Path, func: str) -> set[str]:
    return set(re.findall(rf'{func}\(\s*"([^"]+)"', path.read_text()))


class TestParityWithCallSites:
    """Ningún tool_id/forma usada por los callers reales queda huérfano."""

    def test_generic_shapes_used_by_toolbar_resolve(self) -> None:
        literals = (
            _source_literals(TOOLBAR_SOURCE, "draw_generic_icon")
            | _source_literals(UI_BUILDER_SOURCE, "draw_generic_icon")
        )
        assert literals, "el regex de paridad no encontró nada (¿cambió la fuente?)"
        for literal in sorted(literals):
            assert literal in icons.GENERIC_ICON_MAP, (
                f"draw_generic_icon({literal!r}) sin SVG en el inventario"
            )

    def test_bond_types_used_by_toolbar_resolve(self) -> None:
        for literal in _source_literals(TOOLBAR_SOURCE, "draw_bond_icon"):
            assert literal in icons.BOND_ICON_MAP, (
                f"draw_bond_icon({literal!r}) sin SVG en el inventario"
            )

    def test_arrow_kinds_used_by_toolbar_resolve(self) -> None:
        for literal in _source_literals(TOOLBAR_SOURCE, "draw_arrow_icon"):
            assert literal in icons.ARROW_ICON_MAP, (
                f"draw_arrow_icon({literal!r}) sin SVG en el inventario"
            )


# ---------------------------------------------------------------------------
# 2. IconProvider: API dinámica aditiva
# ---------------------------------------------------------------------------

class TestIconProviderDynamic:
    def test_icon_dynamic_non_null_and_cached(self) -> None:
        provider = IconProvider(dpr=1.0)
        first = provider.icon_dynamic("ring", "#111111", 24, sides=6, aromatic=True)
        second = provider.icon_dynamic("ring", "#111111", 24, sides=6, aromatic=True)
        _require_icon(first, "ring dinámico")
        assert first is second, "misma clave (key, params, color, size) → misma instancia"

    def test_icon_dynamic_params_produce_distinct_icons(self) -> None:
        provider = IconProvider(dpr=1.0)
        hex6 = provider.icon_dynamic("ring", "#111111", 24, sides=6, aromatic=True)
        hex7 = provider.icon_dynamic("ring", "#111111", 24, sides=7, aromatic=False)
        assert hex6 is not hex7
        assert _pixel_signature(hex6) != _pixel_signature(hex7)

    def test_icon_dynamic_unknown_key_fails_visibly(self) -> None:
        provider = IconProvider(dpr=1.0)
        assert provider.icon_dynamic("no-existe", "#111111", 24).isNull()
        assert provider.pixmap_dynamic("no-existe", "#111111", 24).isNull()

    def test_pixmap_dynamic_hidpi(self) -> None:
        provider = IconProvider(dpr=2.0)
        pixmap = provider.pixmap_dynamic("ring", "#111111", 24, sides=6, aromatic=True)
        assert not pixmap.isNull()
        assert pixmap.devicePixelRatio() == 2.0
        assert pixmap.width() == 48

    def test_static_api_unchanged_for_glyph_names(self) -> None:
        """La convención ``glyph:`` de la Fase 1 sigue funcionando igual."""
        provider = IconProvider(dpr=1.0)
        pixmap = provider.pixmap("glyph:C|circle|14", "#0F172A", 20)
        assert not pixmap.isNull() and pixmap.width() == 20
        first = provider.icon("glyph:C|none|12", "#0F172A", 20)
        assert first is provider.icon("glyph:C|none|12", "#0F172A", 20)
        assert first is not provider.icon("glyph:C|none|12", "#F1F5F9", 20)


# ---------------------------------------------------------------------------
# 3. Builders dinámicos (tests puros, sin Qt)
# ---------------------------------------------------------------------------

class TestDynamicBuilders:
    def test_ring_polygon_side_count_and_aromatic(self) -> None:
        for sides in (3, 4, 5, 6, 7, 8):
            svg = icon_svg.ring_svg(sides)
            points = re.search(r'<polygon points="([^"]+)"', svg).group(1)
            assert len(points.split()) == sides, sides
        assert "<circle" in icon_svg.ring_svg(6, aromatic=True)
        assert "<circle" not in icon_svg.ring_svg(6, aromatic=False)

    def test_ring_invalid_size_clamps_to_minimum(self) -> None:
        svg = icon_svg.ring_svg(2)
        points = re.search(r'<polygon points="([^"]+)"', svg).group(1)
        assert len(points.split()) == 3

    def test_atom_uses_cpk_fill_and_contrast_text(self) -> None:
        dark_bg = icon_svg.atom_svg("N", fill="#3050F8")
        assert 'fill="#3050F8"' in dark_bg
        assert '#FFFFFF' in dark_bg  # texto blanco sobre fondo oscuro
        light_bg = icon_svg.atom_svg("S", fill="#FFFF30")
        assert '#000000' in light_bg  # texto negro sobre fondo claro
        assert 'fill="#333333"' in icon_svg.atom_svg("X", fill="#333333")

    def test_atom_escapes_special_chars(self) -> None:
        svg = icon_svg.atom_svg("<A&B>", fill="#333333")
        assert "<A&B>" not in svg
        assert "&lt;A&amp;B&gt;" in svg

    def test_charge_plus_differs_from_minus(self) -> None:
        plus = icon_svg.charge_svg("+")
        minus = icon_svg.charge_svg("-")
        assert plus != minus
        assert plus.count("<line") == 2 and minus.count("<line") == 1

    def test_electrons_dot_count(self) -> None:
        assert icon_svg.electrons_svg(1).count("<circle") == 1
        assert icon_svg.electrons_svg(2).count("<circle") == 2
        assert icon_svg.electrons_svg(0).count("<circle") == 1  # clamp a 1

    def test_radical_dot_and_sign(self) -> None:
        svg = icon_svg.radical_svg("+")
        assert "<circle" in svg and svg.count("<line") == 2
        assert icon_svg.radical_svg("-").count("<line") == 1

    def test_ring_template_bold_edge_and_label(self) -> None:
        svg = icon_svg.ring_template_svg("OH", sides=6)
        assert 'stroke-width="3.2"' in svg and ">OH</text>" in svg

    def test_energy_boxes_count_and_arrow(self) -> None:
        svg = icon_svg.energy_boxes_svg(3, label="p", side="left")
        assert svg.count("<rect") == 3
        assert "<line" in svg and "<path" in svg  # flecha presente
        assert ">p</text>" in svg

    def test_energy_boxes_side_right_and_fill(self) -> None:
        svg = icon_svg.energy_boxes_svg(2, label="sp", side="right", fill="#ECFDF5")
        assert svg.count("<rect") == 2
        assert 'fill="#ECFDF5"' in svg and '>sp</text>' in svg

    def test_energy_boxes_many_fuse_into_band(self) -> None:
        svg = icon_svg.energy_boxes_svg(56)
        # Con N=56 las cajas se fusionan en una banda (1 rect), sin flecha.
        assert svg.count("<rect") == 1
        assert "<line" not in svg

    def test_energy_boxes_stroke_optional(self) -> None:
        svg = icon_svg.energy_boxes_svg(2, stroke=False)
        assert 'stroke="none"' in svg

    def test_sphere_gradient_and_fallback_color(self) -> None:
        svg = icon_svg.sphere_svg("#2E86AB")
        assert "radialGradient" in svg and '#2E86AB' in svg
        fallback = icon_svg.sphere_svg("no-es-hex")
        assert "#8D99A6" in fallback

    def test_glyph_builder_matches_provider_convention(self) -> None:
        svg = icon_svg.glyph_svg("C", font_size=12.5, shape="none")
        assert 'font-size="12.5"' in svg and ">C</text>" in svg
        circled = icon_svg.glyph_svg("C", shape="circle")
        assert "<circle" in circled


# ---------------------------------------------------------------------------
# 4. Fachada de compatibilidad
# ---------------------------------------------------------------------------

class TestFacadeCompatibility:
    def test_historical_public_api_exists(self) -> None:
        for name in HISTORICAL_PUBLIC_API:
            assert hasattr(icons, name), f"falta la API histórica {name}"

    def test_signatures_preserved(self) -> None:
        sig = inspect.signature
        assert list(sig(icons.draw_electron_icon).parameters) == ["count", "spread"]
        assert sig(icons.draw_electron_icon).parameters["count"].default == 1
        assert sig(icons.draw_electron_icon).parameters["spread"].default == 6.0
        assert sig(icons.draw_ring_icon).parameters["size"].default == 6
        assert sig(icons.draw_ring_icon).parameters["aromatic"].default is True
        assert sig(icons.draw_bond_icon).parameters["bond_type"].default == "single"
        assert sig(icons.draw_arrow_icon).parameters["kind"].default == "forward"
        assert sig(icons.draw_coordination_sphere_icon).parameters["color"].default == "#8D99A6"
        params = list(sig(icons.draw_atom_icon).parameters)
        assert params == ["text", "color"]
        assert sig(icons.draw_atom_icon).parameters["color"].default is None
        params = list(sig(icons.draw_energy_diagram_icon).parameters)
        assert params[0] == "box_count"
        for keyword in ("label_text", "label_side", "fill_color", "stroke_visible"):
            assert sig(icons.draw_energy_diagram_icon).parameters[keyword].kind is (
                inspect.Parameter.KEYWORD_ONLY
            )

    def test_constants_preserved(self) -> None:
        assert icons.ICON_SIZE == 32
        assert icons.ATOM_COLORS["C"] == "#333333"
        assert icons.ATOM_COLORS["N"] == "#3050F8"
        assert icons.ATOM_COLORS["O"] == "#FF0D0D"
        assert set(icons.ATOM_COLORS) == {
            "C", "N", "O", "S", "P", "F", "Cl", "Br", "H",
        }

    def test_all_inventory_non_null(self) -> None:
        _app()
        for label, icon in _all_inventory_icons().items():
            _require_icon(icon, label)
            assert _has_ink(icon), f"{label}: icono sin tinta"

    def test_convenience_getters_non_null(self) -> None:
        _app()
        for getter in (
            icons.get_pointer_icon,
            icons.get_eraser_icon,
            icons.get_single_bond_icon,
            icons.get_double_bond_icon,
            icons.get_benzene_icon,
        ):
            _require_icon(getter(), getter.__name__)

    def test_unknown_generic_returns_blank_non_null(self) -> None:
        _app()
        icon = icons.draw_generic_icon("forma-inexistente")
        _require_icon(icon, "fallback genérico")
        assert not _has_ink(icon), "fallback genérico debe ser en blanco (como antes)"
        bond = icons.draw_bond_icon("tipo-inexistente")
        _require_icon(bond, "fallback enlace")
        assert not _has_ink(bond)

    def test_unknown_arrow_falls_back_to_line(self) -> None:
        _app()
        icon = icons.draw_arrow_icon("kind-inexistente")
        _require_icon(icon, "fallback flecha")
        assert _has_ink(icon), "fallback de flecha debe ser la línea (como antes)"

    def test_set_icon_theme_contract(self) -> None:
        icons.set_icon_theme("dark")
        assert icons.icon_foreground_color() == str(get_tokens("dark")["icon"])
        icons.set_icon_theme("light")
        assert icons.icon_foreground_color() == str(get_tokens("light")["icon"])
        # Contrato histórico: cualquier nombre que no sea "dark" → light.
        icons.set_icon_theme("system")
        assert icons.icon_foreground_color() == str(get_tokens("light")["icon"])
        icons.set_icon_theme("nada")
        assert icons.icon_foreground_color() == str(get_tokens("light")["icon"])

    def test_legacy_color_helpers_map_to_tokens(self) -> None:
        for theme_name in ("light", "dark"):
            icons.set_icon_theme(theme_name)
            tokens = get_tokens(theme_name)
            assert icons.icon_foreground_color() == str(tokens["icon"])
            assert icons.icon_muted_color() == str(tokens["text3"])
            assert icons.icon_fill_color() == str(tokens["surface3"])
            assert icons.icon_paper_color() == str(tokens["surface"])
        icons.set_icon_theme("light")


# ---------------------------------------------------------------------------
# 5. Tinte theme-aware (caché por color)
# ---------------------------------------------------------------------------

class TestThemeAwareTint:
    def test_light_dark_light_static(self) -> None:
        _app()
        icons.set_icon_theme("light")
        light = _pixel_signature(icons.draw_generic_icon("pointer"))
        icons.set_icon_theme("dark")
        dark = _pixel_signature(icons.draw_generic_icon("pointer"))
        icons.set_icon_theme("light")
        light_again = _pixel_signature(icons.draw_generic_icon("pointer"))
        assert light != dark, "el tinte dark debe cambiar el icono"
        assert light == light_again, "volver a light reproduce el pixmap light"

    def test_light_dark_light_dynamic(self) -> None:
        _app()
        icons.set_icon_theme("light")
        light = _pixel_signature(icons.draw_ring_icon(6))
        icons.set_icon_theme("dark")
        dark = _pixel_signature(icons.draw_ring_icon(6))
        icons.set_icon_theme("light")
        assert light != dark
        assert light == _pixel_signature(icons.draw_ring_icon(6))

    def test_glyph_explicit_color_ignores_theme(self) -> None:
        _app()
        icons.set_icon_theme("light")
        a = _pixel_signature(icons.draw_glyph_icon("■", color="#DC2626"))
        icons.set_icon_theme("dark")
        b = _pixel_signature(icons.draw_glyph_icon("■", color="#DC2626"))
        assert a == b, "el color explícito es de dominio: no cambia con el tema"

    def test_atom_icons_theme_independent(self) -> None:
        _app()
        icons.set_icon_theme("light")
        light = _pixel_signature(icons.draw_atom_icon("N"))
        icons.set_icon_theme("dark")
        dark = _pixel_signature(icons.draw_atom_icon("N"))
        assert light == dark, "los colores CPK no siguen el tema UI"

    def test_no_cross_theme_cache_poisoning(self) -> None:
        """La clave de caché incluye el color: imposible servir light en dark."""
        provider = IconProvider(dpr=1.0)
        light_icon = provider.icon("pointer", str(get_tokens("light")["icon"]), 28)
        dark_icon = provider.icon("pointer", str(get_tokens("dark")["icon"]), 28)
        assert light_icon is not dark_icon
        assert _pixel_signature(light_icon) != _pixel_signature(dark_icon)


# ---------------------------------------------------------------------------
# 6. Ventana real
# ---------------------------------------------------------------------------

class TestRealWindow:
    @pytest.fixture()
    def window(self):
        from chemuson.gui.main_window import ChemusonWindow

        win = ChemusonWindow()
        yield win
        win.close()
        win.deleteLater()
        QApplication.instance().processEvents()

    def test_window_creation_and_theme_refresh(self, window) -> None:
        _app()
        for theme_name in ("light", "dark", "light"):
            window.current_theme = theme_name
            window._apply_theme()
            toolbar = getattr(window, "toolbar", None)
            if toolbar is not None and hasattr(toolbar, "refresh_icons"):
                toolbar.refresh_icons()
            text_toolbar = getattr(window, "text_toolbar", None)
            if text_toolbar is not None and hasattr(text_toolbar, "refresh_icons"):
                text_toolbar.refresh_icons()
            # Acciones principales con iconos no nulos.
            for action_name in (
                "action_new",
                "action_open",
                "action_save",
                "action_undo",
                "action_redo",
                "action_zoom_in",
            ):
                action = getattr(window, action_name, None)
                if action is not None:
                    assert not action.icon().isNull(), action_name
        # El puntero de la toolbar de dibujo sigue siendo no nulo tras 3 temas.
        select_action = getattr(window, "toolbar", None)
        if select_action is not None:
            act = getattr(select_action, "select_action", None)
            if act is not None:
                assert not act.icon().isNull()
