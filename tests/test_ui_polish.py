"""Tests de la Fase 7 — Polish final de la UI moderna.

Cubre:
- Contraste AA del token claro ``text3`` (y ``text1``/``text2``) sobre los
  fondos claros.
- Estados ``:disabled`` inequívocos en el QSS.
- Tooltips uniformes del rail / app bar / pestañas (convención
  ``Nombre (Shortcut)``).
- Thumbnails de plantillas HiDPI (DPR-aware).
- Refresco de la ``CommandPalette`` tras mutación de plantillas.
- Onboarding nativo de primera ejecución (3 pasos, persistencia, no-repetir).

No toca química (Clean2D/ChemName), persistencia ``.cmsn`` ni orbital math.
"""
from __future__ import annotations

import pytest
from PyQt6.QtCore import QSettings, QStandardPaths, QSize
from PyQt6.QtWidgets import QApplication

from chemuson.gui.app_bar import _tooltip_from_action
from chemuson.gui.main_window import ChemusonWindow
from chemuson.gui.onboarding import OnboardingOverlay
from chemuson.gui.tool_rail import _tooltip_for_action


# ---------------------------------------------------------------------------
# Fixtures
# ---------------------------------------------------------------------------
@pytest.fixture(scope="module", autouse=True)
def _qapp() -> QApplication:
    return QApplication.instance() or QApplication([])


@pytest.fixture(autouse=True)
def _isolated_config_home(tmp_path, monkeypatch):
    """Isola ``XDG_CONFIG_HOME``/``QSettings`` por test (patrón de patrón)."""
    config_location = QStandardPaths.writableLocation(
        QStandardPaths.StandardLocation.ConfigLocation
    )
    monkeypatch.setenv("XDG_CONFIG_HOME", str(tmp_path))
    QSettings.setPath(
        QSettings.Format.NativeFormat, QSettings.Scope.UserScope, str(tmp_path)
    )
    try:
        yield
    finally:
        QSettings.setPath(
            QSettings.Format.NativeFormat, QSettings.Scope.UserScope, config_location
        )


@pytest.fixture
def win() -> ChemusonWindow:
    window = ChemusonWindow()
    window.resize(1440, 900)
    window.show()
    QApplication.processEvents()
    yield window
    window.close()
    QApplication.processEvents()


# ---------------------------------------------------------------------------
# 1. Contraste AA del tema claro (token ``text3``)
# ---------------------------------------------------------------------------
def _rel_luminance(hex_color: str) -> float:
    hex_color = hex_color.lstrip("#")
    r, g, b = (int(hex_color[i:i + 2], 16) / 255.0 for i in (0, 2, 4))

    def chan(c: float) -> float:
        return c / 12.92 if c <= 0.04045 else ((c + 0.055) / 1.055) ** 2.4

    return 0.2126 * chan(r) + 0.7152 * chan(g) + 0.0722 * chan(b)


def _wcag_contrast(fg: str, bg: str) -> float:
    l1, l2 = _rel_luminance(fg), _rel_luminance(bg)
    hi, lo = max(l1, l2), min(l1, l2)
    return (hi + 0.05) / (lo + 0.05)


def test_light_text3_meets_aa_contrast():
    """``text3`` claro debe ser ≥ 4.5:1 sobre los fondos claros (Fase 7)."""
    from chemuson.gui.theme.tokens import get_tokens

    light = get_tokens("light")
    text3 = str(light["text3"])
    for bg in ("#F1F5F9", "#FFFFFF", "#F8FAFC"):
        ratio = _wcag_contrast(text3, bg)
        assert ratio >= 4.5, f"text3 light {text3} sobre {bg}: {ratio:.2f} < 4.5"


def test_light_text1_and_text2_meet_aa():
    from chemuson.gui.theme.tokens import get_tokens

    light = get_tokens("light")
    for key in ("text1", "text2"):
        for bg in ("#F1F5F9", "#FFFFFF", "#F8FAFC"):
            ratio = _wcag_contrast(str(light[key]), bg)
            assert ratio >= 4.5, f"{key} light {light[key]} sobre {bg}: {ratio:.2f} < 4.5"


def test_dark_text3_meets_aa():
    from chemuson.gui.theme.tokens import get_tokens

    dark = get_tokens("dark")
    text3 = str(dark["text3"])
    for bg in ("#0B1120", "#0F172A", "#1E293B"):
        ratio = _wcag_contrast(text3, bg)
        assert ratio >= 4.5, f"text3 dark {text3} sobre {bg}: {ratio:.2f} < 4.5"


# ---------------------------------------------------------------------------
# 2. Estados :disabled inequívocos en el QSS
# ---------------------------------------------------------------------------
def test_qss_disabled_states_have_opacity():
    from chemuson.gui.theme.qss import get_main_stylesheet

    qss = get_main_stylesheet("light")
    for fragment in (
        "QToolButton:disabled",
        "QPushButton:disabled",
        "QLineEdit:disabled",
        "QCheckBox:disabled",
        "QRadioButton:disabled",
    ):
        idx = qss.find(fragment)
        assert idx != -1, f"falta el estado :disabled para {fragment}"
        block = qss[idx: idx + 400]
        assert "opacity" in block, f"sin opacidad en {fragment}"
        opacity = float(block.split("opacity:")[1].split(";")[0].strip())
        assert 0.0 < opacity <= 0.75, f"opacidad {opacity} fuera de rango en {fragment}"


def test_qss_disabled_applies_to_rail_and_flyout():
    from chemuson.gui.theme.qss import get_main_stylesheet

    qss = get_main_stylesheet("dark")
    for fragment in ("QToolButton#railBtn:disabled", "QFrame[cls=\"flyItem\"]:disabled"):
        assert fragment in qss, f"falta el estado :disabled para {fragment}"


# ---------------------------------------------------------------------------
# 3. Tooltips uniformes (convención ``Nombre (Shortcut)``)
# ---------------------------------------------------------------------------
def test_rail_tooltips_follow_action_convention(win):
    rail = win.tool_rail
    for spec in rail._specs:
        button = rail._buttons[spec.key]
        expected = _tooltip_for_action(spec.tooltip, spec.icon_action)
        assert button.toolTip() == expected, (
            f"rail {spec.key}: {button.toolTip()!r} != esperado {expected!r}"
        )
        # Los atajos de tecla única del rail (V, A, L...) viven en el texto
        # base, no en ``QAction.shortcut()``; si alguna acción del rail sí
        # declara un atajo, la convención debe reflejarlo.
        shortcut = (
            spec.icon_action.shortcut().toString() if spec.icon_action else ""
        )
        if shortcut:
            assert shortcut in button.toolTip()


def test_app_bar_tooltips_follow_action_convention(win):
    bar = win.app_bar
    for button in (bar.undo_button, bar.redo_button, bar.theme_button, bar.preferences_button):
        action = button.defaultAction()
        assert action is not None
        assert button.toolTip() == _tooltip_from_action(action)
    # La píldora de búsqueda conserva su tooltip explícito con el atajo.
    assert bar.search_pill.toolTip() == "Buscar o ejecutar… (Ctrl+P)"


def test_document_tabs_new_button_tooltip(win):
    tabs = win.app_bar.tab_bar
    action = win.action_new
    assert action is not None
    # ``set_new_action`` fija el tooltip descriptivo con el atajo real.
    assert tabs.new_button.toolTip() == "Nuevo documento (Ctrl+N)"
    # El atajo se deriva de la ``QAction`` real (Ctrl+N).
    assert action.shortcut().toString() == "Ctrl+N"


# ---------------------------------------------------------------------------
# 4. Thumbnails de plantillas HiDPI (DPR-aware)
# ---------------------------------------------------------------------------
_WATER_MOLBLOCK = (
    "  Water\n"
    "  ChemUSON\n"
    "\n"
    "  3  2  0  0  0  0  0  0  0  0999 V2000\n"
    "    0.0000    0.0000    0.0000 O   0  0\n"
    "    0.7570    0.5860    0.0000 H   0  0\n"
    "   -0.7570    0.5860    0.0000 H   0  0\n"
    "  1  2  1  0\n"
    "  1  3  1  0\n"
    "M  END\n"
    "$$$$\n"
)


def test_template_preview_is_dpr_aware(tmp_path, monkeypatch):
    from types import SimpleNamespace

    from chemuson.gui.template_browser_service import TemplateBrowserService
    from chemuson.gui.template_library import DEFAULT_CATEGORY_USER, TemplateLibrary

    library = TemplateLibrary(tmp_path / "library.json")
    template = library.add_template("Agua", DEFAULT_CATEGORY_USER, _WATER_MOLBLOCK)
    template_id = template["id"]

    # Dos contextos, cada uno con su propia ``preview_cache``: aíslan el
    # efecto del devicePixelRatio sin interferencia de la caché.
    ctx_dpr2 = SimpleNamespace(
        template_library=library,
        preview_cache={},
        show_status=lambda _s: None,
    )
    ctx_dpr1 = SimpleNamespace(
        template_library=library,
        preview_cache={},
        show_status=lambda _s: None,
    )

    # Fuerza un DPR de 2.0: el pixmap físico almacenado es (88×56) × 2 = 176×112.
    monkeypatch.setattr(
        TemplateBrowserService, "_device_pixel_ratio", staticmethod(lambda: 2.0)
    )
    icon_2 = TemplateBrowserService().template_preview_icon(ctx_dpr2, template_id)
    assert icon_2.availableSizes() == [QSize(176, 112)]

    # Con DPR 1.0 (por defecto) el pixmap físico es 88×56.
    monkeypatch.setattr(
        TemplateBrowserService, "_device_pixel_ratio", staticmethod(lambda: 1.0)
    )
    icon_1 = TemplateBrowserService().template_preview_icon(ctx_dpr1, template_id)
    assert icon_1.availableSizes() == [QSize(88, 56)]


def test_template_preview_device_pixel_ratio_falls_back_to_one():
    from chemuson.gui.template_browser_service import TemplateBrowserService

    # Con QApplication activa (offscreen, DPR 1.0) debe devolver 1.0.
    assert TemplateBrowserService._device_pixel_ratio() == pytest.approx(1.0)


# ---------------------------------------------------------------------------
# 5. Refresco de la CommandPalette tras mutación de plantillas
# ---------------------------------------------------------------------------
def test_palette_registry_resyncs_after_template_addition(win):
    from chemuson.gui.template_library import DEFAULT_CATEGORY_USER

    library = win.template_library
    before_titles = [
        e.title for e in win._command_registry.entries()
        if e.section == "Plantillas"
    ]

    library.add_template("Metano F7", DEFAULT_CATEGORY_USER, _WATER_MOLBLOCK)

    # El refresco de vistas de plantilla debe resincronizar el registro.
    win._refresh_template_views()

    titles = [e.title for e in win._command_registry.entries() if e.section == "Plantillas"]
    assert "Metano F7" in titles, f"la paleta no expone la plantilla nueva: {titles}"
    # La nueva plantilla no estaba presente antes del refresco.
    assert "Metano F7" not in before_titles


def test_refresh_template_views_calls_registry_resync(win, monkeypatch):
    """El wiring: ``_refresh_template_views`` invoca
    ``_refresh_command_registry_templates``."""
    calls = []
    monkeypatch.setattr(
        win, "_refresh_command_registry_templates", lambda: calls.append(1)
    )
    win._refresh_template_views()
    assert calls == [1]


# ---------------------------------------------------------------------------
# 6. Onboarding nativo de primera ejecución
# ---------------------------------------------------------------------------
def _find_onboarding(win: ChemusonWindow) -> OnboardingOverlay | None:
    found = win.findChildren(OnboardingOverlay)
    return found[0] if found else None


def test_onboarding_shows_on_first_run(win):
    overlay = _find_onboarding(win)
    assert overlay is not None, "el onboarding debe mostrarse en la primera ejecución"
    assert overlay.isVisible()
    assert overlay.step == 0


def test_onboarding_not_shown_when_completed(_qapp, _isolated_config_home):
    from chemuson.platform.settings import application_settings

    # Se debe fijar ANTES de construir la ventana (la lectura es en el montaje).
    application_settings().setValue("ui/onboarding/completed", True)
    window = ChemusonWindow()
    window.resize(1440, 900)
    window.show()
    QApplication.processEvents()
    assert _find_onboarding(window) is None, "no debe repetirse si ya se completó"
    window.close()
    QApplication.processEvents()


def test_onboarding_three_step_navigation(win):
    overlay = _find_onboarding(win)
    assert overlay is not None
    assert overlay.step == 0
    overlay.go_back()
    assert overlay.step == 0  # en el primer paso no retrocede
    overlay.advance()
    assert overlay.step == 1
    overlay.advance()
    assert overlay.step == 2
    overlay.go_back()
    assert overlay.step == 1


def test_onboarding_persists_completed_on_finish(win):
    from chemuson.platform.settings import application_settings, setting_bool

    overlay = _find_onboarding(win)
    assert overlay is not None
    finished = []
    overlay.finished.connect(lambda: finished.append(True))

    # Completa los 3 pasos (el último ``advance`` cierra).
    overlay.advance()
    overlay.advance()
    overlay.advance()
    QApplication.processEvents()

    assert finished, "la señal ``finished`` debe emitirse al completar"
    assert overlay.isHidden()
    assert setting_bool(
        application_settings().value("ui/onboarding/completed", False), False
    ), "debe persistir ``ui/onboarding/completed`` = True"


def test_onboarding_no_more_flag_default_false(win):
    overlay = _find_onboarding(win)
    assert overlay is not None
    assert overlay.no_more is False
