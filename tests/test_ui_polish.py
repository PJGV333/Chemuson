"""Tests de la Fase 7 — Polish final de la UI moderna.

Cubre:
- Contraste AA del token claro ``text3`` (y ``text1``/``text2``) sobre los
  fondos claros.
- Estados ``:disabled`` inequívocos en el QSS usando solo propiedades
  soportadas (sin ``opacity``, que Qt solo aplica a ``QToolTip``).
- Tooltips uniformes del rail / app bar / pestañas (convención
  ``Nombre (Shortcut)``).
- Thumbnails de plantillas HiDPI (DPR-aware) y ``iconSize`` lógico del árbol
  de PlantillasDock (88×56, sin recorte en el SidePanel de 340 px), con el
  bounding box relativo de la estructura invariante entre DPR 1 y DPR 2.
- Refresco de la ``CommandPalette`` tras mutación de plantillas.
- Onboarding nativo de primera ejecución (3 pasos, semántica de
  persistencia de "No volver a mostrar", no-repetir).

No toca química (Clean2D/ChemName), persistencia ``.cmsn`` ni orbital math.
"""
from __future__ import annotations

import pytest
import re
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
# 2. Estados :disabled inequívocos en el QSS (solo propiedades soportadas)
# ---------------------------------------------------------------------------
# Qt Style Sheets solo soporta ``opacity`` para ``QToolTip``; no funciona en
# ``QToolButton``, ``QPushButton``, ``QLineEdit``, ``QCheckBox``, etc. Por eso
# el QSS no usa ``opacity`` en widgets normales: el estado deshabilitado se
# pinta solo con propiedades soportadas (color, background-color, border-color
# y equivalentes). Estos tests comprueban (a) que NINGÚN selector disabled usa
# ``opacity`` y (b) que cada selector disabled introduce un cambio visual
# soportado respecto al estado normal.
_SUPPORTED_VISUAL_PROPS = {
    "color",
    "background-color",
    "border-color",
    "background",
    "border",
}


def _parse_qss_blocks(qss: str) -> dict[str, dict[str, str]]:
    """Extrae ``selector -> {prop: valor}`` de una hoja QSS.

    Parser mínimo y determinista: localiza cada regla ``selector { props }``
    (sin llaves anidadas) con una regex y descompone las propiedades por
    ``;``. Los selectores múltiples se indexan por cada parte (separada por
    coma). Los comentarios se eliminan antes.
    """
    text = re.sub(r"/\*.*?\*/", "", qss, flags=re.DOTALL)
    blocks: dict[str, dict[str, str]] = {}
    for match in re.finditer(r"([^{}]+)\{([^{}]*)\}", text):
        selector, body = match.group(1).strip(), match.group(2)
        props: dict[str, str] = {}
        for line in body.split(";"):
            if ":" in line:
                prop, _, value = line.partition(":")
                prop, value = prop.strip(), value.strip()
                if prop and value:
                    props[prop] = value
        if not props:
            continue
        for part in selector.split(","):
            part = part.strip()
            if part:
                blocks[part] = {**blocks.get(part, {}), **props}
    return blocks


def _disabled_visual_change(
    blocks: dict[str, dict[str, str]], base_sel: str, disabled_sel: str
) -> list[str]:
    """Devuelve las props soportadas cuyo valor cambia de normal a disabled."""
    base = blocks.get(base_sel, {})
    disabled = blocks.get(disabled_sel, {})
    changed = []
    for prop in _SUPPORTED_VISUAL_PROPS:
        if prop in disabled:
            base_val = base.get(prop)
            if base_val is None or base_val != disabled[prop]:
                changed.append(f"{prop}: {base_val!r} -> {disabled[prop]!r}")
    return changed


# (selector normal, selector :disabled) a auditar. El par ``QToolButton`` usa
# la regla base del rail de app (el primero en la hoja, sin objeto id).
_DISABLED_SELECTOR_PAIRS = (
    ("QToolButton", "QToolButton:disabled"),
    ("QPushButton", "QPushButton:disabled"),
    ('QPushButton[flat="true"]', 'QPushButton[flat="true"]:disabled'),
    ("QLineEdit", "QLineEdit:disabled"),
    ("QSpinBox", "QSpinBox:disabled"),
    ("QCheckBox", "QCheckBox:disabled"),
    ("QRadioButton", "QRadioButton:disabled"),
    ("QToolButton#railBtn", "QToolButton#railBtn:disabled"),
    ('QFrame[cls="flyItem"]', 'QFrame[cls="flyItem"]:disabled'),
    ('QFrame[cls="paletteItem"]', 'QFrame[cls="paletteItem"]:disabled'),
)


def test_qss_disabled_states_use_no_opacity():
    """Ninguna regla ``:disabled`` de widget usa ``opacity`` (solo soporta Qt
    la propiedad en ``QToolTip``)."""
    from chemuson.gui.theme.qss import get_main_stylesheet, get_tool_palette_stylesheet

    for theme in ("light", "dark"):
        for sheet in (get_main_stylesheet(theme), get_tool_palette_stylesheet(theme)):
            offenders = [
                ln.strip() for ln in sheet.splitlines()
                if re.search(r"^\s*opacity\s*:", ln)
            ]
            assert not offenders, (
                f"QSS del tema {theme} usa ``opacity`` en widgets: {offenders}"
            )


def test_qss_disabled_states_have_supported_visual_change():
    """Cada selector ``:disabled`` cambia visualmente el estado normal usando
    solo propiedades soportadas (color / background-color / border-color)."""
    from chemuson.gui.theme.qss import get_main_stylesheet

    qss = get_main_stylesheet("light")
    blocks = _parse_qss_blocks(qss)
    for base_sel, disabled_sel in _DISABLED_SELECTOR_PAIRS:
        assert disabled_sel in blocks, f"falta el estado :disabled para {disabled_sel}"
        changed = _disabled_visual_change(blocks, base_sel, disabled_sel)
        assert changed, (
            f"{disabled_sel} no introduce ningún cambio visual soportado "
            f"respecto a {base_sel}"
        )


def test_qss_disabled_applies_to_rail_and_flyout():
    from chemuson.gui.theme.qss import get_main_stylesheet

    qss = get_main_stylesheet("dark")
    for fragment in ("QToolButton#railBtn:disabled", "QFrame[cls=\"flyItem\"]:disabled"):
        assert fragment in qss, f"falta el estado :disabled para {fragment}"


# ---------------------------------------------------------------------------
# 4b. Thumbnails de plantilla: iconSize lógico del árbol
# ---------------------------------------------------------------------------
def test_templates_tree_sets_thumbnail_icon_size(win):
    """El árbol de PlantillasDock fija un ``iconSize`` lógico razonable (88×56)
    de modo que los thumbnails DPR-aware no se reduzcan al ~16 px por defecto,
    y cabe en el SidePanel de 340 px."""
    from PyQt6.QtCore import QSize

    from chemuson.gui.docks import _TEMPLATE_THUMB_SIZE, PlantillasDock

    assert isinstance(win.templates_dock, PlantillasDock)
    tree = win.templates_dock.tree
    assert tree.iconSize() == _TEMPLATE_THUMB_SIZE
    # El tamaño lógico debe coincidir con el pixmap de
    # ``template_preview_icon`` (88×56).
    assert _TEMPLATE_THUMB_SIZE == QSize(88, 56)
    # Debe caber en el SidePanel (340 px) con margen para el texto del ítem.
    assert win.side_panel.width() == 340
    assert _TEMPLATE_THUMB_SIZE.width() < win.side_panel.width()


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


_BENZENE_MOLBLOCK = (
    "  Benzene\n"
    "  ChemUSON\n"
    "\n"
    "  6  6  0  0  0  0  0  0  0  0999 V2000\n"
    "    1.0000    0.0000    0.0000 C   0  0\n"
    "    0.5000    0.8660    0.0000 C   0  0\n"
    "   -0.5000    0.8660    0.0000 C   0  0\n"
    "   -1.0000    0.0000    0.0000 C   0  0\n"
    "   -0.5000   -0.8660    0.0000 C   0  0\n"
    "    0.5000   -0.8660    0.0000 C   0  0\n"
    "  1  2  2  0\n"
    "  2  3  1  0\n"
    "  3  4  2  0\n"
    "  4  5  1  0\n"
    "  5  6  2  0\n"
    "  6  1  1  0\n"
    "M  END\n"
    "$$$\n"
)

# Mismo esqueleto que el benceno, con un heteroátomo etiquetado (tinta de texto).
_PYRIDINE_MOLBLOCK = _BENZENE_MOLBLOCK.replace(" C   0  0", " N   0  0", 1)


def _relative_ink_bbox(pixmap) -> tuple[float, float, float, float]:
    """Bounding box de la tinta (píxeles no transparentes) normalizado a 0..1.

    Se expresa en coordenadas relativas al pixmap, de modo que un DPR mayor no
    cambia la composición lógica: solo añade resolución.
    """
    image = pixmap.toImage()
    xs: list[int] = []
    ys: list[int] = []
    for y in range(image.height()):
        for x in range(image.width()):
            if image.pixelColor(x, y).alpha() > 0:
                xs.append(x)
                ys.append(y)
    assert xs, "el thumbnail está vacío (no hay tinta que medir)"
    return (
        min(xs) / image.width(),
        min(ys) / image.height(),
        (max(xs) + 1) / image.width(),
        (max(ys) + 1) / image.height(),
    )


def test_template_preview_relative_bbox_is_dpr_invariant(tmp_path, monkeypatch):
    """La estructura conserva el mismo bounding box relativo a DPR 1 y DPR 2.

    No basta con ``availableSizes()``: el DPR debe aportar solo resolución.
    Si el DPR se aplica dos veces (pixmap etiquetado antes de pintar +
    ``painter.scale(dpr, dpr)``), la estructura sale a ``dpr²`` y toca/recorta
    los bordes del thumbnail. Aquí se mide la tinta real del pixmap.
    """
    from types import SimpleNamespace

    from chemuson.gui.template_browser_service import TemplateBrowserService
    from chemuson.gui.template_library import DEFAULT_CATEGORY_USER, TemplateLibrary

    library = TemplateLibrary(tmp_path / "library.json")
    for molblock in (_BENZENE_MOLBLOCK, _PYRIDINE_MOLBLOCK):
        template = library.add_template("Estructura", DEFAULT_CATEGORY_USER, molblock)
        template_id = template["id"]
        boxes: dict[float, tuple[float, float, float, float]] = {}
        ink_pixels: dict[float, int] = {}

        for dpr in (1.0, 2.0):
            ctx = SimpleNamespace(
                template_library=library,
                preview_cache={},
                show_status=lambda _s: None,
            )
            monkeypatch.setattr(
                TemplateBrowserService,
                "_device_pixel_ratio",
                staticmethod(lambda d=dpr: d),
            )
            icon = TemplateBrowserService().template_preview_icon(ctx, template_id)
            physical = icon.availableSizes()[0]
            # Tamaño físico = lógico (88×56) × dpr.
            assert physical == QSize(int(round(88 * dpr)), int(round(56 * dpr)))
            pixmap = icon.pixmap(physical)
            boxes[dpr] = _relative_ink_bbox(pixmap)
            image = pixmap.toImage()
            ink_pixels[dpr] = sum(
                1
                for y in range(image.height())
                for x in range(image.width())
                if image.pixelColor(x, y).alpha() > 0
            )

        # Composición lógica idéntica (tolerancia de redondeo/antialiasing).
        for dpr1, dpr2 in zip(boxes[1.0], boxes[2.0]):
            assert abs(dpr2 - dpr1) <= 0.02, (
                f"bounding box relativo cambiado entre DPR 1 y DPR 2: "
                f"{boxes[1.0]} vs {boxes[2.0]}"
            )
        # La tinta no toca los bordes: el margen lógico de 8 px debe seguir
        # presente a DPR 2 (>= 0.05 del ancho/alto) y el contenido no puede
        # estar recortado (<= 0.95).
        left, top, right, bottom = boxes[2.0]
        assert left >= 0.05 and top >= 0.05, f"la estructura toca el borde inicial: {boxes[2.0]}"
        assert right <= 0.95 and bottom <= 0.95, f"la estructura está recortada: {boxes[2.0]}"
        # DPR 2 aporta resolución: la tinta crece ~dpr² (se permite margen).
        assert ink_pixels[2.0] >= 3.0 * ink_pixels[1.0], (
            f"DPR 2 no añade resolución: {ink_pixels[1.0]} -> {ink_pixels[2.0]}"
        )


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
    overlay.finished.connect(lambda completed: finished.append(completed))

    # Completa los 3 pasos (el último ``advance`` cierra con ``True``).
    overlay.advance()
    overlay.advance()
    overlay.advance()
    QApplication.processEvents()

    assert finished == [True], "completar los 3 pasos debe emitir ``finished(True)``"
    assert overlay.isHidden()
    assert setting_bool(
        application_settings().value("ui/onboarding/completed", False), False
    ), "debe persistir ``ui/onboarding/completed`` = True"


def test_onboarding_no_more_flag_default_false(win):
    overlay = _find_onboarding(win)
    assert overlay is not None
    assert overlay.no_more is False


def test_onboarding_close_with_no_more_persists(win):
    """Cerrar anticipadamente con "No volver a mostrar" → completado (True)."""
    from chemuson.platform.settings import application_settings, setting_bool

    overlay = _find_onboarding(win)
    assert overlay is not None
    finished = []
    overlay.finished.connect(lambda completed: finished.append(completed))

    # Marca "No volver a mostrar" y cierra antes de completar los 3 pasos.
    overlay.set_no_more(True)
    assert overlay.no_more is True
    overlay.request_close()
    QApplication.processEvents()

    assert finished == [True], "cerrar con la casilla marcada debe emitir ``finished(True)``"
    assert overlay.isHidden()
    assert setting_bool(
        application_settings().value("ui/onboarding/completed", False), False
    ), "con 'No volver a mostrar' debe persistir ``completed`` = True"


def test_onboarding_close_without_no_more_not_persisted(win):
    """Cerrar anticipadamente sin marcar → NO completado (False), no persiste."""
    from chemuson.platform.settings import application_settings, setting_bool

    overlay = _find_onboarding(win)
    assert overlay is not None
    finished = []
    overlay.finished.connect(lambda completed: finished.append(completed))

    # Cierra en el paso 1 sin marcar "No volver a mostrar".
    assert overlay.no_more is False
    overlay.request_close()
    QApplication.processEvents()

    assert finished == [False], "cerrar sin la casilla debe emitir ``finished(False)``"
    assert overlay.isHidden()
    assert not setting_bool(
        application_settings().value("ui/onboarding/completed", False), False
    ), "sin 'No volver a mostrar' NO debe persistir ``completed`` (aparece de nuevo)"


def test_onboarding_reappears_when_not_completed(_qapp, _isolated_config_home):
    """Si se cerró sin completar ni marcar, un siguiente arranque lo muestra de nuevo."""
    from chemuson.platform.settings import application_settings, setting_bool

    settings = application_settings()
    # Asegura que la clave NO está presente (primer arranque).
    settings.remove("ui/onboarding/completed")

    # Primer arranque: muestra, se cierra sin completar ni marcar.
    win1 = ChemusonWindow()
    win1.resize(1280, 800)
    win1.show()
    QApplication.processEvents()
    overlay1 = _find_onboarding(win1)
    assert overlay1 is not None, "debe mostrarse en el primer arranque"
    overlay1.request_close()
    QApplication.processEvents()
    win1.close()
    QApplication.processEvents()

    # La clave sigue ausente: un segundo arranque (nueva ventana) debe
    # mostrarlo otra vez.
    assert not setting_bool(
        application_settings().value("ui/onboarding/completed", False), False
    )
    win2 = ChemusonWindow()
    win2.resize(1280, 800)
    win2.show()
    QApplication.processEvents()
    overlay2 = _find_onboarding(win2)
    assert overlay2 is not None, "debe volver a aparecer en el siguiente arranque"
    win2.close()
    QApplication.processEvents()
