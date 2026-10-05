"""Registro de comandos de la Fase 6: construcción del catálogo a partir de la ventana.

Este módulo construye el :class:`CommandRegistry` de la paleta de comandos a
partir de las ``QAction`` **ya existentes** de la ventana (menús, toolbars
ocultos, app bar, side panel y menú dinámico de plantillas). Es la única pieza
que conoce la ventana; ``command_palette.py`` (presentación/filtro/ejecución)
no conoce la ventana ni importa dominio.

Principio: **la ``QAction`` existente es la fuente de verdad**. Aquí se
*registra* cada ``QAction`` (identidad) con metadatos de presentación
(sección, keywords, icono); no se crea ninguna ``QAction`` nueva ni se
duplica ningún handler.

Fuentes registradas (secciones):
- ``Archivo``: new/open/save/recovery/quit.
- ``Editar``: undo/redo, copiar/cortar/pegar/duplicar/eliminar, copiar como
  (SMILES/Molfile/CML/InChI), diagrama electrónico, rotaciones/fragmentos,
  grosor de enlace, redimensionar.
- ``Ver``: carbonos/hidrógenos/aromáticos, zoom, toolbar aux, reglas/cuadrícula,
  numeración, tamaño de lienzo, visibilidad de símbolos/panel lateral.
- ``Estructura``: Clean2D (quick/1 paso/publicación/conformero), validar,
  polímero/R-group, SMILES, nombre a estructura.
- ``Análisis``: los ocho análisis + nombre a estructura.
- ``Exportar``: PNG/SVG/PDF/CML/SMILES.
- ``Panel lateral``: las siete páginas (reutiliza ``side_panel_actions``).
- ``Plantillas``: las ``QAction`` del menú dinámico ``templates_menu`` (contrato
  existente de ``TemplateBrowserService``); adaptador mínimo si el menú no las
  expone.
- ``Preferencias``: preferencias, tema, guía rápida, actualizaciones, acerca de.
- ``Comandos``: la propia paleta (``action_command_palette``).

Ver OpenSpec ``2026-10-01-modernize-ui-command-palette``.
"""

from __future__ import annotations

from typing import Iterable

from PyQt6.QtGui import QAction

from chemuson.gui.command_palette import CommandRegistry
from chemuson.gui.template_library import TemplateLibrary

__all__ = ["build_command_registry"]


def _icon_for(action: QAction) -> str:
    """Nombre de icono presentacional (si la acción ya trae uno, se prefiere)."""
    icon = action.icon()
    if not icon.isNull():
        name = getattr(icon, "chemusonIconName", "") or ""
        if name:
            return str(name)
    return ""


def _register(reg: CommandRegistry, window, attr: str, section: str,
              keywords: Iterable[str] = (), icon: str = "") -> None:
    """Registra ``window.<attr>`` si existe (silencioso si no)."""
    action = getattr(window, attr, None)
    if isinstance(action, QAction):
        reg.register(action, section, list(keywords), icon or _icon_for(action))


def _walk_menu_actions(menu) -> list[tuple[QAction, str]]:
    """Recorre un ``QMenu`` (recursivo) devolviendo ``(action, parent_label)``."""
    out: list[tuple[QAction, str]] = []
    for item in menu.actions():
        if item.menu() is not None:
            for sub_action, _ in _walk_menu_actions(item.menu()):
                out.append((sub_action, item.text()))
        else:
            out.append((item, menu.title()))
    return out


def _register_templates(reg: CommandRegistry, window) -> None:
    """Registra las plantillas alcanzables por contrato ``QAction``.

    Vía principal (design D5): reutiliza las ``QAction`` ya creadas por
    ``TemplateBrowserService.refresh_templates_menu`` en ``window.templates_menu``
    (una por plantilla, ``triggered → start_template_insert_by_id``). Así la
    paleta dispara la ``QAction`` existente sin inventar un segundo sistema de
    inserción química.

    Adaptador mínimo (fallback): si el menú no expone acciones por plantilla,
    se leen ``template_library.grouped_templates()`` y se registra una entrada
    presentacional por plantilla cuyo ``triggered`` delega en el mismo contrato
    público ``window._start_template_insert_by_id`` (no duplica la inserción).
    """
    menu = getattr(window, "templates_menu", None)
    if menu is not None:
        seen: set[int] = set()
        for action, parent in _walk_menu_actions(menu):
            # Solo acciones de plantilla: las acciones de gestión (guardar/
            # categoría/importar/exportar) ya viven en la sección Estructura.
            if id(action) in seen:
                continue
            if action is getattr(window, "action_save_template", None):
                continue
            if action is getattr(window, "action_template_new_category", None):
                continue
            if action is getattr(window, "action_template_import_library", None):
                continue
            if action is getattr(window, "action_template_export_library", None):
                continue
            if not action.isEnabled() and action.text() in ("(Vacío)", "(Sin plantillas)"):
                continue
            seen.add(id(action))
            reg.register(
                action, "Plantillas",
                [parent, "insertar", "mol", "grupo"], _icon_for(action),
            )
        if seen:
            return

    # Fallback adaptador mínimo (sin duplicar inserción: delega en el contrato
    # público existente).
    library = getattr(window, "template_library", None)
    if not isinstance(library, TemplateLibrary):
        return
    start = getattr(window, "_start_template_insert_by_id", None)
    if not callable(start):
        return
    parent = getattr(window, "templates_dock", window)
    for group in library.grouped_templates():
        category = str(group.get("name", "")).strip() or "Plantillas"
        for template in group.get("templates", []):
            template_id = str(template.get("id", "")).strip()
            label = str(template.get("name", "")).strip() or "Plantilla"
            if not template_id:
                continue
            action = QAction(label, window)
            action.triggered.connect(
                lambda _checked=False, tid=template_id: start(tid)
            )
            reg.register(
                action, "Plantillas",
                [category, "insertar", "mol", "grupo"], "",
            )


def build_command_registry(window) -> CommandRegistry:
    """Construye el registro de comandos de ``window`` (una sola vez).

    Las secciones están ordenadas de forma estable (ver diseño); dentro de
    cada sección, el orden de registro se conserva. La deduplicación por
    identidad de ``QAction`` garantiza que reabrir la paleta no duplique.
    """
    reg = CommandRegistry()

    # Archivo
    _register(reg, window, "action_new", "Archivo", ["crear", "nuevo"], "doc")
    _register(reg, window, "action_open", "Archivo", ["abrir"], "doc")
    _register(reg, window, "action_save", "Archivo", ["guardar"], "doc")
    _register(reg, window, "action_recovery_center", "Archivo", ["recuperar", "autosave"])
    _register(reg, window, "action_quit", "Archivo", ["cerrar", "salida"])

    # Editar
    _register(reg, window, "action_undo", "Editar", ["deshacer"], "undo")
    _register(reg, window, "action_redo", "Editar", ["rehacer"], "redo")
    _register(reg, window, "action_copy", "Editar", ["copia"], "copy")
    _register(reg, window, "action_cut", "Editar", ["cortar"])
    _register(reg, window, "action_paste", "Editar", ["pegar"], "paste")
    _register(reg, window, "action_duplicate", "Editar", ["duplicado"])
    _register(reg, window, "action_delete", "Editar", ["borrar", "eliminar"])
    _register(reg, window, "action_copy_smiles", "Editar", ["copiar", "smiles", "exportar"])
    _register(reg, window, "action_copy_molfile", "Editar", ["copiar", "molfile", "mol"])
    _register(reg, window, "action_copy_cml", "Editar", ["copiar", "cml", "xml"])
    _register(reg, window, "action_copy_inchi", "Editar", ["copiar", "inchi", "identificador"])
    _register(reg, window, "action_edit_electronic_diagram", "Editar", ["diagrama", "electrónico"])
    _register(reg, window, "action_scale_selection", "Editar", ["redimensionar", "escala", "tamaño"])
    _register(reg, window, "action_bond_thickness_up", "Editar", ["grosor", "enlace"])
    _register(reg, window, "action_bond_thickness_down", "Editar", ["grosor", "enlace"])
    _register(reg, window, "action_bond_thickness_reset", "Editar", ["grosor", "reset"])

    # Rotaciones y fragmentos (Editar → Rotar)
    _register(reg, window, "action_rotate_left", "Editar", ["rotar", "girar", "90"], "rotate-left")
    _register(reg, window, "action_rotate_right", "Editar", ["rotar", "girar", "90"], "rotate-right")
    _register(reg, window, "action_flip_horizontal", "Editar", ["voltear", "horizontal", "180"], "flip-horizontal")
    _register(reg, window, "action_flip_vertical", "Editar", ["voltear", "vertical", "180"], "flip-vertical")
    _register(reg, window, "action_branch_rotate_minus", "Editar", ["rama", "rotar", "giro"])
    _register(reg, window, "action_branch_rotate_plus", "Editar", ["rama", "rotar", "giro"])
    _register(reg, window, "action_branch_invert", "Editar", ["rama", "inverter", "180"])
    _register(reg, window, "action_branch_auto_arrange", "Editar", ["rama", "auto", "acomodar"], "clean")
    _register(reg, window, "action_fragment_pivot_set", "Editar", ["fragmento", "pivote", "átomo"])
    _register(reg, window, "action_fragment_pivot_clear", "Editar", ["fragmento", "pivote", "limpiar"])
    _register(reg, window, "action_fragment_rotate_minus", "Editar", ["fragmento", "rotar"])
    _register(reg, window, "action_fragment_rotate_plus", "Editar", ["fragmento", "rotar"])
    _register(reg, window, "action_fragment_invert", "Editar", ["fragmento", "inverter"])

    # Ver
    _register(reg, window, "action_show_carbons", "Ver", ["carbonos", "mostrar", "C"])
    _register(reg, window, "action_show_hydrogens", "Ver", ["hidrógenos", "H", "mostrar"])
    _register(reg, window, "action_aromatic_circles", "Ver", ["aromáticos", "círculos", "benceno"])
    _register(reg, window, "action_style", "Ver", ["dimensiones", "dibujo", "estilo"])
    _register(reg, window, "action_zoom_in", "Ver", ["zoom", "acercar"], "zoom-in")
    _register(reg, window, "action_zoom_out", "Ver", ["zoom", "alejar"], "zoom-out")
    _register(reg, window, "action_zoom_reset", "Ver", ["zoom", "100", "reset"])
    _register(reg, window, "action_show_main_toolbar_aux", "Ver", ["barra", "copiar", "pegar", "zoom"])
    _register(reg, window, "action_rules", "Ver", ["reglas", "rulers", "margen"])
    _register(reg, window, "action_crosshair", "Ver", ["cuadrícula", "grid", "rejilla"])
    _register(reg, window, "action_numbering_enabled", "Ver", ["numeración", "números"])
    _register(reg, window, "action_numbering_mode_atoms", "Ver", ["numerar", "átomos"])
    _register(reg, window, "action_numbering_mode_structures", "Ver", ["numerar", "estructuras"])
    _register(reg, window, "action_numbering_mode_both", "Ver", ["numerar", "ambos"])
    _register(reg, window, "action_numbering_recalculate", "Ver", ["recalcular", "número"], "clean")
    _register(reg, window, "action_numbering_export", "Ver", ["numeración", "exportación", "incluir"])
    _register(reg, window, "action_canvas_size_letter_portrait", "Ver", ["tamaño", "lienzo", "letter", "vertical"])
    _register(reg, window, "action_canvas_size_letter_landscape", "Ver", ["tamaño", "lienzo", "letter", "horizontal"])
    _register(reg, window, "action_canvas_size_a4_portrait", "Ver", ["tamaño", "lienzo", "a4", "vertical"])
    _register(reg, window, "action_canvas_size_a4_landscape", "Ver", ["tamaño", "lienzo", "a4", "horizontal"])
    _register(reg, window, "action_canvas_size_a3_portrait", "Ver", ["tamaño", "lienzo", "a3", "vertical"])
    _register(reg, window, "action_canvas_size_a3_landscape", "Ver", ["tamaño", "lienzo", "a3", "horizontal"])
    _register(reg, window, "action_canvas_size_custom", "Ver", ["tamaño", "lienzo", "personalizado"])
    _register(reg, window, "action_toggle_side_panel", "Ver", ["panel", "lateral", "ocultar"])

    # Estructura
    _register(reg, window, "action_clean_2d", "Estructura", ["limpiar", "2d", "clean"], "clean")
    _register(reg, window, "action_clean_2d_full", "Estructura", ["limpiar", "2d", "1 paso", "quick", "clean"], "clean")
    _register(reg, window, "action_clean_2d_publication", "Estructura", ["limpiar", "2d", "publicación", "publicar"], "clean")
    _register(reg, window, "action_clean_2d_propose", "Estructura", ["conformero", "proponer", "2d"], "clean")
    _register(reg, window, "action_validate_structure", "Estructura", ["validar", "valencias", "errores"])
    _register(reg, window, "action_validation_next", "Estructura", ["validación", "siguiente", "error"])
    _register(reg, window, "action_validation_previous", "Estructura", ["validación", "anterior", "error"])
    _register(reg, window, "action_set_polymer_repeat", "Estructura", ["polímero", "repetición", "monómero"])
    _register(reg, window, "action_set_r_group_substituents", "Estructura", ["sustituyentes", "r", "grupo"])
    _register(reg, window, "action_import_smiles", "Estructura", ["importar", "smiles", "dibujar"], "flask")
    _register(reg, window, "action_export_smiles", "Estructura", ["exportar", "smiles"], "flask")
    _register(reg, window, "action_name_to_structure", "Estructura", ["nombre", "estructura", "n2s", "iupac"])
    _register(
        reg,
        window,
        "action_ai_molecular_assistant",
        "Estructura",
        ["ia", "ai", "generar", "estructura", "molécula", "smiles"],
    )

    # Análisis
    _register(reg, window, "action_analysis_name", "Análisis", ["nombre", "smiles", "nomenclatura"])
    _register(reg, window, "action_analysis_formula", "Análisis", ["fórmula", "química"])
    _register(reg, window, "action_analysis_exact", "Análisis", ["masa", "exacta"])
    _register(reg, window, "action_analysis_weight", "Análisis", ["peso", "molecular"])
    _register(reg, window, "action_analysis_mz", "Análisis", ["m/z", "masa", "carga"])
    _register(reg, window, "action_analysis_elemental", "Análisis", ["elemental", "composición"])
    _register(reg, window, "action_analysis_all", "Análisis", ["todo", "análisis", "completo"])

    # Exportar (QActions reales de exportación)
    _register(reg, window, "action_export_png", "Exportar", ["imagen", "png", "exportar"], "doc")
    _register(reg, window, "action_export_svg", "Exportar", ["vectorial", "svg", "exportar"], "doc")
    _register(reg, window, "action_export_pdf", "Exportar", ["pdf", "documento", "exportar"], "doc")
    _register(reg, window, "action_export_cml", "Exportar", ["cml", "xml", "chemistry", "exportar"], "doc")

    # Panel lateral (reutiliza side_panel_actions; no llama show_page en paralelo)
    side_panel_actions = getattr(window, "side_panel_actions", None)
    if isinstance(side_panel_actions, dict):
        for key, action in side_panel_actions.items():
            if isinstance(action, QAction):
                reg.register(action, "Panel lateral", [key, "panel", "pestaña"], "sliders")

    # Plantillas (menú dinámico existente)
    _register_templates(reg, window)
    _register(reg, window, "action_save_template", "Plantillas", ["guardar", "selección", "plantilla"], "doc")
    _register(reg, window, "action_template_new_category", "Plantillas", ["categoría", "nueva"])
    _register(reg, window, "action_template_import_library", "Plantillas", ["importar", "biblioteca"])
    _register(reg, window, "action_template_export_library", "Plantillas", ["exportar", "biblioteca"])
    _register(reg, window, "action_template_linear_chain", "Plantillas", ["cadena", "lineal", "insertar"])

    # Preferencias / tema / ayuda
    _register(reg, window, "action_preferences", "Preferencias", ["ajustes", "settings"], "sliders")
    _register(reg, window, "action_theme_toggle", "Preferencias", ["tema", "claro", "oscuro", "dark", "light"], "moon")
    _register(reg, window, "action_quick_start", "Preferencias", ["guía", "rápida", "ayuda"])
    _register(reg, window, "action_check_updates_now", "Preferencias", ["actualizaciones", "nueva", "versión"])
    _register(reg, window, "action_about", "Preferencias", ["acerca", "información", "quien"])

    # Comandos (la propia paleta)
    _register(reg, window, "action_command_palette", "Comandos", ["buscar", "comando", "paleta"], "search")

    return reg
