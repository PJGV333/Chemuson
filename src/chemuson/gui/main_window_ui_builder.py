from __future__ import annotations

from PyQt6.QtCore import QSize, Qt
from PyQt6.QtGui import QAction, QActionGroup, QIcon
from PyQt6.QtWidgets import QMenu, QToolBar

from chemuson.gui.energy_diagrams import ENERGY_DIAGRAM_MENU_ORDER, energy_diagram_display_name
from chemuson.gui.icons import draw_atom_icon, draw_generic_icon
from chemuson.gui.orbitals import ORBITAL_MENU_ORDER, orbital_display_name

from .controllers.view_controller import ViewController


class MainWindowUiBuilder:
    """Construcción declarativa de menús y barras principales."""

    def create_local_actions(self, window) -> None:
        self._create_template_actions(window)
        self._create_canvas_size_actions(window)
        self._create_misc_actions(window)
        self._create_analysis_actions(window)
        self._create_text_actions(window)

    def build_menu_bar(self, window) -> None:
        menubar = window.menuBar()

        file_menu = menubar.addMenu("Archivo")
        file_menu.addAction(window.action_new)
        file_menu.addAction(window.action_open)
        file_menu.addAction(window.action_save)
        window.recent_menu = file_menu.addMenu("Archivos recientes")
        window._update_recent_menu()
        file_menu.addAction(window.action_recovery_center)
        file_menu.addSeparator()

        export_menu = file_menu.addMenu("Exportar como")
        self._add_export_actions(window, export_menu)
        file_menu.addSeparator()
        file_menu.addAction(window.action_quit)

        edit_menu = menubar.addMenu("Editar")
        edit_menu.addAction(window.action_undo)
        edit_menu.addAction(window.action_redo)
        edit_menu.addSeparator()
        edit_menu.addAction(window.action_copy)
        edit_menu.addAction(window.action_cut)
        edit_menu.addAction(window.action_paste)
        edit_menu.addAction(window.action_duplicate)
        edit_menu.addAction(window.action_delete)
        self._add_copy_as_menu(window, edit_menu)
        edit_menu.addAction(window.action_paste)
        edit_menu.addAction(window.action_edit_electronic_diagram)
        edit_menu.addSeparator()
        self._build_rotate_menu(window, edit_menu)
        edit_menu.addAction(window.action_scale_selection)
        edit_menu.addSeparator()
        bond_thickness_menu = edit_menu.addMenu("Grosor de enlace/flecha/corchete")
        bond_thickness_menu.addAction(window.action_bond_thickness_up)
        bond_thickness_menu.addAction(window.action_bond_thickness_down)
        bond_thickness_menu.addAction(window.action_bond_thickness_reset)
        edit_menu.addSeparator()
        window.action_preferences = QAction("Preferencias...", window)
        window.action_preferences.triggered.connect(window._on_preferences)
        edit_menu.addAction(window.action_preferences)

        view_menu = menubar.addMenu("Ver")
        window.view_menu = view_menu
        self._build_view_menu(window, view_menu)

        structure_menu = menubar.addMenu("Estructura")
        self._build_structure_menu(window, structure_menu)

        reaction_menu = menubar.addMenu("Reacción")
        placeholder_reaction = QAction("Próximamente", window)
        placeholder_reaction.setEnabled(False)
        reaction_menu.addAction(placeholder_reaction)

        help_menu = menubar.addMenu("Ayuda")
        window.action_quick_start = QAction("Guía rápida...", window)
        window.action_quick_start.triggered.connect(window._on_quick_start)
        help_menu.addAction(window.action_quick_start)
        help_menu.addAction(window.action_check_updates_now)
        help_menu.addSeparator()
        window.action_about = QAction("Acerca de Chemuson...", window)
        window.action_about.triggered.connect(window._on_about)
        help_menu.addAction(window.action_about)

    def build_main_toolbar(self, window) -> None:
        window.main_toolbar = QToolBar("Principal")
        window.main_toolbar.setMovable(False)
        window.main_toolbar.setFloatable(False)
        window.main_toolbar.setIconSize(QSize(24, 24))
        window.main_toolbar.setToolButtonStyle(Qt.ToolButtonStyle.ToolButtonIconOnly)
        window.addToolBar(Qt.ToolBarArea.TopToolBarArea, window.main_toolbar)

        self.refresh_main_toolbar_icons(window)

        for action in (
            window.action_new,
            window.action_open,
            window.action_save,
        ):
            window.main_toolbar.addAction(action)
        window.main_toolbar.addSeparator()
        for action in (
            window.action_undo,
            window.action_redo,
        ):
            window.main_toolbar.addAction(action)
        window.main_toolbar.addSeparator()
        for action in (
            window.action_rotate_left,
            window.action_rotate_right,
            window.action_flip_horizontal,
            window.action_flip_vertical,
        ):
            window.main_toolbar.addAction(action)
        window.main_toolbar.addSeparator()
        window.main_toolbar.addAction(window.action_clean_2d)
        window.main_toolbar.addSeparator()
        window.main_toolbar.addAction(window.action_draw_smiles)

        window._main_toolbar_aux_separator_1 = QAction(window)
        window._main_toolbar_aux_separator_1.setSeparator(True)
        window._main_toolbar_aux_separator_2 = QAction(window)
        window._main_toolbar_aux_separator_2.setSeparator(True)
        self.set_main_toolbar_aux_visible(
            window,
            window.action_show_main_toolbar_aux.isChecked(),
        )

    def refresh_main_toolbar_icons(self, window) -> None:
        """Regenera iconos de la barra principal según el tema activo."""
        window.action_new.setIcon(draw_generic_icon("document_new"))
        window.action_open.setIcon(draw_generic_icon("document_open"))
        window.action_save.setIcon(draw_generic_icon("document_save"))
        window.action_undo.setIcon(draw_generic_icon("undo"))
        window.action_redo.setIcon(draw_generic_icon("redo"))
        window.action_copy.setIcon(draw_generic_icon("copy"))
        window.action_paste.setIcon(draw_generic_icon("paste"))
        window.action_zoom_in.setIcon(draw_generic_icon("zoom_in"))
        window.action_zoom_out.setIcon(draw_generic_icon("zoom_out"))
        window.action_rotate_left.setIcon(draw_generic_icon("rotate_left"))
        window.action_rotate_right.setIcon(draw_generic_icon("rotate_right"))
        window.action_flip_horizontal.setIcon(draw_generic_icon("flip_horizontal"))
        window.action_flip_vertical.setIcon(draw_generic_icon("flip_vertical"))
        window.action_branch_rotate_minus.setIcon(draw_generic_icon("rotate_left"))
        window.action_branch_rotate_plus.setIcon(draw_generic_icon("rotate_right"))
        window.action_branch_invert.setIcon(draw_generic_icon("flip_horizontal"))
        window.action_branch_auto_arrange.setIcon(draw_generic_icon("clean"))
        window.action_clean_2d.setIcon(draw_generic_icon("clean"))
        window.action_draw_smiles.setIcon(draw_atom_icon("SMI"))

    def set_main_toolbar_aux_visible(self, window, visible: bool) -> None:
        if not hasattr(window, "main_toolbar"):
            return
        aux_actions = (
            window.action_copy,
            window.action_paste,
            window._main_toolbar_aux_separator_1,
            window.action_zoom_in,
            window.action_zoom_out,
            window._main_toolbar_aux_separator_2,
        )
        current = window.main_toolbar.actions()
        if visible:
            anchor = window.action_rotate_left if window.action_rotate_left in current else None
            for action in aux_actions:
                if action in window.main_toolbar.actions():
                    continue
                if anchor is not None:
                    window.main_toolbar.insertAction(anchor, action)
                else:
                    window.main_toolbar.addAction(action)
            return
        for action in aux_actions:
            if action in window.main_toolbar.actions():
                window.main_toolbar.removeAction(action)

    def tool_status_label(self, window, tool_id: str) -> str:
        ring_label = f"Anillo {window.canvas.state.active_ring_size}"
        energy_diagram_label = energy_diagram_display_name(
            window.canvas.state.active_energy_diagram_kind
        )
        orbital_label = orbital_display_name(window.canvas.state.active_orbital_kind)
        if window.canvas.state.active_ring_template:
            template_name = {
                "haworth": "Haworth",
                "chair": "Silla",
            }.get(window.canvas.state.active_ring_template, "Anillo")
            anomeric = window.canvas.state.active_ring_anomeric
            suffix = ""
            if anomeric == "alpha":
                suffix = " α"
            elif anomeric == "beta":
                suffix = " β"
            ring_label = f"{template_name}{suffix}".strip()
        tool_names = {
            "tool_none": "Ninguna",
            "tool_select": "Seleccionar",
            "tool_select_lasso": "Seleccion (lazo)",
            "tool_bond": "Enlace",
            "tool_coordination_bond": "Enlace coordinativo",
            "tool_rotate_3d_precise": "Rotación 3D precisa",
            "tool_ring": ring_label,
            "tool_atom": f"Elemento {window.canvas.state.default_element}",
            "tool_energy_diagram": energy_diagram_label,
            "tool_orbital": orbital_label,
            "tool_coordination_center": "Centro de coordinación (esfera)",
            "tool_chain": "Cadena",
            "tool_arrow_line": "Linea",
            "tool_arrow_line_dashed": "Linea discontinua",
            "tool_arrow_forward": "Flecha directa",
            "tool_arrow_forward_open": "Flecha directa abierta",
            "tool_arrow_forward_dashed": "Flecha directa discontinua",
            "tool_arrow_retro": "Flecha retro",
            "tool_arrow_retro_open": "Flecha retro abierta",
            "tool_arrow_retro_dashed": "Flecha retro discontinua",
            "tool_arrow_both": "Flecha doble",
            "tool_arrow_both_open": "Flecha doble abierta",
            "tool_arrow_both_dashed": "Flecha doble discontinua",
            "tool_arrow_equilibrium": "Equilibrio",
            "tool_arrow_equilibrium_dashed": "Equilibrio discontinuo",
            "tool_arrow_retrosynthetic": "Flecha retrosintesis",
            "tool_arrow_curved": "Flecha curva",
            "tool_arrow_curved_fishhook": "Flecha curva (1 e-)",
            "tool_brackets_round": "Parentesis",
            "tool_brackets_square": "Corchetes",
            "tool_brackets_square_left": "Corchete izquierdo",
            "tool_brackets_square_right": "Corchete derecho",
            "tool_brackets_corner": "Esquinas",
            "tool_brackets_curly": "Llaves",
            "tool_brackets_curly_left": "Llave izquierda",
            "tool_brackets_curly_right": "Llave derecha",
            "tool_brackets_frame": "Marco",
            "tool_brackets_frame_rounded": "Marco redondeado",
            "tool_charge_plus": "Carga positiva",
            "tool_charge_minus": "Carga negativa",
            "tool_charge": "Carga alterna",
            "tool_symbol_plus": "Signo más",
            "tool_symbol_minus": "Signo menos",
            "tool_symbol_radical": "Electrón desapareado",
            "tool_symbol_lone_pair": "Par solitario",
            "tool_symbol_wavy_anchor": "Ancla ondulada",
            "tool_symbol_radical_cation": "Radical catión",
            "tool_symbol_radical_anion": "Radical anión",
            "tool_symbol_partial_plus": "Carga parcial (+)",
            "tool_symbol_partial_minus": "Carga parcial (-)",
        }
        tool_names.update(
            {
                f"tool_energy_diagram_{kind}": energy_diagram_display_name(kind)
                for kind in ENERGY_DIAGRAM_MENU_ORDER
            }
        )
        tool_names.update(
            {
                f"tool_orbital_{kind}": orbital_display_name(kind)
                for kind in ORBITAL_MENU_ORDER
            }
        )
        return tool_names.get(tool_id, tool_id)

    def _create_template_actions(self, window) -> None:
        window.action_template_linear_chain = QAction("Cadena lineal", window)
        window.action_template_linear_chain.triggered.connect(window._on_insert_linear_chain)
        window.action_template_new_category = QAction("Nueva categoría...", window)
        window.action_template_new_category.triggered.connect(window._on_new_template_category)
        window.action_template_import_library = QAction("Importar biblioteca...", window)
        window.action_template_import_library.triggered.connect(window._on_import_template_library)
        window.action_template_export_library = QAction("Exportar biblioteca...", window)
        window.action_template_export_library.triggered.connect(window._on_export_template_library)

    def _create_canvas_size_actions(self, window) -> None:
        presets = (
            ("action_canvas_size_letter_portrait", self.LETTER_PORTRAIT),
            ("action_canvas_size_letter_landscape", self.LETTER_LANDSCAPE),
            ("action_canvas_size_a4_portrait", self.A4_PORTRAIT),
            ("action_canvas_size_a4_landscape", self.A4_LANDSCAPE),
            ("action_canvas_size_a3_portrait", self.A3_PORTRAIT),
            ("action_canvas_size_a3_landscape", self.A3_LANDSCAPE),
        )
        for attr_name, preset in presets:
            action = QAction(preset.label, window)
            action.triggered.connect(
                lambda checked=False, width=preset.width, height=preset.height: window._set_canvas_size(width, height)
            )
            setattr(window, attr_name, action)
        window.action_canvas_size_custom = QAction("Personalizado...", window)
        window.action_canvas_size_custom.triggered.connect(window._on_canvas_custom_size)

    LETTER_PORTRAIT = ViewController.LETTER_PORTRAIT
    LETTER_LANDSCAPE = ViewController.LETTER_LANDSCAPE
    A4_PORTRAIT = ViewController.A4_PORTRAIT
    A4_LANDSCAPE = ViewController.A4_LANDSCAPE
    A3_PORTRAIT = ViewController.A3_PORTRAIT
    A3_LANDSCAPE = ViewController.A3_LANDSCAPE

    def _create_misc_actions(self, window) -> None:
        window.action_style = QAction("Dimensiones del dibujo...", window)
        window.action_style.triggered.connect(window._on_style_dialog)
        window.action_import_smiles = QAction("Importar SMILES...", window)
        window.action_import_smiles.triggered.connect(window._on_import_smiles)
        window.action_export_smiles = QAction("Exportar SMILES...", window)
        window.action_export_smiles.triggered.connect(window._on_export_smiles)
        window.action_draw_smiles = QAction("Dibujar desde SMILES...", window)
        window.action_draw_smiles.triggered.connect(window._on_import_smiles)
        window.action_scale_selection = QAction("Redimensionar selección...", window)
        window.action_scale_selection.setShortcut("Ctrl+Alt+S")
        window.action_scale_selection.triggered.connect(
            lambda checked=False: window.canvas.open_selection_scale_dialog()
        )
        window.action_bond_thickness_up = QAction(
            "Aumentar grosor de enlace/flecha/corchete",
            window,
        )
        window.action_bond_thickness_up.setShortcut("Ctrl+Shift+Up")
        window.action_bond_thickness_up.triggered.connect(window._on_bond_thickness_up)
        window.action_bond_thickness_down = QAction(
            "Reducir grosor de enlace/flecha/corchete",
            window,
        )
        window.action_bond_thickness_down.setShortcut("Ctrl+Shift+Down")
        window.action_bond_thickness_down.triggered.connect(window._on_bond_thickness_down)
        window.action_bond_thickness_reset = QAction(
            "Restablecer grosor de enlace/flecha/corchete",
            window,
        )
        window.action_bond_thickness_reset.setShortcut("Ctrl+Shift+0")
        window.action_bond_thickness_reset.triggered.connect(window._on_bond_thickness_reset)

    def _create_analysis_actions(self, window) -> None:
        window.action_analysis_name = QAction("Nombre (SMILES)", window)
        window.action_analysis_name.triggered.connect(lambda: window.canvas.run_analysis("name"))
        window.action_analysis_formula = QAction("Fórmula química", window)
        window.action_analysis_formula.triggered.connect(lambda: window.canvas.run_analysis("formula"))
        window.action_analysis_exact = QAction("Masa exacta", window)
        window.action_analysis_exact.triggered.connect(lambda: window.canvas.run_analysis("exact"))
        window.action_analysis_weight = QAction("Peso molecular", window)
        window.action_analysis_weight.triggered.connect(lambda: window.canvas.run_analysis("weight"))
        window.action_analysis_mz = QAction("m/z", window)
        window.action_analysis_mz.triggered.connect(lambda: window.canvas.run_analysis("mz"))
        window.action_analysis_elemental = QAction("Análisis elemental...", window)
        window.action_analysis_elemental.triggered.connect(window._open_elemental_analysis_dialog)
        window.action_analysis_all = QAction("Todo", window)
        window.action_analysis_all.triggered.connect(lambda: window.canvas.run_analysis("all"))

    def _create_text_actions(self, window) -> None:
        window.action_label_font = QAction("Fuente de etiquetas...", window)
        window.action_label_font.triggered.connect(window._on_label_font)
        window.action_label_size_set = QAction("Tamaño...", window)
        window.action_label_size_set.triggered.connect(window._on_label_font_size_dialog)
        window.action_label_bold = QAction("Negrita", window)
        window.action_label_bold.setCheckable(True)
        window.action_label_bold.triggered.connect(window._on_label_bold)
        window.action_label_italic = QAction("Cursiva", window)
        window.action_label_italic.setCheckable(True)
        window.action_label_italic.triggered.connect(window._on_label_italic)
        window.action_label_underline = QAction("Subrayado", window)
        window.action_label_underline.setCheckable(True)
        window.action_label_underline.triggered.connect(window._on_label_underline)
        window.action_label_subscript = QAction("Subíndice...", window)
        window.action_label_subscript.triggered.connect(window._on_label_subscript)
        window.action_label_superscript = QAction("Superíndice...", window)
        window.action_label_superscript.triggered.connect(window._on_label_superscript)
        window.action_label_size_up = QAction("Aumentar tamaño", window)
        window.action_label_size_up.triggered.connect(lambda: window._on_label_font_size(1.0))
        window.action_label_size_down = QAction("Reducir tamaño", window)
        window.action_label_size_down.triggered.connect(lambda: window._on_label_font_size(-1.0))
        window.action_text_align_left = QAction("Alinear a la izquierda", window)
        window.action_text_align_left.triggered.connect(
            lambda: window.canvas.update_text_alignment(Qt.AlignmentFlag.AlignLeft)
        )
        window.action_text_align_center = QAction("Centrar", window)
        window.action_text_align_center.triggered.connect(
            lambda: window.canvas.update_text_alignment(Qt.AlignmentFlag.AlignHCenter)
        )
        window.action_text_align_justify = QAction("Justificar", window)
        window.action_text_align_justify.triggered.connect(
            lambda: window.canvas.update_text_alignment(Qt.AlignmentFlag.AlignJustify)
        )
        window.action_label_color_element = QAction("Por elemento", window)
        window.action_label_color_element.setCheckable(True)
        window.action_label_color_element.triggered.connect(
            lambda checked=False: window._on_label_color_mode(True)
        )
        window.action_label_color_black = QAction("Negro", window)
        window.action_label_color_black.setCheckable(True)
        window.action_label_color_black.triggered.connect(
            lambda checked=False: window._on_label_color_mode(False)
        )
        window._label_color_group = QActionGroup(window)
        window._label_color_group.setExclusive(True)
        window._label_color_group.addAction(window.action_label_color_element)
        window._label_color_group.addAction(window.action_label_color_black)

    def _add_export_actions(self, window, file_menu: QMenu) -> None:
        window.action_export_png = QAction("PNG...", window)
        window.action_export_png.triggered.connect(lambda: window._on_export("png"))
        file_menu.addAction(window.action_export_png)
        window.action_export_svg = QAction("SVG...", window)
        window.action_export_svg.triggered.connect(lambda: window._on_export("svg"))
        file_menu.addAction(window.action_export_svg)
        window.action_export_pdf = QAction("PDF...", window)
        window.action_export_pdf.triggered.connect(lambda: window._on_export("pdf"))
        file_menu.addAction(window.action_export_pdf)

    def _add_copy_as_menu(self, window, edit_menu: QMenu) -> None:
        copy_as_menu = edit_menu.addMenu("Copiar como")
        window.action_copy_smiles = QAction("SMILES", window)
        window.action_copy_smiles.triggered.connect(lambda: window._on_copy_as("smiles"))
        copy_as_menu.addAction(window.action_copy_smiles)
        window.action_copy_molfile = QAction("Molfile", window)
        window.action_copy_molfile.triggered.connect(lambda: window._on_copy_as("molfile"))
        copy_as_menu.addAction(window.action_copy_molfile)
        window.action_copy_inchi = QAction("InChI", window)
        window.action_copy_inchi.triggered.connect(lambda: window._on_copy_as("inchi"))
        copy_as_menu.addAction(window.action_copy_inchi)

    def _build_rotate_menu(self, window, edit_menu: QMenu) -> None:
        rotate_menu = edit_menu.addMenu("Rotar")
        rotate_menu.addAction(window.action_rotate_left)
        rotate_menu.addAction(window.action_rotate_right)
        rotate_menu.addSeparator()
        rotate_menu.addAction(window.action_flip_horizontal)
        rotate_menu.addAction(window.action_flip_vertical)
        rotate_menu.addSeparator()
        rotate_menu.addAction(window.action_branch_rotate_minus)
        rotate_menu.addAction(window.action_branch_rotate_plus)
        rotate_menu.addAction(window.action_branch_invert)
        rotate_menu.addAction(window.action_branch_auto_arrange)
        rotate_menu.addSeparator()
        window.fragment_rotate_menu = rotate_menu.addMenu("Fragmento con pivote")
        window.fragment_rotate_menu.addAction(window.action_fragment_pivot_set)
        window.fragment_rotate_menu.addAction(window.action_fragment_pivot_clear)
        window.fragment_rotate_menu.addSeparator()
        window.fragment_rotate_menu.addAction(window.action_fragment_rotate_minus)
        window.fragment_rotate_menu.addAction(window.action_fragment_rotate_plus)
        window.fragment_rotate_menu.addAction(window.action_fragment_invert)

    def _build_view_menu(self, window, view_menu: QMenu) -> None:
        view_menu.addAction(window.action_show_carbons)
        view_menu.addAction(window.action_show_hydrogens)
        view_menu.addSeparator()
        view_menu.addAction(window.action_aromatic_circles)
        view_menu.addSeparator()
        view_menu.addAction(window.action_style)
        view_menu.addSeparator()
        view_menu.addAction(window.action_zoom_in)
        view_menu.addAction(window.action_zoom_out)
        view_menu.addAction(window.action_zoom_reset)
        view_menu.addSeparator()
        view_menu.addAction(window.action_show_main_toolbar_aux)
        view_menu.addSeparator()
        view_menu.addAction(window.action_rules)
        view_menu.addAction(window.action_crosshair)
        view_menu.addSeparator()
        numbering_menu = view_menu.addMenu("Numeración")
        numbering_menu.addAction(window.action_numbering_enabled)
        numbering_menu.addSeparator()
        numbering_menu.addAction(window.action_numbering_mode_atoms)
        numbering_menu.addAction(window.action_numbering_mode_structures)
        numbering_menu.addAction(window.action_numbering_mode_both)
        numbering_menu.addSeparator()
        numbering_menu.addAction(window.action_numbering_recalculate)
        numbering_menu.addAction(window.action_numbering_export)
        view_menu.addSeparator()
        canvas_size_menu = view_menu.addMenu("Tamaño de lienzo")
        for action in (
            window.action_canvas_size_letter_portrait,
            window.action_canvas_size_letter_landscape,
            None,
            window.action_canvas_size_a4_portrait,
            window.action_canvas_size_a4_landscape,
            None,
            window.action_canvas_size_a3_portrait,
            window.action_canvas_size_a3_landscape,
            None,
            window.action_canvas_size_custom,
        ):
            if action is None:
                canvas_size_menu.addSeparator()
            else:
                canvas_size_menu.addAction(action)
        view_menu.addSeparator()
        view_menu.addAction(window.symbols_toolbar.toggleViewAction())
        view_menu.addAction(window.templates_dock.toggleViewAction())
        view_menu.addAction(window.inspector_dock.toggleViewAction())
        view_menu.addAction(window.chemical_properties_dock.toggleViewAction())
        view_menu.addAction(window.appearance_dock.toggleViewAction())

    def _build_structure_menu(self, window, structure_menu: QMenu) -> None:
        structure_menu.addAction(window.action_clean_2d)
        structure_menu.addAction(window.action_clean_2d_full)
        structure_menu.addSeparator()
        structure_menu.addAction(window.action_validate_structure)
        structure_menu.addAction(window.action_validation_next)
        structure_menu.addAction(window.action_validation_previous)
        structure_menu.addSeparator()
        window.templates_menu = structure_menu.addMenu("Plantillas")
        window.action_save_template = QAction("Guardar selección como plantilla...", window)
        window.action_save_template.triggered.connect(window._on_save_template)
        window._refresh_templates_menu()
        structure_menu.addAction(window.action_save_template)
        structure_menu.addAction(window.action_template_import_library)
        structure_menu.addAction(window.action_template_export_library)
        structure_menu.addSeparator()
        structure_menu.addAction(window.action_import_smiles)
        structure_menu.addAction(window.action_export_smiles)
        structure_menu.addSeparator()
        analysis_menu = structure_menu.addMenu("Análisis")
        analysis_menu.addAction(window.action_analysis_name)
        analysis_menu.addAction(window.action_analysis_formula)
        analysis_menu.addAction(window.action_analysis_exact)
        analysis_menu.addAction(window.action_analysis_weight)
        analysis_menu.addAction(window.action_analysis_mz)
        analysis_menu.addAction(window.action_analysis_elemental)
        analysis_menu.addSeparator()
        analysis_menu.addAction(window.action_analysis_all)
