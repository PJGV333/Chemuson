"""Generador de capturas autÃ©nticas de ChemUSON 0.3.0-beta.4 para el manual de usuario.

Corre sin pantalla (QT_QPA_PLATFORM=offscreen) sobre la ventana real de la
aplicaciÃ³n, es decir sobre el cÃ³digo del tag v0.3.0-beta.4. No dibuja ni retoca
ilustraciones: captura el widget real (ventana, menÃº, paleta del rail, lienzo o
panel) y recorta la regiÃ³n de interÃ©s.

Uso (desde la raÃ­z del repositorio, con un entorno con PyQt6):

    QT_QPA_PLATFORM=offscreen python docs/manual-figures/capturar_figuras.py ventana menu flyout contenido contextual paneles

Grupos: ventana, menu, flyout, contenido, contextual, paneles.
Las imÃ¡genes se escriben en este mismo directorio.
"""

from __future__ import annotations

import math
import os
import sys
import tempfile

os.environ.setdefault("QT_QPA_PLATFORM", "offscreen")
REPO = os.environ.get(
    "CHEMUSON_REPO",
    os.path.abspath(os.path.join(os.path.dirname(os.path.abspath(__file__)), os.pardir, os.pardir)),
)
sys.path.insert(0, os.path.join(REPO, "src"))

OUT = os.environ.get("CHEMUSON_OUT", os.path.dirname(os.path.abspath(__file__)))
os.makedirs(OUT, exist_ok=True)

from PyQt6.QtCore import QPoint, QPointF, QRectF, QSettings, Qt  # noqa: E402
from PyQt6.QtWidgets import QApplication  # noqa: E402

from chemuson.gui.energy_diagrams import ENERGY_DIAGRAM_MENU_ORDER  # noqa: E402
from chemuson.gui.orbitals import ORBITAL_MENU_ORDER  # noqa: E402

REPORT: list[str] = []
APP: QApplication = None
W = None
CANVAS = None

FIGS = {}


def fig(name):
    def deco(fn):
        FIGS[name] = fn
        return fn

    return deco


# ------------------------------------------------------------------- utilices
def isolated_settings():
    path = os.path.join(tempfile.mkdtemp(), "chemuson-docs.ini")

    def factory():
        store = QSettings(path, QSettings.Format.IniFormat)
        store.clear()
        store.setValue("ui/onboarding/completed", True)
        store.sync()
        return store

    return factory


def pump(n=12):
    for _ in range(n):
        APP.processEvents()


def save(pm, name):
    if pm is None or pm.isNull():
        line = f"FALLO  {name}: pixmap nulo"
        REPORT.append(line)
        print(line, flush=True)
        return
    pm.save(os.path.join(OUT, name))
    line = f"OK    {name}  {pm.width()}x{pm.height()}"
    REPORT.append(line)
    print(line, flush=True)


def new_window(dark=False, width=1440, height=900):
    global W, CANVAS
    print(f"-- nueva ventana {width}x{height}", flush=True)
    from chemuson.gui.main_window import ChemusonWindow

    if W is not None:
        W.close()
        pump(6)
    W = ChemusonWindow()
    W.resize(width, height)
    W.show()
    pump(20)
    if dark:
        W.toggle_theme(True)
        pump(16)
    CANVAS = W.canvas
    return W


def draw_ring(cx=520.0, cy=420.0, r=42.0, aromatic=True, element="C"):
    model = CANVAS.model
    ids = []
    for i in range(6):
        ang = math.radians(60 * i - 90)
        ids.append(
            model.add_atom(element, cx + r * math.cos(ang), cy + r * math.sin(ang))
        )
    for i in range(6):
        model.add_bond(ids[i].id, ids[(i + 1) % 6].id, order=1, is_aromatic=aromatic)
    CANVAS._rebuild_items_from_model()
    pump(6)
    return ids


def fit_and_grab_canvas(name, margin=28):
    scene = CANVAS.scene
    if not scene.items():
        REPORT.append(f"FALLO {name}: escena vacia")
        return
    rect = scene.itemsBoundingRect().adjusted(-margin, -margin, margin, margin)
    CANVAS.fitInView(rect, Qt.AspectRatioMode.KeepAspectRatio)
    pump(8)
    tl = CANVAS.mapFromScene(rect.topLeft())
    br = CANVAS.mapFromScene(rect.bottomRight())
    box = (
        QRectF(QPointF(tl), QPointF(br)).normalized().toRect().intersected(CANVAS.rect())
    )
    save(CANVAS.grab(box), name)


def reset_canvas():
    for item in list(CANVAS.scene.items()):
        CANVAS.scene.removeItem(item)
    CANVAS.model.clear()
    CANVAS.atom_items.clear()
    CANVAS.bond_items.clear()
    CANVAS.undo_stack.clear()
    pump(6)


def grab_region(name, x, y, w, h, scale=1.4):
    """Encuadre y captura de una region concreta del papel del lienzo."""
    rect = QRectF(x, y, w, h)
    CANVAS.setSceneRect(rect)
    CANVAS.fitInView(rect, Qt.AspectRatioMode.KeepAspectRatio)
    pump(8)
    tl = CANVAS.mapFromScene(rect.topLeft())
    br = CANVAS.mapFromScene(rect.bottomRight())
    box = (
        QRectF(QPointF(tl), QPointF(br)).normalized().toRect().intersected(CANVAS.rect())
    )
    save(CANVAS.grab(box), name)


def menubar_menu(title):
    for act in W.menuBar().actions():
        if act.text() == title:
            return act.menu()
    return None


def grab_menu(name, title, open_sub=None):
    """Captura un menu de la barra; ``open_sub`` admite ruta anidada."""
    menu = menubar_menu(title)
    menu.popup(W.mapToGlobal(QPoint(140, 140)))
    pump(16)
    target = menu
    if open_sub:
        path = open_sub if isinstance(open_sub, (tuple, list)) else (open_sub,)
        for label in path:
            found = None
            for act in target.actions():
                if act.text() == label and act.menu() is not None:
                    found = act.menu()
                    break
            if found is None:
                REPORT.append(f"FALLO {name}: submenu {label!r} ausente")
                break
            parent = target
            target = found
            target.popup(parent.mapToGlobal(QPoint(parent.width() - 10, 30)))
            pump(16)
    save(target.grab(), name)
    if target is not menu:
        target.close()
    menu.close()
    pump(4)


def grab_flyout(key, name):
    fly = W.tool_rail.open_flyout(key)
    pump(14)
    if fly is None:
        REPORT.append(f"FALLO {name}: sin flyout {key!r}")
        return
    save(fly.grab(), name)
    W.tool_rail.open_flyout(key)
    pump(4)


# --------------------------------------------------------------------- grupos
def grupo_ventana():
    new_window()
    draw_ring(560.0, 430.0, 70.0)
    CANVAS.fitInView(
        CANVAS.scene.itemsBoundingRect().adjusted(-60, -60, 60, 60),
        Qt.AspectRatioMode.KeepAspectRatio,
    )
    W.side_panel.show_page("inspector")
    W.action_toggle_side_panel.setChecked(True)
    pump(20)
    save(W.grab(), "01-ventana-tema-claro.png")
    W.toggle_theme(True)
    pump(20)
    save(W.grab(), "02-ventana-tema-oscuro.png")


def grupo_menu():
    new_window()
    grab_menu("03-menu-estructura.png", "Estructura")
    grab_menu("04-menu-ayuda.png", "Ayuda")
    grab_menu("05-menu-editar.png", "Editar")
    grab_menu("06-menu-editar-rotar.png", "Editar", open_sub="Rotar")
    grab_menu(
        "07-menu-rotar-fragmento-pivote.png",
        "Editar",
        open_sub=("Rotar", "Fragmento con pivote"),
    )
    grab_menu("08-menu-ver-paneles.png", "Ver", open_sub="Paneles")
    grab_menu("09-menu-estructura-plantillas.png", "Estructura", open_sub="Plantillas")


def grupo_flyout():
    new_window()
    for key, name in (
        ("select", "10-flyout-seleccion.png"),
        ("atom", "11-flyout-atomos.png"),
        ("arrows", "12-flyout-flechas.png"),
        ("brackets", "13-flyout-corchetes.png"),
        ("symbols", "14-flyout-simbolos.png"),
        ("plates", "15-flyout-placas.png"),
        ("energy", "16-flyout-energia.png"),
        ("orbitals", "17-flyout-orbitales.png"),
    ):
        grab_flyout(key, name)


def grupo_paneles():
    new_window()
    draw_ring(560.0, 430.0, 60.0)
    CANVAS.atom_items[next(iter(CANVAS.atom_items))].setSelected(True)
    CANVAS._sync_selection_from_scene()
    W.action_toggle_side_panel.setChecked(True)
    pump(16)
    W._refresh_spectroscopy_dock(CANVAS.model)
    rows, _ = W._chemical_properties_rows()
    W.chemical_properties_dock.update_properties(rows)
    W._refresh_validation_dock()
    pump(20)
    for page, name in (
        ("inspector", "40-panel-inspector.png"),
        ("properties", "41-panel-propiedades.png"),
        ("validation", "42-panel-validacion.png"),
        ("templates", "43-panel-plantillas.png"),
        ("appearance", "44-panel-apariencia.png"),
        ("spectroscopy", "45-panel-espectroscopia.png"),
        ("compchem", "46-panel-compchem.png"),
    ):
        W.side_panel.show_page(page)
        pump(14)
        save(W.side_panel.grab(), name)
    # Paleta de comandos abierta.
    W.action_command_palette.trigger()
    pump(14)
    palette = getattr(W, "command_palette", None)
    if palette is not None:
        save(palette.grab(), "47-paleta-comandos.png")
    else:
        REPORT.append("FALLO 47-paleta-comandos.png: sin objeto command_palette")


def grupo_fig():
    wanted = sys.argv[2:]
    for name, fn in FIGS.items():
        if not wanted or name in wanted:
            print(f"-- figura {name}", flush=True)
            fn()


def grupo_contenido():
    grupo_fig()


@fig("flechas")
def _fig_flechas():
    # Flechas y corchetes junto a una estructura pequena.
    new_window(width=1100, height=760)
    draw_ring(170.0, 210.0, 44.0)
    CANVAS.add_arrow_item(QPointF(300, 140), QPointF(450, 140), "forward")
    CANVAS.add_arrow_item(QPointF(300, 190), QPointF(450, 190), "equilibrium")
    CANVAS.add_arrow_item(QPointF(300, 240), QPointF(450, 240), "retrosynthetic")
    CANVAS.add_arrow_item(QPointF(300, 290), QPointF(450, 290), "line")
    CANVAS.add_arrow_item(QPointF(300, 340), QPointF(450, 340), "line_dashed")
    CANVAS.add_arrow_item(QPointF(300, 390), QPointF(450, 390), "both")
    CANVAS.add_arrow_item(QPointF(540, 130), QPointF(660, 260), "curved", 0.6)
    CANVAS.add_arrow_item(QPointF(700, 130), QPointF(820, 260), "curved_fishhook", -0.6)
    CANVAS.add_bracket_item(QRectF(120, 155, 105, 110), "[]")
    grab_region("20-canvas-flechas-corchetes.png", 60, 60, 800, 340)

    # Placas: TLC y gel de electroforesis (valores por defecto del dialogo).
    reset_canvas()
    from chemuson.gui.dialogs import GelInsertDialog, TLCInsertDialog

    TLCInsertDialog.exec = lambda self, *a: 1
    GelInsertDialog.exec = lambda self, *a: 1
    CANVAS._insert_tlc_plate_item(QPointF(300, 250))
    CANVAS._insert_gel_electrophoresis_item(QPointF(620, 250))
    pump(12)
    # Distribuir manchas y bandas a distinto Rf (equivale a "Set Rf..." del
    # menu contextual) para una figura representativa.
    rf_values = [0.25, 0.55, 0.8, 0.4, 0.65]
    n = 0
    for item in list(CANVAS.scene.items()):
        lane = getattr(item, "lane_ref", None)
        if lane is None or not hasattr(lane, "baseline_y"):
            continue
        dist_total = lane.baseline_y() - lane.solvent_front_y()
        if dist_total == 0:
            continue
        rf = rf_values[n % len(rf_values)]
        n += 1
        item.setY(lane.scenePos().y() + lane.baseline_y() - rf * dist_total)
        if hasattr(item, "_show_rf_label"):
            item._show_rf_label = True
        lane.update_rf_labels()
    CANVAS.scene.clearSelection()
    pump(8)
    CANVAS.centerOn(0, 0)
    grab_region("21-canvas-placas.png", 0, 0, 1100, 700)

    print("-- simbolos", flush=True)
    # Simbolos quimicos sobre una estructura sencilla (morfolino).
    reset_canvas()
    cx, cy = 300.0, 300.0
    ids = draw_ring(cx, cy, 52.0, aromatic=False)
    model = CANVAS.model
    model.get_atom(ids[0].id).element = "N"
    model.get_atom(ids[3].id).element = "O"
    CANVAS._rebuild_items_from_model()
    pump(6)

    def outward(aid, dist=46.0):
        atom = model.get_atom(aid)
        dx, dy = atom.x - cx, atom.y - cy
        norm = max(1e-6, (dx * dx + dy * dy) ** 0.5)
        return QPointF(atom.x + dx / norm * dist, atom.y + dy / norm * dist)

    # Carga formal +1 sobre el nitrogeno (se dibuja como superindice propio).
    model.get_atom(ids[0].id).formal_charge = 1
    CANVAS._rebuild_items_from_model()
    pump(4)
    # Carga parcial negativa junto al oxigeno y su par solitario.
    o_id = ids[3].id
    o = model.get_atom(o_id)
    CANVAS._insert_electron_dots(QPointF(o.x, o.y - 4), o_id, 2, 1.0)
    CANVAS._insert_symbol_item(
        "δ-", QPointF(o.x + 46, o.y + 20), None, 1.0, False
    )
    # Carga parcial positiva sobre el carbono adyacente al oxigeno.
    c_id = ids[2].id
    c = model.get_atom(c_id)
    CANVAS._insert_symbol_item(
        "δ+", QPointF(c.x + 34, c.y + 26), c_id, 1.0, True
    )
    pump(10)
    grab_region("22-canvas-simbolos-cargas.png", 150, 180, 320, 260)

    print("-- energia", flush=True)
    # Diagramas de energia: presets de la paleta.
    reset_canvas()
    for i, kind in enumerate(ENERGY_DIAGRAM_MENU_ORDER[:8]):
        col, row = i % 4, i // 4
        CANVAS._insert_energy_diagram_item(
            QPointF(180 + col * 200, 180 + row * 230), kind
        )
    pump(10)
    grab_region("23-canvas-diagramas-energia.png", 40, 40, 840, 660)

    print("-- orbitales", flush=True)
    # Orbitales: muestra representativa de la paleta.
    reset_canvas()
    for i, kind in enumerate(ORBITAL_MENU_ORDER[:12]):
        col, row = i % 4, i // 4
        CANVAS._insert_orbital_item(QPointF(180 + col * 190, 200 + row * 210), kind)
    pump(10)
    grab_region("24-canvas-orbitales.png", 40, 40, 780, 660)

    print("-- enlaces", flush=True)
    # Ordenes y estilos de enlace (seccion 15.4): dos capturas distintas.
    reset_canvas()
    from chemuson.core.model import BondStyle

    model = CANVAS.model
    b = [
        model.add_atom("C", 180 + 60 * i, 220 if i % 2 == 0 else 265)
        for i in range(4)
    ]
    model.add_bond(b[0].id, b[1].id, order=1)
    model.add_bond(b[1].id, b[2].id, order=2)
    model.add_bond(b[2].id, b[3].id, order=3)
    CANVAS._rebuild_items_from_model()
    pump(8)
    grab_region("25-canvas-ordenes-enlace.png", 80, 120, 520, 240)

    reset_canvas()
    model = CANVAS.model
    a = [
        model.add_atom("C", 150 + 60 * i, 220 if i % 2 == 0 else 265)
        for i in range(7)
    ]
    for i, st in enumerate(
        [
            BondStyle.PLAIN,
            BondStyle.BOLD,
            BondStyle.WEDGE,
            BondStyle.HASHED,
            BondStyle.WAVY,
            BondStyle.FLEX,
        ]
    ):
        model.add_bond(a[i].id, a[i + 1].id, order=1, style=st)
    CANVAS._rebuild_items_from_model()
    pump(8)
    grab_region("26-canvas-estilos-enlace.png", 60, 120, 720, 240)


def grupo_contextual():
    new_window(width=1100, height=760)
    ids = draw_ring(300.0, 260.0, 48.0)
    atom_item = CANVAS.atom_items[ids[0].id]
    bond_item = CANVAS.bond_items[next(iter(CANVAS.bond_items))]

    def context_shot(name, scene_pos, clicked_item=None, with_sub=False):
        captured = {}

        def fake_exec(self, pos=None):
            pump(16)
            if with_sub:
                for act in self.actions():
                    sub = act.menu()
                    if sub is not None and (sub.title() or "").lower().startswith("an"):
                        sub.popup(self.mapToGlobal(QPoint(self.width() - 14, 14)))
                        pump(14)
                        captured["sub"] = sub.grab()
                        sub.close()
                        break
            captured["menu"] = self.grab()
            return None

        from PyQt6.QtWidgets import QMenu

        original = QMenu.exec
        QMenu.exec = fake_exec
        try:
            CANVAS._show_context_menu(
                scene_pos, CANVAS.mapToGlobal(QPointF(scene_pos)), clicked_item
            )
        finally:
            QMenu.exec = original
        if "menu" in captured:
            save(captured["menu"], name)
        else:
            REPORT.append(f"FALLO {name}: no se construyo el menu")
        if captured.get("sub") is not None:
            save(captured["sub"], "33-contextual-analisis-submenu.png")

    context_shot("30-contextual-atomo.png", atom_item.pos(), atom_item)
    context_shot("31-contextual-enlace.png", bond_item.pos(), bond_item)
    context_shot("32-contextual-lienzo-vacio.png", QPointF(120, 120), None)
    for aid in ids[:3]:
        CANVAS.atom_items[aid.id].setSelected(True)
    CANVAS._sync_selection_from_scene()
    pump(6)
    context_shot("34-contextual-seleccion.png", atom_item.pos(), None)
    CANVAS.scene.clearSelection()
    CANVAS._sync_selection_from_scene()
    pump(4)


def main():
    global APP
    import chemuson.gui.shell.assembly as asm
    import chemuson.platform.settings as ps

    factory = isolated_settings()
    ps.application_settings = factory
    asm.application_settings = factory

    APP = QApplication([])
    groups = set(sys.argv[1:])
    all_groups = {
        "ventana": grupo_ventana,
        "menu": grupo_menu,
        "flyout": grupo_flyout,
        "contenido": grupo_contenido,
        "contextual": grupo_contextual,
        "paneles": grupo_paneles,
        "fig": grupo_fig,
    }
    for key, fn in all_groups.items():
        if not groups or key in groups:
            print(f"### {key}", flush=True)
            fn()
    print("\n".join(REPORT))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
