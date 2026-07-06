from __future__ import annotations

import builtins
from types import SimpleNamespace

import pytest
from PyQt6.QtGui import QAction
from PyQt6.QtWidgets import QApplication, QMenu, QWidget

from chemuson.core.model import MolGraph
from chemuson.gui.template_browser_service import TemplateBrowserService
from chemuson.gui.template_library import DEFAULT_CATEGORY_USER, TemplateLibrary


@pytest.fixture(scope="module", autouse=True)
def _qapp():
    return QApplication.instance() or QApplication([])


def _migration_context(library: TemplateLibrary, statuses: list[str]) -> SimpleNamespace:
    return SimpleNamespace(
        template_library=library,
        default_category_user=DEFAULT_CATEGORY_USER,
        show_status=statuses.append,
    )


def test_migrate_legacy_templates_imports_new_mol_file(tmp_path, monkeypatch):
    templates_dir = tmp_path / "legacy"
    templates_dir.mkdir()
    (templates_dir / "mi_plantilla.mol").write_text("MOL LEGADO\n", encoding="utf-8")
    monkeypatch.setattr(
        TemplateBrowserService,
        "_get_templates_dir",
        staticmethod(lambda: str(templates_dir)),
    )
    library = TemplateLibrary(tmp_path / "library.json")
    statuses: list[str] = []

    TemplateBrowserService().migrate_legacy_templates(_migration_context(library, statuses))

    imported = [template for template in library.list_templates() if template["name"] == "mi plantilla"]
    assert len(imported) == 1
    assert imported[0]["category"] == DEFAULT_CATEGORY_USER
    assert imported[0]["molblock"] == "MOL LEGADO"
    assert statuses == ["Plantillas legadas importadas a la biblioteca."]


def test_migrate_legacy_templates_skips_duplicates(tmp_path, monkeypatch):
    templates_dir = tmp_path / "legacy"
    templates_dir.mkdir()
    (templates_dir / "Duplicada.mol").write_text("MOL DUP\n", encoding="utf-8")
    monkeypatch.setattr(
        TemplateBrowserService,
        "_get_templates_dir",
        staticmethod(lambda: str(templates_dir)),
    )
    library = TemplateLibrary(tmp_path / "library.json")
    library.add_template("Duplicada", DEFAULT_CATEGORY_USER, "MOL DUP")
    statuses: list[str] = []

    TemplateBrowserService().migrate_legacy_templates(_migration_context(library, statuses))

    matches = [template for template in library.list_templates() if template["name"] == "Duplicada"]
    assert len(matches) == 1
    assert statuses == []


class _FakeTemplateLibrary:
    def __init__(self, grouped: list[dict], graph_error: Exception | None = None) -> None:
        self._grouped = grouped
        self._graph_error = graph_error

    def grouped_templates(self) -> list[dict]:
        return self._grouped

    def graph_from_template(self, _template_id: str) -> MolGraph:
        if self._graph_error is not None:
            raise self._graph_error
        graph = MolGraph()
        a1 = graph.add_atom("C", 0.0, 0.0)
        a2 = graph.add_atom("O", 40.0, 0.0)
        graph.add_bond(a1.id, a2.id, order=1)
        return graph


class _FakeTemplatesDock:
    def __init__(self) -> None:
        self.calls: list[list[dict]] = []

    def set_templates(self, grouped_templates: list[dict]) -> None:
        self.calls.append(grouped_templates)


def _browser_context(library: object, menu: QMenu | None = None) -> SimpleNamespace:
    parent = QWidget()
    return SimpleNamespace(
        parent=parent,
        template_library=library,
        templates_menu=menu if menu is not None else QMenu(),
        templates_dock=_FakeTemplatesDock(),
        action_save_template=QAction("Guardar", parent),
        action_template_new_category=QAction("Nueva categoría", parent),
        action_template_import_library=QAction("Importar", parent),
        action_template_export_library=QAction("Exportar", parent),
        preview_cache={},
        show_status=lambda _message: None,
        start_template_insert_by_id=lambda _template_id: None,
        insert_template=lambda _name, _graph: None,
        default_category_user=DEFAULT_CATEGORY_USER,
    )


def test_refresh_templates_menu_renders_groups_empty_entries_and_fixed_actions() -> None:
    menu = QMenu()
    library = _FakeTemplateLibrary(
        [
            {"name": "Personal", "templates": [{"id": "tpl_1", "name": "Uno"}]},
            {"name": "Vacía", "templates": []},
        ]
    )
    context = _browser_context(library, menu)

    TemplateBrowserService().refresh_templates_menu(context)

    actions = menu.actions()
    submenus = [action.menu() for action in actions if action.menu() is not None]
    assert [submenu.title() for submenu in submenus] == ["Personal", "Vacía"]
    assert submenus[0].actions()[0].text() == "Uno"
    assert not submenus[1].actions()[0].isEnabled()
    assert submenus[1].actions()[0].text() == "(Vacío)"
    assert [action.text() for action in actions[-4:]] == [
        "Guardar",
        "Nueva categoría",
        "Importar",
        "Exportar",
    ]


def test_refresh_templates_menu_renders_empty_placeholder() -> None:
    menu = QMenu()
    context = _browser_context(_FakeTemplateLibrary([]), menu)

    TemplateBrowserService().refresh_templates_menu(context)

    actions = menu.actions()
    assert actions[0].text() == "(Sin plantillas)"
    assert not actions[0].isEnabled()
    assert [action.text() for action in actions[-4:]] == [
        "Guardar",
        "Nueva categoría",
        "Importar",
        "Exportar",
    ]


def test_refresh_template_views_adds_icons_and_clears_cache() -> None:
    library = _FakeTemplateLibrary(
        [{"name": "Personal", "templates": [{"id": "tpl_1", "name": "Uno"}]}]
    )
    context = _browser_context(library)
    context.preview_cache["stale"] = object()

    TemplateBrowserService().refresh_template_views(context)

    sent = context.templates_dock.calls[-1]
    assert "stale" not in context.preview_cache
    assert sent[0]["templates"][0]["id"] == "tpl_1"
    assert sent[0]["templates"][0]["icon"] is context.preview_cache["tpl_1"]


def test_refresh_template_views_keeps_template_visible_when_preview_fails() -> None:
    library = _FakeTemplateLibrary(
        [{"name": "Personal", "templates": [{"id": "tpl_1", "name": "Uno"}]}],
        graph_error=ValueError("preview rota"),
    )
    context = _browser_context(library)

    TemplateBrowserService().refresh_template_views(context)

    sent = context.templates_dock.calls[-1]
    assert sent[0]["templates"][0]["id"] == "tpl_1"
    assert "icon" in sent[0]["templates"][0]


def test_migrate_legacy_templates_ignores_unreadable_or_empty_files(tmp_path, monkeypatch):
    templates_dir = tmp_path / "legacy"
    templates_dir.mkdir()
    unreadable = templates_dir / "rota.mol"
    unreadable.write_text("MOL ROTO", encoding="utf-8")
    (templates_dir / "vacia.mol").write_text("", encoding="utf-8")
    monkeypatch.setattr(
        TemplateBrowserService,
        "_get_templates_dir",
        staticmethod(lambda: str(templates_dir)),
    )
    real_open = builtins.open

    def _open(path, *args, **kwargs):
        if str(path) == str(unreadable):
            raise OSError("sin permisos")
        return real_open(path, *args, **kwargs)

    monkeypatch.setattr(builtins, "open", _open)
    library = TemplateLibrary(tmp_path / "library.json")
    statuses: list[str] = []

    TemplateBrowserService().migrate_legacy_templates(_migration_context(library, statuses))

    assert not any(template["name"] in {"rota", "vacia"} for template in library.list_templates())
    assert statuses == []


def test_migrate_legacy_templates_ignores_discovery_errors(tmp_path, monkeypatch):
    monkeypatch.setattr(
        TemplateBrowserService,
        "_get_templates_dir",
        staticmethod(lambda: (_ for _ in ()).throw(OSError("sin directorio"))),
    )
    library = TemplateLibrary(tmp_path / "library.json")
    statuses: list[str] = []

    TemplateBrowserService().migrate_legacy_templates(_migration_context(library, statuses))

    assert statuses == []
