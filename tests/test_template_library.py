"""Pruebas para biblioteca de plantillas de usuario."""

import os
import sys

import pytest

sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), "..", "src")))

from core.model import MolGraph
from gui.template_library import DEFAULT_CATEGORY_USER, TemplateLibrary


def _simple_cc_graph() -> MolGraph:
    graph = MolGraph()
    a1 = graph.add_atom("C", 0.0, 0.0, is_explicit=False)
    a2 = graph.add_atom("C", 40.0, 0.0, is_explicit=False)
    graph.add_bond(a1.id, a2.id, order=1)
    return graph


def test_template_library_bootstraps_defaults(tmp_path):
    """Debe crear categorías/plantillas iniciales cuando no existe archivo."""
    lib = TemplateLibrary(tmp_path / "library.json")
    categories = lib.categories()
    assert "Aromáticos" in categories
    assert DEFAULT_CATEGORY_USER in categories
    names = {tpl["name"] for tpl in lib.list_templates()}
    assert "Benceno" in names


def test_add_rename_and_delete_template_and_category(tmp_path):
    """Operaciones básicas de gestión deben persistir sin inconsistencias."""
    lib = TemplateLibrary(tmp_path / "library.json")
    tpl = lib.add_template_from_graph("Mi plantilla", "Personal", _simple_cc_graph())
    tpl_id = tpl["id"]

    renamed = lib.rename_template(tpl_id, "Plantilla renombrada")
    assert renamed == "Plantilla renombrada"

    new_category = lib.rename_category("Personal", "Orgánicos")
    assert new_category == "Orgánicos"
    updated = lib.get_template(tpl_id)
    assert updated["category"] == "Orgánicos"
    assert updated["name"] == "Plantilla renombrada"

    lib.delete_template(tpl_id)
    with pytest.raises(ValueError):
        lib.get_template(tpl_id)


def test_import_export_template_library_merge(tmp_path):
    """Importar en modo merge debe incorporar nuevas plantillas sin borrar las existentes."""
    source = TemplateLibrary(tmp_path / "source.json")
    source.add_template_from_graph("Cadena prueba", "Compartidas", _simple_cc_graph())
    export_path = tmp_path / "exported.json"
    source.export_to_file(export_path)

    target = TemplateLibrary(tmp_path / "target.json")
    before = len(target.list_templates())
    added = target.import_from_file(export_path, merge=True)
    after = len(target.list_templates())

    assert added >= 1
    assert after >= before + 1
    assert any(tpl["name"] == "Cadena prueba" for tpl in target.list_templates())
