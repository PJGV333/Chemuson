"""Pruebas para biblioteca de plantillas de usuario."""

import json
import os
import sys
from pathlib import Path

import pytest

sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), "..", "src")))

import gui.template_library as template_library
from chemio.rdkit_io import molgraph_to_molfile
from core.model import MolGraph
from gui.template_library import DEFAULT_CATEGORY_USER, TemplateLibrary
from gui.templates import build_haworth_template


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


def test_builtin_chair_is_simple_cyclohexane_chair(tmp_path):
    """La plantilla Silla β debe ser un anillo silla sin heteroátomos ni sustituyentes."""
    lib = TemplateLibrary(tmp_path / "library.json")
    chairs = [
        tpl
        for tpl in lib.list_templates()
        if tpl.get("name") == "Silla β" and tpl.get("category") == "Bioquímicos"
    ]
    assert len(chairs) == 1
    chair_graph = lib.graph_from_template(chairs[0]["id"])
    assert len(chair_graph.atoms) == 6
    assert len(chair_graph.bonds) == 6
    assert all(atom.element == "C" for atom in chair_graph.atoms.values())


def test_load_replaces_legacy_builtin_chair_template(tmp_path):
    """Al cargar biblioteca existente, Silla β debe sincronizarse a la versión corregida."""
    library_path = tmp_path / "library.json"
    legacy_graph = build_haworth_template(40.0, anomeric_up=True, bold_front=True)
    raw = {
        "version": 1,
        "categories": ["Bioquímicos", DEFAULT_CATEGORY_USER],
        "templates": [
            {
                "id": "tpl_legacy_chair",
                "name": "Silla β",
                "category": "Bioquímicos",
                "molblock": molgraph_to_molfile(legacy_graph),
                "smiles": "",
            }
        ],
    }
    library_path.write_text(json.dumps(raw, ensure_ascii=False), encoding="utf-8")

    lib = TemplateLibrary(library_path)
    chairs = [
        tpl
        for tpl in lib.list_templates()
        if tpl.get("name") == "Silla β" and tpl.get("category") == "Bioquímicos"
    ]
    assert len(chairs) == 1
    assert chairs[0]["id"] == "tpl_legacy_chair"
    chair_graph = lib.graph_from_template(chairs[0]["id"])
    assert len(chair_graph.atoms) == 6
    assert len(chair_graph.bonds) == 6
    assert all(atom.element == "C" for atom in chair_graph.atoms.values())


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


def test_normalize_repair_molblock_header_on_load(tmp_path):
    """Debe reparar MOL con línea de comentario CTAB faltante."""
    library_path = tmp_path / "library.json"
    broken_molblock = (
        "RDKit          3D\n\n"
        "  2  1  0  0  0  0  0  0  0  0999 V2000\n"
        "   34.6410   20.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n"
        "    0.0000   40.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n"
        "  1  2  1  0\n"
        "M  END\n"
    )
    raw = {
        "version": 1,
        "categories": [DEFAULT_CATEGORY_USER],
        "templates": [
            {
                "id": "tpl_header_fix",
                "name": "HeaderFix",
                "category": DEFAULT_CATEGORY_USER,
                "molblock": broken_molblock,
                "smiles": "CC",
            }
        ],
    }
    library_path.write_text(json.dumps(raw, ensure_ascii=False), encoding="utf-8")
    lib = TemplateLibrary(library_path)

    molblock = lib.get_template("tpl_header_fix")["molblock"]
    lines = molblock.splitlines()
    assert lines[2] == ""
    assert "V2000" in lines[3]


def test_graph_from_template_falls_back_to_smiles_when_molblock_fails(tmp_path, monkeypatch):
    """Si el MOL falla, debe usar SMILES como respaldo."""
    lib = TemplateLibrary(tmp_path / "library.json")
    tpl = lib.add_template(
        "Fallback",
        DEFAULT_CATEGORY_USER,
        "MOL INVALIDO",
        smiles="C",
    )

    def _fail_mol(_molblock: str):
        raise ValueError("Mol inválido")

    def _from_smiles(_smiles: str):
        graph = MolGraph()
        graph.add_atom("C", 0.0, 0.0)
        return graph

    monkeypatch.setattr("gui.template_library.molfile_to_molgraph", _fail_mol)
    monkeypatch.setattr("gui.template_library.smiles_to_molgraph", _from_smiles)

    graph = lib.graph_from_template(tpl["id"])
    assert len(graph.atoms) == 1
    assert len(graph.bonds) == 0


def test_default_library_path_linux_uses_xdg(monkeypatch):
    """En Linux debe preferir XDG_CONFIG_HOME."""
    monkeypatch.setattr(template_library, "_is_windows_platform", lambda: False)
    monkeypatch.delenv("CHEMUSON_CONFIG_HOME", raising=False)
    monkeypatch.setenv("XDG_CONFIG_HOME", "/tmp/chemuson-xdg")
    path = template_library._default_library_path()
    assert path == Path("/tmp/chemuson-xdg") / "Chemuson" / "template_library.json"


def test_default_library_path_windows_uses_appdata(monkeypatch):
    """En Windows debe guardar en APPDATA cuando está disponible."""
    monkeypatch.setattr(template_library, "_is_windows_platform", lambda: True)
    monkeypatch.delenv("CHEMUSON_CONFIG_HOME", raising=False)
    monkeypatch.setenv("APPDATA", r"C:\Users\Ana\AppData\Roaming")
    monkeypatch.delenv("LOCALAPPDATA", raising=False)
    monkeypatch.delenv("XDG_CONFIG_HOME", raising=False)
    path = template_library._default_library_path()
    assert str(path).startswith(r"C:\Users\Ana\AppData\Roaming")
    assert str(path).endswith("Chemuson/template_library.json")


def test_default_library_path_honors_chemuson_config_home(monkeypatch):
    """CHEMUSON_CONFIG_HOME debe tener prioridad en cualquier plataforma."""
    monkeypatch.setattr(template_library, "_is_windows_platform", lambda: True)
    monkeypatch.setenv("CHEMUSON_CONFIG_HOME", "/tmp/chemuson-portable")
    monkeypatch.setenv("APPDATA", r"C:\Users\Ana\AppData\Roaming")
    path = template_library._default_library_path()
    assert path == Path("/tmp/chemuson-portable") / "Chemuson" / "template_library.json"
