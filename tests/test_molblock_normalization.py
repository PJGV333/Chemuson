"""Pruebas unitarias para normalización de cabecera MOL."""

import math
from types import SimpleNamespace

import pytest

from chemuson.chemio import rdkit_io
from chemuson.chemio.rdkit_io import (
    _should_use_molfile_fallback,
    molfile_to_molgraph,
    molgraph_to_molfile,
    normalize_molblock_header,
)
from chemuson.core.model import MolGraph


MINIMAL_MOLFILE = (
    "Chemuson\n"
    "Chemuson\n"
    "\n"
    "  2  1  0  0  0  0  0  0  0  0  0  0  0  0 V2000\n"
    "    0.0000    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n"
    "   40.0000    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n"
    "  1  2  1  0  0  0  0\n"
    "M  END\n"
)


def _bond_length(graph: MolGraph) -> float:
    bond = next(iter(graph.bonds.values()))
    first = graph.atoms[bond.a1_id]
    second = graph.atoms[bond.a2_id]
    return math.hypot(second.x - first.x, second.y - first.y)


def test_normalize_molblock_header_inserts_missing_comment_line():
    """Inserta línea de comentario si la línea de conteo quedó adelantada."""
    broken = (
        "RDKit          3D\n\n"
        "  2  1  0  0  0  0  0  0  0  0999 V2000\n"
        "   34.6410   20.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n"
        "    0.0000   40.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n"
        "  1  2  1  0\n"
        "M  END\n"
    )
    fixed = normalize_molblock_header(broken)
    lines = fixed.splitlines()
    assert lines[2] == ""
    assert "V2000" in lines[3]


def test_normalize_molblock_header_keeps_valid_header():
    """No altera bloques que ya tienen encabezado CTAB correcto."""
    valid = (
        "Chemuson\n"
        "Chemuson\n"
        "\n"
        "  2  1  0  0  0  0  0  0  0  0  0  0  0  0 V2000\n"
        "    0.0000    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n"
        "   40.0000    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n"
        "  1  2  1  0  0  0  0\n"
        "M  END\n"
    )
    assert normalize_molblock_header(valid) == valid


def test_molfile_to_molgraph_fallback_supports_pseudoatoms():
    """Parsea pseudoátomos sin depender del parser estricto de RDKit."""
    pseudo = (
        "Chemuson\n"
        "Chemuson\n"
        "\n"
        "  2  1  0  0  0  0  0  0  0  0999 V2000\n"
        "    0.0000    0.0000    0.0000 CHO 0  0  0  0  0  0  0  0  0  0  0  0\n"
        "   40.0000    0.0000    0.0000 OCH3 0  0  0  0  0  0  0  0  0  0  0  0\n"
        "  1  2  1  0  0  0  0\n"
        "M  END\n"
    )
    graph = molfile_to_molgraph(pseudo)
    elements = [atom.element for atom in graph.atoms.values()]
    assert elements == ["CHO", "OCH3"]
    assert len(graph.bonds) == 1


def test_molgraph_to_molfile_fallback_writes_alias_and_roundtrips(monkeypatch):
    """El escritor fallback debe emitir alias y poder recuperarlos al importar."""
    graph = MolGraph()
    a1 = graph.add_atom("CHO", 0.0, 0.0, is_explicit=True)
    a2 = graph.add_atom("OCH3", 40.0, 0.0, is_explicit=True)
    graph.add_bond(a1.id, a2.id, order=1)

    monkeypatch.setattr("chemuson.chemio.rdkit_io._rdkit_available", lambda: False)
    molblock = molgraph_to_molfile(graph)
    assert "V2000" in molblock.splitlines()[3]
    assert "\nA  " in molblock
    assert "\nCHO\n" in molblock
    assert "\nOCH3\n" in molblock

    restored = molfile_to_molgraph(molblock)
    restored_elements = [atom.element for atom in restored.atoms.values()]
    assert restored_elements == ["CHO", "OCH3"]


def test_should_use_molfile_fallback_when_alias_lines_present():
    """Si un MOL contiene alias `A  n`, debe usarse el parser fallback."""
    molblock = (
        "Chemuson\n"
        "Chemuson\n"
        "\n"
        "  2  1  0  0  0  0  0  0  0  0999 V2000\n"
        "    0.0000    0.0000    0.0000 *   0  0  0  0  0  0  0  0  0  0  0  0\n"
        "   40.0000    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n"
        "  1  2  1  0  0  0  0\n"
        "A    1\n"
        "OH\n"
        "M  END\n"
    )
    assert _should_use_molfile_fallback(molblock)


def test_molfile_internal_parser_applies_requested_target_bond_length():
    graph = molfile_to_molgraph(MINIMAL_MOLFILE, target_bond_length=55.0)

    assert _bond_length(graph) == pytest.approx(55.0)


def test_molfile_simulated_rdkit_parser_applies_requested_target_bond_length(monkeypatch):
    fallback_error = ValueError("internal parser failed")
    graph = MolGraph()
    first = graph.add_atom("C", 0.0, 0.0)
    second = graph.add_atom("C", 20.0, 0.0)
    graph.add_bond(first.id, second.id)
    fake_mol = object()
    fake_chem = SimpleNamespace(
        MolFromMolBlock=lambda _text, sanitize: fake_mol,
        SanitizeMol=lambda _mol, _flags: None,
        SanitizeFlags=SimpleNamespace(SANITIZE_ALL=3, SANITIZE_KEKULIZE=1),
    )

    monkeypatch.setattr(rdkit_io, "_molfile_to_molgraph_fallback", lambda _text: (_ for _ in ()).throw(fallback_error))
    monkeypatch.setattr(rdkit_io, "_rdkit_available", lambda: True)
    monkeypatch.setattr(rdkit_io, "Chem", fake_chem)
    monkeypatch.setattr(rdkit_io, "rdkit_to_molgraph", lambda mol: graph if mol is fake_mol else None)

    result = molfile_to_molgraph(MINIMAL_MOLFILE, target_bond_length=55.0)

    assert result is graph
    assert _bond_length(graph) == pytest.approx(55.0)


def test_molfile_preserves_internal_parser_error_when_both_parsers_fail(monkeypatch):
    fallback_error = ValueError("internal parser failed")

    monkeypatch.setattr(rdkit_io, "_molfile_to_molgraph_fallback", lambda _text: (_ for _ in ()).throw(fallback_error))
    monkeypatch.setattr(rdkit_io, "_rdkit_available", lambda: False)

    with pytest.raises(ValueError, match="internal parser failed") as exc_info:
        molfile_to_molgraph(MINIMAL_MOLFILE)

    assert str(exc_info.value) == "internal parser failed"
