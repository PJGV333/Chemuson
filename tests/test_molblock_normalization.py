"""Pruebas unitarias para normalización de cabecera MOL."""

import os
import sys

sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), "..", "src")))

from chemuson.chemio.rdkit_io import (
    _should_use_molfile_fallback,
    molfile_to_molgraph,
    molgraph_to_molfile,
    normalize_molblock_header,
)
from chemuson.core.model import MolGraph


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
