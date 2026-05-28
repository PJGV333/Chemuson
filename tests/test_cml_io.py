"""Pruebas de interoperabilidad CML."""

from __future__ import annotations

import os
import sys

from PyQt6.QtWidgets import QApplication

sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), "..", "src")))

from chemuson.chemio.cml_io import cml_to_molgraph, molgraph_to_cml
from chemuson.core.model import BondStyle, MolGraph
from chemuson.gui.canvas import ChemusonCanvas
from chemuson.gui.controllers.file_controller import FileController


def _sample_graph() -> MolGraph:
    graph = MolGraph()
    n = graph.add_atom("N", 0.0, 0.0, charge=1, is_explicit=True, r_group_substituents=("Me", "Et"))
    c = graph.add_atom("C", 42.0, 0.0)
    o = graph.add_atom("O", 84.0, 0.0, isotope=18, radical_electrons=1, is_explicit=True)
    graph.add_bond(n.id, c.id, order=1, style=BondStyle.PLAIN)
    graph.add_bond(c.id, o.id, order=2)
    return graph


def test_cml_roundtrip_preserves_atoms_bonds_and_metadata() -> None:
    graph = _sample_graph()

    cml = molgraph_to_cml(graph)
    restored = cml_to_molgraph(cml)

    assert "<cml:cml" in cml
    assert len(restored.atoms) == 3
    assert len(restored.bonds) == 2
    restored_n = restored.get_atom(1)
    restored_o = restored.get_atom(3)
    assert restored_n.element == "N"
    assert restored_n.charge == 1
    assert restored_n.r_group_substituents == ("Me", "Et")
    assert restored_o.isotope == 18
    assert restored_o.radical_electrons == 1
    assert restored.find_bond_between(2, 3).order == 2


def test_file_controller_loads_cml_into_canvas(tmp_path) -> None:
    QApplication.instance() or QApplication([])
    path = tmp_path / "sample.cml"
    path.write_text(molgraph_to_cml(_sample_graph()), encoding="utf-8")
    canvas = ChemusonCanvas()

    FileController().load_file_into_canvas(str(path), canvas)

    assert len(canvas.model.atoms) == 3
    assert len(canvas.model.bonds) == 2


def test_cml_import_accepts_minimal_external_cml() -> None:
    cml = """<?xml version='1.0'?>
    <cml xmlns='http://www.xml-cml.org/schema'>
      <molecule>
        <atomArray>
          <atom id='a1' elementType='C' x2='0' y2='0'/>
          <atom id='a2' elementType='O' x2='1' y2='0' formalCharge='-1'/>
        </atomArray>
        <bondArray><bond id='b1' atomRefs2='a1 a2' order='2'/></bondArray>
      </molecule>
    </cml>
    """

    graph = cml_to_molgraph(cml)

    assert len(graph.atoms) == 2
    assert graph.get_atom(2).charge == -1
    assert graph.find_bond_between(1, 2).order == 2


def test_cml_import_accepts_minimal_external_cml_without_namespace() -> None:
    cml = """
    <cml>
      <molecule>
        <atomArray><atom id='a1' elementType='C'/><atom id='a2' elementType='N'/></atomArray>
        <bondArray><bond atomRefs2='a1 a2' order='3'/></bondArray>
      </molecule>
    </cml>
    """

    graph = cml_to_molgraph(cml)

    assert len(graph.atoms) == 2
    assert graph.find_bond_between(1, 2).order == 3
