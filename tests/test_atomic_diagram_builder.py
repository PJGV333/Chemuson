"""Pruebas para el builder de diagramas atómicos."""

from __future__ import annotations



from chemuson.gui.energy_diagrams import build_atomic_subshell_diagram


def test_atomic_diagram_two_electrons_is_1s2() -> None:
    diagram = build_atomic_subshell_diagram(2)
    level_map = {level.id: level for level in diagram.levels}

    assert [level.id for level in diagram.levels] == ["1s"]
    assert level_map["1s"].occupancies == [2]


def test_atomic_diagram_seven_electrons_applies_hund_rule() -> None:
    diagram = build_atomic_subshell_diagram(7)
    level_map = {level.id: level for level in diagram.levels}

    assert [level.id for level in diagram.levels] == ["1s", "2s", "2p"]
    assert level_map["1s"].occupancies == [2]
    assert level_map["2s"].occupancies == [2]
    assert level_map["2p"].occupancies == [1, 1, 1]
    assert diagram.metadata["electron_count"] == 7
    assert diagram.metadata["filling_rule"] == "aufbau_hund"


def test_atomic_diagram_twenty_six_electrons_fills_expected_subshells() -> None:
    diagram = build_atomic_subshell_diagram(26)

    assert [level.id for level in diagram.levels] == ["1s", "2s", "2p", "3s", "3p", "4s", "3d"]
    assert len(diagram.levels) == 7
    assert sum(sum(level.occupancies) for level in diagram.levels) == 26
