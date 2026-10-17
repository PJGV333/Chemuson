from __future__ import annotations

import pytest

from chemuson.core.model import BondStyle, BondStereo
from chemuson.gui.projections import (
    ChairParameters,
    ChairSubstituent,
    FischerParameters,
    HaworthParameters,
    HaworthSubstituent,
    ProjectionAtom,
    ProjectionBond,
    ProjectionGeometry,
    build_chair_geometry,
    build_fischer_geometry,
    build_haworth_geometry,
)
from chemuson.gui.templates import (
    build_chair_template,
    build_haworth_template,
    build_linear_chain_template,
)


def _stereo_values(graph):
    return {bond.stereo for bond in graph.bonds.values()}


def _substituent_vertical_directions(graph, ring_size: int) -> dict[int, str]:
    ring_atoms = list(graph.atoms.values())[:ring_size]
    result = {}
    for bond in graph.bonds.values():
        if bond.a1_id not in {atom.id for atom in ring_atoms}:
            continue
        group = graph.get_atom(bond.a2_id)
        is_hydroxyl = group.element == "O" and group.group_h_cap == 1
        if not is_hydroxyl and group.element != "CH2OH":
            continue
        carbon_position = next(
            index for index, atom in enumerate(ring_atoms) if atom.id == bond.a1_id
        )
        result[carbon_position] = "up" if group.y < graph.get_atom(bond.a1_id).y else "down"
    return result


def test_haworth_beta_d_glucose_geometry_is_deterministic_and_explicit() -> None:
    parameters = HaworthParameters()
    first = build_haworth_geometry(parameters)
    second = build_haworth_geometry(parameters)
    graph = first.to_molgraph()
    repeat = second.to_molgraph()

    assert first == second
    assert tuple(vars(atom) for atom in graph.atoms.values()) == tuple(
        vars(atom) for atom in repeat.atoms.values()
    )
    assert len(graph.atoms) == 11
    assert len(graph.bonds) == 11
    assert graph.get_atom(1).element == "O"
    assert _stereo_values(graph) == {BondStereo.NONE}
    assert _substituent_vertical_directions(graph, 6) == {
        1: "up",
        2: "down",
        3: "up",
        4: "down",
        5: "up",
    }
    ring_styles = [bond.style for bond in graph.bonds.values()][:6]
    assert ring_styles == [
        BondStyle.PLAIN,
        BondStyle.BOLD,
        BondStyle.BOLD,
        BondStyle.BOLD,
        BondStyle.PLAIN,
        BondStyle.PLAIN,
    ]


def test_haworth_alpha_beta_and_d_l_are_explicit_projection_inputs() -> None:
    alpha_d = build_haworth_geometry(HaworthParameters(anomer="alpha")).to_molgraph()
    beta_l = build_haworth_geometry(
        HaworthParameters(anomer="beta", configuration="L")
    ).to_molgraph()
    alpha_l = build_haworth_geometry(
        HaworthParameters(anomer="alpha", configuration="L")
    ).to_molgraph()

    assert _substituent_vertical_directions(alpha_d, 6)[1] == "down"
    assert _substituent_vertical_directions(beta_l, 6)[1] == "down"
    assert _substituent_vertical_directions(alpha_l, 6)[1] == "up"
    assert _stereo_values(alpha_d) == _stereo_values(beta_l) == {BondStereo.NONE}


def test_generic_furanose_and_custom_substituent_policy_do_not_infer_stereo() -> None:
    parameters = HaworthParameters(
        ring_size=5,
        sugar_family="generic",
        front_edge_policy="none",
        substituents=(HaworthSubstituent(2, "OH", "up"),),
    )
    geometry = build_haworth_geometry(parameters)
    graph = geometry.to_molgraph()

    assert geometry.family == "haworth"
    assert len(graph.atoms) == 6
    assert len(graph.bonds) == 6
    assert all(bond.style == BondStyle.PLAIN for bond in graph.bonds.values())
    assert _substituent_vertical_directions(graph, 5)[2] == "up"
    assert _stereo_values(graph) == {BondStereo.NONE}


def test_haworth_parameters_reject_invalid_projection_geometry() -> None:
    with pytest.raises(ValueError, match="ring_size"):
        HaworthParameters(ring_size=7)
    with pytest.raises(ValueError, match="aspect_ratio"):
        HaworthParameters(aspect_ratio=3.0)
    with pytest.raises(ValueError, match="ring carbon"):
        HaworthParameters(substituents=(HaworthSubstituent(6, "OH", "up"),))


def test_fischer_constructor_has_canonical_backbone_and_editable_crosses() -> None:
    geometry = build_fischer_geometry()
    graph = geometry.to_molgraph()
    customized = build_fischer_geometry(
        FischerParameters(
            vertical_backbone_length=72.0,
            horizontal_cross_length=58.0,
            substituent_sides=("left", "right"),
            top_terminal_group="CHO",
            bottom_terminal_group="CH3",
        )
    )

    assert geometry == build_fischer_geometry(FischerParameters())
    assert len(graph.atoms) == 14
    assert len(graph.bonds) == 13
    assert graph.atoms[1].element == "CHO"
    assert graph.atoms[6].element == "CH2OH"
    assert _stereo_values(graph) == {BondStereo.NONE}
    assert len(customized.atoms) == 8
    assert len(customized.bonds) == 7


def test_chair_constructor_exposes_conformation_vectors_and_substituent_policy() -> None:
    parameters = ChairParameters(
        substituents=(ChairSubstituent(0, "OH", "up", "axial"),),
    )
    geometry = build_chair_geometry(parameters)
    graph = geometry.to_molgraph()
    repeated = build_chair_geometry(parameters)
    flipped = build_chair_geometry(ChairParameters(conformation="chair_b"))

    assert geometry == repeated
    assert len(graph.atoms) == 7
    assert len(graph.bonds) == 7
    substituent = graph.atoms[max(graph.atoms)]
    assert substituent.element == "O" and substituent.group_h_cap == 1
    assert substituent.y < graph.atoms[1].y
    assert len(flipped.atoms) == 6
    assert geometry.family == "chair"
    assert _stereo_values(graph) == {BondStereo.NONE}
    assert len(ChairParameters().axial_vectors) == len(ChairParameters().equatorial_vectors) == 6


def test_chair_rejects_orientation_placement_mismatches_for_the_declared_conformation() -> None:
    with pytest.raises(ValueError, match="inconsistent"):
        ChairParameters(
            conformation="chair_a",
            substituents=(ChairSubstituent(0, "OH", "down", "axial"),),
        )


def test_legacy_template_entry_points_remain_callable_and_deterministic() -> None:
    haworth = build_haworth_template(40.0, anomeric_up=True, bold_front=True)
    repeated = build_haworth_template(40.0, anomeric_up=True, bold_front=True)
    assert tuple(vars(atom) for atom in haworth.atoms.values()) == tuple(
        vars(atom) for atom in repeated.atoms.values()
    )
    assert (len(haworth.atoms), len(haworth.bonds)) == (11, 11)
    assert (len(build_linear_chain_template(40.0).atoms), len(build_linear_chain_template(40.0).bonds)) == (14, 13)
    chair = build_chair_template(40.0, anomeric_up=False, bold_front=False)
    assert (len(chair.atoms), len(chair.bonds)) == (6, 6)
    assert all(bond.style == BondStyle.PLAIN for bond in chair.bonds.values())


def test_projection_geometry_rejects_invalid_endpoints_and_nonfinite_coordinates() -> None:
    with pytest.raises(ValueError, match="finite"):
        ProjectionGeometry(
            "invalid",
            40.0,
            (ProjectionAtom("C", float("nan"), 0.0, False),),
            (),
        )
    with pytest.raises(ValueError, match="out of range"):
        ProjectionGeometry(
            "invalid",
            40.0,
            (ProjectionAtom("C", 0.0, 0.0, False),),
            (ProjectionBond(0, 1),),
        )
    with pytest.raises(ValueError, match="positive integer"):
        ProjectionGeometry(
            "invalid",
            40.0,
            (ProjectionAtom("C", 0.0, 0.0, False), ProjectionAtom("C", 1.0, 0.0, False)),
            (ProjectionBond(0, 1, order=1.5),),
        )
