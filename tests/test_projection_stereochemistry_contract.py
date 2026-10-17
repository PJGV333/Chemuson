from __future__ import annotations

import math

from chemuson.core.model import BondStereo, BondStyle
from chemuson.gui.projections import (
    ChairParameters,
    ChairSubstituent,
    FischerParameters,
    HaworthParameters,
    HaworthProjectionSemantics,
    build_chair_geometry,
    build_fischer_geometry,
    build_haworth_geometry,
)


def _orientation_by_position(geometry: HaworthProjectionSemantics) -> dict[int, str]:
    return {
        item.ring_position: item.orientation
        for item in geometry.substituents
    }


def _assert_no_bond_stereo(graph) -> None:
    assert {bond.stereo for bond in graph.bonds.values()} == {BondStereo.NONE}


def test_haworth_semantics_preserve_reference_alpha_beta_and_d_l_inversions() -> None:
    beta_d = build_haworth_geometry()
    alpha_d = build_haworth_geometry(HaworthParameters(anomer="alpha"))
    beta_l = build_haworth_geometry(HaworthParameters(configuration="L"))
    alpha_l = build_haworth_geometry(
        HaworthParameters(anomer="alpha", configuration="L")
    )

    assert isinstance(beta_d.semantics, HaworthProjectionSemantics)
    assert beta_d.semantics.anomer == "beta" and beta_d.semantics.configuration == "D"
    assert _orientation_by_position(beta_d.semantics) == {
        1: "up", 2: "down", 3: "up", 4: "down", 5: "up"
    }
    assert _orientation_by_position(alpha_d.semantics) == {
        1: "down", 2: "down", 3: "up", 4: "down", 5: "up"
    }
    assert _orientation_by_position(beta_l.semantics) == {
        1: "down", 2: "up", 3: "down", 4: "up", 5: "down"
    }
    assert _orientation_by_position(alpha_l.semantics) == {
        1: "up", 2: "up", 3: "down", 4: "up", 5: "down"
    }
    assert all(geometry.semantics.front_edge_indices == (1, 2, 3) for geometry in (beta_d, alpha_d, beta_l, alpha_l))


def test_projection_instance_retains_semantics_while_graph_adapter_stays_render_only() -> None:
    geometry = build_haworth_geometry()
    instance = geometry.to_projection_instance()
    legacy_graph = geometry.to_molgraph()

    assert instance.geometry is geometry
    assert instance.semantics is geometry.semantics
    assert instance.render_graph is not legacy_graph
    _assert_no_bond_stereo(instance.render_graph)
    _assert_no_bond_stereo(legacy_graph)
    assert geometry.semantics is not None


def test_fischer_contract_records_crossing_order_without_inference() -> None:
    geometry = build_fischer_geometry(
        FischerParameters(
            substituent_sides=("left", "right"),
            center_configurations=("R", None),
        )
    )
    semantics = geometry.semantics

    assert semantics.horizontal_toward_observer is True
    assert semantics.vertical_away_from_observer is True
    first, second = semantics.centers
    assert first.absolute_configuration == "R"
    assert second.absolute_configuration is None
    assert geometry.atoms[first.left_atom_index].element == "OH"
    assert geometry.atoms[first.right_atom_index].element == "H"
    assert geometry.atoms[second.left_atom_index].element == "H"
    assert geometry.atoms[second.right_atom_index].element == "OH"
    assert first.upper_atom_index == 0 and first.lower_atom_index == second.center_atom_index
    _assert_no_bond_stereo(geometry.to_molgraph())


def test_chair_flip_changes_conformation_without_losing_declared_relations() -> None:
    substituents = (
        ChairSubstituent(0, "OH", "up", "axial"),
        ChairSubstituent(3, "F", "up", "equatorial"),
    )
    parameters = ChairParameters(conformation="chair_a", substituents=substituents)
    chair_a = build_chair_geometry(parameters)
    chair_b = build_chair_geometry(parameters.flipped())

    assert chair_a.semantics.conformation == "chair_a"
    assert chair_b.semantics.conformation == "chair_b"
    assert [(item.orientation, item.placement) for item in chair_a.semantics.substituents] == [
        ("up", "axial"), ("up", "equatorial")
    ]
    assert [(item.orientation, item.placement) for item in chair_b.semantics.substituents] == [
        ("up", "equatorial"), ("up", "axial")
    ]
    assert chair_a.atoms[:6] != chair_b.atoms[:6]
    for geometry in (chair_a, chair_b):
        for item in geometry.semantics.substituents:
            group = geometry.atoms[item.group_atom_index]
            ring = geometry.atoms[item.ring_atom_index]
            assert (group.y < ring.y) == (item.orientation == "up")
        _assert_no_bond_stereo(geometry.to_molgraph())


def test_haworth_canonical_rings_are_distinct_nonregular_and_front_edges_are_explicit() -> None:
    for ring_size in (5, 6):
        geometry = build_haworth_geometry(
            HaworthParameters(ring_size=ring_size, sugar_family="generic")
        )
        points = geometry.atoms[:ring_size]
        headings = []
        for index, point in enumerate(points):
            other = points[(index + 1) % ring_size]
            headings.append(math.atan2(other.y - point.y, other.x - point.x))
        turns = [
            (headings[(index + 1) % ring_size] - headings[index]) % (2.0 * math.pi)
            for index in range(ring_size)
        ]
        assert max(turns) - min(turns) > math.radians(3.0)
        assert geometry.semantics.front_edge_indices == (1, 2, 3)

    custom = build_haworth_geometry(
        HaworthParameters(
            sugar_family="generic",
            front_edge_policy="none",
            front_edge_indices=(0, 4),
            perspective=0.18,
            tilt_degrees=12.0,
        )
    )
    assert custom.semantics.front_edge_indices == (0, 4)
    assert [bond.style for bond in custom.bonds[:6]] == [
        BondStyle.BOLD if index in (0, 4) else BondStyle.PLAIN
        for index in range(6)
    ]
