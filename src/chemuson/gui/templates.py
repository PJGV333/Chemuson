"""Backward-compatible adapters for special molecular projection templates.

The projection package owns geometry construction; the editor renderer owns global
DrawingStyle policy. These wrappers keep the established template API stable.
"""
from __future__ import annotations

from chemuson.core.model import MolGraph
from chemuson.gui.projections import (
    ChairParameters,
    FischerParameters,
    HaworthParameters,
    build_chair_geometry,
    build_fischer_geometry,
    build_haworth_geometry,
)


def build_linear_chain_template(bond_length: float) -> MolGraph:
    """Build the historical render-only Fischer-style glucose graph.

    Use ``build_fischer_geometry(...).to_projection_instance()`` when the
    horizontal/vertical observer convention must remain available in memory.
    """
    length = float(bond_length)
    return build_fischer_geometry(
        FischerParameters(
            vertical_backbone_length=1.3 * length,
            horizontal_cross_length=1.2 * length,
        )
    ).to_molgraph()


def build_haworth_template(
    bond_length: float,
    *,
    anomeric_up: bool = True,
    bold_front: bool = True,
) -> MolGraph:
    """Build the compatible beta-D-glucose Haworth template.

    `anomeric_up` and `bold_front` retain their historical meaning; the returned
    graph is a rendering-only adapter: it carries neither projection semantics
    nor inferred ``BondStereo`` metadata.
    """
    return build_haworth_geometry(
        HaworthParameters(
            bond_length=float(bond_length),
            front_edge_policy="lower_arc" if bold_front else "none",
            anomeric_up_override=anomeric_up,
        )
    ).to_molgraph()


def build_chair_template(
    bond_length: float,
    *,
    anomeric_up: bool = True,
    bold_front: bool = True,
) -> MolGraph:
    """Build the compatible render-only six-carbon chair template.

    Axial/equatorial and up/down semantics are available from the typed
    projection API, not from this legacy ``MolGraph`` adapter.
    """
    # This argument remains for API compatibility; alpha/beta does not apply to
    # an unsubstituted cyclohexane projection.
    _ = anomeric_up
    return build_chair_geometry(
        ChairParameters(
            ring_bond_length=float(bond_length),
            front_edge_policy="lower_arc" if bold_front else "none",
        )
    ).to_molgraph()
