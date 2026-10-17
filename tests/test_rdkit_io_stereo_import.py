from __future__ import annotations

from chemuson.chemio.rdkit_io import molfile_to_molgraph
from chemuson.core.model import BondStereo, BondStyle


def _v2000_molfile(
    atoms: list[tuple[float, float, str]],
    bonds: list[tuple[int, int, int, int]],
) -> str:
    lines = [
        "ChemUSON stereo import",
        "  ChemUSON          2D",
        "",
        f"{len(atoms):3d}{len(bonds):3d}  0  0  0  0            999 V2000",
    ]
    for x, y, symbol in atoms:
        lines.append(
            f"{x:10.4f}{y:10.4f}{0.0:10.4f} {symbol:<3}"
            " 0  0  0  0  0  0  0  0  0  0  0  0"
        )
    lines.extend(f"{first:3d}{second:3d}{order:3d}{stereo:3d}" for first, second, order, stereo in bonds)
    lines.append("M  END")
    return "\n".join(lines) + "\n"


def test_molfile_import_preserves_wedge_direction_and_cip_identity():
    molfile = _v2000_molfile(
        [
            (1.5, 0.0, "F"),
            (0.0, 0.0, "C"),
            (-1.5, 0.0, "Cl"),
            (0.0, 1.5, "Br"),
            (0.0, -1.5, "I"),
        ],
        [(2, 1, 1, 6), (2, 3, 1, 0), (2, 4, 1, 0), (2, 5, 1, 0)],
    )

    graph = molfile_to_molgraph(molfile)

    hashed = next(bond for bond in graph.bonds.values() if bond.style == BondStyle.HASHED)
    center = graph.atoms[2]
    assert (hashed.a1_id, hashed.a2_id) == (2, 1)
    assert hashed.stereo == BondStereo.DOWN
    assert center.stereo_cip == "S"


def test_molfile_import_keeps_ambiguous_double_bond_mark_explicit():
    molfile = _v2000_molfile(
        [(-1.5, 0.0, "F"), (-0.5, 0.5, "C"), (0.5, -0.5, "C"), (1.5, 0.0, "Cl")],
        [(1, 2, 1, 0), (2, 3, 2, 3), (3, 4, 1, 0)],
    )

    graph = molfile_to_molgraph(molfile)

    double_bond = next(bond for bond in graph.bonds.values() if bond.order == 2)
    assert double_bond.style == BondStyle.WAVY
    assert double_bond.stereo == BondStereo.EITHER
    assert double_bond.stereo_ez is None
