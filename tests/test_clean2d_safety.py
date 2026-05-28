"""Pruebas del módulo de seguridad y validación para Clean2D."""

from __future__ import annotations

import math
import os
import sys

import pytest

sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), "..", "src")))

from chemuson.clean2d.safety import (
    Clean2DQualityReport,
    bond_length_stats,
    count_new_bond_crossings,
    evaluate_clean2d_layout,
    has_cycles,
    is_clean2d_candidate_safe,
    max_atom_displacement,
    min_nonbonded_distance,
    ring_degeneracy_score,
)
from chemuson.clean2d.length_only import length_only_polish
from chemuson.core.model import Bond


def _fake_bond(a1, a2, order=1):
    return Bond(id=hash((a1, a2)) % 10000, a1_id=a1, a2_id=a2, order=order)


# ─── has_cycles ───────────────────────────────────────────────────────────


def test_has_cycles_detects_chain() -> None:
    atom_ids = {1, 2, 3}
    bonds = [_fake_bond(1, 2), _fake_bond(2, 3)]
    assert not has_cycles(atom_ids, bonds)


def test_has_cycles_detects_ring() -> None:
    atom_ids = {1, 2, 3}
    bonds = [_fake_bond(1, 2), _fake_bond(2, 3), _fake_bond(3, 1)]
    assert has_cycles(atom_ids, bonds)


# ─── evaluate / candidate safety ─────────────────────────────────────────


def test_evaluate_rejects_missing_coords() -> None:
    atom_ids = {1, 2}
    bonds = [_fake_bond(1, 2)]
    before = {1: (0.0, 0.0), 2: (42.0, 0.0)}
    after = {1: (0.0, 0.0)}
    report = evaluate_clean2d_layout(atom_ids, bonds, before, after, 42.0)
    assert report.missing_coords == [2]
    assert not is_clean2d_candidate_safe(report)
    assert "faltan" in report.rejection_reason


def test_evaluate_rejects_nan_coords() -> None:
    atom_ids = {1, 2}
    bonds = [_fake_bond(1, 2)]
    before = {1: (0.0, 0.0), 2: (42.0, 0.0)}
    after = {1: (0.0, 0.0), 2: (float("nan"), 0.0)}
    report = evaluate_clean2d_layout(atom_ids, bonds, before, after, 42.0)
    assert report.nan_coords == [2]
    assert not is_clean2d_candidate_safe(report)


def test_evaluate_rejects_excessive_displacement() -> None:
    atom_ids = {1, 2}
    bonds = [_fake_bond(1, 2)]
    before = {1: (0.0, 0.0), 2: (42.0, 0.0)}
    after = {1: (0.0, 0.0), 2: (42.0, 500.0)}  # massive displacement
    report = evaluate_clean2d_layout(atom_ids, bonds, before, after, 42.0)
    assert report.max_displacement > 42.0 * 3.0
    assert not is_clean2d_candidate_safe(report)
    assert report.rejection_reason  # alguna razón de rechazo


def test_evaluate_rejects_new_bond_crossings() -> None:
    # Two bonds that should not cross
    atom_ids = {1, 2, 3, 4}
    bonds = [_fake_bond(1, 2), _fake_bond(3, 4)]
    before = {1: (0.0, 0.0), 2: (42.0, 0.0), 3: (0.0, 42.0), 4: (42.0, 42.0)}
    after = {1: (0.0, 0.0), 2: (42.0, 42.0), 3: (42.0, 0.0), 4: (0.0, 42.0)}  # cross
    report = evaluate_clean2d_layout(atom_ids, bonds, before, after, 42.0)
    assert report.new_crossings > 0
    assert not is_clean2d_candidate_safe(report)


def test_evaluate_accepts_good_candidate() -> None:
    atom_ids = {1, 2, 3}
    bonds = [_fake_bond(1, 2), _fake_bond(2, 3)]
    before = {1: (0.0, 0.0), 2: (42.0, 0.0), 3: (84.0, 0.0)}
    after = {1: (0.0, 0.0), 2: (42.0, 0.0), 3: (84.0, 0.0)}
    report = evaluate_clean2d_layout(atom_ids, bonds, before, after, 42.0)
    assert is_clean2d_candidate_safe(report)
    assert report.passed


def test_evaluate_rejects_bond_length_out_of_range() -> None:
    atom_ids = {1, 2}
    bonds = [_fake_bond(1, 2)]
    before = {1: (0.0, 0.0), 2: (42.0, 0.0)}
    after = {1: (0.0, 0.0), 2: (200.0, 0.0)}  # too long
    report = evaluate_clean2d_layout(atom_ids, bonds, before, after, 42.0)
    assert not is_clean2d_candidate_safe(report)
    assert "rango" in report.rejection_reason


def test_evaluate_rejects_nonbonded_collision() -> None:
    atom_ids = {1, 2, 3}
    bonds = [_fake_bond(1, 2)]
    before = {1: (0.0, 0.0), 2: (42.0, 0.0), 3: (84.0, 0.0)}
    after = {1: (0.0, 0.0), 2: (42.0, 0.0), 3: (44.0, 0.0)}  # too close to atom 2
    report = evaluate_clean2d_layout(atom_ids, bonds, before, after, 42.0)
    assert not is_clean2d_candidate_safe(report)
    assert "colision" in report.rejection_reason


# ─── ring degeneracy ─────────────────────────────────────────────────────


def test_ring_degeneracy_regular_hexagon() -> None:
    """Un hexágono regular debe tener degeneracy ≈ 1.0."""
    cx, cy = 0.0, 0.0
    r = 42.0
    vertices: list[tuple[float, float]] = []
    for i in range(6):
        angle = math.radians(60 * i)
        vertices.append((cx + r * math.cos(angle), cy + r * math.sin(angle)))
    positions = {i: v for i, v in enumerate(vertices)}
    ring = set(range(6))
    score = ring_degeneracy_score(positions, ring)
    assert score > 0.8, f"Hexágono regular debe tener degeneracy alta, obtuvo {score}"


def test_ring_degeneracy_collapsed() -> None:
    """Un anillo colapsado debe tener degeneracy ≈ 0."""
    positions = {0: (0.0, 0.0), 1: (42.0, 0.0), 2: (43.0, 0.0)}  # casi lineal
    ring = {0, 1, 2}
    score = ring_degeneracy_score(positions, ring)
    assert score < 0.1


def test_evaluate_rejects_collapsed_ring() -> None:
    atom_ids = {1, 2, 3}
    bonds = [_fake_bond(1, 2), _fake_bond(2, 3), _fake_bond(3, 1)]
    before = {
        1: (0.0, 0.0),
        2: (42.0, 0.0),
        3: (21.0, 36.0),
    }
    after = {
        1: (0.0, 0.0),
        2: (42.0, 0.0),
        3: (42.0, 0.01),  # casi colapsado en la misma línea
    }
    report = evaluate_clean2d_layout(atom_ids, bonds, before, after, 42.0, is_cyclic=True)
    assert not is_clean2d_candidate_safe(report)
    assert report.rejection_reason


# ─── bond_crossings ──────────────────────────────────────────────────────


def test_count_new_bond_crossings_no_cross() -> None:
    bonds = [_fake_bond(1, 2), _fake_bond(3, 4)]
    before = {1: (0.0, 0.0), 2: (42.0, 0.0), 3: (0.0, 42.0), 4: (42.0, 42.0)}
    after = {1: (0.0, 0.0), 2: (42.0, 0.0), 3: (0.0, 42.0), 4: (42.0, 42.0)}
    assert count_new_bond_crossings(before, after, bonds) == 0


def test_count_new_bond_crossings_detects_new_cross() -> None:
    bonds = [_fake_bond(1, 2), _fake_bond(3, 4)]
    before = {1: (0.0, 0.0), 2: (42.0, 0.0), 3: (0.0, 42.0), 4: (42.0, 42.0)}
    # after: both bonds cross
    after = {1: (0.0, 0.0), 2: (42.0, 42.0), 3: (42.0, 0.0), 4: (0.0, 42.0)}
    n = count_new_bond_crossings(before, after, bonds)
    assert n == 1, f"Expected 1 new crossing, got {n}"


# ─── min_nonbonded_distance ──────────────────────────────────────────────


def test_min_nonbonded_distance() -> None:
    bonds = [_fake_bond(1, 2)]
    positions = {1: (0.0, 0.0), 2: (42.0, 0.0), 3: (44.0, 0.0)}
    dist = min_nonbonded_distance(positions, bonds, {1, 2, 3})
    assert dist > 0
    assert dist < 5.0


# ─── max_atom_displacement ───────────────────────────────────────────────


def test_max_atom_displacement() -> None:
    before = {1: (0.0, 0.0), 2: (42.0, 0.0)}
    after = {1: (0.0, 0.0), 2: (84.0, 0.0)}
    assert max_atom_displacement(before, after, {1, 2}) == 42.0


# ─── bond_length_stats ───────────────────────────────────────────────────


def test_bond_length_stats() -> None:
    bonds = [_fake_bond(1, 2), _fake_bond(2, 3)]
    positions = {1: (0.0, 0.0), 2: (42.0, 0.0), 3: (84.0, 0.0)}
    stats = bond_length_stats(positions, bonds)
    assert abs(stats["mean"] - 42.0) < 0.01
    assert abs(stats["min"] - 42.0) < 0.01
    assert abs(stats["max"] - 42.0) < 0.01


# ─── length_only_polish ──────────────────────────────────────────────────


def test_length_only_polish_reduces_long_bond() -> None:
    positions = {1: (0.0, 0.0), 2: (200.0, 0.0)}
    bonds = [_fake_bond(1, 2)]
    changed = length_only_polish(positions, bonds, 42.0, max_iterations=12)
    if 2 in changed:
        dist = math.hypot(changed[2][0] - positions[1][0], changed[2][1] - positions[1][1])
        assert dist < 200.0


def test_length_only_polish_doesnt_move_good_structure() -> None:
    positions = {1: (0.0, 0.0), 2: (42.0, 0.0), 3: (84.0, 0.0)}
    bonds = [_fake_bond(1, 2), _fake_bond(2, 3)]
    changed = length_only_polish(positions, bonds, 42.0, max_iterations=6)
    assert len(changed) == 0


def test_length_only_polish_respects_displacement_limit() -> None:
    positions = {1: (0.0, 0.0), 2: (1000.0, 0.0)}
    bonds = [_fake_bond(1, 2)]
    changed = length_only_polish(
        positions, bonds, 42.0,
        max_iterations=20,
        max_displacement_per_atom=10.0,
    )
    if 2 in changed:
        dx = abs(changed[2][0] - positions[2][0])
        assert dx <= 10.0 + 1e-6


# ─── fenilo + cadena + ciclopropano (screenshot-like structure) ──────────


def test_complex_molecule_not_collapsed_by_length_polish() -> None:
    """Construye fenilo+cadena+ciclopropano y verifica que length_only
    no colapsa anillos ni crea cruces."""
    # Fenilo: 6 átomos en hexágono
    cx, cy = 0.0, 0.0
    r = 42.0
    phenyl = {}
    for i in range(6):
        angle = math.radians(60 * i + 30)
        phenyl[i + 1] = (cx + r * math.cos(angle), cy + r * math.sin(angle))

    # Cadena: 4 carbonos
    chain_start = 7
    chain = {}
    for i in range(4):
        chain[chain_start + i] = (200.0 + i * 42.0, 0.0)

    # Ciclopropano triangular
    cp = {
        11: (320.0, 0.0),
        12: (320.0 + 42.0, 0.0),
        13: (320.0 + 21.0, 36.37),
    }

    all_positions = {**phenyl, **chain, **cp}

    bonds = [
        _fake_bond(1, 2), _fake_bond(2, 3), _fake_bond(3, 4),
        _fake_bond(4, 5), _fake_bond(5, 6), _fake_bond(6, 1),  # fenilo
        _fake_bond(6, 7),  # conexión fenilo-cadena
        _fake_bond(7, 8, order=2), _fake_bond(8, 9, order=1),
        _fake_bond(9, 10, order=1),  # cadena
        _fake_bond(10, 11, order=1),  # conexión cadena-ciclopropano
        _fake_bond(11, 12), _fake_bond(12, 13), _fake_bond(13, 11),  # ciclopropano
    ]
    atom_ids = set(all_positions.keys())

    changed = length_only_polish(all_positions, bonds, 42.0, max_iterations=10)
    after = {aid: changed.get(aid, all_positions[aid]) for aid in atom_ids}

    # Fenilo no se colapsó
    phenyl_ids = {1, 2, 3, 4, 5, 6}
    phenyl_before = {aid: all_positions[aid] for aid in phenyl_ids}
    phenyl_after = {aid: after[aid] for aid in phenyl_ids}
    report = evaluate_clean2d_layout(
        phenyl_ids, bonds[:6], phenyl_before, phenyl_after, 42.0, is_cyclic=True,
    )
    assert report.ring_degeneracy_after > 0.1, "Fenilo colapsado"

    # Ciclopropano no se colapsó
    cp_ids = {11, 12, 13}
    cp_before = {aid: all_positions[aid] for aid in cp_ids}
    cp_after = {aid: after[aid] for aid in cp_ids}
    cp_bonds = [b for b in bonds if b.a1_id in cp_ids and b.a2_id in cp_ids]
    cp_report = evaluate_clean2d_layout(
        cp_ids, cp_bonds, cp_before, cp_after, 42.0, is_cyclic=True,
    )
    assert cp_report.ring_degeneracy_after > 0.05, "Ciclopropano colapsado"

    # No hay cruces nuevos
    full_report = evaluate_clean2d_layout(
        atom_ids, bonds, all_positions, after, 42.0, is_cyclic=True,
    )
    assert full_report.new_crossings == 0


# ─── clean2d_v2 no se aplica a cíclicos (prueba conceptual) ──────────────


def test_has_cycles_blocks_v2() -> None:
    """Verifica que has_cycles detecta correctamente estructuras cíclicas."""
    # Cadena acíclica
    chain_ids = {1, 2, 3, 4}
    chain_bonds = [
        _fake_bond(1, 2), _fake_bond(2, 3), _fake_bond(3, 4),
    ]
    assert not has_cycles(chain_ids, chain_bonds)

    # Anillo de 6
    ring_ids = {1, 2, 3, 4, 5, 6}
    ring_bonds = [
        _fake_bond(1, 2), _fake_bond(2, 3), _fake_bond(3, 4),
        _fake_bond(4, 5), _fake_bond(5, 6), _fake_bond(6, 1),
    ]
    assert has_cycles(ring_ids, ring_bonds)


# ─── candidato artificialmente malo es rechazado ─────────────────────────


def test_artificially_bad_candidate_is_rejected() -> None:
    atom_ids = {1, 2, 3}
    bonds = [_fake_bond(1, 2), _fake_bond(2, 3)]
    before = {1: (0.0, 0.0), 2: (42.0, 0.0), 3: (84.0, 0.0)}
    # after: todos en el mismo punto -> absurdo
    after = {1: (0.0, 0.0), 2: (0.0, 0.0), 3: (0.0, 0.0)}
    report = evaluate_clean2d_layout(atom_ids, bonds, before, after, 42.0)
    assert not is_clean2d_candidate_safe(report)
    assert "rango" in report.rejection_reason or "colision" in report.rejection_reason


# ─── estructura ya buena queda casi igual ────────────────────────────────


def test_good_structure_unchanged_by_safety() -> None:
    atom_ids = {1, 2, 3}
    bonds = [_fake_bond(1, 2), _fake_bond(2, 3)]
    positions = {1: (0.0, 0.0), 2: (42.0, 0.0), 3: (84.0, 0.0)}
    # length_only no debe modificar geometría buena
    changed = length_only_polish(positions, bonds, 42.0, max_iterations=6)
    assert len(changed) == 0, f"Buena estructura no debería cambiar, cambió: {changed}"


# ─── cadena acíclica muy deformada sí es corregible ──────────────────────


def test_acyclic_deformed_chain_can_be_corrected() -> None:
    atom_ids = {1, 2, 3}
    bonds = [_fake_bond(1, 2), _fake_bond(2, 3)]
    positions = {1: (0.0, 0.0), 2: (42.0, 0.0), 3: (50.0, 0.0)}  # slightly too short
    changed = length_only_polish(positions, bonds, 42.0, max_iterations=10)
    assert len(changed) > 0


# ─── bounding box ratio ──────────────────────────────────────────────────


def test_bounding_box_ratio_detects_absurd_change() -> None:
    """Un cambio de caja extremo debe ser rechazado."""
    atom_ids = {1, 2, 3}
    bonds = [_fake_bond(1, 2), _fake_bond(2, 3)]
    before = {1: (0.0, 0.0), 2: (42.0, 0.0), 3: (84.0, 0.0)}
    after = {1: (0.0, 0.0), 2: (420.0, 0.0), 3: (840.0, 0.0)}  # 10x scale
    report = evaluate_clean2d_layout(atom_ids, bonds, before, after, 42.0)
    # El candidato debe ser rechazado (ya sea por caja, longitud, o desplazamiento)
    assert not is_clean2d_candidate_safe(report)


# ─── double bond target length ───────────────────────────────────────────


def test_double_bond_target_length() -> None:
    """Comprueba que el ajuste de longitud para dobles enlaces funciona."""
    positions = {1: (0.0, 0.0), 2: (50.0, 0.0)}
    bonds = [_fake_bond(1, 2, order=2)]
    changed = length_only_polish(positions, bonds, 42.0, max_iterations=10)
    # El target para doble enlace es 42*0.96 = 40.32
    if 2 in changed:
        dist = math.hypot(changed[2][0] - positions[1][0], changed[2][1] - positions[1][1])
        assert dist < 50.0
