"""Pruebas del módulo de seguridad y validación para Clean2D."""

from __future__ import annotations

import math

import pytest


from chemuson.clean2d.safety import (
    bond_length_stats,
    count_new_bond_crossings,
    evaluate_clean2d_layout,
    has_cycles,
    is_clean2d_candidate_safe,
    max_atom_displacement,
    min_nonbonded_distance,
    ring_degeneracy_score,
)
from chemuson.clean2d.length_only import (
    length_only_polish,
    structure_preserving_geometry_polish,
    structure_preserving_length_polish,
)
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


def test_structure_preserving_length_polish_moves_only_distorted_branch() -> None:
    """Normaliza una rama larga sin arrastrar el anillo ni la cadena principal."""
    positions = {
        1: (0.0, 0.0),
        2: (36.0, -21.0),
        3: (72.0, 0.0),
        4: (72.0, 42.0),
        5: (36.0, 63.0),
        6: (0.0, 42.0),
        7: (114.0, 0.0),
        8: (156.0, 0.0),
        9: (198.0, 0.0),
        10: (198.0, 190.0),
    }
    bonds = [
        _fake_bond(1, 2), _fake_bond(2, 3), _fake_bond(3, 4),
        _fake_bond(4, 5), _fake_bond(5, 6), _fake_bond(6, 1),
        _fake_bond(3, 7), _fake_bond(7, 8), _fake_bond(8, 9),
        _fake_bond(9, 10),
    ]

    changed = structure_preserving_length_polish(
        positions,
        bonds,
        42.0,
        max_iterations=4,
    )
    after = {aid: changed.get(aid, positions[aid]) for aid in positions}

    assert math.dist(after[9], after[10]) == pytest.approx(42.0, rel=0.02)
    for atom_id in range(1, 10):
        assert after[atom_id] == pytest.approx(positions[atom_id], abs=1e-6)


def test_structure_preserving_geometry_polish_rebuilds_cyclic_linker_with_double_bond() -> None:
    """Corrige una cadena con doble enlace sin deformar los anillos terminales."""
    positions = {
        1: (0.0, 0.0),
        2: (36.0, -21.0),
        3: (72.0, 0.0),
        4: (72.0, 42.0),
        5: (36.0, 63.0),
        6: (0.0, 42.0),
        7: (104.0, 6.0),
        8: (122.0, -42.0),
        9: (140.0, 10.0),
        10: (168.0, -32.0),
        11: (193.0, 8.0),
        12: (214.0, -6.0),
        13: (244.0, 0.0),
        14: (284.0, 0.0),
        15: (264.0, 35.0),
    }
    bonds = [
        _fake_bond(1, 2), _fake_bond(2, 3), _fake_bond(3, 4),
        _fake_bond(4, 5), _fake_bond(5, 6), _fake_bond(6, 1),
        _fake_bond(3, 7), _fake_bond(7, 8),
        _fake_bond(8, 9, order=2), _fake_bond(9, 10),
        _fake_bond(10, 11), _fake_bond(11, 12), _fake_bond(12, 13),
        _fake_bond(13, 14), _fake_bond(14, 15), _fake_bond(15, 13),
    ]

    changed = structure_preserving_geometry_polish(
        positions,
        bonds,
        42.0,
    )
    after = {aid: changed.get(aid, positions[aid]) for aid in positions}

    expected_lengths = {
        (3, 7): 42.0,
        (7, 8): 42.0,
        (8, 9): 42.0 * 0.96,
        (9, 10): 42.0,
        (10, 11): 42.0,
        (11, 12): 42.0,
        (12, 13): 42.0,
    }
    for (a_id, b_id), expected in expected_lengths.items():
        assert math.dist(after[a_id], after[b_id]) == pytest.approx(expected, abs=1.0)

    assert _angle_between(after[8], after[7], after[9]) == pytest.approx(120.0, abs=8.0)
    assert _angle_between(after[9], after[8], after[10]) == pytest.approx(120.0, abs=8.0)
    for center_id, left_id, right_id in (
        (7, 3, 8), (10, 9, 11), (11, 10, 12), (12, 11, 13),
    ):
        angle = _angle_between(after[center_id], after[left_id], after[right_id])
        assert 95.0 <= angle <= 135.0

    for atom_id in range(1, 7):
        assert after[atom_id] == pytest.approx(positions[atom_id], abs=1e-6)
    for a_id, b_id in ((13, 14), (14, 15), (15, 13)):
        assert math.dist(after[a_id], after[b_id]) == pytest.approx(
            math.dist(positions[a_id], positions[b_id]),
            abs=1e-6,
        )


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


# ─── Clean2DAttempt flow control ─────────────────────────────────────────


def test_clean2d_attempt_unavailable_does_not_stop_flow() -> None:
    """Unavailable no debe interpretarse como 'aplicado'."""
    from chemuson.gui.controllers.clean2d_controller import Clean2DAttempt

    a = Clean2DAttempt(unavailable=True, message="error")
    assert not a.applied
    assert not a.rejected
    assert a.unavailable
    assert a.result_state == "failed-controlled"


def test_clean2d_attempt_rejected_does_not_stop_flow() -> None:
    """Rejected no debe interpretarse como 'aplicado'."""
    from chemuson.gui.controllers.clean2d_controller import Clean2DAttempt

    a = Clean2DAttempt(rejected=True, message="malo")
    assert not a.applied
    assert a.rejected
    assert not a.unavailable
    assert a.result_state == "rejected"


# ─── _apply_clean2d_candidate helper tests ────────────────────────────────


def test_apply_candidate_already_clean_returns_no_change() -> None:
    """Si candidato = before, no debe aplicar movimiento."""
    from unittest.mock import MagicMock
    from chemuson.gui.controllers.clean2d_controller import Clean2DController

    ctrl = Clean2DController()
    window = MagicMock()
    window.canvas.model.atoms = {}
    window.canvas.state.bond_length = 42.0
    window.canvas.undo_stack = MagicMock()
    window.canvas._update_selection_overlay = MagicMock()

    target_ids = {1, 2}
    bonds = [_fake_bond(1, 2)]
    before = {1: (0.0, 0.0), 2: (42.0, 0.0)}
    candidate = {1: (0.0, 0.0), 2: (42.0, 0.0)}

    attempt = ctrl._apply_clean2d_candidate(
        window, target_ids, bonds, before, candidate,
        42.0, "quick", "test_label", False,
    )
    assert not attempt.applied
    assert not attempt.rejected
    assert attempt.result_state == "no-op"
    assert "ya estaba limpia" in attempt.message
    window.canvas.undo_stack.push.assert_not_called()


def test_apply_candidate_rejected_on_bad_geometry() -> None:
    """Un candidato que deforma debe ser rechazado."""
    from unittest.mock import MagicMock
    from chemuson.gui.controllers.clean2d_controller import Clean2DController

    ctrl = Clean2DController()
    window = MagicMock()
    window.canvas.model.atoms = {}
    window.canvas.state.bond_length = 42.0
    window.canvas.undo_stack = MagicMock()
    window.canvas._update_selection_overlay = MagicMock()

    target_ids = {1, 2}
    bonds = [_fake_bond(1, 2)]
    before = {1: (0.0, 0.0), 2: (42.0, 0.0)}
    candidate = {1: (0.0, 0.0), 2: (999.0, 0.0)}  # absurdo

    attempt = ctrl._apply_clean2d_candidate(
        window, target_ids, bonds, before, candidate,
        42.0, "quick", "test_label", False,
    )
    assert not attempt.applied
    assert attempt.rejected
    assert attempt.result_state == "rejected"
    assert attempt.stable_reason == "excessive-displacement"
    window.canvas.undo_stack.push.assert_not_called()


def test_apply_candidate_without_motion_does_not_apply_empty_command() -> None:
    """Un candidato idéntico pero con longitudes malas debe dejar seguir el fallback."""
    from unittest.mock import MagicMock
    from chemuson.gui.controllers.clean2d_controller import Clean2DController

    ctrl = Clean2DController()
    window = MagicMock()
    window.canvas.model.atoms = {}
    window.canvas.state.bond_length = 42.0
    window.canvas.undo_stack = MagicMock()
    window.canvas._update_selection_overlay = MagicMock()

    target_ids = {1, 2}
    bonds = [_fake_bond(1, 2)]
    before = {1: (0.0, 0.0), 2: (64.0, 0.0)}
    candidate = dict(before)

    attempt = ctrl._apply_clean2d_candidate(
        window, target_ids, bonds, before, candidate,
        42.0, "quick", "test_label", False,
    )

    assert not attempt.applied
    assert attempt.unavailable
    assert attempt.result_state == "failed-controlled"
    window.canvas.undo_stack.push.assert_not_called()


def test_rotatable_reflection_pose_proposes_alternative_for_clean_chain() -> None:
    """Una cadena ya limpia puede recibir un conformero 2D alternativo seguro."""
    from unittest.mock import MagicMock
    from chemuson.gui.controllers.clean2d_controller import Clean2DController

    ctrl = Clean2DController()
    window = MagicMock()
    window.canvas.model.atoms = {1: object(), 2: object(), 3: object(), 4: object()}
    window.canvas.undo_stack = MagicMock()
    window.canvas._update_selection_overlay = MagicMock()

    target_ids = {1, 2, 3, 4}
    bonds = [_fake_bond(1, 2), _fake_bond(2, 3), _fake_bond(3, 4)]
    before = {
        1: (0.0, 0.0),
        2: (42.0, 0.0),
        3: (56.0, 39.6),
        4: (98.0, 39.6),
    }

    attempt = ctrl._try_rotatable_reflection_pose(
        window,
        target_ids,
        bonds,
        before,
        42.0,
        "quick",
        False,
    )

    assert attempt.applied
    assert attempt.result_state == "applied"
    window.canvas.undo_stack.push.assert_called_once()
    command = window.canvas.undo_stack.push.call_args.args[0]
    moved = getattr(command, "_after", {})
    assert moved
    assert max(
        math.hypot(moved[aid][0] - before[aid][0], moved[aid][1] - before[aid][1])
        for aid in moved
    ) > 1.0


def test_apply_candidate_already_clean_uses_rotatable_conformer_fallback() -> None:
    """Si RDKit no aporta cambio, clean2d debe proponer una pose alternativa."""
    from unittest.mock import MagicMock, patch
    from chemuson.gui.controllers.clean2d_controller import Clean2DController

    ctrl = Clean2DController()
    window = MagicMock()
    window.canvas.model.atoms = {1: object(), 2: object(), 3: object(), 4: object()}
    window.canvas.state.bond_length = 42.0
    window.canvas.graph = MagicMock()
    window.canvas.undo_stack = MagicMock()
    window.canvas._update_selection_overlay = MagicMock()

    target_ids = {1, 2, 3, 4}
    bonds = [_fake_bond(1, 2), _fake_bond(2, 3), _fake_bond(3, 4)]
    before = {
        1: (0.0, 0.0),
        2: (42.0, 0.0),
        3: (56.0, 39.6),
        4: (98.0, 39.6),
    }

    with patch(
        "chemuson.gui.controllers.clean2d_controller.conformer_3d_for_graph",
        side_effect=RuntimeError("rdkit unavailable"),
    ):
        attempt = ctrl._apply_clean2d_candidate(
            window, target_ids, bonds, before, dict(before),
            42.0, "quick", "test_label", False,
        )

    assert attempt.applied
    window.canvas.undo_stack.push.assert_called_once()


def test_alternative_pose_for_cyclic_structure_skips_3d_projection() -> None:
    """En ciclos, la alternativa debe preservar geometría con reflexión local."""
    from unittest.mock import MagicMock, patch
    from chemuson.gui.controllers.clean2d_controller import Clean2DController

    ctrl = Clean2DController()
    window = MagicMock()
    window.canvas.model.atoms = {
        aid: object()
        for aid in range(1, 7)
    }
    window.canvas.undo_stack = MagicMock()
    window.canvas._update_selection_overlay = MagicMock()

    target_ids = {1, 2, 3, 4, 5, 6}
    bonds = [
        _fake_bond(1, 2), _fake_bond(2, 3), _fake_bond(3, 1),
        _fake_bond(3, 4), _fake_bond(4, 5), _fake_bond(5, 6),
    ]
    before = {
        1: (0.0, 0.0),
        2: (40.0, 0.0),
        3: (20.0, 34.6),
        4: (60.0, 34.6),
        5: (80.0, 69.2),
        6: (120.0, 69.2),
    }

    with patch(
        "chemuson.gui.controllers.clean2d_controller.conformer_3d_for_graph",
    ) as mock_conformer:
        attempt = ctrl._try_alternative_conformer_pose(
            window,
            target_ids,
            bonds,
            before,
            40.0,
            "quick",
            True,
        )

    assert attempt.applied
    mock_conformer.assert_not_called()


def test_apply_candidate_conformer_mode_can_propose_cyclic_alternative() -> None:
    """Una estructura cíclica ya limpia aún puede proponer conformero local."""
    from unittest.mock import MagicMock, patch
    from chemuson.gui.controllers.clean2d_controller import Clean2DController

    ctrl = Clean2DController()
    window = MagicMock()
    window.canvas.model.atoms = {aid: object() for aid in range(1, 7)}
    window.canvas.state.bond_length = 40.0
    window.canvas.undo_stack = MagicMock()
    window.canvas._update_selection_overlay = MagicMock()

    target_ids = {1, 2, 3, 4, 5, 6}
    bonds = [
        _fake_bond(1, 2), _fake_bond(2, 3), _fake_bond(3, 1),
        _fake_bond(3, 4), _fake_bond(4, 5), _fake_bond(5, 6),
    ]
    before = {
        1: (0.0, 0.0),
        2: (40.0, 0.0),
        3: (20.0, 34.6),
        4: (60.0, 34.6),
        5: (80.0, 69.2),
        6: (120.0, 69.2),
    }

    with patch(
        "chemuson.gui.controllers.clean2d_controller.conformer_3d_for_graph",
    ) as mock_conformer:
        attempt = ctrl._apply_clean2d_candidate(
            window,
            target_ids,
            bonds,
            before,
            dict(before),
            40.0,
            "conformer",
            "test_label",
            True,
        )

    assert attempt.applied
    mock_conformer.assert_not_called()
    window.canvas.undo_stack.push.assert_called_once()


def test_cyclic_rotatable_reflection_keeps_multiring_structure_recognizable() -> None:
    """No debe aceptar una reflexión que mueva media molécula cíclica."""
    from unittest.mock import MagicMock
    from chemuson.gui.controllers.clean2d_controller import Clean2DController

    ctrl = Clean2DController()
    window = MagicMock()
    coords = {
        1: (0.0, 0.0), 2: (36.0, -21.0), 3: (72.0, 0.0),
        4: (72.0, 42.0), 5: (36.0, 63.0), 6: (0.0, 42.0),
        7: (112.0, 18.0), 8: (152.0, -5.0), 9: (192.0, 20.0),
        10: (232.0, -2.0), 11: (272.0, 20.0), 12: (312.0, 0.0),
        13: (352.0, 23.0), 14: (392.0, 0.0), 15: (432.0, 23.0),
        16: (472.0, 5.0), 17: (510.0, 25.0), 18: (510.0, 65.0),
        19: (550.0, 45.0), 20: (392.0, -42.0), 21: (312.0, -42.0),
        22: (232.0, -42.0),
    }
    edges = [
        (1, 2), (2, 3), (3, 4), (4, 5), (5, 6), (6, 1),
        (3, 7), (7, 8), (8, 9), (9, 10), (10, 11), (11, 12),
        (12, 13), (13, 14), (14, 15), (15, 16),
        (16, 17), (17, 18), (18, 19), (19, 17),
        (14, 20), (12, 21), (10, 22),
    ]
    bonds = [_fake_bond(a_id, b_id) for a_id, b_id in edges]
    window.canvas.model.atoms = {atom_id: object() for atom_id in coords}
    window.canvas.undo_stack = MagicMock()
    window.canvas._update_selection_overlay = MagicMock()

    attempt = ctrl._try_rotatable_reflection_pose(
        window,
        set(coords),
        bonds,
        coords,
        40.0,
        "quick",
        True,
    )

    assert attempt.applied
    assert attempt.report is not None
    assert attempt.report.max_displacement <= 100.0
    assert attempt.report.mean_displacement <= 40.0


# ─── Bridge test: run_clean_2d never silent ────────────────────────────────


def test_run_clean_2d_never_silent_on_no_target() -> None:
    """Sin target debe mostrar 'Nada que limpiar'."""
    from unittest.mock import MagicMock
    from chemuson.gui.controllers.clean2d_controller import Clean2DController

    ctrl = Clean2DController()
    window = MagicMock()
    window.canvas._selected_structure_ids.return_value = (set(), [])
    window.canvas.model.atoms = {}
    window.canvas.model.bonds = {}

    ctrl.run_clean_2d(window, 1.0, 200, "(test)")
    status_calls = [args[0][0] for args in window.statusBar().showMessage.call_args_list]
    assert any("Nada que limpiar" in str(msg) for msg in status_calls)


def test_run_clean_2d_shows_message_when_isolated_unavailable() -> None:
    """El controlador aplica el candidato elegido por el motor unificado."""
    from unittest.mock import MagicMock, patch
    from chemuson.clean2d import Clean2DCandidate, Clean2DMode, Clean2DResult
    from chemuson.core.model import MolGraph
    from chemuson.gui.controllers.clean2d_controller import Clean2DController

    ctrl = Clean2DController()
    graph = MolGraph()
    a1 = graph.add_atom("C", 0.0, 0.0, atom_id=1)
    a2 = graph.add_atom("C", 120.0, 0.0, atom_id=2)
    bond = graph.add_bond(a1.id, a2.id, order=1)
    window = MagicMock()
    window.canvas._selected_structure_ids.return_value = ({1, 2}, [bond])
    window.canvas.model = graph
    window.canvas.state.bond_length = 42.0
    window.canvas.undo_stack = MagicMock()
    window.canvas._update_selection_overlay = MagicMock()
    window.canvas._build_selection_graph = MagicMock(return_value=graph)
    window.canvas.clean_2d_fallback = MagicMock()

    candidate = Clean2DCandidate(
        source="internal_templates",
        coords={1: (0.0, 0.0), 2: (42.0, 0.0)},
        message="Estructura 2D limpiada (motor interno)",
    )
    result = Clean2DResult(
        mode=Clean2DMode.QUICK,
        atom_ids={1, 2},
        before_coords={1: (0.0, 0.0), 2: (120.0, 0.0)},
        candidates=(candidate,),
        selected=candidate,
    )

    with patch("chemuson.gui.controllers.clean2d_controller.run_clean2d_engine", return_value=result) as mock_engine:
        ctrl.run_clean_2d(window, 1.0, 200, "(test)")

    mock_engine.assert_called_once()
    window.canvas.undo_stack.push.assert_called_once()
    window.canvas.clean_2d_fallback.assert_not_called()
    status_calls = window.statusBar().showMessage.call_args_list
    assert any("motor interno" in str(args[0][0]) for args in status_calls)


def test_run_clean_2d_cyclic_current_tries_smart_propose_instead_of_length_early_return() -> None:
    """Un current limpio pasa por Smart Clean2D; no hay early return length-only."""
    from unittest.mock import MagicMock, patch
    from chemuson.clean2d import Clean2DCandidate, Clean2DMode, Clean2DResult
    from chemuson.core.model import MolGraph
    from chemuson.gui.controllers.clean2d_controller import Clean2DController

    ctrl = Clean2DController()
    graph = MolGraph()
    for atom_id, coord in enumerate([(0.0, 0.0), (42.0, 0.0), (21.0, 36.0)], 1):
        graph.add_atom("C", coord[0], coord[1], atom_id=atom_id)
    cyclic_bonds = [
        graph.add_bond(1, 2, order=1),
        graph.add_bond(2, 3, order=1),
        graph.add_bond(3, 1, order=1),
    ]
    window = MagicMock()
    window.canvas._selected_structure_ids.return_value = ({1, 2, 3}, cyclic_bonds)
    window.canvas.model = graph
    window.canvas.state.bond_length = 42.0
    window.canvas.undo_stack = MagicMock()
    window.canvas._update_selection_overlay = MagicMock()
    window.canvas._build_selection_graph = MagicMock(return_value=graph)
    window.canvas.clean_2d_fallback = MagicMock()

    candidate = Clean2DCandidate(
        source="current",
        coords={1: (0.0, 0.0), 2: (42.0, 0.0), 3: (21.0, 36.0)},
        message="Estructura 2D ya estaba limpia",
        metadata={"current_canonical_enough": True},
    )
    result = Clean2DResult(
        mode=Clean2DMode.QUICK,
        atom_ids={1, 2, 3},
        before_coords=dict(candidate.coords),
        candidates=(candidate,),
        selected=candidate,
    )

    with (
        patch("chemuson.gui.controllers.clean2d_controller.run_clean2d_engine", return_value=result) as mock_engine,
        patch.object(ctrl, "_safe_length_polish_only") as mock_length,
    ):
        ctrl.run_clean_2d(window, 1.0, 200, "(test)")

    assert mock_engine.call_count == 2
    assert mock_engine.call_args_list[0].kwargs["mode"] == "quick"
    assert mock_engine.call_args_list[0].kwargs["avoid_hashes"] == set()
    assert mock_engine.call_args_list[1].kwargs["mode"] == "propose"
    assert mock_engine.call_args_list[1].kwargs["avoid_hashes"]
    mock_length.assert_not_called()
    window.canvas.clean_2d_fallback.assert_not_called()


def test_run_clean_2d_cyclic_engine_candidate_can_move_atoms() -> None:
    """Un candidato cíclico validado se aplica con un solo MoveAtomsCommand."""
    from unittest.mock import MagicMock, patch
    from chemuson.clean2d import Clean2DCandidate, Clean2DMode, Clean2DResult
    from chemuson.core.model import MolGraph
    from chemuson.gui.controllers.clean2d_controller import Clean2DController

    ctrl = Clean2DController()
    graph = MolGraph()
    for atom_id, coord in enumerate([(0.0, 0.0), (4.0, 0.0), (2.0, 1.0)], 1):
        graph.add_atom("C", coord[0], coord[1], atom_id=atom_id)
    cyclic_bonds = [
        graph.add_bond(1, 2, order=1),
        graph.add_bond(2, 3, order=1),
        graph.add_bond(3, 1, order=1),
    ]
    window = MagicMock()
    window.canvas._selected_structure_ids.return_value = ({1, 2, 3}, cyclic_bonds)
    window.canvas.model = graph
    window.canvas.state.bond_length = 42.0
    window.canvas.undo_stack = MagicMock()
    window.canvas._update_selection_overlay = MagicMock()
    window.canvas._build_selection_graph = MagicMock(return_value=graph)
    window.canvas.clean_2d_fallback = MagicMock()

    candidate = Clean2DCandidate(
        source="internal_templates",
        coords={1: (0.0, 24.0), 2: (36.0, -12.0), 3: (-36.0, -12.0)},
        message="Estructura 2D limpiada (motor interno)",
    )
    result = Clean2DResult(
        mode=Clean2DMode.PUBLICATION,
        atom_ids={1, 2, 3},
        before_coords={1: (0.0, 0.0), 2: (4.0, 0.0), 3: (2.0, 1.0)},
        candidates=(candidate,),
        selected=candidate,
    )

    with patch("chemuson.gui.controllers.clean2d_controller.run_clean2d_engine", return_value=result):
        ctrl.run_clean_2d(window, 1.0, 200, "(test)")

    window.canvas.undo_stack.push.assert_called_once()
    window.canvas.clean_2d_fallback.assert_not_called()


def test_run_clean_2d_quick_does_not_pass_history_as_avoid_hashes() -> None:
    """Ctrl+K/quick debe limpiar, no usar la memoria anti-repetición."""
    from unittest.mock import MagicMock, patch
    from chemuson.clean2d import Clean2DCandidate, Clean2DMode, Clean2DResult
    from chemuson.core.model import MolGraph
    from chemuson.gui.controllers.clean2d_controller import Clean2DController

    ctrl = Clean2DController()
    graph = MolGraph()
    a1 = graph.add_atom("C", 0.0, 0.0, atom_id=1)
    a2 = graph.add_atom("C", 120.0, 0.0, atom_id=2)
    bond = graph.add_bond(a1.id, a2.id, order=1)
    ctrl._recent_geometry_hashes[ctrl._history_key(graph, {1, 2})] = ["previous-clean-hash"]

    window = MagicMock()
    window.canvas._selected_structure_ids.return_value = ({1, 2}, [bond])
    window.canvas.model = graph
    window.canvas.state.bond_length = 42.0
    window.canvas.undo_stack = MagicMock()
    window.canvas._update_selection_overlay = MagicMock()
    window.canvas._build_selection_graph = MagicMock(return_value=graph)

    candidate = Clean2DCandidate(
        source="internal_templates",
        coords={1: (0.0, 0.0), 2: (42.0, 0.0)},
        message="Estructura 2D limpiada (motor interno)",
    )
    result = Clean2DResult(
        mode=Clean2DMode.QUICK,
        atom_ids={1, 2},
        before_coords={1: (0.0, 0.0), 2: (120.0, 0.0)},
        candidates=(candidate,),
        selected=candidate,
    )

    with patch("chemuson.gui.controllers.clean2d_controller.run_clean2d_engine", return_value=result) as mock_engine:
        ctrl.run_clean_2d(window, 1.0, 200, "(test)", mode="quick")

    assert mock_engine.call_args.kwargs["mode"] == "quick"
    assert mock_engine.call_args.kwargs["avoid_hashes"] == set()


def test_run_clean_2d_cyclic_rejection_safe_polish_moves_distorted_branch() -> None:
    """El ajuste seguro debe corregir una rama larga sin re-layout global."""
    from unittest.mock import MagicMock, patch
    from PyQt6.QtWidgets import QApplication
    from chemuson.gui.canvas import ChemusonCanvas
    from chemuson.gui.controllers.clean2d_controller import Clean2DController

    QApplication.instance() or QApplication([])
    canvas = ChemusonCanvas()
    coords = {
        1: (0.0, 0.0),
        2: (36.0, -21.0),
        3: (72.0, 0.0),
        4: (72.0, 42.0),
        5: (36.0, 63.0),
        6: (0.0, 42.0),
        7: (114.0, 0.0),
        8: (156.0, 0.0),
        9: (198.0, 0.0),
        10: (240.0, 0.0),
        11: (282.0, 0.0),
        12: (324.0, 0.0),
        13: (198.0, 190.0),
    }
    for atom_id, (x, y) in coords.items():
        canvas.model.add_atom("C", x, y, atom_id=atom_id)
    for a_id, b_id in (
        (1, 2), (2, 3), (3, 4), (4, 5), (5, 6), (6, 1),
        (3, 7), (7, 8), (8, 9), (9, 10), (10, 11), (11, 12),
        (9, 13),
    ):
        canvas.model.add_bond(a_id, b_id, order=1)
    canvas._sync_scene_with_model()
    canvas.state.selected_atoms = set(coords)
    canvas.clean_2d_fallback = MagicMock(wraps=canvas.clean_2d_fallback)

    window = MagicMock()
    window.canvas = canvas
    ctrl = Clean2DController()

    before_a = canvas.model.get_atom(9)
    before_b = canvas.model.get_atom(13)
    before_len = math.hypot(before_b.x - before_a.x, before_b.y - before_a.y)

    with (
        patch.object(ctrl, "_try_isolated_rdkit2d") as mock_isolated,
        patch.object(ctrl, "_try_direct_rdkit") as mock_direct,
        patch.object(ctrl, "_polish_with_v2") as mock_polish,
    ):
        mock_isolated.return_value = _make_attempt(rejected=True, message="rechazado")
        mock_direct.return_value = _make_attempt(unavailable=True, message="sin direct")
        mock_polish.return_value = _make_attempt(unavailable=True, message="sin polish")

        ctrl.run_clean_2d(window, 1.0, 200, "(test)")

    after_a = canvas.model.get_atom(9)
    after_b = canvas.model.get_atom(13)
    after_len = math.hypot(after_b.x - after_a.x, after_b.y - after_a.y)

    assert before_len > float(canvas.state.bond_length) * 3.0
    assert after_len == pytest.approx(float(canvas.state.bond_length), rel=0.15)
    mock_isolated.assert_not_called()
    mock_direct.assert_not_called()
    mock_polish.assert_not_called()
    canvas.clean_2d_fallback.assert_not_called()


def test_run_clean_2d_cyclic_geometry_polish_repairs_double_bond_linker() -> None:
    """El flujo real debe corregir cadena deformada con C=C sin RDKit global."""
    from unittest.mock import MagicMock, patch
    from PyQt6.QtWidgets import QApplication
    from chemuson.gui.canvas import ChemusonCanvas
    from chemuson.gui.controllers.clean2d_controller import Clean2DController

    QApplication.instance() or QApplication([])
    canvas = ChemusonCanvas()
    coords = {
        1: (0.0, 0.0),
        2: (36.0, -21.0),
        3: (72.0, 0.0),
        4: (72.0, 42.0),
        5: (36.0, 63.0),
        6: (0.0, 42.0),
        7: (104.0, 6.0),
        8: (122.0, -42.0),
        9: (140.0, 10.0),
        10: (168.0, -32.0),
        11: (193.0, 8.0),
        12: (214.0, -6.0),
        13: (244.0, 0.0),
        14: (284.0, 0.0),
        15: (264.0, 35.0),
    }
    for atom_id, (x, y) in coords.items():
        canvas.model.add_atom("C", x, y, atom_id=atom_id)
    for a_id, b_id, order in (
        (1, 2, 1), (2, 3, 1), (3, 4, 1),
        (4, 5, 1), (5, 6, 1), (6, 1, 1),
        (3, 7, 1), (7, 8, 1), (8, 9, 2),
        (9, 10, 1), (10, 11, 1), (11, 12, 1), (12, 13, 1),
        (13, 14, 1), (14, 15, 1), (15, 13, 1),
    ):
        canvas.model.add_bond(a_id, b_id, order=order)
    canvas._sync_scene_with_model()
    canvas.state.selected_atoms = set(coords)
    canvas.clean_2d_fallback = MagicMock(wraps=canvas.clean_2d_fallback)

    window = MagicMock()
    window.canvas = canvas
    ctrl = Clean2DController()

    with (
        patch.object(ctrl, "_try_isolated_rdkit2d") as mock_isolated,
        patch.object(ctrl, "_try_direct_rdkit") as mock_direct,
    ):
        ctrl.run_clean_2d(window, 1.0, 200, "(test)")

    after = {
        atom_id: (canvas.model.get_atom(atom_id).x, canvas.model.get_atom(atom_id).y)
        for atom_id in coords
    }

    assert math.dist(after[8], after[9]) == pytest.approx(
        float(canvas.state.bond_length) * 0.96,
        abs=1.0,
    )
    assert _angle_between(after[8], after[7], after[9]) == pytest.approx(120.0, abs=8.0)
    assert _angle_between(after[9], after[8], after[10]) == pytest.approx(120.0, abs=8.0)
    mock_isolated.assert_not_called()
    mock_direct.assert_not_called()
    canvas.clean_2d_fallback.assert_not_called()


def _make_attempt(*, applied=False, rejected=False, unavailable=False, message="", report=None):
    from chemuson.gui.controllers.clean2d_controller import Clean2DAttempt
    return Clean2DAttempt(
        applied=applied, rejected=rejected, unavailable=unavailable,
        message=message, report=report,
    )


def _angle_between(
    center: tuple[float, float],
    left: tuple[float, float],
    right: tuple[float, float],
) -> float:
    a1 = math.atan2(left[1] - center[1], left[0] - center[0])
    a2 = math.atan2(right[1] - center[1], right[0] - center[0])
    return abs((math.degrees(a2 - a1) + 180.0) % 360.0 - 180.0)
