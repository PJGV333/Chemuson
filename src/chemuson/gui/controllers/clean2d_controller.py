from __future__ import annotations

from dataclasses import dataclass, field
import math
from typing import Any

from chemuson.clean2d import (
    Clean2DCandidate,
    Clean2DParameters,
    Clean2DQualityReport,
    Clean2DResult,
    assert_clean2d_invariants,
    clean2d_geometry_hash,
    evaluate_clean2d_layout,
    has_cycles,
    is_clean2d_candidate_safe,
    length_only_polish,
    optimize_clean2d_positions,
    run_clean2d_engine,
    structure_preserving_geometry_polish,
)
from chemuson.core.model import bond_affects_valence
from chemuson.geometry3d import Rotation3D, conformer_3d_for_graph, project_conformer_to_2d


@dataclass(frozen=True)
class Clean2DAttempt:
    """Resultado explícito de un intento de limpieza 2D."""
    applied: bool = False
    rejected: bool = False
    unavailable: bool = False
    message: str = ""
    report: Clean2DQualityReport | None = None


class Clean2DController:
    """Pipeline de clean-2D desacoplado de ChemusonWindow."""

    def __init__(self) -> None:
        self._recent_geometry_hashes: dict[str, list[str]] = {}

    @staticmethod
    def _is_acyclic_structure(atom_ids: set[int], bonds) -> bool:
        if not atom_ids:
            return True
        adjacency: dict[int, list[int]] = {}
        for bond in bonds:
            if bond.a1_id not in atom_ids or bond.a2_id not in atom_ids:
                continue
            adjacency.setdefault(bond.a1_id, []).append(bond.a2_id)
            adjacency.setdefault(bond.a2_id, []).append(bond.a1_id)

        visited: set[int] = set()
        for start in atom_ids:
            if start in visited:
                continue
            stack: list[tuple[int, int]] = [(start, -1)]
            while stack:
                node, parent = stack.pop()
                if node in visited:
                    continue
                visited.add(node)
                for neigh in adjacency.get(node, []):
                    if neigh == parent:
                        continue
                    if neigh in visited:
                        return False
                    stack.append((neigh, node))
        return True

    def _apply_clean2d_candidate(
        self,
        window: Any,
        target_ids: set[int],
        bonds: list,
        before: dict[int, tuple[float, float]],
        candidate: dict[int, tuple[float, float]],
        target_len: float,
        mode: str,
        label: str,
        cyclic: bool,
    ) -> Clean2DAttempt:
        """Helper único que completa, alínea, evalúa safety y aplica un candidato.

        Returns:
            Clean2DAttempt con applied/rejected según corresponda.
        """
        canvas = window.canvas

        after: dict[int, tuple[float, float]] = {}
        for aid in target_ids:
            if aid in candidate:
                after[aid] = candidate[aid]
            elif aid in before:
                after[aid] = before[aid]

        if not after:
            return Clean2DAttempt(unavailable=True, message="sin_coordenadas_candidato")

        before_cx, before_cy = _coords_center(before)
        after_cx, after_cy = _coords_center(after)
        shift_x = before_cx - after_cx
        shift_y = before_cy - after_cy
        if abs(shift_x) > 1e-6 or abs(shift_y) > 1e-6:
            after = {
                aid: (x + shift_x, y + shift_y)
                for aid, (x, y) in after.items()
            }

        changed: dict[int, tuple[float, float]] = {}
        for aid in target_ids:
            if aid not in before or aid not in after:
                continue
            d = _position_delta(before[aid], after[aid])
            if d > 0.5:
                changed[aid] = after[aid]

        if not changed:
            if mode in {"quick", "publication", "conformer"}:
                is_bad = False
                for b in bonds:
                    try:
                        x1, y1 = before[b.a1_id]
                        x2, y2 = before[b.a2_id]
                        curr_len = math.hypot(x2 - x1, y2 - y1)
                        target = float(getattr(window.canvas.state, "bond_length", 42.0))
                        if curr_len < target * 0.8 or curr_len > target * 1.2:
                            is_bad = True
                            break
                    except Exception:
                        continue

                if not is_bad:
                    alternative = self._try_alternative_conformer_pose(
                        window,
                        target_ids,
                        bonds,
                        before,
                        target_len,
                        mode,
                        cyclic,
                    )
                    if alternative.applied:
                        return alternative
                    return Clean2DAttempt(message="Estructura 2D ya estaba limpia")
                return Clean2DAttempt(
                    unavailable=True,
                    message="candidato_clean2d_sin_movimiento",
                )
            return Clean2DAttempt(
                unavailable=True,
                message="candidato_clean2d_sin_movimiento",
            )

        report = evaluate_clean2d_layout(
            target_ids,
            bonds,
            before,
            {aid: changed.get(aid, before[aid]) for aid in target_ids},
            target_len,
            is_cyclic=cyclic,
        )
        if not is_clean2d_candidate_safe(report, mode):
            return Clean2DAttempt(
                rejected=True,
                message="Limpieza 2D omitida: el candidato deformaba la estructura",
                report=report,
            )

        from chemuson.gui.commands import MoveAtomsCommand

        canvas.undo_stack.push(MoveAtomsCommand(canvas.model, canvas, before, changed))
        canvas._update_selection_overlay()
        return Clean2DAttempt(applied=True, message=label)

    def _try_alternative_conformer_pose(
        self,
        window: Any,
        target_ids: set[int],
        bonds: list,
        before: dict[int, tuple[float, float]],
        target_len: float,
        mode: str,
        cyclic: bool,
    ) -> Clean2DAttempt:
        """Propone una geometría 2D alternativa basada en conformero 3D."""
        canvas = window.canvas
        if len(target_ids) < 3 or not bonds:
            return Clean2DAttempt(unavailable=True, message="sin_pose_alternativa")
        if cyclic:
            return self._try_rotatable_reflection_pose(
                window,
                target_ids,
                bonds,
                before,
                target_len,
                mode,
                cyclic,
            )

        conf3d_positions: dict[int, tuple[float, float, float]] | None = None
        try:
            graph = (
                canvas._build_selection_graph(target_ids, bonds)
                if target_ids != set(canvas.model.atoms.keys())
                else canvas.graph
            )
            conf3d = conformer_3d_for_graph(graph, timeout_s=3.5)
            if conf3d.ok:
                conf3d_positions = conf3d.atom_positions
        except Exception:
            conf3d_positions = None

        center = _coords_center(before)
        rotations = (
            Rotation3D(pitch=math.radians(24.0), yaw=math.radians(32.0)),
            Rotation3D(pitch=math.radians(-18.0), yaw=math.radians(57.0)),
            Rotation3D(pitch=math.radians(36.0), yaw=math.radians(-41.0)),
        )

        if conf3d_positions:
            for rotation in rotations:
                projected = project_conformer_to_2d(
                    conf3d_positions,
                    rotation=rotation,
                    center=center,
                    scale=max(1.0, target_len),
                )
                after = {
                    atom.atom_id: (float(atom.x), float(atom.y))
                    for atom in projected
                    if atom.atom_id in target_ids
                }
                if len(after) != len(target_ids):
                    continue

                after = window._rescale_coords_to_bond_length(after, bonds, target_len)
                after = window._align_coords_to_reference(before, after)
                changed = {
                    aid: after[aid]
                    for aid in target_ids
                    if aid in before and _position_delta(before[aid], after[aid]) > 0.9
                }
                if not changed:
                    continue

                report = evaluate_clean2d_layout(
                    target_ids,
                    bonds,
                    before,
                    {aid: changed.get(aid, before[aid]) for aid in target_ids},
                    target_len,
                    is_cyclic=cyclic,
                )
                if not is_clean2d_candidate_safe(report, mode):
                    continue
                if report.mean_displacement < target_len * 0.20:
                    continue

                from chemuson.gui.commands import MoveAtomsCommand

                canvas.undo_stack.push(MoveAtomsCommand(canvas.model, canvas, before, changed))
                canvas._update_selection_overlay()
                label = "Estructura 2D alternativa propuesta"
                if len(target_ids) < len(canvas.model.atoms):
                    label = "Selección 2D alternativa propuesta"
                return Clean2DAttempt(applied=True, message=label, report=report)

        return self._try_rotatable_reflection_pose(
            window,
            target_ids,
            bonds,
            before,
            target_len,
            mode,
            cyclic,
        )

    def _try_rotatable_reflection_pose(
        self,
        window: Any,
        target_ids: set[int],
        bonds: list,
        before: dict[int, tuple[float, float]],
        target_len: float,
        mode: str,
        cyclic: bool,
    ) -> Clean2DAttempt:
        """Genera un conformero 2D reflejando un lado de un enlace rotatable."""
        canvas = window.canvas
        if len(target_ids) < 3 or not bonds:
            return Clean2DAttempt(unavailable=True, message="sin_pose_alternativa")

        for bond, side in _rotatable_reflection_sides(target_ids, bonds):
            candidate = _reflect_side_across_bond(before, bond, side)
            if not candidate:
                continue
            after = dict(before)
            after.update(candidate)
            changed = {
                aid: after[aid]
                for aid in target_ids
                if aid in before and _position_delta(before[aid], after[aid]) > 0.9
            }
            if not changed:
                continue

            report = evaluate_clean2d_layout(
                target_ids,
                bonds,
                before,
                {aid: changed.get(aid, before[aid]) for aid in target_ids},
                target_len,
                is_cyclic=cyclic,
            )
            if cyclic:
                safety_mode = "publication" if mode == "publication" else "quick"
            else:
                safety_mode = "conformer" if mode in {"quick", "publication"} else mode
            if not is_clean2d_candidate_safe(report, safety_mode):
                continue
            if cyclic and report.max_displacement > target_len * 2.5:
                continue
            if cyclic and report.mean_displacement > target_len * 1.0:
                continue
            if report.mean_displacement < target_len * 0.05:
                continue

            from chemuson.gui.commands import MoveAtomsCommand

            canvas.undo_stack.push(MoveAtomsCommand(canvas.model, canvas, before, changed))
            canvas._update_selection_overlay()
            label = "Estructura 2D alternativa propuesta"
            if len(target_ids) < len(canvas.model.atoms):
                label = "Selección 2D alternativa propuesta"
            return Clean2DAttempt(applied=True, message=label, report=report)

        return Clean2DAttempt(unavailable=True, message="sin_pose_alternativa")

    def _polish_with_v2(
        self,
        window: Any,
        atom_ids: set[int],
        *,
        mode: str,
        status_suffix: str,
    ) -> Clean2DAttempt:
        """Aplica pulido clean2d_v2 como un único comando undoable."""
        canvas = window.canvas
        if not atom_ids:
            return Clean2DAttempt(unavailable=True, message="sin_seleccion")
        before = {
            aid: (canvas.model.get_atom(aid).x, canvas.model.get_atom(aid).y)
            for aid in atom_ids
            if aid in canvas.model.atoms
        }
        if not before:
            return Clean2DAttempt(unavailable=True, message="sin_coordenadas")

        bonds = _get_bonds_from_model(canvas.model, atom_ids)
        target_len = float(getattr(canvas.state, "bond_length", 42.0) or 42.0)
        cyclic = has_cycles(atom_ids, bonds)

        if cyclic:
            after = structure_preserving_geometry_polish(
                before, bonds, target_len,
                max_displacement_per_atom=target_len * 8.0,
            )
            label = f"{status_suffix} (ajuste)".strip()
            return self._apply_clean2d_candidate(
                window, atom_ids, bonds, before, after,
                target_len, "conformer", label, cyclic=True,
            )

        params = (
            Clean2DParameters.publication(target_len)
            if mode == "publication"
            else Clean2DParameters.quick(target_len)
        )
        after_v2 = optimize_clean2d_positions(canvas.model, atom_ids, params)
        label = f"{status_suffix} (clean2d_v2)".strip()
        return self._apply_clean2d_candidate(
            window, atom_ids, bonds, before, after_v2,
            target_len, mode, label, cyclic=False,
        )

    def run_clean_2d(
        self,
        window: Any,
        step_ratio: float,
        fallback_iterations: int,
        status_suffix: str,
        mode: str = "quick",
    ) -> None:
        canvas = window.canvas
        atom_ids, bonds = canvas._selected_structure_ids()
        target_ids = atom_ids if atom_ids else set(canvas.model.atoms.keys())
        if not target_ids:
            window.statusBar().showMessage("Nada que limpiar")
            return
        scale_bonds = bonds if atom_ids else list(canvas.model.bonds.values())
        graph = (
            canvas._build_selection_graph(target_ids, scale_bonds)
            if atom_ids
            else canvas.graph
        )
        target_bond_len = float(getattr(canvas.state, "bond_length", 42.0) or 42.0)
        before = {
            aid: (canvas.model.get_atom(aid).x, canvas.model.get_atom(aid).y)
            for aid in target_ids
            if aid in canvas.model.atoms
        }
        alternative_mode = self._is_alternative_mode(mode)
        avoid_hashes = self._avoid_hashes_for_proposal(graph, target_ids, before) if alternative_mode else set()
        result = run_clean2d_engine(
            graph,
            target_ids,
            mode=mode,
            target_bond_length=target_bond_len,
            avoid_hashes=avoid_hashes,
        )
        attempt = self._apply_engine_result(
            window,
            graph,
            target_ids,
            scale_bonds,
            before,
            result,
            target_bond_len,
            require_alternative=alternative_mode,
            remember_hashes=alternative_mode,
        )

        if attempt.applied:
            window.statusBar().showMessage(
                self._selection_status_message(attempt.message, bool(atom_ids))
                or "Estructura 2D limpiada"
            )
            return

        if self._smart_propose_enabled(mode) and self._result_can_try_proposal(result, attempt):
            propose_attempt = self._run_propose_pass(
                window,
                graph,
                target_ids,
                scale_bonds,
                before,
                target_bond_len,
                remember_hashes=True,
            )
            if propose_attempt.applied:
                window.statusBar().showMessage(
                    self._selection_status_message(propose_attempt.message, bool(atom_ids))
                    or "Estructura 2D alternativa propuesta"
                )
                return
            window.statusBar().showMessage(
                self._selection_status_message(self._no_safe_alternative_message(), bool(atom_ids))
            )
            return

        if alternative_mode:
            window.statusBar().showMessage(
                self._selection_status_message(
                    self._no_safe_alternative_message()
                    if attempt.unavailable or attempt.rejected
                    else attempt.message or self._no_safe_alternative_message(),
                    bool(atom_ids),
                )
            )
            return

        message = attempt.message or result.message or "Limpieza 2D omitida: sin candidato seguro"
        window.statusBar().showMessage(self._selection_status_message(message, bool(atom_ids)))

    def _apply_engine_result(
        self,
        window: Any,
        graph: Any,
        target_ids: set[int],
        bonds: list,
        before: dict[int, tuple[float, float]],
        result: Clean2DResult,
        target_bond_len: float,
        *,
        require_alternative: bool = False,
        remember_hashes: bool = False,
    ) -> Clean2DAttempt:
        candidate = result.selected
        if candidate is None:
            return Clean2DAttempt(
                unavailable=True,
                message=result.message or "Limpieza 2D omitida: sin candidato seguro",
            )
        return self._apply_engine_candidate(
            window,
            graph,
            target_ids,
            bonds,
            before,
            candidate,
            target_bond_len,
            require_alternative=require_alternative,
            remember_hashes=remember_hashes,
        )

    def _apply_engine_candidate(
        self,
        window: Any,
        graph: Any,
        target_ids: set[int],
        bonds: list,
        before: dict[int, tuple[float, float]],
        candidate: Clean2DCandidate,
        target_bond_len: float,
        *,
        require_alternative: bool = False,
        remember_hashes: bool = False,
    ) -> Clean2DAttempt:
        if require_alternative and candidate.source == "current":
            return Clean2DAttempt(unavailable=True, message="sin_pose_alternativa")

        after = {
            aid: (
                (float(candidate.coords[aid][0]), float(candidate.coords[aid][1]))
                if aid in candidate.coords
                else before[aid]
            )
            for aid in target_ids
            if aid in before
        }
        changed = {
            aid: coord
            for aid, coord in after.items()
            if _position_delta(before[aid], coord) > 0.5
        }

        if not changed:
            message = (
                "sin_pose_alternativa"
                if require_alternative
                else candidate.message or "Estructura 2D ya estaba limpia"
            )
            return Clean2DAttempt(unavailable=require_alternative, message=message)

        report = candidate.report
        if require_alternative:
            try:
                assert_clean2d_invariants(graph, graph, before, after, atom_ids=target_ids)
            except Exception as exc:
                return Clean2DAttempt(
                    rejected=True,
                    message=f"Limpieza 2D omitida: {exc}",
                    report=report,
                )

            report = evaluate_clean2d_layout(
                target_ids,
                bonds,
                before,
                after,
                target_bond_len,
                is_cyclic=has_cycles(target_ids, bonds),
            )
            if not is_clean2d_candidate_safe(report, "conformer"):
                return Clean2DAttempt(
                    rejected=True,
                    message="Limpieza 2D omitida: el candidato deformaba la estructura",
                    report=report,
                )
            if report.mean_displacement < target_bond_len * 0.05:
                return Clean2DAttempt(
                    unavailable=True,
                    message="sin_pose_alternativa",
                    report=report,
                )

        from chemuson.gui.commands import MoveAtomsCommand

        window.canvas.undo_stack.push(MoveAtomsCommand(window.canvas.model, window.canvas, before, changed))
        window.canvas._update_selection_overlay()
        if remember_hashes:
            self._remember_geometry_hashes(graph, target_ids, before, after)
        return Clean2DAttempt(applied=True, message=candidate.message or "Estructura 2D limpiada", report=report)

    def _run_propose_pass(
        self,
        window: Any,
        graph: Any,
        target_ids: set[int],
        bonds: list,
        before: dict[int, tuple[float, float]],
        target_bond_len: float,
        *,
        remember_hashes: bool,
    ) -> Clean2DAttempt:
        result = run_clean2d_engine(
            graph,
            target_ids,
            mode="propose",
            target_bond_length=target_bond_len,
            avoid_hashes=self._avoid_hashes_for_proposal(graph, target_ids, before),
        )
        return self._apply_engine_result(
            window,
            graph,
            target_ids,
            bonds,
            before,
            result,
            target_bond_len,
            require_alternative=True,
            remember_hashes=remember_hashes,
        )

    @staticmethod
    def _is_alternative_mode(mode: str) -> bool:
        return str(mode or "").strip().lower() in {"propose", "conformer", "alternative"}

    @staticmethod
    def _smart_propose_enabled(mode: str) -> bool:
        return str(mode or "").strip().lower() in {"quick", "publication", "publish"}

    @staticmethod
    def _result_can_try_proposal(result: Clean2DResult, attempt: Clean2DAttempt) -> bool:
        candidate = result.selected
        return (
            not attempt.applied
            and not attempt.rejected
            and candidate is not None
            and (
                candidate.source == "current"
                or attempt.message in {"Estructura 2D ya estaba limpia", "candidato_clean2d_sin_movimiento"}
                or not attempt.message
            )
        )

    def _avoid_hashes_for_proposal(
        self,
        graph: Any,
        atom_ids: set[int],
        before: dict[int, tuple[float, float]],
    ) -> set[str]:
        history_key = self._history_key(graph, atom_ids)
        avoid_hashes = set(self._recent_geometry_hashes.get(history_key, []))
        try:
            avoid_hashes.add(clean2d_geometry_hash(graph, before, atom_ids))
        except Exception:
            pass
        return avoid_hashes

    @staticmethod
    def _no_safe_alternative_message() -> str:
        return "Estructura 2D ya estaba limpia; no se encontró alternativa segura"

    def _history_key(self, graph: Any, atom_ids: set[int]) -> str:
        bonds = _get_bonds_from_model(graph, atom_ids) if hasattr(graph, "bonds") else []
        atom_part = tuple(
            (atom_id, getattr(graph.atoms[atom_id], "element", ""))
            for atom_id in sorted(atom_ids)
            if atom_id in getattr(graph, "atoms", {})
        )
        bond_part = tuple(
            (
                min(int(bond.a1_id), int(bond.a2_id)),
                max(int(bond.a1_id), int(bond.a2_id)),
                int(getattr(bond, "order", 1) or 1),
                bool(getattr(bond, "is_aromatic", False)),
            )
            for bond in sorted(bonds, key=lambda item: item.id)
        )
        return repr((atom_part, bond_part))

    def _remember_geometry_hashes(
        self,
        graph: Any,
        atom_ids: set[int],
        before: dict[int, tuple[float, float]],
        after: dict[int, tuple[float, float]],
    ) -> None:
        key = self._history_key(graph, atom_ids)
        hashes = self._recent_geometry_hashes.setdefault(key, [])
        for coords in (before, after):
            try:
                value = clean2d_geometry_hash(graph, coords, atom_ids)
            except Exception:
                continue
            if value in hashes:
                hashes.remove(value)
            hashes.append(value)
        del hashes[:-8]

    @staticmethod
    def _selection_status_message(message: str, is_selection: bool) -> str:
        if not message:
            return ""
        if not is_selection:
            return message
        return message.replace("Estructura 2D", "Selección 2D", 1)

    def _try_isolated_rdkit2d(
        self,
        window: Any,
        target_ids: set[int],
        scale_bonds: list,
        cyclic: bool,
        mode: str,
    ) -> Clean2DAttempt:
        """Intenta RDKit en worker aislado. Retorna Clean2DAttempt explícito."""
        canvas = window.canvas
        try:
            from chemuson.chemio.rdkit_safe import clean2d_isolated

            graph = canvas._build_selection_graph(
                target_ids, scale_bonds
            ) if target_ids != set(canvas.model.atoms.keys()) else canvas.graph
            positions_2d, error = clean2d_isolated(graph, timeout_s=8.0)
            if error or not positions_2d:
                return Clean2DAttempt(
                    unavailable=True,
                    message=f"RDKit aislado no disponible: {error or 'sin_coordenadas'}",
                )

            before = {
                aid: (canvas.model.get_atom(aid).x, canvas.model.get_atom(aid).y)
                for aid in target_ids
            }
            target_bond_len = float(getattr(canvas.state, "bond_length", 42.0) or 42.0)

            after = window._rescale_coords_to_bond_length(
                positions_2d, scale_bonds, target_bond_len
            )
            after = window._align_coords_to_reference(before, after)

            label = "Estructura 2D limpiada (RDKit aislado)"
            if len(target_ids) < len(canvas.model.atoms):
                label = "Selección 2D limpiada (RDKit aislado)"
            return self._apply_clean2d_candidate(
                window, target_ids, scale_bonds, before, after,
                target_bond_len, mode, label, cyclic,
            )
        except Exception as exc:
            return Clean2DAttempt(
                unavailable=True,
                message=f"RDKit aislado no disponible: {exc}",
            )

    def _try_direct_rdkit(
        self,
        window: Any,
        atom_ids: set[int],
        target_ids: set[int],
        scale_bonds: list,
        step_ratio: float,
        mode: str,
        status_suffix: str,
        cyclic: bool,
    ) -> Clean2DAttempt:
        """Ruta legacy deshabilitada: RDKit se usa solo vía worker aislado."""
        return Clean2DAttempt(
            unavailable=True,
            message="RDKit directo deshabilitado; use RDKit aislado",
        )

    def _safe_length_polish_only(
        self,
        window: Any,
        target_ids: set[int],
        scale_bonds: list,
        status_suffix: str,
    ) -> Clean2DAttempt:
        """Aplica solo ajuste de longitudes de enlace sin deformar topología."""
        canvas = window.canvas
        before = {
            aid: (canvas.model.get_atom(aid).x, canvas.model.get_atom(aid).y)
            for aid in target_ids
            if aid in canvas.model.atoms
        }
        if not before:
            return Clean2DAttempt(unavailable=True, message="sin_coordenadas")

        target_len = float(getattr(canvas.state, "bond_length", 42.0) or 42.0)
        cyclic = has_cycles(target_ids, scale_bonds)
        if cyclic:
            after = structure_preserving_geometry_polish(
                before, scale_bonds, target_len,
                max_displacement_per_atom=target_len * 8.0,
            )
            mode = "conformer"
        else:
            after = length_only_polish(
                before, scale_bonds, target_len,
                max_iterations=36, damping=0.6,
                max_displacement_per_atom=target_len * 5.0,
            )
            mode = "quick"
        label = "Estructura 2D ajustada con fallback seguro"
        return self._apply_clean2d_candidate(
            window, target_ids, scale_bonds, before, after,
            target_len, mode, label, cyclic=cyclic,
        )


def _get_bonds_from_model(model, atom_ids: set[int]) -> list:
    bonds = []
    for bond in model.bonds.values():
        if bond.a1_id in atom_ids and bond.a2_id in atom_ids:
            bonds.append(bond)
    return bonds


def _position_delta(
    before: tuple[float, float],
    after: tuple[float, float],
) -> float:
    dx = float(after[0]) - float(before[0])
    dy = float(after[1]) - float(before[1])
    return (dx * dx + dy * dy) ** 0.5


def _coords_center(coords: dict[int, tuple[float, float]]) -> tuple[float, float]:
    if not coords:
        return 0.0, 0.0
    return (
        sum(x for x, _y in coords.values()) / len(coords),
        sum(y for _x, y in coords.values()) / len(coords),
    )


def _rotatable_reflection_sides(
    atom_ids: set[int],
    bonds: list,
) -> list[tuple[Any, set[int]]]:
    """Devuelve lados candidatos para reflejar alrededor de enlaces rotatables."""
    adjacency = _adjacency_for_clean_bonds(atom_ids, bonds)
    ring_core = _cycle_core_atoms(atom_ids, adjacency)
    candidates: list[tuple[tuple[int, int, int, int], Any, set[int]]] = []

    for bond in bonds:
        if not _is_rotatable_clean2d_bond(bond, atom_ids):
            continue
        left = _component_without_edge(
            int(bond.a1_id),
            int(bond.a2_id),
            adjacency,
            atom_ids,
        )
        right = _component_without_edge(
            int(bond.a2_id),
            int(bond.a1_id),
            adjacency,
            atom_ids,
        )
        if not left or not right or left & right:
            continue

        side = _preferred_reflection_side(left, right, ring_core)
        if not side:
            continue
        ring_count = len(side & ring_core)
        other_count = len((right if side is left else left) & ring_core)
        side_size = len(side)
        balance = min(len(left), len(right))
        score = (
            0 if ring_count == 0 and other_count > 0 else 1,
            0 if side_size >= 2 else 1,
            -balance,
            side_size,
        )
        candidates.append((score, bond, side))

    candidates.sort(key=lambda item: item[0])
    return [(bond, side) for _score, bond, side in candidates]


def _is_rotatable_clean2d_bond(bond: Any, atom_ids: set[int]) -> bool:
    try:
        a1 = int(getattr(bond, "a1_id", -1))
        a2 = int(getattr(bond, "a2_id", -1))
    except Exception:
        return False
    if a1 not in atom_ids or a2 not in atom_ids:
        return False
    try:
        if not bond_affects_valence(bond):
            return False
    except Exception:
        return False
    if bool(getattr(bond, "is_aromatic", False)):
        return False
    try:
        order = int(getattr(bond, "order", 1) or 1)
    except Exception:
        order = 1
    return order == 1


def _preferred_reflection_side(
    left: set[int],
    right: set[int],
    ring_core: set[int],
) -> set[int]:
    left_ring = len(left & ring_core)
    right_ring = len(right & ring_core)
    if left_ring != right_ring:
        return left if left_ring < right_ring else right
    if len(left) != len(right):
        return left if len(left) < len(right) else right
    return left if min(left) <= min(right) else right


def _adjacency_for_clean_bonds(atom_ids: set[int], bonds: list) -> dict[int, set[int]]:
    adjacency: dict[int, set[int]] = {int(atom_id): set() for atom_id in atom_ids}
    for bond in bonds:
        try:
            a1 = int(bond.a1_id)
            a2 = int(bond.a2_id)
        except Exception:
            continue
        if a1 not in atom_ids or a2 not in atom_ids:
            continue
        try:
            if not bond_affects_valence(bond):
                continue
        except Exception:
            continue
        adjacency.setdefault(a1, set()).add(a2)
        adjacency.setdefault(a2, set()).add(a1)
    return adjacency


def _component_without_edge(
    start: int,
    blocked_neighbor: int,
    adjacency: dict[int, set[int]],
    atom_ids: set[int],
) -> set[int]:
    if start not in atom_ids:
        return set()
    visited: set[int] = set()
    stack = [start]
    while stack:
        node = stack.pop()
        if node in visited:
            continue
        visited.add(node)
        for neigh in adjacency.get(node, set()):
            if {node, neigh} == {start, blocked_neighbor}:
                continue
            if neigh in atom_ids and neigh not in visited:
                stack.append(neigh)
    return visited


def _cycle_core_atoms(
    atom_ids: set[int],
    adjacency: dict[int, set[int]],
) -> set[int]:
    core = {int(atom_id) for atom_id in atom_ids}
    degrees = {
        atom_id: len([neigh for neigh in adjacency.get(atom_id, set()) if neigh in core])
        for atom_id in core
    }
    queue = [atom_id for atom_id, degree in degrees.items() if degree <= 1]
    while queue:
        atom_id = queue.pop()
        if atom_id not in core:
            continue
        core.remove(atom_id)
        for neigh in adjacency.get(atom_id, set()):
            if neigh not in core:
                continue
            degrees[neigh] = degrees.get(neigh, 0) - 1
            if degrees[neigh] <= 1:
                queue.append(neigh)
    return core


def _reflect_side_across_bond(
    before: dict[int, tuple[float, float]],
    bond: Any,
    side: set[int],
) -> dict[int, tuple[float, float]]:
    a1 = int(getattr(bond, "a1_id", -1))
    a2 = int(getattr(bond, "a2_id", -1))
    if a1 not in before or a2 not in before:
        return {}

    x1, y1 = before[a1]
    x2, y2 = before[a2]
    dx = x2 - x1
    dy = y2 - y1
    denom = dx * dx + dy * dy
    if denom <= 1e-9:
        return {}

    reflected: dict[int, tuple[float, float]] = {}
    for atom_id in side:
        if atom_id not in before:
            continue
        px, py = before[atom_id]
        rel_x = px - x1
        rel_y = py - y1
        projection = (rel_x * dx + rel_y * dy) / denom
        foot_x = x1 + projection * dx
        foot_y = y1 + projection * dy
        reflected[atom_id] = (2.0 * foot_x - px, 2.0 * foot_y - py)
    return reflected
