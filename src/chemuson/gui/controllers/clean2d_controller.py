from __future__ import annotations

from dataclasses import dataclass, field
import math
from typing import Any

from chemuson.clean2d import (
    Clean2DParameters,
    Clean2DQualityReport,
    evaluate_clean2d_layout,
    has_cycles,
    is_clean2d_candidate_safe,
    length_only_polish,
    optimize_clean2d_positions,
)


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
            if mode in {"quick", "publication"}:
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
                    return Clean2DAttempt(message="Estructura 2D ya estaba limpia")

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
            after = length_only_polish(
                before, bonds, target_len,
                max_iterations=6, damping=0.4,
                max_displacement_per_atom=target_len * 1.5,
            )
            label = f"{status_suffix} (ajuste)".strip()
            return self._apply_clean2d_candidate(
                window, atom_ids, bonds, before, after,
                target_len, mode, label, cyclic=True,
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
        cyclic = has_cycles(target_ids, scale_bonds)

        # 1) RDKit aislado es el backend preferido para quick/publication
        if mode in {"quick", "publication"}:
            attempt = self._try_isolated_rdkit2d(
                window, target_ids, scale_bonds, cyclic, mode,
            )
            if attempt.applied:
                window.statusBar().showMessage(attempt.message)
                return
            if attempt.rejected:
                window.statusBar().showMessage(attempt.message)
                if cyclic:
                    self._safe_length_polish_only(
                        window, target_ids, scale_bonds, status_suffix
                    )
                    return
                if self._polish_with_v2(
                    window, target_ids, mode=mode, status_suffix=status_suffix
                ).applied:
                    return
                self._safe_length_polish_only(
                    window, target_ids, scale_bonds, status_suffix
                )
                return
            if attempt.message:
                # unavailable con mensaje — mostrar para debug, pero continuar
                pass

        # 2) Fallback a importación directa de RDKit
        if mode in {"quick", "publication"}:
            direct_attempt = self._try_direct_rdkit(
                window, atom_ids, target_ids, scale_bonds,
                step_ratio, mode, status_suffix, cyclic,
            )
            if direct_attempt.applied:
                window.statusBar().showMessage(direct_attempt.message)
                return
            if direct_attempt.rejected:
                window.statusBar().showMessage(direct_attempt.message)
                if cyclic:
                    self._safe_length_polish_only(
                        window, target_ids, scale_bonds, status_suffix
                    )
                    return
                if self._polish_with_v2(
                    window, target_ids, mode=mode, status_suffix=status_suffix
                ).applied:
                    return
                self._safe_length_polish_only(
                    window, target_ids, scale_bonds, status_suffix
                )
                return

        # 3) Fallback interno del canvas (force-layout básico)
        canvas.clean_2d_fallback(target_ids, iterations=fallback_iterations)
        polish = self._polish_with_v2(
            window, target_ids, mode=mode, status_suffix=status_suffix
        )
        if polish.applied:
            window.statusBar().showMessage(polish.message)
            return
        if polish.rejected:
            window.statusBar().showMessage(polish.message)
            return
        msg = "Selección 2D limpiada" if atom_ids else "Estructura 2D limpiada"
        window.statusBar().showMessage(f"{msg} {status_suffix} (básico)")

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
        """Ruta de RDKit directo en el proceso principal.

        Returns Clean2DAttempt — nunca levanta excepción por flujo de
        control; los errores de importación se reportan como unavailable.
        """
        canvas = window.canvas
        try:
            from chemuson.chemio.rdkit_io import molgraph_to_rdkit_with_map
            from rdkit import Chem
            from rdkit.Chem import AllChem
        except Exception as exc:
            return Clean2DAttempt(
                unavailable=True,
                message=f"RDKit directo no disponible: {exc}",
            )

        try:
            graph = canvas._build_selection_graph(atom_ids, scale_bonds) if atom_ids else canvas.graph

            mol, id_map = molgraph_to_rdkit_with_map(graph)
            for aid, rd_idx in id_map.items():
                try:
                    mol.GetAtomWithIdx(rd_idx).SetIntProp("_chemuson_id", int(aid))
                except Exception:
                    continue

            before = {
                aid: (canvas.model.get_atom(aid).x, canvas.model.get_atom(aid).y)
                for aid in target_ids
            }
            before_avg_len = window._average_bond_length(before, scale_bonds)

            raw_after: dict[int, tuple[float, float]] = {}
            used_no_h_layout = False
            try:
                mol_no_h = Chem.RemoveHs(Chem.Mol(mol), sanitize=True)
                if 0 < mol_no_h.GetNumAtoms() < mol.GetNumAtoms():
                    AllChem.Compute2DCoords(mol_no_h)
                    conf_no_h = mol_no_h.GetConformer()
                    for atom in mol_no_h.GetAtoms():
                        if not atom.HasProp("_chemuson_id"):
                            continue
                        aid = int(atom.GetIntProp("_chemuson_id"))
                        if aid not in target_ids:
                            continue
                        pos = conf_no_h.GetAtomPosition(atom.GetIdx())
                        raw_after[aid] = (pos.x, pos.y)
                    if raw_after:
                        atom_elements = {
                            aid: canvas.model.get_atom(aid).element for aid in target_ids
                        }
                        raw_after = window._project_missing_hydrogen_coords(
                            before=before,
                            after=raw_after,
                            bonds=scale_bonds,
                            atom_elements=atom_elements,
                        )
                        used_no_h_layout = True
            except Exception:
                used_no_h_layout = False

            if not used_no_h_layout:
                AllChem.Compute2DCoords(mol)
                conf = mol.GetConformer()
                for aid, rd_idx in id_map.items():
                    if aid not in target_ids:
                        continue
                    pos = conf.GetAtomPosition(rd_idx)
                    raw_after[aid] = (pos.x, pos.y)

            target_bond_len = before_avg_len if before_avg_len > 1e-6 else float(canvas.state.bond_length)
            raw_after = window._rescale_coords_to_bond_length(raw_after, scale_bonds, target_bond_len)
            raw_after = window._align_coords_to_reference(before, raw_after)
            after_cx, after_cy = _coords_center(raw_after)
            before_cx, before_cy = _coords_center(before)
            after = {}
            for aid, (x, y) in raw_after.items():
                x = x - after_cx + before_cx
                y = y - after_cy + before_cy
                if step_ratio < 1.0:
                    bx, by = before[aid]
                    x = bx + (x - bx) * step_ratio
                    y = by + (y - by) * step_ratio
                after[aid] = (x, y)

            after = window._rescale_coords_to_bond_length(after, scale_bonds, target_bond_len)
            final_cx, final_cy = _coords_center(after)
            if abs(final_cx - before_cx) > 1e-9 or abs(final_cy - before_cy) > 1e-9:
                after = {
                    aid: (x - final_cx + before_cx, y - final_cy + before_cy)
                    for aid, (x, y) in after.items()
                }

            label = f"Estructura 2D limpiada {status_suffix}".strip()
            if len(target_ids) < len(canvas.model.atoms):
                label = f"Selección 2D limpiada {status_suffix}".strip()
            return self._apply_clean2d_candidate(
                window, target_ids, scale_bonds, before, after,
                target_bond_len, mode, label, cyclic,
            )
        except Exception as exc:
            return Clean2DAttempt(
                unavailable=True,
                message=f"RDKit directo falló: {exc}",
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
        after = length_only_polish(
            before, scale_bonds, target_len,
            max_iterations=8, damping=0.4,
            max_displacement_per_atom=target_len * 1.5,
        )
        label = "Estructura 2D ajustada con fallback seguro"
        return self._apply_clean2d_candidate(
            window, target_ids, scale_bonds, before, after,
            target_len, "quick", label, cyclic=True,
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
