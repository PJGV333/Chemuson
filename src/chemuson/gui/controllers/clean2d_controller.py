from __future__ import annotations

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

    def _polish_with_v2(
        self,
        window: Any,
        atom_ids: set[int],
        *,
        mode: str,
        status_suffix: str,
    ) -> bool:
        """Aplica pulido clean2d_v2 como un único comando undoable.

        Para moléculas con ciclos NO aplica clean2d_v2, solo un
        reescalado ligero de longitudes de enlace si es necesario.
        """
        canvas = window.canvas
        if not atom_ids:
            return False
        before = {
            aid: (canvas.model.get_atom(aid).x, canvas.model.get_atom(aid).y)
            for aid in atom_ids
            if aid in canvas.model.atoms
        }
        if not before:
            return False

        bonds = _get_bonds_from_model(canvas.model, atom_ids)
        target_len = float(getattr(canvas.state, "bond_length", 42.0) or 42.0)
        cyclic = has_cycles(atom_ids, bonds)

        if cyclic:
            # Para estructuras cíclicas: solo reescalado ligero + safety gate
            after = length_only_polish(
                before, bonds, target_len,
                max_iterations=6,
                damping=0.4,
                max_displacement_per_atom=target_len * 1.5,
            )
            after = {
                aid: pos
                for aid, pos in after.items()
                if aid in before and _position_delta(before[aid], pos) > 0.5
            }
            if not after:
                return False

            report = evaluate_clean2d_layout(
                atom_ids, bonds, before,
                {aid: after.get(aid, before[aid]) for aid in atom_ids},
                target_len, is_cyclic=True,
            )
            if not is_clean2d_candidate_safe(report, mode):
                window.statusBar().showMessage(
                    "Limpieza 2D omitida: el candidato deformaba la estructura"
                )
                return False

            from chemuson.gui.commands import MoveAtomsCommand

            canvas.undo_stack.push(MoveAtomsCommand(canvas.model, canvas, before, after))
            canvas._update_selection_overlay()
            msg = (
                "Selección 2D limpiada"
                if len(atom_ids) < len(canvas.model.atoms)
                else "Estructura 2D limpiada"
            )
            window.statusBar().showMessage(f"{msg} {status_suffix} (ajuste)".strip())
            return True

        # Para estructuras acíclicas: clean2d_v2 completo + safety gate
        params = (
            Clean2DParameters.publication(target_len)
            if mode == "publication"
            else Clean2DParameters.quick(target_len)
        )
        after_v2 = optimize_clean2d_positions(canvas.model, atom_ids, params)
        after = {
            aid: pos
            for aid, pos in after_v2.items()
            if aid in before and _position_delta(before[aid], pos) > 0.01
        }
        if not after:
            return False

        report = evaluate_clean2d_layout(
            atom_ids, bonds, before,
            {aid: after.get(aid, before[aid]) for aid in atom_ids},
            target_len, is_cyclic=False,
        )
        if not is_clean2d_candidate_safe(report, mode):
            window.statusBar().showMessage(
                "Limpieza 2D omitida: el candidato deformaba la estructura"
            )
            return False

        from chemuson.gui.commands import MoveAtomsCommand

        canvas.undo_stack.push(MoveAtomsCommand(canvas.model, canvas, before, after))
        canvas._update_selection_overlay()
        msg = (
            "Selección 2D limpiada"
            if len(atom_ids) < len(canvas.model.atoms)
            else "Estructura 2D pulida"
        )
        window.statusBar().showMessage(f"{msg} {status_suffix} (clean2d_v2)".strip())
        return True

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
            return
        scale_bonds = bonds if atom_ids else list(canvas.model.bonds.values())
        cyclic = has_cycles(target_ids, scale_bonds)

        # Usar worker RDKit aislado (preferido sobre importación directa)
        attempted_isolated = False
        if mode in {"quick", "publication"}:
            result = self._try_isolated_rdkit2d(window, target_ids, scale_bonds, cyclic)
            if result is not None:
                return result

        # Fallback a importación directa de RDKit (solo si habilitado)
        try:
            self._run_direct_rdkit(
                window, atom_ids, target_ids, scale_bonds,
                step_ratio, mode, status_suffix, cyclic,
            )
            return
        except Exception as exc:
            message = str(exc)
            if "No module named" in message and "rdkit" in message:
                pass
            # Intentar polish seguro si falló RDKit
            if cyclic:
                self._safe_length_polish_only(
                    window, target_ids, scale_bonds, status_suffix
                )
                return
            canvas.clean_2d_fallback(target_ids, iterations=fallback_iterations)
            if self._polish_with_v2(window, target_ids, mode=mode, status_suffix=status_suffix):
                return
            msg = "Selección 2D limpiada" if atom_ids else "Estructura 2D limpiada"
            window.statusBar().showMessage(f"{msg} {status_suffix} (básico)")
            return

    def _try_isolated_rdkit2d(
        self,
        window: Any,
        target_ids: set[int],
        scale_bonds: list,
        cyclic: bool,
    ) -> None | bool:
        """Intenta usar RDKit en worker aislado para Clean2D. Retorna True/None."""
        canvas = window.canvas
        try:
            from chemuson.chemio.rdkit_safe import clean2d_isolated

            graph = canvas._build_selection_graph(
                target_ids, scale_bonds
            ) if target_ids != set(canvas.model.atoms.keys()) else canvas.graph
            positions_2d, error = clean2d_isolated(graph, timeout_s=8.0)
            if error or not positions_2d:
                return None

            before = {
                aid: (canvas.model.get_atom(aid).x, canvas.model.get_atom(aid).y)
                for aid in target_ids
            }
            target_bond_len = float(getattr(canvas.state, "bond_length", 42.0) or 42.0)

            after = window._rescale_coords_to_bond_length(
                positions_2d, scale_bonds, target_bond_len
            )
            after = window._align_coords_to_reference(before, after)
            before_cx, before_cy = window._coords_center(before)
            after_cx, after_cy = window._coords_center(after)
            final = {}
            for aid, (x, y) in after.items():
                adjusted_x = x - after_cx + before_cx
                adjusted_y = y - after_cy + before_cy
                final[aid] = (adjusted_x, adjusted_y)

            report = evaluate_clean2d_layout(
                target_ids, scale_bonds, before, final,
                target_bond_len, is_cyclic=cyclic,
            )
            if not is_clean2d_candidate_safe(report):
                window.statusBar().showMessage(
                    "Limpieza 2D omitida: el candidato deformaba la estructura"
                )
                return False

            from chemuson.gui.commands import MoveAtomsCommand

            cmd = MoveAtomsCommand(canvas.model, canvas, before, final)
            canvas.undo_stack.push(cmd)
            canvas._update_selection_overlay()
            msg = "Selección 2D limpiada" if len(target_ids) < len(canvas.model.atoms) else "Estructura 2D limpiada"
            window.statusBar().showMessage(f"{msg} (RDKit aislado)")
            return True
        except Exception:
            return None

    def _safe_length_polish_only(
        self,
        window: Any,
        target_ids: set[int],
        scale_bonds: list,
        status_suffix: str,
    ) -> None:
        """Aplica solo ajuste de longitudes de enlace sin deformar topología."""
        canvas = window.canvas
        before = {
            aid: (canvas.model.get_atom(aid).x, canvas.model.get_atom(aid).y)
            for aid in target_ids
            if aid in canvas.model.atoms
        }
        target_len = float(getattr(canvas.state, "bond_length", 42.0) or 42.0)
        after = length_only_polish(
            before, scale_bonds, target_len,
            max_iterations=8, damping=0.4,
            max_displacement_per_atom=target_len * 1.5,
        )
        after = {
            aid: pos
            for aid, pos in after.items()
            if aid in before and _position_delta(before[aid], pos) > 0.5
        }
        if not after:
            window.statusBar().showMessage(
                "Estructura 2D ya estaba limpia"
            )
            return

        report = evaluate_clean2d_layout(
            target_ids, scale_bonds, before,
            {aid: after.get(aid, before[aid]) for aid in target_ids},
            target_len, is_cyclic=True,
        )
        if not is_clean2d_candidate_safe(report):
            window.statusBar().showMessage(
                "Limpieza 2D omitida: el candidato deformaba la estructura"
            )
            return

        from chemuson.gui.commands import MoveAtomsCommand

        canvas.undo_stack.push(MoveAtomsCommand(canvas.model, canvas, before, after))
        canvas._update_selection_overlay()
        window.statusBar().showMessage(
            "Estructura 2D ajustada con fallback seguro"
        )

    def _run_direct_rdkit(
        self,
        window: Any,
        atom_ids: set[int],
        target_ids: set[int],
        scale_bonds: list,
        step_ratio: float,
        mode: str,
        status_suffix: str,
        cyclic: bool,
    ) -> None:
        """Ruta de RDKit directo (importación directa, solo si habilitado)."""
        from chemuson.chemio.rdkit_io import molgraph_to_rdkit_with_map
        from rdkit import Chem
        from rdkit.Chem import AllChem
        from chemuson.gui.commands import MoveAtomsCommand

        canvas = window.canvas
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
        before_cx, before_cy = window._coords_center(before)
        before_avg_len = window._average_bond_length(before, scale_bonds)

        raw_after = {}
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
        after_cx, after_cy = window._coords_center(raw_after)
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
        final_cx, final_cy = window._coords_center(after)
        if abs(final_cx - before_cx) > 1e-9 or abs(final_cy - before_cy) > 1e-9:
            after = {
                aid: (x - final_cx + before_cx, y - final_cy + before_cy)
                for aid, (x, y) in after.items()
            }

        # Safety gate para RDKit directo
        report = evaluate_clean2d_layout(
            target_ids, bonds, before, after,
            target_bond_len, is_cyclic=cyclic,
        )
        if not is_clean2d_candidate_safe(report, mode):
            if cyclic:
                self._safe_length_polish_only(
                    window, target_ids, scale_bonds, status_suffix
                )
                return
            canvas.clean_2d_fallback(target_ids, iterations=fallback_iterations)
            if self._polish_with_v2(window, target_ids, mode=mode, status_suffix=status_suffix):
                return
            window.statusBar().showMessage(
                "Limpieza 2D omitida: el candidato deformaba la estructura"
            )
            return

        cmd = MoveAtomsCommand(canvas.model, canvas, before, after)
        canvas.undo_stack.push(cmd)
        canvas._update_selection_overlay()
        if mode in {"quick", "publication"} and not cyclic:
            self._polish_with_v2(window, target_ids, mode=mode, status_suffix=status_suffix)
            return
        msg = "Selección 2D limpiada" if atom_ids else "Estructura 2D limpiada"
        window.statusBar().showMessage(f"{msg} {status_suffix}".strip())
        return


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
