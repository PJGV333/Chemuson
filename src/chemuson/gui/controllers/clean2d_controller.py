from __future__ import annotations

from typing import Any


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

    def run_clean_2d(self, window: Any, step_ratio: float, fallback_iterations: int, status_suffix: str) -> None:
        canvas = window.canvas
        atom_ids, bonds = canvas._selected_structure_ids()
        target_ids = atom_ids if atom_ids else set(canvas.model.atoms.keys())
        if not target_ids:
            return
        scale_bonds = bonds if atom_ids else list(canvas.model.bonds.values())

        if self._is_acyclic_structure(target_ids, scale_bonds):
            canvas.clean_2d_fallback(target_ids, iterations=max(40, fallback_iterations))
            msg = "Selección 2D limpiada" if atom_ids else "Estructura 2D limpiada"
            window.statusBar().showMessage(f"{msg} {status_suffix} (acíclico)")
            return

        try:
            from chemuson.chemio.rdkit_io import molgraph_to_rdkit_with_map
            from rdkit import Chem
            from rdkit.Chem import AllChem
            from chemuson.gui.commands import MoveAtomsCommand

            graph = canvas._build_selection_graph(atom_ids, bonds) if atom_ids else canvas.graph
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

            cmd = MoveAtomsCommand(canvas.model, canvas, before, after)
            canvas.undo_stack.push(cmd)
            canvas._update_selection_overlay()
            msg = "Selección 2D limpiada" if atom_ids else "Estructura 2D limpiada"
            window.statusBar().showMessage(f"{msg} {status_suffix}".strip())
            return
        except Exception as exc:
            message = str(exc)
            if "No module named" in message and "rdkit" in message:
                canvas.clean_2d_fallback(target_ids, iterations=fallback_iterations)
                msg = "Selección 2D limpiada" if atom_ids else "Estructura 2D limpiada"
                window.statusBar().showMessage(f"{msg} {status_suffix} (básico)")
                return
            window.statusBar().showMessage(f"Error: {exc}")
