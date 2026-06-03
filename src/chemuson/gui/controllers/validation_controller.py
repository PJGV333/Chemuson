from __future__ import annotations

"""Controller for interactive, undoable chemical validation fixes."""

import copy
from dataclasses import dataclass
from typing import Any, Callable

from PyQt6.QtGui import QUndoCommand

from chemuson.core.model import MolGraph, ValidationCorrectionAction, ValidationIssue


@dataclass(frozen=True)
class ValidationFixResult:
    """Outcome of an attempted validation fix."""

    applied: bool = False
    message: str = ""
    action_id: str = ""


class _ClearAssignedHydrogensCommand(QUndoCommand):
    """Undo command that clears label-assigned hydrogens on one atom."""

    def __init__(self, model: MolGraph, view: Any, atom_id: int) -> None:
        super().__init__("Clear assigned hydrogens")
        self._model = model
        self._view = view
        self._atom_id = int(atom_id)
        atom = model.get_atom(atom_id)
        self._old_explicit_h = getattr(atom, "explicit_h", None)
        self._old_group_h_cap = getattr(atom, "group_h_cap", None)
        self._new_explicit_h = 0 if self._old_explicit_h is not None else None
        self._new_group_h_cap = 0 if self._old_group_h_cap is not None else None

    def redo(self) -> None:
        self._apply(self._new_explicit_h, self._new_group_h_cap)

    def undo(self) -> None:
        self._apply(self._old_explicit_h, self._old_group_h_cap)

    def _apply(self, explicit_h: int | None, group_h_cap: int | None) -> None:
        if self._atom_id not in self._model.atoms:
            return
        atom = self._model.get_atom(self._atom_id)
        atom.explicit_h = explicit_h
        atom.group_h_cap = group_h_cap
        refresher = getattr(self._view, "refresh_atom_labels", None)
        if callable(refresher):
            refresher([self._atom_id])
        validator = getattr(self._view, "validate_structure", None)
        if callable(validator):
            validator()


class ValidationController:
    """Calculates available validation fixes and applies safe ones via undo."""

    @staticmethod
    def issues_for_canvas(canvas: Any) -> dict[int, ValidationIssue]:
        """Run validation on a canvas and attach context-sensitive actions."""
        validator = getattr(canvas, "validate_structure", None)
        if callable(validator):
            validator()
        issues = dict(getattr(canvas, "current_validation_issues", lambda: {})() or {})
        return {
            atom_id: ValidationController.with_available_actions(canvas, issue)
            for atom_id, issue in issues.items()
        }

    @staticmethod
    def with_available_actions(canvas: Any, issue: Any) -> Any:
        """Return an issue copy with only actions that are currently safe."""
        if not isinstance(issue, ValidationIssue):
            return issue
        actions = ValidationController.available_actions(canvas, issue)
        return issue.with_correction_actions(
            ValidationCorrectionAction(
                str(action["id"]),
                str(action["label"]),
                str(action.get("description", "")),
            )
            for action in actions
        )

    @staticmethod
    def available_actions(canvas: Any, issue: ValidationIssue) -> list[dict[str, object]]:
        """Return deterministic actions that pass validation simulation."""
        graph = getattr(canvas, "model", None)
        if graph is None or int(issue.atom_id) not in getattr(graph, "atoms", {}):
            return []
        actions: list[dict[str, object]] = []
        selected_bond_id = ValidationController._selected_single_bond_id(canvas)
        if (
            issue.code == "VALENCE_EXCEEDED"
            and selected_bond_id is not None
            and selected_bond_id in graph.bonds
        ):
            bond = graph.get_bond(selected_bond_id)
            if (
                int(getattr(bond, "order", 1) or 1) > 1
                and (bond.a1_id == issue.atom_id or bond.a2_id == issue.atom_id)
            ):
                new_order = int(bond.order) - 1
                if ValidationController._mutation_is_safe(
                    graph,
                    issue.atom_id,
                    lambda trial: trial.update_bond(selected_bond_id, order=new_order),
                ):
                    actions.append(
                        {
                            "id": "reduce_selected_bond",
                            "label": "Reducir enlace seleccionado",
                            "description": f"Orden {bond.order} -> {new_order}.",
                            "bond_id": selected_bond_id,
                            "new_order": new_order,
                        }
                    )
        charge_action = ValidationController._charge_action(graph, issue.atom_id)
        if charge_action is not None:
            actions.append(charge_action)
        hydrogen_action = ValidationController._clear_hydrogen_action(graph, issue.atom_id)
        if hydrogen_action is not None:
            actions.append(hydrogen_action)
        return actions

    @staticmethod
    def apply_correction(canvas: Any, atom_id: int, action_id: str) -> ValidationFixResult:
        """Apply one selected safe correction through the canvas undo stack."""
        graph = getattr(canvas, "model", None)
        if graph is None or int(atom_id) not in getattr(graph, "atoms", {}):
            return ValidationFixResult(False, "Átomo no disponible", str(action_id))
        issue = dict(getattr(canvas, "current_validation_issues", lambda: {})() or {}).get(int(atom_id))
        if not isinstance(issue, ValidationIssue):
            issue = graph.validate_detailed().get(int(atom_id))
        if not isinstance(issue, ValidationIssue):
            return ValidationFixResult(False, "No hay issue de validación activo", str(action_id))
        actions = ValidationController.available_actions(canvas, issue)
        selected = next((action for action in actions if action.get("id") == action_id), None)
        if selected is None:
            return ValidationFixResult(False, "La corrección no es inequívoca o ya no es segura", str(action_id))

        from chemuson.gui.commands import ChangeBondCommand, ChangeChargeCommand

        command: QUndoCommand
        if action_id == "reduce_selected_bond":
            command = ChangeBondCommand(
                graph,
                canvas,
                int(selected["bond_id"]),
                new_order=int(selected["new_order"]),
            )
        elif action_id == "adjust_charge":
            command = ChangeChargeCommand(graph, canvas, int(atom_id), int(selected["new_charge"]))
        elif action_id == "clear_assigned_h":
            command = _ClearAssignedHydrogensCommand(graph, canvas, int(atom_id))
        else:
            return ValidationFixResult(False, "Acción de corrección desconocida", str(action_id))

        undo_stack = getattr(canvas, "undo_stack", None)
        if undo_stack is None:
            return ValidationFixResult(False, "Undo no disponible", str(action_id))
        undo_stack.push(command)
        return ValidationFixResult(True, str(selected.get("label", "Corrección aplicada")), str(action_id))

    @staticmethod
    def _charge_action(graph: MolGraph, atom_id: int) -> dict[str, object] | None:
        atom = graph.get_atom(atom_id)
        current = int(getattr(atom, "charge", 0) or 0)
        candidates: list[dict[str, object]] = []
        for delta in (1, -1):
            new_charge = current + delta
            if ValidationController._mutation_is_safe(
                graph,
                atom_id,
                lambda trial, charge=new_charge: trial.update_atom_charge(atom_id, charge),
            ):
                candidates.append(
                    {
                        "id": "adjust_charge",
                        "label": f"Ajustar carga formal a {new_charge:+d}",
                        "description": f"Carga {current:+d} -> {new_charge:+d}.",
                        "new_charge": new_charge,
                    }
                )
        return candidates[0] if len(candidates) == 1 else None

    @staticmethod
    def _clear_hydrogen_action(graph: MolGraph, atom_id: int) -> dict[str, object] | None:
        atom = graph.get_atom(atom_id)
        if getattr(atom, "explicit_h", None) is None and getattr(atom, "group_h_cap", None) is None:
            return None

        def mutate(trial: MolGraph) -> None:
            trial_atom = trial.get_atom(atom_id)
            if getattr(trial_atom, "explicit_h", None) is not None:
                trial_atom.explicit_h = 0
            if getattr(trial_atom, "group_h_cap", None) is not None:
                trial_atom.group_h_cap = 0

        if not ValidationController._mutation_is_safe(graph, atom_id, mutate):
            return None
        return {
            "id": "clear_assigned_h",
            "label": "Limpiar H asignados",
            "description": "Quita H declarados en la etiqueta del átomo.",
        }

    @staticmethod
    def _mutation_is_safe(graph: MolGraph, atom_id: int, mutate: Callable[[MolGraph], None]) -> bool:
        before = graph.validate_detailed()
        before_target = before.get(int(atom_id))
        if before_target is None:
            return False
        before_neighbors = ValidationController._neighbor_ids(graph, atom_id)
        trial = copy.deepcopy(graph)
        try:
            mutate(trial)
        except Exception:
            return False
        after = trial.validate_detailed()
        after_target = after.get(int(atom_id))
        if after_target is None:
            target_improved = True
        else:
            target_improved = float(after_target.observed_valence) < float(before_target.observed_valence)
        if not target_improved:
            return False

        before_neighbor_issues = before_neighbors & set(before)
        after_neighbor_issues = before_neighbors & set(after)
        if after_neighbor_issues - before_neighbor_issues:
            return False
        if set(after) - set(before) - {int(atom_id)}:
            return False
        return True

    @staticmethod
    def _neighbor_ids(graph: MolGraph, atom_id: int) -> set[int]:
        neighbors: set[int] = set()
        for bond in graph.bonds.values():
            if bond.a1_id == atom_id:
                neighbors.add(int(bond.a2_id))
            elif bond.a2_id == atom_id:
                neighbors.add(int(bond.a1_id))
        return neighbors

    @staticmethod
    def _selected_single_bond_id(canvas: Any) -> int | None:
        method = getattr(canvas, "_selected_single_bond_id", None)
        if callable(method):
            return method()
        state = getattr(canvas, "state", None)
        selected = getattr(state, "selected_bonds", set()) if state is not None else set()
        valid = [int(bond_id) for bond_id in selected if bond_id in getattr(canvas.model, "bonds", {})]
        return valid[0] if len(valid) == 1 else None
