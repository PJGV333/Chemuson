from __future__ import annotations

"""Backend RDKit aislado para conformaciones y optimización 3D."""

from collections.abc import Callable, Iterator, Mapping
import math
import threading
import time
from typing import Any

from chemuson.core.model import MolGraph
from chemuson.geometry3d.model import CoordinateSet3D, OptimizationFrame, OptimizationResult, OptimizationSettings


def generate_conformer(
    graph: MolGraph,
    settings: OptimizationSettings | None = None,
    *,
    cancel_event: threading.Event | None = None,
    on_progress: Callable[[Mapping[str, object]], None] | None = None,
) -> OptimizationResult:
    """Genera geometría con worker RDKit aislado, etapas y errores trazables."""
    settings = settings or OptimizationSettings()
    from chemuson.chemio.rdkit_safe import conformer_3d_stream_isolated

    started = time.monotonic()
    terminal: dict[str, Any] | None = None
    last_stage = "worker_start"
    try:
        events = conformer_3d_stream_isolated(
            graph,
            forcefield=settings.forcefield.value,
            max_iters=settings.max_iters,
            seed=settings.seed,
            timeout_s=settings.timeout_s,
            cancel_event=cancel_event,
        )
        for event in events:
            if event.get("event") == "progress":
                last_stage = str(event.get("stage") or last_stage)
                if on_progress is not None:
                    try:
                        on_progress(dict(event))
                    except Exception:
                        pass
            if event.get("event") in {"done", "error", "cancelled"}:
                terminal = event
    except Exception as exc:
        diagnostics = {
            "error_code": "worker_failed",
            "stage": last_stage,
            "detail": str(exc),
            "graph": {
                "atom_count": len(getattr(graph, "atoms", {})),
                "bond_count": len(getattr(graph, "bonds", {})),
            },
            "timings_s": {"total": max(0.0, time.monotonic() - started)},
        }
        return OptimizationResult(
            None,
            method="rdkit",
            message=f"worker_failed: {exc}",
            metadata={"diagnostics": diagnostics},
        )

    if terminal is None:
        diagnostics = {
            "error_code": "missing_terminal_event",
            "stage": last_stage,
            "graph": {
                "atom_count": len(getattr(graph, "atoms", {})),
                "bond_count": len(getattr(graph, "bonds", {})),
            },
            "timings_s": {"total": max(0.0, time.monotonic() - started)},
        }
        return OptimizationResult(
            None,
            method="rdkit",
            message="missing_terminal_event",
            metadata={"diagnostics": diagnostics},
        )

    raw_diagnostics = terminal.get("diagnostics", {})
    diagnostics = (
        dict(raw_diagnostics) if isinstance(raw_diagnostics, Mapping) else {}
    )
    error_code = str(terminal.get("error", "") or "")
    detail = str(terminal.get("detail", "") or "")
    if terminal.get("event") != "done" or not terminal.get("ok"):
        diagnostics.setdefault("error_code", error_code or terminal.get("event", "worker_failed"))
        diagnostics.setdefault("stage", terminal.get("stage") or last_stage)
        diagnostics.setdefault("detail", detail)
        message = ": ".join(part for part in (error_code, detail) if part)
        if not message:
            message = str(terminal.get("event", "worker_failed"))
        return OptimizationResult(
            None,
            method="rdkit",
            message=message,
            metadata={"diagnostics": diagnostics},
        )

    positions = _event_positions(terminal.get("positions"))
    if not positions:
        diagnostics.update(
            {
                "error_code": "invalid_positions",
                "stage": "publication",
            }
        )
        return OptimizationResult(
            None,
            method="rdkit",
            message="invalid_positions",
            metadata={"diagnostics": diagnostics},
        )
    raw_metadata = terminal.get("metadata", {})
    metadata = dict(raw_metadata) if isinstance(raw_metadata, Mapping) else {}
    metadata["diagnostics"] = diagnostics
    energy = _float_or_none(terminal.get("energy", metadata.get("energy")))
    method = str(terminal.get("method", metadata.get("method", "rdkit")))
    coordinates = CoordinateSet3D(
        positions=positions,
        source=str(metadata.get("source", "rdkit")),
        method=method,
        energy=energy,
        metadata=metadata,
    )
    message = str(terminal.get("message", "") or "")
    converged = bool(metadata.get("converged", False))
    return OptimizationResult(
        coordinates,
        converged=converged,
        energy=energy,
        method=method,
        message=message,
        metadata={**metadata, "diagnostics": diagnostics},
    )


def optimize(
    graph: MolGraph,
    coordset: CoordinateSet3D | None = None,
    settings: OptimizationSettings | None = None,
) -> OptimizationResult:
    """Optimiza geometría con UFF/MMFF usando RDKit aislado."""
    settings = settings or OptimizationSettings()
    try:
        from chemuson.chemio.rdkit_safe import optimize_3d_isolated

        positions, metadata, error = optimize_3d_isolated(
            graph,
            coordinates=coordset.normalized_positions() if coordset is not None else None,
            forcefield=settings.forcefield.value,
            max_iters=settings.max_iters,
            steps_per_update=settings.steps_per_update,
            seed=settings.seed,
            timeout_s=settings.timeout_s,
        )
    except Exception as exc:
        return OptimizationResult(None, method="rdkit", message=str(exc))
    if error or not positions:
        return OptimizationResult(None, method="rdkit", message=str(error or "empty_optimization"))
    coordset = CoordinateSet3D(
        positions=positions,
        source=str(metadata.get("source", "rdkit")),
        method=str(metadata.get("method", settings.forcefield.value)),
        energy=_float_or_none(metadata.get("energy")),
        metadata=metadata,
    )
    converged = bool(metadata.get("converged", False))
    frames = _frames_from_metadata(metadata, coordset.source if coordset is not None else "rdkit")
    return OptimizationResult(
        coordset,
        converged=converged,
        energy=coordset.energy,
        method=coordset.method,
        message=str(metadata.get("message", "")),
        frames=frames,
    )


def optimize_iter(
    graph: MolGraph,
    coordset: CoordinateSet3D | None = None,
    settings: OptimizationSettings | None = None,
    *,
    fixed_atom_ids: set[int] | tuple[int, ...] = (),
    cancel_event: threading.Event | None = None,
) -> Iterator[OptimizationFrame]:
    """Emite cada frame NDJSON recibido mientras RDKit sigue optimizando."""
    settings = settings or OptimizationSettings()
    from chemuson.chemio.rdkit_safe import optimize_3d_stream_isolated

    events = optimize_3d_stream_isolated(
        graph,
        coordinates=coordset.normalized_positions() if coordset is not None else None,
        forcefield=settings.forcefield.value,
        max_iters=settings.max_iters,
        steps_per_update=settings.steps_per_update,
        seed=settings.seed,
        fixed_atom_ids=fixed_atom_ids,
        timeout_s=settings.timeout_s,
        cancel_event=cancel_event,
    )
    for event in events:
        frame = _frame_from_event(event, source="rdkit")
        if frame is not None:
            yield frame
        elif event.get("event") in {"error", "cancelled"}:
            error = str(event.get("error", event.get("event", "worker_error")))
            detail = str(event.get("detail", "") or "")
            raise RuntimeError(f"{error}: {detail}".rstrip(": "))


def optimize_stream(
    graph: MolGraph,
    coordset: CoordinateSet3D | None = None,
    settings: OptimizationSettings | None = None,
    *,
    fixed_atom_ids: set[int] | tuple[int, ...] = (),
    cancel_event: threading.Event | None = None,
    on_frame: Callable[[OptimizationFrame], None] | None = None,
) -> OptimizationResult:
    """Consume el flujo RDKit, notificando frames antes del resultado terminal."""
    settings = settings or OptimizationSettings()
    from chemuson.chemio.rdkit_safe import optimize_3d_stream_isolated

    frames: list[OptimizationFrame] = []
    terminal: dict | None = None
    try:
        events = optimize_3d_stream_isolated(
            graph,
            coordinates=coordset.normalized_positions() if coordset is not None else None,
            forcefield=settings.forcefield.value,
            max_iters=settings.max_iters,
            steps_per_update=settings.steps_per_update,
            seed=settings.seed,
            fixed_atom_ids=fixed_atom_ids,
            timeout_s=settings.timeout_s,
            cancel_event=cancel_event,
        )
        for event in events:
            frame = _frame_from_event(event, source="rdkit")
            if frame is not None:
                frames.append(frame)
                if on_frame is not None:
                    on_frame(frame)
            if event.get("event") in {"done", "error", "cancelled"}:
                terminal = event
    except Exception as exc:
        return OptimizationResult(
            None,
            method=settings.forcefield.value,
            message=f"stream_failed:{exc}",
            frames=tuple(frames),
        )

    if terminal is None:
        return OptimizationResult(
            None,
            method=settings.forcefield.value,
            message="missing_terminal_event",
            frames=tuple(frames),
        )
    if terminal.get("event") != "done" or not terminal.get("ok"):
        error = str(terminal.get("error", terminal.get("event", "worker_error")))
        detail = str(terminal.get("detail", "") or "")
        return OptimizationResult(
            None,
            method=str(terminal.get("method", settings.forcefield.value)),
            message=f"{error}: {detail}".rstrip(": "),
            frames=tuple(frames),
        )

    positions = _event_positions(terminal.get("positions"))
    if not positions:
        return OptimizationResult(
            None,
            method=str(terminal.get("method", settings.forcefield.value)),
            message="invalid_positions",
            frames=tuple(frames),
        )
    metadata = terminal.get("metadata", {})
    metadata = dict(metadata) if isinstance(metadata, dict) else {}
    metadata.setdefault("energy_unit", str(terminal.get("energy_unit", "kcal/mol")))
    metadata.setdefault("forcefield", str(terminal.get("forcefield", terminal.get("method", ""))))
    energy = _float_or_none(terminal.get("energy"))
    coordinates = CoordinateSet3D(
        positions=positions,
        source=str(metadata.get("source", "rdkit")),
        method=str(terminal.get("method", metadata.get("method", settings.forcefield.value))),
        energy=energy,
        metadata=metadata,
    )
    return OptimizationResult(
        coordinates,
        converged=bool(terminal.get("converged", False)),
        energy=energy,
        method=coordinates.method,
        frames=tuple(frames),
        metadata=metadata,
    )


def _derived_hydrogen_metadata(payload: dict) -> dict[str, object]:
    """Copy scene-only hydrogen geometry carried alongside documentary positions."""
    return {
        key: payload[key]
        for key in ("generated_atoms", "generated_bonds")
        if isinstance(payload.get(key), list)
    }


def _frame_from_event(event: dict, *, source: str) -> OptimizationFrame | None:
    if event.get("event") != "frame":
        return None
    positions = _event_positions(event.get("positions"))
    energy = _float_or_none(event.get("energy"))
    if not positions or energy is None or not math.isfinite(energy):
        return None
    method = str(event.get("forcefield", "rdkit"))
    step = int(event.get("step", event.get("sequence", 0)) or 0)
    return OptimizationFrame(
        step=step,
        coordinates=CoordinateSet3D(
            positions=positions,
            source=source,
            method=method,
            energy=energy,
            metadata={
                "forcefield": method,
                "energy_unit": str(event.get("energy_unit", "kcal/mol")),
                "frame": True,
                **_derived_hydrogen_metadata(event),
            },
        ),
        energy=energy,
        converged=bool(event.get("converged", False)),
    )


def _event_positions(raw_positions: object) -> dict[int, tuple[float, float, float]]:
    if not isinstance(raw_positions, dict):
        return {}
    positions: dict[int, tuple[float, float, float]] = {}
    try:
        for atom_id, raw_point in raw_positions.items():
            if not isinstance(raw_point, (list, tuple)) or len(raw_point) != 3:
                continue
            point = tuple(float(component) for component in raw_point)
            if not all(math.isfinite(component) for component in point):
                continue
            positions[int(atom_id)] = (point[0], point[1], point[2])
    except (TypeError, ValueError, OverflowError):
        return {}
    return positions


def _float_or_none(value: object) -> float | None:
    try:
        return None if value is None else float(value)
    except Exception:
        return None


def _frames_from_metadata(metadata: dict, source: str) -> tuple[OptimizationFrame, ...]:
    raw_frames = metadata.get("frames", [])
    if not isinstance(raw_frames, list):
        return ()
    frames: list[OptimizationFrame] = []
    method = str(metadata.get("method", "rdkit"))
    for item in raw_frames:
        if not isinstance(item, dict):
            continue
        raw_positions = item.get("positions", {})
        if not isinstance(raw_positions, dict):
            continue
        positions: dict[int, tuple[float, float, float]] = {}
        for atom_id, coords in raw_positions.items():
            if not isinstance(coords, (list, tuple)) or len(coords) != 3:
                continue
            positions[int(atom_id)] = (float(coords[0]), float(coords[1]), float(coords[2]))
        if not positions:
            continue
        energy = _float_or_none(item.get("energy"))
        frames.append(
            OptimizationFrame(
                step=int(item.get("step", len(frames) + 1) or len(frames) + 1),
                coordinates=CoordinateSet3D(
                    positions=positions,
                    source=source or "rdkit",
                    method=method,
                    energy=energy,
                    metadata={"frame": True, **_derived_hydrogen_metadata(item)},
                ),
                energy=energy,
                converged=bool(item.get("converged", False)),
            )
        )
    return tuple(frames)
