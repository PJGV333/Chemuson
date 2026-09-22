# Design: Tighten Campaign 3 global placement safety

## Explicit global budget

Compute `displacement_budget` in `engine.py` from the requested target bond length and the topology block counts. The budget is finite, recorded in candidate metadata, and applies only to `global_block_placement`. A candidate above it is rejected with a stable reason. The existing local safety helper is not globally relaxed; the only allowed difference is an explicit global displacement budget.

## Explicit safety record

After placement, evaluate finite coordinates, selected-graph invariants, stereo signature, new crossings, collision distance, ring degeneracy, bounding-box ratio, bond-length ratio, and the global displacement budget independently. Store a JSON-safe `hard_gate_checks` map, `max_displacement`, `displacement_budget`, `bounding_box_ratio`, and `rejection_reason` in metadata.

## Complex-preserve competition

Generate global, scaffold, and block-unwrap candidates before selecting. Keep only candidates already accepted by their constructors, compare them using existing quality classification and visual-quality metrics, and break exact ties by source name. Annotate the selected candidate with `competing_safe_sources` and `selected_source`. If no safe alternative exists, preserve the existing `preserve-only` fallback.

## Evidence

Extend the Campaign 3 corrective evidence record with the safety fields and competing safe sources. Do not modify prior Campaign 3 baseline files.
