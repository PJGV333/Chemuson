# Add Clean 2D Debug Snapshots

## Summary

Add an opt-in Clean 2D debug snapshot system that captures a stable, serializable JSON record of a Clean 2D run for reproducing complex algorithmic bugs.

## Motivation

Clean 2D now has stable states, stable reasons, and `quality_diagnostic` reporting from the previous `document-clean2d-decision-contract` and `unify-clean2d-quality-reporting` changes. However, complex failures are still hard to reproduce because the information needed to compare a run across commits is spread across geometry, candidate generation, candidate evaluation, final decision state, and selection context.

The broad Clean 2D suite exposed algorithmic failures that should not be fixed by this change:

- `test_naphthalene_fused_system_does_not_collapse`
- `test_tetrandrine_like_hierarchical_blocks_do_not_select_local_graph`
- `test_smart_clean_repairs_distorted_structure_before_proposing`
- `test_complex_engine_does_not_call_global_redraw_candidates`

This change provides debugging infrastructure so future work can capture reproducible evidence for those cases without altering Clean 2D behavior.

## Capabilities

- `clean-2d-debug-snapshots`: capture opt-in JSON snapshots for Clean 2D runs, including molecule identity data, geometry inputs and outputs, candidate provenance, candidate diagnostics, and final decision state.

## Scope

- Define a stable JSON snapshot schema for Clean 2D debug snapshots.
- Capture snapshot schema/version information.
- Capture Clean 2D mode.
- Capture target atom IDs or whole-structure intent.
- Capture initial coordinates.
- Capture final coordinates when a final geometry exists.
- Capture atom IDs, bond IDs, connectivity, and bond orders.
- Capture initial selection.
- Capture evaluated candidate sources.
- Capture `quality_diagnostic` per candidate when available.
- Capture final state.
- Capture final stable reason when available.
- Allow optional internal metadata for troubleshooting.
- Make snapshot capture opt-in through an environment variable, an explicit debug parameter, or test-only helpers.
- Define helpers for writing and reading snapshots in tests.
- Add snapshot-focused tests that do not change Clean 2D algorithms.

## Out of Scope

- Implementing a new Clean 2D algorithm.
- Changing expected geometry.
- Changing candidate ranking or candidate selection.
- Performing a global refactor of `engine.py`.
- Changing RDKit or CoordGen integration.
- Changing MolGraph or canvas behavior.
- Fixing the four known algorithmic failures in this change.
- Enabling snapshots by default in normal use.

## Success Criteria

- Clean 2D behavior remains unchanged when snapshots are not enabled.
- Snapshot capture can be enabled explicitly in tests or debugging sessions.
- Snapshot JSON is stable, serializable, and suitable for storing as test fixtures.
- Tests prove that snapshot capture records the required fields without altering algorithmic decisions.
