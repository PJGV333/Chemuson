# Proposal: Complete Clean2D Campaign 4 target improvement

## Why

The archived Campaign 4 implementation emits a topology-derived `rigid_multiring_layout` candidate, but its evidence does not demonstrate the required target-family promotion: none of the twelve fixtures selects that candidate. The candidate metadata also collapses hard-gate safety with post-ranking acceptance and stores gate names without boolean values.

This corrective change keeps the Campaign 4 implementation history intact and makes promotion evidence truthful. It adds a real geometry improvement for a reusable rigid-system family, with spiro as the preferred target, while preserving safe fallback for bridged systems and leaving Campaign 3 global placement ownership unchanged.

## Scope

### Allowed

- `src/chemuson/clean2d/complex_policy.py`
- `src/chemuson/clean2d/engine.py`
- `src/chemuson/clean2d/__init__.py` only if a public contract must be exported
- Campaign 4 corrective tests and topology-built fixtures
- This OpenSpec and its deterministic evidence

### Explicitly out of scope

- Campaign 5 or any later campaign
- Source-priority or ranking bias introduced solely to select Campaign 4
- Relaxing hard gates or quality thresholds
- Branch routing, macrocycle layout, or global block-placement redesign
- GUI, persistence, `architecture/modules.yml`, `src/architecture`, and historical failures
- Fixture-name or molecule-name routing in production

## Promotion contract

The corrective campaign remains active until at least one legitimate target satisfies all of:

1. `rigid_multiring_layout` passes its hard gates.
2. Its geometry produces measurable target-family improvement over the declared baseline.
3. The general engine accepts it through normal evaluation and ranking.
4. `run_clean2d_engine(...).selected.source == "rigid_multiring_layout"`.
5. Selection is caused by geometry/metrics, not source-priority bias.

If these conditions cannot be demonstrated, the OpenSpec must remain unarchived.

## Baseline

The previous implementation commit is `23b7046`; the archived Campaign 4 state is `9a60968`. The corrective evidence compares baseline `9a60968` with the current implementation and retains `346b229` as the historical Campaign 3 safety reference.
