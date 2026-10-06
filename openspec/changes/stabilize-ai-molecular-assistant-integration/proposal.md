# Stabilize AI Molecular Assistant Integration

## Why

The Phase 4.5 changes are implemented across separate OpenSpecs, but integration hygiene remains incomplete: oversized baseline logs are tracked, runtime requirements omit Pillow, development dependencies are not separated, Clean2D evaluation reports insufficient candidate provenance, and the evaluator implicitly reads `OPENAI_API_KEY` even for local endpoints. The repository also retains an obsolete Phase 5 design-block section after the user selected whole-molecule transformation.

## What Changes

- Consolidate the useful campaign history in `docs/history/CAMPAIGNS.md`, remove the obsolete Phase 5 block from `AGENT_REPORT.md`, and delete redundant tracked baseline-output logs while retaining concise baseline/validation records.
- Align `requirements.txt` with the runtime package dependencies and add a documented development dependency path.
- Extend the evaluator report with selected/rejected candidate summaries and selected source/outcome/reason/JSON-safe score, without changing Clean2D.
- Make API-key environment lookup opt-in; the default CLI invocation never reads `OPENAI_API_KEY`.
- Record bounded, shard-based verification; do not run a known 19-minute suite as one command.

## Boundaries

No Clean2D production code, chemistry behavior, provider payload behavior other than the explicit key source, or external dependencies are changed. Open Babel remains an optional system executable.
