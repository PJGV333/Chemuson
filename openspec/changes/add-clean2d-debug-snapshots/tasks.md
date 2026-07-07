## 1. Snapshot schema and activation

- [x] 1.1 Identify the smallest Clean 2D integration point that can observe run context, candidates, diagnostics, and final decision without changing algorithm flow.
- [x] 1.2 Define the JSON snapshot data shape with schema identifier, version, Clean 2D mode, target, topology, coordinates, selection, candidates, final state, final reason, and optional metadata.
- [x] 1.3 Add opt-in activation through a documented environment variable.
- [x] 1.4 Add opt-in activation through an explicit debug parameter where Clean 2D is invoked in tests or developer-facing APIs.
- [x] 1.5 Add a test-only helper for scoped snapshot capture.

## 2. Snapshot capture and serialization

- [x] 2.1 Capture atom IDs, bond IDs, connectivity, and bond orders from existing molecule state.
- [x] 2.2 Capture initial coordinates before Clean 2D changes geometry.
- [x] 2.3 Capture final coordinates when final coordinates exist.
- [x] 2.4 Capture initial selection separately from target atom IDs or whole-structure targeting.
- [x] 2.5 Capture evaluated candidate sources without adding, removing, reordering, rescoring, or selecting candidates.
- [x] 2.6 Capture per-candidate `quality_diagnostic` when available.
- [x] 2.7 Capture final Clean 2D state and final stable reason when available.
- [x] 2.8 Ensure snapshot output uses only stable JSON-serializable primitives.

## 3. Test helpers

- [x] 3.1 Add helpers for writing Clean 2D debug snapshots to test-controlled paths.
- [x] 3.2 Add helpers for reading snapshot JSON in tests.
- [x] 3.3 Ensure helpers do not write snapshots outside explicit test/debug locations.

## 4. Tests

- [x] 4.1 Add a test proving normal Clean 2D execution does not capture snapshots by default.
- [x] 4.2 Add tests proving each opt-in path enables snapshot capture.
- [x] 4.3 Add tests asserting required schema/version fields and JSON serializability.
- [x] 4.4 Add tests asserting topology, initial coordinates, target, and selection are captured.
- [x] 4.5 Add tests asserting final coordinates are captured when available and represented safely when unavailable.
- [x] 4.6 Add tests asserting candidate sources and available `quality_diagnostic` data are captured.
- [x] 4.7 Add tests asserting final state and final stable reason are captured when available.
- [x] 4.8 Add snapshot-focused coverage for the known complex Clean 2D failure shapes without asserting algorithm, geometry, ranking, or selection fixes.

## 5. Verification

- [x] 5.1 Run the targeted Clean 2D snapshot tests.
- [x] 5.2 Run relevant existing Clean 2D tests to confirm behavior is unchanged when snapshots are disabled.
- [x] 5.3 Confirm the four known algorithmic failures are not fixed or reclassified by this change unless separately addressed by a future proposal.
