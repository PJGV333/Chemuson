# PACKAGED-RDKIT-WORKER-001 — RDKit isolated worker in distributed executables

## Purpose

Restore RDKit-backed descriptors and other existing isolated operations in PyInstaller Windows portable, Linux portable, and AppImage packages without requiring an external Python installation or changing chemical algorithms.

## ADDED Requirements

### Requirement: A frozen application can launch its RDKit worker without bootstrapping the GUI
A PyInstaller executable MUST dispatch a private worker invocation before normal CLI/GUI startup. The worker MUST execute in a separate process using the bundled Python runtime, MUST NOT open a GUI or import RDKit into the parent application process, and MUST remain isolated from native RDKit crashes. Frozen worker transport MUST work with Windows `console=False` and MUST preserve the existing JSON request/response schema and caller timeouts. The source-Python worker path remains supported.

#### Scenario: Frozen worker dispatch
- **WHEN** a frozen ChemUSON process requests an isolated RDKit operation
- **THEN** it starts the embedded worker mode of the same executable rather than passing a `.py` worker script to that executable or requiring system Python
- **AND** the child handles the request without starting the main window or recursively launching another worker.

#### Scenario: Windows without console streams
- **WHEN** the parent is the Windows GUI executable built with `console=False`
- **THEN** worker communication succeeds without relying on child `sys.stdin`/`sys.stdout` console handles
- **AND** a worker timeout terminates and waits for the child and leaves no temporary request/response files or orphaned process.

### Requirement: Packaged RDKit and native extensions are proven from the bundle
Both preview and official release packaging MUST run a fail-closed test against the actual frozen executable. It MUST separately prove that RDKit imports and that required compiled RDKit extensions load from the executable's own extracted bundle, not from an external Python/RDKit installation. The test MUST calculate known ethanol descriptors (logP, TPSA, HBD, HBA) and exercise bounded isolated 3D conformer and SMILES operations. The packaged smoke MUST NOT skip when RDKit or a worker is unavailable.

#### Scenario: Windows and Linux package smoke
- **WHEN** the Windows portable or Linux PyInstaller executable is built for preview or official release
- **THEN** the packaging gate invokes that exact executable and validates a structured smoke report, bundled RDKit extension paths, known descriptor values, 3D coordinates, and canonical SMILES
- **AND** any import, startup, timeout, malformed response, missing extension, or descriptor mismatch fails the job rather than being skipped.

#### Scenario: AppImage package smoke
- **WHEN** the AppImage Type 2 artifact is extracted for validation
- **THEN** the same RDKit smoke runs against the executable inside the extracted AppDir in addition to the standalone Linux executable test
- **AND** AppImage validation fails if either frozen process cannot use bundled RDKit.

### Requirement: Worker failures have accurate, nonfatal diagnostics
Worker import/native-extension failures, worker start/exit failures, timeouts, invalid JSON/payload, and chemical calculation errors MUST remain distinguishable. The Properties pane MUST label RDKit unavailable only for an actual RDKit import failure; a worker/protocol/calculation error MUST identify its own cause and MUST leave the GUI responsive with partial properties intact.

#### Scenario: Worker failure while updating Properties
- **WHEN** an isolated descriptor request fails
- **THEN** the UI reports the specific failure category and retains available formula, mass, and estimated spectra
- **AND** the application does not block, close unexpectedly, or misreport every failure as a missing RDKit installation.

## Compatibility

The existing worker JSON payload fields, descriptor definitions, timeout values, chemical algorithms, Clean2D behavior, package/import identities, and `.cmsn` persistence remain unchanged. The known Qt teardown/SIGSEGV issue remains a separate unresolved debt.

## Acceptance

Source-level RDKit import and unit tests are necessary but not sufficient. `WINDOWS DESCRIPTORS` and `APPIMAGE DESCRIPTORS` remain NOT VERIFIED until the corrected frozen artifacts pass their real Build Preview smoke gates. Beta acceptance remains pending the owner's manual retest of the new Windows/Linux packages.
