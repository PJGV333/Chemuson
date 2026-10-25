## ADDED Requirements

### Requirement: Release preflight path resolution is correct and regression-protected

The official tag preflight SHALL resolve its repository root to the root of the checked-out repository, and that resolution SHALL be verified by automated tests against the real Git metadata of a checkout.

#### Scenario: Preflight constant matches the Git checkout root
- **GIVEN** the release preflight script in any checkout of this repository
- **WHEN** the test suite compares its resolved repository root with `git rev-parse --show-toplevel`
- **THEN** both paths are identical
- **AND** Git commands issued by the preflight execute inside the checkout, never outside it

#### Scenario: Real tag checkout logic is exercised, not only helpers
- **GIVEN** a throwaway Git checkout containing an annotated tag whose commit equals the checked-out HEAD
- **WHEN** the test suite executes the preflight's tag/checkout verification logic and its full entry point under a simulated tag event
- **THEN** the annotated tag, HEAD and event SHA identity checks pass and the expected `tag`, `version`, `channel` and `sha` outputs are produced
- **AND** mismatched event SHAs, non-tag refs and missing tags still fail closed

### Requirement: The official publication preflight is rehearsed in GitHub Actions before the official tag is created

The preview workflow SHALL contain a non-publishing rehearsal job that executes the exact official release preflight code on a GitHub-hosted runner under a fully simulated tag event for the candidate version, and the overall preview result SHALL require that rehearsal to succeed.

#### Scenario: Rehearsal replays the official preflight faithfully
- **GIVEN** a candidate commit on a `release/**-prep` branch
- **WHEN** the preview workflow runs the rehearsal job
- **THEN** the runner performs a full checkout (`fetch-depth: 0`) of the exact candidate SHA
- **AND** it first fails closed if the candidate tag already exists on the remote
- **AND** it creates a temporary annotated tag for the candidate version only inside the runner, pointing at the exact candidate SHA
- **AND** it simulates the tag event context: ref name, ref type `tag`, exact event SHA, event payload file and dedicated output file
- **AND** it runs the official `packaging/release/validate_release_tag.py` unchanged, authorized by a read-only repository token for the GitHub release-absence lookup
- **AND** it verifies the validator exits with code 0 and emits the expected tag, version, channel and SHA outputs

#### Scenario: Rehearsal never publishes and always cleans up
- **WHEN** the rehearsal job finishes, regardless of outcome
- **THEN** the temporary local tag has been deleted and its absence is re-checked locally and on the remote
- **AND** no push, GitHub Release, remote or Pages operation was executed

#### Scenario: Preview is not publish-fit without a green rehearsal
- **WHEN** the preview build report is written
- **THEN** the rehearsal result is one of the statuses that must all be `success` for the preview to be complete
- **AND** any failure or inability to run the rehearsal faithfully marks the whole preview as incomplete

### Requirement: Bounded release gate checks are rehearsed under gate-equivalent runner conditions

The preview workflow SHALL rehearse the bounded verification steps of the official `release_gate` job after its preflight, using the same runner image, Python version, runtime dependencies and byte-for-byte identical commands, and the rehearsal SHALL be kept synchronized with the official workflow by static contracts.

#### Scenario: Gate rehearsal mirrors the official bounded commands
- **GIVEN** the official `release_gate` job definition
- **WHEN** the static contract compares the gate rehearsal job
- **THEN** every mirrored bounded step has an identical normalized command
- **AND** the rehearsal result gates the overall preview report exactly like the preflight rehearsal
