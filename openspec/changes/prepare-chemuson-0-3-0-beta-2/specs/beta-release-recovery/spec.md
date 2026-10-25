## ADDED Requirements

### Requirement: A blocked beta release recovers with a new immutable candidate identity

When an official beta workflow did not run for an existing tag, the project SHALL preserve that tag and its target commit, prepare any recovery under the next monotonic beta number, and SHALL NOT create the successor tag until the candidate SHA and version change have been reported to and confirmed by the owner.

#### Scenario: Recover a tag whose push workflow was suppressed
- **WHEN** a beta tag exists but no official release workflow completed
- **THEN** the existing tag object and target commit remain unchanged
- **AND** the successor candidate is based on the designated preparation branch and uses the next beta number
- **AND** no stable release or manual asset publication is substituted for the tag workflow
- **AND** the successor tag is withheld pending owner confirmation of the exact candidate SHA.

### Requirement: Release notes match the candidate version

The official release workflow SHALL consume public notes for the exact candidate version being built, and those notes SHALL distinguish verified evidence from unverified or historical acceptance observations.

#### Scenario: Prepare version-specific beta notes
- **WHEN** a beta candidate is prepared for the official tag-only workflow
- **THEN** the canonical application version and current AppStream version agree
- **AND** the release workflow's `body_path` resolves to the candidate's public notes
- **AND** the notes do not attribute prior package observations to candidate artifacts without matching SHA or artifact evidence.
