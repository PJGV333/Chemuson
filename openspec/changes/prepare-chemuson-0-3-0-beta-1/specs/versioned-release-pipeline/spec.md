# Spec Delta — versioned release pipeline

## Purpose

Define repeatable SemVer releases whose application version, source tag, package metadata and update manifests agree, and whose builds are gated on the exact source commit before reaching an isolated distribution channel.

## ADDED Requirements

### Requirement: Application version has one canonical source
The application SHALL obtain its version from `src/chemuson/_version.py`; package metadata SHALL derive from that value. Release metadata and generated installers/manifests SHALL be verified against the version encoded in the release tag. Release CI MUST NOT rewrite source version files while building a tagged commit.

#### Scenario: Build a tagged version
- **WHEN** a valid release tag is processed
- **THEN** the canonical application version, dynamic package version, AppStream release entry, installer metadata, artifact names and updater manifest all identify the tag's version
- **AND** the checked-out source tree is unchanged by version synchronization during the build.

### Requirement: Release versions follow the project SemVer sequence
The project SHALL use `MAJOR.MINOR.PATCH[-prerelease]` without reusing published versions or tags. While major is zero, compatible fixes increment PATCH, additive feature lines increment MINOR, and breaking/API-generation changes also increment MINOR until the owner explicitly approves 1.0.0. Prerelease order SHALL be `-dev`, `-beta.N`, optional `-rc.N`, then the stable version, with monotonic prerelease numbers.

#### Scenario: Prepare a beta
- **WHEN** the owner prepares the next 0.3.0 beta
- **THEN** the canonical source is `0.3.0-beta.1`
- **AND** a subsequent beta uses `0.3.0-beta.2` rather than reusing or replacing beta.1.

#### Scenario: Promote or correct a release
- **WHEN** acceptance approves stable promotion
- **THEN** the stable tag is `v0.3.0` and no beta artifact is relabeled as stable
- **AND** a correction after stable `0.3.0` uses a new patch version such as `0.3.1`.

### Requirement: Release channels derive only from a protected SemVer tag
Release publication SHALL start from a `v*` Git tag whose strict SemVer version matches the canonical source version. Stable versions without prerelease identifiers SHALL map only to `stable`; supported `beta.N` and `rc.N` prereleases SHALL map only to the beta distribution channel and be marked prerelease. Development, unknown prerelease labels, mismatches, deleted refs, and malformed tags MUST fail closed. Manual publication MUST NOT accept independently selectable version and channel values.

#### Scenario: Beta tag routes only to beta
- **WHEN** `v0.3.0-beta.1` is pushed after review
- **THEN** its GitHub Release is prerelease and its Flatpak/AppImage update metadata targets beta
- **AND** stable manifests/remotes are untouched.

#### Scenario: Stable tag routes only to stable
- **WHEN** `v0.3.0` is pushed after explicit owner approval
- **THEN** its release is stable and its update metadata targets stable
- **AND** beta manifests/remotes are untouched.

#### Scenario: Invalid or mismatched release ref
- **WHEN** a tag is malformed, unsupported, already released, or disagrees with `_version.py`
- **THEN** preflight fails before package builds or publication.

### Requirement: Release tests and artifacts use one validated commit
A release SHALL run an authoritative bounded release-gate on the exact commit peeled from the tag before any platform build. Every build and publication job MUST check out and identify that same SHA, depend on the gate, and fail on any new/unrecognized test, static-check, version, or packaging error. Independent branch/PR CI results MUST NOT be treated as evidence for a different tag SHA.

#### Scenario: Same-SHA gate passes
- **WHEN** the release-gate verifies the tag commit and required focused tests
- **THEN** only dependent artifact jobs may proceed, each recording the validated SHA.

#### Scenario: Gate has a new failure or times out
- **WHEN** any required release check fails, exceeds its timeout, or its result is unavailable
- **THEN** all package and publication jobs are blocked.

### Requirement: Published version identifiers and artifacts are immutable
A published tag or release SHALL NOT be moved, deleted and reused, or have existing assets silently overwritten. Preflight SHALL fail closed if a release for that tag already exists or GitHub cannot determine that it does not exist. Repository policy SHALL restrict creation/update/deletion of `v*` tags to authorized maintainers. Corrections SHALL use a new higher version and new artifact names.

#### Scenario: Existing release collision
- **WHEN** preflight finds an existing GitHub Release for the requested tag
- **THEN** the workflow stops without replacing release assets or updating a channel.

### Requirement: Artifact integrity and provenance are auditable
Every release artifact SHALL have a SHA-256 checksum. The release payload SHALL identify its tag, canonical version, channel and validated source SHA. HMAC and Flatpak GPG signatures SHALL be reported accurately as optional unless required secrets are present; an unsigned artifact MUST NOT be described as cryptographically signed. Publication SHALL preserve separate beta/stable targets.

#### Scenario: Verify release payload
- **WHEN** package artifacts are assembled
- **THEN** checksums and provenance are generated from the assembled files and validated SHA
- **AND** channel manifests reference the same version and tag.

### Requirement: Rollback does not rewrite published history
A defective release SHALL be recovered by stopping/withdrawing distribution where possible and publishing a corrected higher version. The prior tag and release assets SHALL remain immutable; beta corrections advance the beta number and stable corrections advance PATCH.

#### Scenario: Defective beta
- **WHEN** a beta is rejected during acceptance
- **THEN** it is not promoted or relabeled stable
- **AND** any correction is published as the next beta version.

### Requirement: Baseline exceptions are exact and traceable
Known baseline exceptions SHALL be enumerated by exact test/node or narrowly described reproducer, linked to recorded evidence and a tracking owner/reference, and given a removal condition. The release gate MUST NOT use blanket ignores, wildcard deselection, or convert unknown failures to skips. A new failure or an application-level P0/P1 regression blocks release.

#### Scenario: Previously recorded test failure is observed
- **WHEN** a test matches a catalogued baseline exception exactly
- **THEN** its evidence and owner approval are recorded separately from passing gates
- **AND** no other failing test is accepted by that exception.

#### Scenario: Unknown failure or packaged-app crash
- **WHEN** a failure is not an exact catalog entry, or a packaged application has a reproducible crash/data-loss/security issue
- **THEN** beta/stable promotion is blocked until reviewed and resolved or explicitly reclassified with evidence.
