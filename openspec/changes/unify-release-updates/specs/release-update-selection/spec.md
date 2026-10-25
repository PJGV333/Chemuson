# Release update selection Specification

## Purpose

Define one SemVer-based selection rule for public GitHub Releases while retaining distribution-channel contracts already used by installed ChemUSON clients.

## ADDED Requirements

### Requirement: Official release identity and precedence SHALL be SemVer-based

The updater SHALL identify an official candidate from its exact `v<SemVer>` tag and SHALL order candidates by SemVer precedence, not by GitHub publication date or the GitHub `prerelease` boolean. Invalid tags or malformed SemVer identifiers SHALL NOT be offered.

#### Scenario: Publication date differs from semantic order

- **GIVEN** two valid stable release tags where the older SemVer tag has a later `published_at`
- **WHEN** the provider selects the stable release
- **THEN** it SHALL select the greater SemVer version

#### Scenario: GitHub prerelease metadata contradicts SemVer

- **GIVEN** a `vX.Y.Z-beta.N` release marked `prerelease=false` and a `vX.Y.Z` release marked `prerelease=true`
- **WHEN** the provider filters by channel
- **THEN** stable SHALL exclude the beta and may accept the final version; beta SHALL recognize the beta by its SemVer tag

### Requirement: Existing channel preferences SHALL remain compatible

The updater SHALL preserve the persisted `stable` and `beta` preference. Stable SHALL accept only final SemVer versions. Beta SHALL accept newer beta/RC prereleases and newer final versions so an opted-in beta installation can graduate to the reviewed version. Channel selection SHALL NOT be inferred solely from GitHub metadata.

#### Scenario: Beta installation moves to a newer beta

- **GIVEN** an installed beta at `0.3.0-beta.2` and an authorized newer beta `0.3.0-beta.3`
- **WHEN** the user has selected beta updates
- **THEN** the updater SHALL offer beta.3 if its platform asset is available and integrity checks pass

#### Scenario: Beta installation moves to the reviewed final

- **GIVEN** an installed beta at `0.3.0-beta.3` and a final `0.3.0`
- **WHEN** the user has selected beta updates
- **THEN** the updater SHALL offer `0.3.0` if its platform asset is available and integrity checks pass

### Requirement: Updates SHALL be monotonic and fail closed

The updater SHALL NOT offer an equal or older version. It SHALL preserve platform-specific asset selection and SHALL refuse to apply a download with a missing required checksum or mismatched SHA-256. Network, malformed metadata, and missing compatible assets SHALL NOT cause an older package to be offered.

#### Scenario: No newer compatible release exists

- **GIVEN** the highest compatible release is equal to or older than the installed version
- **WHEN** the updater checks the catalog
- **THEN** it SHALL report no available update

#### Scenario: Latest release has no compatible platform asset

- **GIVEN** the newest compatible SemVer release has no asset for the current platform/flavor
- **WHEN** the updater checks the release
- **THEN** it SHALL report no update rather than select an older package

### Requirement: One protected release workflow SHALL serve all SemVer release notes

The official tag-only workflow SHALL build and publish from one validated version/SHA and SHALL load notes from the version-specific `docs/release/<version>-public.md` path. Preview and diagnostic builds SHALL remain non-publishing.

#### Scenario: Versioned release notes

- **GIVEN** a protected tag and matching canonical version
- **WHEN** the official workflow prepares the GitHub Release
- **THEN** it SHALL use the public notes file corresponding to that SemVer without a version-specific workflow edit
- **AND** it SHALL preserve existing gates, permissions, package checks, checksums, provenance and platform update mechanisms

#### Scenario: Diagnostic build on the primary or feature branch

- **GIVEN** a clean checkout whose exact source SHA and canonical version match
- **WHEN** a diagnostic build is requested from `main`, `feature/*` or `fix/*`
- **THEN** diagnostic identity validation SHALL accept the branch without publishing a tag, Release or update metadata
- **AND** release-preparation branches SHALL continue to use the separate Preview path
