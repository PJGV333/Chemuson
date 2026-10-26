## ADDED Requirements

### Requirement: Official project logo is preserved and used consistently

The preparation SHALL use only the owner-approved ChemUSON PNG as visual source. Generated variants SHALL preserve its artwork, proportions, colors, and lettering; no substitute logo or lossy vector tracing is allowed. The logo identifies the ChemUSON project and SHALL NOT be described as an official University of Sonora seal.

#### Scenario: Approved master and variants are verifiable
- **GIVEN** the approved source PNG
- **WHEN** it is integrated into the repository and packaged resources
- **THEN** the master retains the approved image bytes and expected dimensions
- **AND** raster/ICO sizes and the Linux SVG container are derived from that master
- **AND** the SVG embeds the source raster rather than tracing or changing the design

#### Scenario: Branded desktop packages display the same project identity
- **GIVEN** a ChemUSON desktop package on Windows or Linux
- **WHEN** its window/taskbar/launcher icon is rendered
- **THEN** it uses the approved logo through the existing stable application ID
- **AND** the Qt resource is available in both light and dark themes
- **AND** Flatpak/AppStream, AppImage, Windows executable and installer reference matching variants

### Requirement: Local build cleanup preserves historical evidence and redirects outputs

Cleanup SHALL remove only inspected, ignored generated build artifacts that are no longer needed. Tracked diagnostic provenance/checksum manifests SHALL remain unchanged. Local builds SHALL support `CHEMUSON_OUTPUT_DIR` where a corresponding builder currently writes into the checkout, while official workflow defaults and package names remain compatible.

#### Scenario: Historical diagnostic manifests are preserved
- **GIVEN** `dist-appimage-diagnostic/checksums.sha256` and `diagnostic-provenance.json` are tracked
- **WHEN** obsolete local build outputs are cleaned
- **THEN** both tracked files remain byte-identical
- **AND** source, tests, scripts, docs and historical tags/releases are untouched

#### Scenario: Local Windows installer output is relocated
- **GIVEN** `CHEMUSON_OUTPUT_DIR` is set before Inno Setup runs
- **WHEN** the installer is built
- **THEN** its output goes under that external directory
- **AND** with no override the existing workflow path `dist-installer` remains valid

### Requirement: Release documentation accurately describes beta.3's one-off GitHub classification

Documentation SHALL state that the owner manually set beta.3's GitHub Release `prerelease` metadata to false as an exceptional decision. It SHALL NOT establish that change as a general beta policy or allow it to override SemVer-based update selection or user settings.

#### Scenario: Beta.3 appears as GitHub Latest
- **GIVEN** tag `v0.3.0-beta.3` and the API metadata `prerelease=false`
- **WHEN** the README and versioning policy explain GitHub's **Latest** label
- **THEN** they describe the owner-approved exception and current display accurately
- **AND** explain that beta versions continue to be identified and ordered by SemVer
- **AND** future beta releases retain the official workflow's normal prerelease behavior

#### Scenario: Update preferences and integrity remain protected
- **GIVEN** an installation on stable or beta with persisted settings
- **WHEN** update candidates are selected or verified
- **THEN** stable excludes SemVer prereleases, beta may receive newer beta/RC or final versions, downgrades are blocked, and integrity checks remain active
- **AND** no persisted preference is changed by this campaign

### Requirement: Beta.4 metadata and preview packages identify one exact source SHA

The canonical version SHALL be `0.3.0-beta.4`, with matching AppStream release metadata and Spanish-then-English public notes limited to branding and maintenance changes. Preview builds SHALL be non-publishing and use one exact candidate SHA for four package groups.

#### Scenario: Version and release notes match
- **GIVEN** the beta.4 preparation commit
- **WHEN** version metadata is validated
- **THEN** `_version.py`, AppStream and packaging identity all agree on `0.3.0-beta.4`
- **AND** the notes do not present beta.3 scientific capabilities as new beta.4 features

#### Scenario: Preview is complete but not published
- **GIVEN** a candidate commit on `release/**-prep`
- **WHEN** the official Build Preview completes
- **THEN** Windows portable, Windows installer, Linux AppImage Type 2 and Linux Flatpak share its exact SHA and matching checksums/provenance
- **AND** no tag, Release or public update channel is created
- **AND** tag creation remains pending explicit owner acceptance

### Requirement: AppImageUpdate release tracks use selectors supported by AppImageUpdate

Official AppImage metadata SHALL use `latest` for the stable track and `latest-all` for the beta track. It SHALL NOT use `prerelease` as a track value: AppImageUpdate treats unrecognized values as literal GitHub tags. `latest-all` preserves the established beta-to-final graduation behavior and includes releases independently of the GitHub `prerelease` flag. Existing published artifacts SHALL NOT be rewritten.

#### Scenario: Release track mapping is validated
- **GIVEN** a stable or beta AppImage release build
- **WHEN** the builder and release validators produce/check update information
- **THEN** stable uses `gh-releases-zsync|owner|repo|latest|asset.zsync`
- **AND** beta uses `gh-releases-zsync|owner|repo|latest-all|asset.zsync`
- **AND** the update JSON, sidecar, and embedded AppImage value match exactly

#### Scenario: Historical beta.3 remains immutable
- **GIVEN** beta.3 embeds the unsupported literal `prerelease` track
- **WHEN** beta.4 fixes the builder contract
- **THEN** beta.3 assets and tag remain unchanged
- **AND** documentation tells beta.3 AppImage users to install the new package manually
- **AND** beta.4 and later use the supported selector without altering user preferences.

### Requirement: The non-Qt startup harness models application icon initialization

The AST-isolated `run_app()` tests SHALL provide mocks for each runtime name used by the startup flow, including the icon path resolver, `QIcon`, application/window icon assignment, and crash-log reporting. Success tests SHALL keep the startup-order assertions strict and verify the same non-null icon reaches both the application and main window. A missing icon SHALL be reported through the existing startup error contract without constructing the main window. Production startup and error-handling behavior SHALL NOT be changed to accommodate the test harness.

#### Scenario: Startup loads and assigns the bundled logo
- **GIVEN** a temporary icon resource exists and the fake `QIcon` reports non-null
- **WHEN** `run_app()` starts in the AST-isolated harness
- **THEN** the icon path is resolved, loaded, and assigned to QApplication before window construction
- **AND** the identical icon object is assigned to the main window before autosave checks/show/event loop
- **AND** the existing startup event order remains asserted

#### Scenario: Missing icon follows existing startup error handling
- **GIVEN** the resolved icon path does not exist
- **WHEN** `run_app()` reaches logo initialization
- **THEN** it raises the expected `FileNotFoundError` internally and does not construct QIcon or the main window
- **AND** it writes the crash log and shows the existing critical startup message when QApplication exists
- **AND** existing generic startup failure tests and their reporting contracts remain intact

### Requirement: Author listing recognizes the additional author without access changes

The Spanish and English project author lists SHALL add Jesús David Albino-Espíndola after the existing authors, with no corresponding-author asterisk. The existing author names, order and asterisks SHALL remain intact.

#### Scenario: Documented authorship only
- **GIVEN** the README author lists
- **WHEN** the new name is incorporated
- **THEN** both language versions contain the same author order and corresponding-author markers
- **AND** no GitHub collaborator, member or permission is invited or changed

### Requirement: Flatpak host toolchain preparation is bounded, cached and origin-verified

The Flatpak packaging jobs SHALL prepare their host toolchain through a shared, documented
procedure that selects a healthy official Ubuntu package origin, bounds every network phase
with an explicit wall-clock limit, reuses a package cache only after apt verifies each
archive against the current signed index, and records the measured preparation and packaging
times. Increasing a job `timeout-minutes` SHALL NOT be used as the fix, and third-party
repositories or unsigned package sources SHALL NOT be introduced.

#### Scenario: Degraded assigned mirror no longer consumes the job
- **GIVEN** a GitHub-hosted runner whose priority-1 Ubuntu mirror is throttled
- **WHEN** the Flatpak job prepares `flatpak`, `flatpak-builder` and `ostree`
- **THEN** the preparation probes the official origins listed in the lock file, uses the fastest healthy one and fails over to the next if it fails
- **AND** each phase ends within its configured limit or the job fails closed with the phase and limit in the error

#### Scenario: Warm package cache cannot introduce stale or tampered archives
- **GIVEN** a previous run published its downloaded `.deb` archives to the Actions cache
- **WHEN** a later Flatpak job restores that cache and installs
- **THEN** apt installs with `--no-download` so every archive is verified against the signed index hashes
- **AND** any missing, outdated or mismatched archive is re-fetched from an official origin instead of being installed

#### Scenario: Packaging phases report real timings and stay bounded
- **GIVEN** the Flatpak bundle build downloads the pinned runtime and wheel sources
- **WHEN** `flatpak-builder`, `build-update-repo` or `build-bundle` runs
- **THEN** each phase is limited by a wall-clock timeout and its elapsed seconds are reported in the step summary
- **AND** an exhausted phase reports the phase name and its limit rather than letting the job be cancelled without a cause

#### Scenario: Scientific packaging and release safety stay unchanged
- **GIVEN** the preparation repair touches workflows and host tooling only
- **WHEN** the packaging contracts are checked
- **THEN** the Flatpak manifest keeps its runtime pin, the five sha256-pinned wheel sources, the packaged ChemName smoke, desktop/metainfo/icon installation and sandbox permissions
- **AND** exact-SHA verification, checksums, provenance, GPG/OSTree metadata checks, Pages verification and fail-closed preview reporting remain active
