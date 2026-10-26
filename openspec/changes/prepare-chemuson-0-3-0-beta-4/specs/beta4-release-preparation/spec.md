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

### Requirement: Author listing recognizes the additional author without access changes

The Spanish and English project author lists SHALL add Jesús David Albino-Espíndola after the existing authors, with no corresponding-author asterisk. The existing author names, order and asterisks SHALL remain intact.

#### Scenario: Documented authorship only
- **GIVEN** the README author lists
- **WHEN** the new name is incorporated
- **THEN** both language versions contain the same author order and corresponding-author markers
- **AND** no GitHub collaborator, member or permission is invited or changed
