# Spec Delta — isolated Actions preview builds

## Purpose

Generate downloadable Windows and Linux/Flatpak test packages from an exact preparation-branch commit, without publishing a release or changing any public update channel.

## ADDED Requirements

### Requirement: Preview builds are restricted to preparation branches and exact source identity
The preview workflow SHALL support `workflow_dispatch` on a selected preparation branch and MAY run on pushes matching `release/**-prep`. Before any build, it MUST require a branch ref matching that preparation pattern, resolve the event's full Git SHA and canonical `_version.py` version, and make every platform job check out and verify that exact SHA. It MUST NOT create, move, delete, or push tags.

#### Scenario: Owner dispatches a preview
- **WHEN** the owner selects `release/v0.3.0-beta.1-prep` in Actions
- **THEN** the workflow records the selected branch, full event SHA, canonical version and UTC build time
- **AND** every package job verifies its checkout equals the recorded SHA and uses the same version without rewriting source files.

#### Scenario: A non-preparation ref is selected
- **WHEN** dispatch or a push does not identify a branch matching `release/**-prep`
- **THEN** the preparation gate fails before platform builds.

### Requirement: Preview artifacts are complete, distinct, and verifiable
A successful preview run SHALL publish four separate Actions artifacts named `chemuson-preview-windows-portable`, `chemuson-preview-windows-installer`, `chemuson-preview-linux-appimage`, and `chemuson-preview-linux-flatpak`. Each group MUST include its nonempty package, SHA-256 checksum(s), and a provenance manifest recording version, `Build type: preview`, exact Git SHA, source branch, UTC build time, operating system, successful build status and `Publication: false`. Package names MUST be distinguishable from official releases, and checks MUST fail before upload if expected output is missing or empty. A summary MUST report every build job status; any failed/missing platform means the overall run is not a successful complete preview.

#### Scenario: All preview packages build
- **WHEN** Windows portable, Windows Inno installer, Linux portable executable and Flatpak bundle are produced from the validated commit
- **THEN** each is uploaded under its dedicated artifact group with checksums and matching provenance.

#### Scenario: One package is missing or a build fails
- **WHEN** any expected file is absent/empty or a build job fails
- **THEN** that artifact is not reported as built, the run summary records the failure, and the overall workflow conclusion is failure rather than partial success.

### Requirement: Preview workflow cannot publish or affect public update channels
The preview workflow SHALL use no more than `contents: read`, SHALL NOT receive release/signing/publishing secrets, and SHALL NOT invoke GitHub Release actions/commands, tag operations, `git push`, `gh-pages`, public Flatpak remote publication, updater-channel manifest generation, or public beta/stable URLs. The preview Linux portable build MUST omit `.updateinfo`, `.update.json`, and `.zsync` publication metadata. The Flatpak preview MUST build a local bundle without configuring a public Chemuson remote URL, signing key, or `gh-pages` deployment. Preview artifacts SHALL be available only as Actions run artifacts and MUST NOT become visible as an application update.

#### Scenario: Preview artifacts are downloaded
- **WHEN** a user downloads the preview AppImage-named executable or Flatpak from Actions
- **THEN** no public updater manifest, remote repository, stable/beta channel or GitHub Release has been modified.

#### Scenario: Static isolation contract is evaluated
- **WHEN** the preview workflow contract tests run
- **THEN** they assert triggers, permissions, SHA/version identity, artifact groups, existence checks and absence of publication operations
- **AND** their PASS result is described only as static isolation evidence, not as proof that platform packages compiled.

### Requirement: Preview use and execution status are documented
Maintainer documentation SHALL explain GitHub UI dispatch, artifact download and manual testing, plus authorized GitHub CLI/API execution. A real preview run MAY be initiated only after the workflow is available on GitHub and its isolation has been reviewed; lack of authorization MUST leave real artifacts explicitly `NOT BUILT` without escalating credentials or inventing another publication route.

#### Scenario: Workflow is not yet on GitHub or agent is unauthenticated
- **WHEN** no authorized Actions run can be started
- **THEN** the report provides the owner steps and marks preview infrastructure and real artifact state separately.
