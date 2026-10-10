# Tasks — ChemUSON 0.3.0-beta.1 preparation

## 1. Scope, baseline and contracts
- [x] 1.1 Fetch origin, record clean `origin/main` SHA, create the prep branch from it, and verify ancestry/status.
- [x] 1.2 Audit published GitHub Releases/tags, live beta/stable Flatpak endpoints, canonical version sources, existing CI and package builders; record actual facts in `baseline.md`.
- [x] 1.3 Define strict tag/channel policy, same-SHA gate, baseline exception policy, manual acceptance contract, and preview-build isolation; validate this OpenSpec before implementation.

## 2. Version policy and official release preflight
- [x] 2.1 Prepare `0.3.0-beta.1` with `_version.py` canonical and AppStream synchronized via the existing script; retain prior release entries and verify `pyproject` remains dynamic.
- [x] 2.2 Add a standard-library release-tag validator for strict stable/beta/rc grammar, canonical/AppStream equality, exact checked-out/tag/GitHub SHA, channel derivation and fail-closed existing-release lookup; cover valid, invalid, mismatch and API-error cases.
- [x] 2.3 Add package artifact/provenance validation and ensure channel manifests cannot pair beta versions with stable or vice versa; verify SHA-256 covers release payload/provenance.
- [x] 2.4 Make Windows installer compilation fail without `CHEMUSON_VERSION`; validate the existing Inno smoke still provides its explicit CI version.

## 3. Official CI and channel safeguards
- [x] 3.1 Remove separately selectable `workflow_dispatch` version/channel inputs; make release tag-only and fail closed on malformed/deleted/mismatched/already-published refs.
- [x] 3.2 Add bounded same-SHA release gate; require every Windows/Linux/Flatpak build and publication job to depend on it and check out its validated SHA.
- [x] 3.3 Narrow `GITHUB_TOKEN` permissions, prevent asset-name collisions/overwrites, preserve separate beta/stable Flatpak targets, and ensure stable/beta/rc channel mapping is exact.
- [x] 3.4 Review `test.yml` independently; document why its branch/PR result is not substituted for the authoritative same-tag-SHA release gate. Do not weaken or disable existing tests.
- [x] 3.5 Add exact Known Baseline Exceptions register and schema test; release checks must not deselect broad patterns, skip unknowns, or treat a suite crash as green.

## 4. Release documents and manual acceptance
- [x] 4.1 Write `docs/release/VERSIONING_POLICY.md` covering SemVer while MAJOR=0, dev→beta→optional RC→stable, tags/channels, hotfixes, ownership, immutability and recovery.
- [x] 4.2 Prepare draft notes in `docs/release/0.3.0-beta.1.md` with candidate scope, limitations, experimental AI, Qt debt, Clean2D scope and accurate Linux artifact label; avoid claiming unverified acceptance.
- [x] 4.3 Create `docs/release/manual-acceptance-0.3.0.md` with at least 50 stable-ID cases across startup/close, drawing, persistence/import, export, Clean2D, Assistant, UI and distribution; every case has steps, data, expected result and PASS/FAIL/BLOCKED/NOT TESTED fields.
- [x] 4.4 State P0/P1/P2/P3 triage and beta/stable promotion criteria, including package installed-app crash/data-loss blocks and required owner sign-off.
- [x] 4.5 Update README and the obsolete hotfix guide so they describe tag-based, single-version publication and do not instruct manual version/channel dispatch.

## 5. Actions preview builds (addendum)
- [x] 5.1 Add the preview-build capability and safety requirements to this OpenSpec before implementation.
- [x] 5.2 Add explicit preview modes to shared Linux builders so they validate the canonical version, produce distinguishable artifacts, and emit no public updater metadata or remote Flatpak publication configuration.
- [x] 5.3 Add `.github/workflows/build-preview.yml` for manual dispatch and only `release/**-prep` pushes; use `contents: read`, the exact event SHA, four separately named artifact groups, per-artifact checksums/provenance, and fail-closed file-existence checks.
- [x] 5.4 Add static workflow-contract tests for isolation, minimal permissions, exact SHA/version, artifact names/existence, and separation from `release.yml`; distinguish those from real package build evidence.
- [x] 5.5 Document GitHub UI and authorized `gh`/API invocation, artifact download/manual testing, no-publication guarantees, and the rule not to improvise with elevated credentials.
- [x] 5.6 Apply the execution gate: workflow is not yet on GitHub and `gh auth status` is unauthenticated, so no run was started; all four artifacts are recorded NOT BUILT and owner instructions are documented.

## 6. Bounded validation and delivery
- [x] 6.1 Run strict OpenSpec for this change, architecture, version/release/updater/packaging tests, compileall, focused Ruff and `git diff --check`; cap every automated test command at 10 minutes and do not rerun the monolithic suite.
- [x] 6.2 Validate AppStream XML, Flatpak YAML, release/preview workflow policy and manifest/checksum/provenance smoke with available tools; explicitly report unavailable Windows/Flatpak-builder/AppImage Type 2 validation.
- [x] 6.3 Verify no tag/release/gh-pages/channel publication occurred and no Clean2D/chemistry/`.cmsn` changes exist.
- [x] 6.4 Commit focused changes on this prep branch and push only if normal authentication is available; verify remote branch SHA. Never merge/rebase/force-push or modify `main`, tags, GitHub Releases or `gh-pages`.

## 7. Genuine AppImage Type 2 addendum
- [x] 7.1 Add a reproducible AppDir with `AppRun`, validated desktop entry, ChemUSON SVG icon, AppStream metainfo and the existing PyInstaller binary/resources.
- [x] 7.2 Pin the official AppImageKit appimagetool asset by upstream URL/asset/version and SHA-256; verify the tool before execution and document why linuxdeploy is unnecessary.
- [x] 7.3 Replace executable-copy/rename behavior with appimagetool Type 2 generation in the shared Linux builder; pin preview and official Linux build jobs to a compatible runner.
- [x] 7.4 Validate ELF/x86_64 + `AI\\x02`, extract with `--appimage-extract` without FUSE, validate AppDir/desktop/icon/AppStream, inspect bundled PyQt6/ChemUSON resources, verify internal version and perform a bounded headless launch.
- [x] 7.5 Preserve official artifact names and update-information string; embed the exact existing AppImageUpdate value and retain/validate `.updateinfo`, `.update.json`, `.zsync`, channel, tag and source SHA. Previews remain without update metadata.
- [x] 7.6 Update both preview and official release workflow steps and add regression tests for false AppImages, Type 2 signature, extraction/AppDir/launch validation, updater compatibility, SHA/provenance and workflow integration.
- [x] 7.7 Update current release/preview/Linux packaging documentation and manual acceptance criteria. Historical campaign records remain historical facts.

## 8. Python CI dependency addendum
- [x] 8.1 Install `requirements.txt`, `requirements-dev.txt`, and the editable project in `.github/workflows/test.yml`; do not duplicate PyYAML or add runtime dependencies.
- [x] 8.2 Add a static workflow contract proving dev dependencies are installed before full pytest collection and that the real suite is not skipped/ignored or masked.

## 9. Addendum validation
- [x] 9.1 Run the AppImage packaging, release workflow, preview helper/workflow, CI workflow, architecture and OpenSpec tests with a 10-minute cap per test command; run compileall, focused Ruff and diff checks.
- [x] 9.2 Perform real Linux PyInstaller→AppImageTool preview/release builds with bounded time, full Type 2/resource/icon/updater validation; do not claim graphical acceptance from a headless runner.
- [x] 9.3 Authorization recorded: `gh auth status` is unauthenticated. Packaging and CI fixes are in separate commits (`ca7f89e`, `ec20a70`, `46af127`); the noninteractive normal push failed before remote update (`could not read Username`). No password was used. Do not create tag/release or merge; push remains pending.

## 10. P1 packaged UI icon defect addendum
- [x] 10.1 Record the owner's confirmed missing-icon failures for Windows and Linux portable packages from preview run `37826597134`; mark manual acceptance `FAILED — P1 blocks beta acceptance` and preserve owner retest as a separate gate.
- [x] 10.2 Reproduce the shared cause from a real Linux PyInstaller binary: `collect_all("chemuson")` skips because it is not installed as a package in the build environment; binary archive has zero SVGs while QtSvg bindings are present.
- [x] 10.3 Add deterministic PyInstaller `datas` for exactly 69 package-relative static SVGs and fail if any are missing; preserve source/frozen lookup and themes.
- [x] 10.4 Add environment-gated frozen binary diagnostics for path resolution, all 69 SVG resources, QtSvg and visible raster pixels for essential icon categories in light/dark themes and DPR 2.
- [x] 10.5 Linux PyInstaller executable, extracted preview AppImage and Windows portable pass frozen-process path/resource/QtSvg/raster, QToolButton and RDKit checks in Build Preview `37854485440` on `d339b71a6ea8b7085572e59b37492c27b35759c8`; owner manual retest remains pending.
- [x] 10.6 Update tests/docs without altering the SVG inventory or chemical/UI behavior; owner manual retest remains NOT PASSED and blocks beta publication.

## 11. UI-ONBOARDING-001 and BRANDING-001
- [x] 11.1 Record the second Build Preview's Windows onboarding alignment report; preserve its prior manual result and mark the corrected-package retest pending in the acceptance matrix.
- [x] 11.2 Defer automatic onboarding until the main window's first shown layout; track parent/target geometry, keep the three steps and existing QSettings semantics, and test overlay cleanup, alignment, card bounds, target-geometry stability, resize and move.
- [x] 11.3 Normalize application-owned display strings and Windows/Linux presentation metadata to `ChemUSON`, including Help/About, update messages, installer display metadata, desktop entries and AppStream.
- [x] 11.4 Add static/focused regression tests for displayed branding, package metadata and technical identity/installer/updater compatibility; preserve all package IDs, commands, paths, update routes and artifact names.
- [x] 11.5 Update the active beta notes and manual acceptance matrix with `UI-ONBOARDING-001` and `BRANDING-001`; keep new package/manual results pending owner retest.
- [ ] 11.6 Owner manually retests the new preview packages on Windows and Linux, including three window sizes, common DPI scales, branding surfaces, installer upgrade/uninstall identity and updater compatibility.
- [x] 11.7 Replace the About description with the owner-approved ChemUSON project description and use the non-repeating em-dash main-window title; audit visible surfaces without changing historical credits, licenses, attributions or internal identifiers.

## 12. P1 — RDKit isolated backend unavailable in packaged executable
- [x] 12.1 Record the owner's Windows portable failure (formula/mass/estimated spectra work; RDKit descriptors report unavailable), mark this P1 as failed/blocking, and stop manual package verification until a corrected preview is available.
- [x] 12.2 Reproduce the pre-fix invocation failure from an actual frozen executable; separately verify RDKit/native imports in the installed runtime and make the corrected frozen smoke prove imports resolve inside each packaged bundle. Corrected Windows, Linux PyInstaller and extracted AppImage gates passed in Build Preview `37854485440` on SHA `d339b71a6ea8b7085572e59b37492c27b35759c8`.
- [x] 12.3 Implement a frozen-only worker dispatch in the existing executable; keep source mode compatible, retain process isolation/JSON API/timeouts, support Windows `console=False` without stdin/stdout, prevent GUI bootstrap/recursion, and ensure timeout cleanup/no orphaned process.
- [x] 12.4 Add required diagnostics for import/native-extension, worker start/exit, timeout, malformed response and chemistry errors; make the Properties pane label only actual RDKit import failure as unavailable while retaining partial results.
- [x] 12.5 Add a non-skipping frozen executable validator/smoke: require `sys.executable` to be the tested executable, extensions to load from its extracted bundle, ethanol logP/TPSA/HBD/HBA values, and bounded isolated 3D plus SMILES calls; assert the parent stays RDKit/GUI-free.
- [x] 12.6 Wire fail-closed frozen smoke gates into Windows portable, Linux PyInstaller and extracted AppImage jobs in both Build Preview and official release workflows; run the real Build Preview package jobs. Preview `37854485440` passed all package jobs; official release workflow was not triggered because no tag/release is authorized.
- [x] 12.7 Add focused unit/integration tests for source worker errors/protocol, frozen dispatch, app-UI messages, smoke validator, workflows and architecture catalog; do not alter chemistry algorithms, Clean2D, `.cmsn`, or the independent Qt teardown debt.
- [x] 12.8 Run bounded RDKit/app packaging tests, architecture, strict OpenSpec, scoped Ruff, compileall and diff check; push a normal commit only to `release/v0.3.0-beta.1-prep` and wait for the matching Preview/CI jobs, without tags/releases/publication. Preview `37854485440` passed; the separate full pytest CI job still aborts on the documented Qt teardown SIGSEGV and is not represented as passing.
- [ ] 12.9 Owner manually retests corrected Windows portable and AppImage/Linux packages; beta acceptance remains awaiting manual retest.

## 13. P1 — ChemName templates omitted from frozen packages
- [x] 13.1 Inspect actual Windows portable and extracted AppImage CArchive tables and the prior Flatpak bundle before changing code; compare against the nine source `.mol` templates.
- [x] 13.2 Record the source naming baseline, diagnose the masked `FileNotFoundError` with `return_nd_on_fail=False`, and inspect the status-bar and analysis-annotation call paths; do not modify naming rules.
- [x] 13.3 Add an explicit nine-template PyInstaller inventory and package-relative `datas` entries; keep existing setuptools/Flatpak package-data behavior and catalog the private M19→M04 smoke dependency.
- [x] 13.4 Extend the ChemName acceptance corpus with ethanol, acetamide, ethane and cyclohexane; add source/frozen smoke coverage for every template and representative molecule names.
- [x] 13.5 Add fail-closed Windows/Linux executable validators that compare actual frozen names/resources against Python source, and add a Flatpak `/app` package-data smoke.
- [x] 13.6 Gate Build Preview and official package workflows, record the historical missing-resource P1 and manual retest cases, and preserve the owner-pending state.
- [x] 13.7 Run bounded focused tests, strict OpenSpec and a new Build Preview for the pushed prep-branch SHA; report per-format evidence. Build Preview `37866732971` passed all four artifact groups on `1ad1db00aeb4f76e06cbffc5e23312b19e1dfa8f`; per-format hashes/automated smoke are in `validation.md`. Owner manual retest remains required before beta acceptance.
- [ ] 13.8 Owner repeats ChemName naming/resource checks on the corrected Windows portable/setup, extracted AppImage and Flatpak, recording the exact artifact SHA and both status-bar/analysis names.

## 14. ChemName live invalidation and independent Qt/CI diagnosis
- [x] 14.1 Separate the historical frozen-template omission from source algorithm correctness; rerun the bounded source acceptance corpus without changing nomenclature rules.
- [x] 14.2 Verify the actual status-bar signal path on a real molecular edit plus Undo/Redo; ensure each state is recomputed from the current graph before returning control.
- [x] 14.3 Reproduce the Qt teardown SIGSEGV in the ordered drag/lifecycle test shard and identify the retained `indexChanged` callback from `test_drag_move_undo.py`.
- [x] 14.4 Disconnect the test callback in `finally`; verify the same ordered shard completes without a crash. Do not refactor production Qt lifecycle code.
- [x] 14.5 Map the two ordinary failures in test run `37867439026` to exact nodes, reproduce boundedly where possible, and preserve Clean2D/CompChem code and tests without suppressing failures.
- [x] 14.6 Record the deferred `chemname/iupac-robustness` campaign, including authoritative reference names and a separate verdict for GUI staleness versus algorithm output. Do not begin presentation enhancements.
- [ ] 14.7 Owner supplies any disputed structure/name and completes the artifact-backed ChemName update/Undo/Redo and manual platform retests; no native Windows verification is inferred from CI.

## 15. Integración de las campañas ChemIO stereo y CI aprobadas
- [x] 15.1 Revalidar refs/upstream: origen `fix/ci-pytest-stabilization` en `571e012`, descendiente estricto del destino; destino base `6aeef19`, divergencia 0/10, CI de origen #379990 green.
- [x] 15.2 Integrar sólo con `git merge --ff-only`, preservando los commits y sin incorporar `chemname/iupac-robustness` ni cambios de subsistemas protegidos.
- [x] 15.3 Verificar localmente ChemIO stereo, panel, smokes/contratos de empaquetado, arquitectura, OpenSpec estricto, plan de ocho shards, versión, Ruff focal, compileall y diff.
- [x] 15.4 Publicar sólo `release/v0.3.0-beta.1-prep` (SHA `57d481d21fee2db6ecacf6a862839caa63a1c56d`) y verificar CI #38006731537 + Build Preview #38006731375 para los cuatro formatos.
- [x] 15.5 Descargar/verificar cada artifact, checksum/provenance y `publication=false`; registrar run/IDs/SHA en `validation.md` y notas beta.
- [ ] 15.6 Propietario completa aceptación manual en Windows y Linux; automation no cambia estados manuales.
