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
- [ ] 6.4 Commit focused changes on this prep branch and push only if normal authentication is available; verify remote branch SHA. Never merge/rebase/force-push or modify `main`, tags, GitHub Releases or `gh-pages`.

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
- [ ] 9.1 Run the AppImage packaging, release workflow, preview helper/workflow, CI workflow, architecture and OpenSpec tests with a 10-minute cap per test command; run compileall, focused Ruff and diff checks.
- [x] 9.2 Perform real Linux PyInstaller→AppImageTool preview/release builds with bounded time, full Type 2/resource/icon/updater validation; do not claim graphical acceptance from a headless runner.
- [ ] 9.3 Record whether normal push authorization exists. Commit AppImage and CI fixes separately; push only this prep branch if possible. Do not create tag/release or merge.

## 10. P1 packaged UI icon defect addendum
- [x] 10.1 Record the owner's confirmed missing-icon failures for Windows and Linux portable packages from preview run `37826597134`; mark manual acceptance `FAILED — P1 blocks beta acceptance` and preserve owner retest as a separate gate.
- [x] 10.2 Reproduce the shared cause from a real Linux PyInstaller binary: `collect_all("chemuson")` skips because it is not installed as a package in the build environment; binary archive has zero SVGs while QtSvg bindings are present.
- [x] 10.3 Add deterministic PyInstaller `datas` for exactly 69 package-relative static SVGs and fail if any are missing; preserve source/frozen lookup and themes.
- [x] 10.4 Add environment-gated frozen binary diagnostics for path resolution, all 69 SVG resources, QtSvg and visible raster pixels for essential icon categories in light/dark themes and DPR 2.
- [ ] 10.5 Linux actual PyInstaller executable and authentic preview/release AppImages pass frozen-process path/resource/QtSvg/raster checks. Windows preview/release executable checks are wired fail-closed but remain unexecuted locally; verify in the next Actions run.
- [x] 10.6 Update tests/docs without altering the SVG inventory or chemical/UI behavior; owner manual retest remains NOT PASSED and blocks beta publication.
