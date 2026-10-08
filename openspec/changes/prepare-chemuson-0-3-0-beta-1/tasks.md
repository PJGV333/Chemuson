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
