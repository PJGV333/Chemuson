# Design — release preparation and channel-safe publishing

## Context

See `proposal.md` for motivation. Current state and baseline are in `baseline.md`. The project already has a SemVer parser, `_version.py` as the dynamic package source, a version-sync script, updater manifests/checksums, Windows/Flatpak builds and a Linux portable executable called `.AppImage`. The current release workflow independently accepts version/channel, rewrites source during each build, has no same-SHA test gate, and has write permission at workflow scope.

## Goals / Non-Goals

**Goals:** prepare `0.3.0-beta.1` on this branch; enforce one source/tag/channel; ensure release gate and official platform builds use the exact same SHA; provide a separately downloadable, no-publication preview build from preparation branches; preserve artifact integrity/provenance; provide human-reviewable baseline exceptions and acceptance materials.

**Non-Goals:** create/push any tag, run a release workflow, publish GitHub Release assets, change `gh-pages`, install public channels, expose preview builds to the updater, add a real AppImage technology, or change chemical, updater runtime, GUI, or `.cmsn` behavior.

## Decisions

1. **Canonical version and immutable build input.** Set `_version.py` to `0.3.0-beta.1` using the existing `set_version.py`, with a fixed preparation date in AppStream. Keep `pyproject.toml` dynamic. Release CI validates, but never edits, the canonical version; tag, canonical version, AppStream's current release entry and build outputs must agree. Inno receives the already-validated version as `CHEMUSON_VERSION`. Alternative rejected: rewriting `_version.py` independently in each build, because the artifact then differs from the tag's source tree.

2. **Tag-only release trigger.** Keep `push.tags: v*`; remove `workflow_dispatch` and its independently selectable, stale version/channel. A stdlib preflight accepts only strict stable `X.Y.Z`, `X.Y.Z-beta.N`, or `X.Y.Z-rc.N`; `beta` and `rc` map to the beta update channel, stable maps only to stable. It checks tag type/name, source version, AppStream metadata, checked-out/tag/GITHUB SHA equality, and GitHub API absence of an existing release. Any API error is fail-closed. A published tag/release is never overwritten; maintainers must configure an immutable `v*` tag ruleset and protected stable environment in GitHub settings.

3. **Same-SHA release gate inside `release.yml`.** Add a read-only `release_gate` job that checks out the peeled tag SHA and runs bounded architecture, release/version/updater tests, compileall, and scoped Ruff. All platform build jobs depend on it and check out its SHA output; the final assembly/publication jobs depend on every build. `test.yml` remains normal branch/PR CI and is not treated as evidence for a different tag SHA. Alternative rejected: gating on a separate workflow's latest green status, which can refer to another commit or be absent.

4. **Narrow permissions and no overwrite.** Default workflow permission is read-only. Only GitHub Release and Flatpak `gh-pages` publication jobs receive `contents: write`; they run after build/gate and use the validated channel. Existing channel-specific Flatpak directories and updater manifest layout are retained. Release assets are created only if absent; collisions/unmatched expected assets fail.

5. **Integrity and provenance.** Keep SHA-256 checksums mandatory. Add a small provenance JSON containing tag, version, channel and full source SHA; include it in the checksum list and release payload. Keep HMAC and Flatpak GPG as optional existing features and label their status accurately. No cryptographic or build dependency is added.

6. **Bounded exact baseline register, not pytest suppression.** Add a version-controlled catalog with exact node IDs/reproducers, evidence, owner role/tracking reference and retirement condition. Validate its schema in release-gate tests. It does not deselect tests, convert failures to skips, or mark an aborted suite green. Release gate uses focused groups and fails on every failure in those groups. A packaged-application P0/P1 or any unlisted release-gate failure blocks promotion.

7. **Packaging claims remain factual.** Preserve the official PyInstaller/Flatpak/Inno mechanisms. Document that the current `.AppImage`-named file is a PyInstaller portable executable, not a Type 2 AppImage; do not install or introduce appimagetool in this campaign. Windows interactive runtime, Flatpak bundle build, and AppImage Type 2 cannot be verified on this host; existing GitHub runner jobs and manual matrix are the remaining evidence path.

8. **Manual acceptance and recovery.** Prepare curated beta notes and a 50+ case matrix with stable IDs/status fields. Beta may be published for installation/acceptance before every manual case is complete. Stable requires explicit owner approval, priority cases, zero open candidate-attributable P0/P1, and documented acceptance of any P2/P3. Never move/reuse tags or assets; correct a rejected beta with the next beta number, and a stable regression with a new patch version. A partially failed publication is not repaired by overwriting the published tag; pause the affected channel and use a new higher version.

9. **Preview builds are an isolated Actions-only path.** Add `build-preview.yml` with manual `workflow_dispatch` and a narrow push filter for `release/**-prep`. The selected ref must be a branch matching the preparation pattern; the workflow resolves the canonical `_version.py` value, records the event SHA, checks out that exact SHA in every platform job, and never changes it. It reuses the existing PyInstaller spec, Inno script, and Flatpak build script. `build_appimage.sh` gains an explicit preview mode that emits only a uniquely named portable binary—no public `.updateinfo`, `.update.json`, `.zsync`, channel manifest, or release URL. Flatpak preview builds omit remote URL and signing credentials and upload only the local bundle. Windows/Flatpak preview outputs use preview/SHA names. Every package group contains SHA-256 checksums and a provenance record (version, build type, source SHA/branch, UTC build time, OS, status, publication=false); missing/empty files fail before upload. No release/publish action or write permission is present.

10. **Preview outputs are not official releases.** Four named Actions artifacts separate Windows portable, Windows installer, Linux portable `.AppImage`-named executable, and Linux Flatpak. A final report records each platform job result, including failures. Downloaded preview artifacts are manual QA inputs only; they do not update an updater manifest or any public beta/stable channel. If the workflow code is not present on GitHub or authorization is unavailable, no API/CLI workaround is attempted; document the GitHub UI/authorized CLI steps and leave real builds pending.

11. **Static isolation contracts.** Tests parse the preview workflow and assert read-only permissions, exact-SHA checkout, expected artifact groups/existence checks, no publication/tag/`gh-pages` operations, and separation from the official workflow. The tests establish workflow structure only; they are not evidence that Windows, AppImage-named portable, or Flatpak packaging actually built. Real preview validation requires a successful GitHub Actions run with all four outputs.

## Risks / Trade-offs

- [GitHub environment/tag rulesets are repository settings, not versioned workflow files] → document exact required settings; do not claim stable protection is active unless the owner verifies it.
- [GitHub API is unavailable or rate-limited during preflight] → fail closed before builds; retry the same not-yet-published tag only after confirming no release exists.
- [Release build or channel publication partially succeeds] → preserve all published bytes/tags; stop channel promotion and recover with a higher version, never overwrite.
- [Cross-platform tools are absent locally] → validate source/config/tests here; rely on the release gate's Windows/Linux CI builds and the manual matrix before stable.
- [Preview workflow is not yet available on GitHub or lacks Actions authorization] → do not launch it through an alternate privileged path; document the pending real-run status and use the owner-facing guide.
- [Known monolithic Qt teardown abort] → do not rerun the monolith; keep it explicit in the exception register, run bounded release-critical families, and require packaged-app close tests. Do not claim the debt resolved.
- [Linux file is named `.AppImage` but is not Type 2] → disclose in notes/docs and test its actual PyInstaller executable behavior; no claim of Type 2 support.

## Migration Plan and rollback

1. Use `set_version.py` on this branch to set the canonical beta version and AppStream entry; preserve old metadata/tags/assets.
2. Validate all release policy/unit tests, strict OpenSpec, package config and bounded CI-equivalent gates.
3. Push this preparation branch only. Owner reviews and later decides whether to merge/tag; no release operation occurs here.
4. For future publication, owner creates a protected tag on the reviewed commit; tag-only workflow validates and builds the exact SHA. A failure before publication leaves no version released. A failure after partial publication blocks promotion; use a higher version for recovery.
5. Stable promotion is a separate reviewed stable tag with separately protected environment/channel. Beta artifacts are never relabeled stable.

## Open Questions

None. GitHub environment/ruleset configuration is an explicit external prerequisite and is not changed by this branch.
