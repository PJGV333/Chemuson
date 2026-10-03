# Baseline — repository hygiene closure

**Date:** 2026-10-03 (local environment)  
**Branch:** `maintenance/repository-hygiene-closure`  
**HEAD / origin/main:** `1db4f63b52af79247745b3a8a220fb728348218c`  
**Initial status:** clean.

## Test baseline

- `python -m compileall src tests tools packaging`: exit 0 (captured at start of this campaign; `/tmp/chemuson-before-compileall.log`).
- `pytest --collect-only -q`: exit 0; 1816 collected (`/tmp/chemuson-before-collect.log`).
- `pytest -q`: 1760 passed, 55 skipped, 1 failed; the exact same source HEAD had just been exercised during Fase 8 and no code paths changed in the intervening documentation-only commit. The known failure is `tests/test_compchem3d_dock.py::test_compchem_controller_generates_async_with_fake_backend` (`/tmp/f8_final_full_suite.log`; Fase 8 report records 1152.05 s). This is the current baseline; post-clean must retain its identity and introduce no new failure.
- Scoped Ruff (`F401,F811,F821,E722,E741`): one F401, unused `math` import in `tests/test_clean2d_para_disubstituted_aromatic_layout_v1.py` (`/tmp/chemuson-before-ruff.log`). No other selected diagnostics.
- Architecture at the same source baseline: 269 passed in 9.42 s (`/tmp/f8_final_architecture.log`).

## Size baseline

Exact measurements and helper outputs are in `/tmp/chemuson-repository-before.txt`, `/tmp/chemuson-largest-blobs.txt`, `/tmp/chemuson-head-largest-files.txt`.

| Measure | Before |
|---|---:|
| Tracked HEAD logical file contents | 22,926,388 bytes / 1244 files (~21.87 MiB). |
| Checkout disk, excluding `.git` but including ignored local environment | 639 MiB. |
| Checkout disk, excluding `.git`, `.venv`, caches and build/dist directories | 25 MiB allocated. |
| `.git` directory | 169 MiB. |
| Loose objects | 784 objects / 8.05 MiB. |
| Packed objects | 14,142 objects; two packfiles / 160.27 MiB. |
| Requested directory allocation | `src` 30 MiB; `tests` 22 MiB; `docs` 3.0 MiB; `openspec` 3.9 MiB; `packaging` 288 KiB. |
| Tracked PNGs | 222 files / 4,597,068 bytes (4.38 MiB). |
| `src/sys` | 11,708,416 bytes; `file` identifies DSC Level 3 PostScript; raw header starts `%!PS-Adobe-3.0` and `%%Creator: (ImageMagick)`. |

`.venv` accounts for roughly the difference between the raw 639 MiB checkout and source/docs allocation; it is an ignored local Python environment and is not included in tracked-tree size.

## Blob/top-file measurement method

Historical blobs were enumerated with `git rev-list --objects --all` piped to `git cat-file --batch-check`, filtered to blobs and sorted by uncompressed object size. The complete list is `/tmp/chemuson-largest-blobs.txt`; only its top 30 was reviewed here. No blob was decoded as UTF-8. Current HEAD files were measured separately using `git ls-tree -rl -z HEAD`; top entries are in `/tmp/chemuson-head-largest-files.txt`.

The top historical entries include `flatpak/beta/repo/objects/*.filez` (largest 12,372,162 bytes) and `src/sys` (11,708,416 bytes). These historical blobs will remain in Git history; this campaign does not rewrite history.

## Remote branches

The complete pre-clean inventory, branch SHA, ancestry result, ahead/behind counts and tip subject are in `/tmp/chemuson-remote-branches-before.txt`. Of 33 remote refs (including `origin` and `origin/main`), 25 refs were ancestry-contained and 8 tips were not. The user-specified twelve `architecture/*` tips were all confirmed ancestors of `origin/main`; so were `origin/fix/ui-openspec-post-integration` and `origin/release/ui-modernization-qa`. Non-ancestor branches require content review before any decision.
