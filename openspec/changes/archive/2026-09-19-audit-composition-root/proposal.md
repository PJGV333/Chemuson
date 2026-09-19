# Audit composition root

## Why

M19 is the architectural composition root and must not be split merely to fill
the next module slot. Its current files are small and their dependency direction
is intentional.

## What Changes

- Document the M19 composition-root audit and consolidated ownership.
- Preserve `__main__.py` and `app/bootstrap.py` as the single startup boundary.
- Reserve M24 without creating a placeholder module.

## Non-goals

- No startup behavior, CLI behavior or GUI composition changes.
- No new module ID, package, shim or dependency edge.
