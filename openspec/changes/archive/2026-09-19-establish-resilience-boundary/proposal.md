# Establish resilience boundary

## Why

Autosave and crash reporting are generic runtime-resilience services but are
currently owned by `utils`, while GUI recovery consumes them through historical
paths. The ownership is broad and obscures the failure-isolation boundary.

## What Changes

- Add M22 `resilience` with canonical autosave and crash-reporting modules.
- Preserve `utils.autosave` and `utils.crash_reporter` as import-only shims.
- Make bootstrap, tab management and recovery consume canonical M22 imports.
- Register M22 and update only the direct dependency edges required by those
  imports.

## Non-goals

- Do not move RecoveryController, PersistenceManager, tabs or widgets.
- Do not change autosave JSON, recovery metadata, crash-log format or timers.
- Do not add GUI imports to the canonical resilience package.
