# Reconcile final module boundaries

## Why

The module catalog and archived specifications still describe historical ownership:
M19 omits M22, M15 describes canonical autosave/crash ownership, M21 uses a broad
package path, and M20/M09 wording does not distinguish canonical helpers from
legacy shims. Recovery filesystem policy also remains embedded in the GUI
controller.

## What Changes

- Reconcile the catalog and OpenSpec wording for M15, M19, M20, M21 and M22.
- Make M21's owned platform files explicit.
- Move the three deterministic autosave filesystem-policy functions to M22.
- Keep RecoveryController as the Qt/dialog/document orchestration layer.
- Remove the remaining global OpenSpec validation failure.

## Non-goals

- Do not change autosave JSON, recovery metadata, archive naming or UI behavior.
- Do not move document loading, tabs, dialogs, widgets or controller orchestration.
- Do not create M23 or M24, add dependencies, or alter public historical imports.
- No credentials or secrets are stored; sensitive values are represented as `[REDACTED]`.
