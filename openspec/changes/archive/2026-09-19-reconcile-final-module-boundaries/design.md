# Design

M22 receives a GUI-free `resilience/recovery.py` containing the existing
`read_autosave_metadata`, `list_autosave_entries` and `archive_autosave`
operations. RecoveryController retains thin static wrappers for historical
class access and keeps document loading, Qt dialogs and window coordination.

The catalog uses explicit file ownership for M21's package files, records M22's
recovery policy API, and states that M15 owns only compatibility shims. M19's
composition-root dependency set is corrected to M08/M18/M22. M20 owns only its
selection namespace; M09 owns the legacy canvas shim paths and depends on M20.

The archived application-composition-root specification receives a normative
Purpose so `openspec validate --all --strict` has no unrelated failure.
