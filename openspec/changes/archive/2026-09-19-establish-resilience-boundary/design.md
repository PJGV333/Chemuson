# Design

M22 owns the existing implementations moved verbatim to `src/chemuson/resilience/`.
The package contains no ChemUSON imports; crash notification keeps its existing
external QtWidgets dependency. GUI and bootstrap callers consume M22 directly.

The old `utils` modules re-export the exact public symbols and remain available
to historical consumers. This leaves M15 as a compatibility owner with direct
M21/M22 shim dependencies, while the canonical services remain outside `utils`.

RecoveryController keeps document loading and dialogs in M10; it uses only the
M22 autosave manager for the shared autosave directory policy.
