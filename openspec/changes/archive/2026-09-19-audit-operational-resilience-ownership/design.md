# Design

The audit treats resilience as a cross-cutting concern with explicit owners:
M22 for generic failure persistence and isolation, M08/M10 for Qt task
containment, M14 for update telemetry, and M19 for bootstrap installation.

The document and architecture test are the contract. No code moves because each
owner already has the dependencies needed for its operational responsibility.
