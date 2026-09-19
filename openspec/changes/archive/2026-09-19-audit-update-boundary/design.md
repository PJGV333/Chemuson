# Design

Use the existing M14 catalog entry as the source of truth and record each
cohesive update cluster in a maintenance audit. The audit contract checks that
M14 still owns `src/chemuson/update/` and that M23 is only a reserved slot.

No code extraction is justified: moving individual update files would create
cross-file coupling without reducing policy, provider, security or platform
responsibilities.
