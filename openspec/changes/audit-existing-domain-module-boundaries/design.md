# Design

The catalog remains authoritative. The audit report mirrors the catalog entry
for each selected module and records the conclusion `audited / no structural
change required` when current and target ownership already agree.

A small architecture test checks that every audited ID appears in the report and
that its catalog entry has equal current/target dependencies with no temporary
exceptions or circular dependencies. No production import or path changes are
needed.
