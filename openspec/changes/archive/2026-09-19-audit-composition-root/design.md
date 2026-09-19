# Design

Use M19's existing catalog entry and startup tests as the source of truth. The
new audit records that CLI parsing and GUI composition are already separate
responsibilities inside one small composition boundary.

No extraction is justified: moving the QApplication bootstrap or CLI parser
would add indirection without improving ownership or testability.
