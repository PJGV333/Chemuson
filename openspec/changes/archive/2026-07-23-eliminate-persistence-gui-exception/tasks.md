## 0. Baseline

- [x] 0.1 Capture architecture, persistence, collection, suite and Ruff baselines; acceptance: expected one TYPE_CHECKING exception and no unexplained failures.

## 1. Reconnaissance

- [x] 1.1 Inventory persistence consumers, model operations, canvas methods, autosave compatibility and type-checker availability; acceptance: temporary reconnaissance records the minimal structural surface and risks.

## 2. Characterization

- [x] 2.1 Add pure-fake persistence characterization; acceptance: payload, validation, load ordering, rebuild count, file I/O and import isolation are covered without GUI.

## 3. PersistenceDocument Contract

- [x] 3.1 Define the internal `PersistenceDocument` Protocol beside `PersistenceManager`; acceptance: it contains only model and persistence operations and no GUI type.

## 4. Public Canvas Rebuild Hook

- [x] 4.1 Add the public canvas reconstruction hook; acceptance: it delegates once to the existing private rebuild implementation without changing it.

## 5. PersistenceManager Migration

- [x] 5.1 Migrate persistence annotations and load reconstruction call to the structural contract; acceptance: no GUI import, string annotation or private canvas call remains in `persistence.py`.

## 6. Call-Site Compatibility

- [x] 6.1 Validate FileController, RecoveryController and autosave call sites without production changes; acceptance: real canvas remains structurally accepted and no adapter is introduced.

## 7. Catalog and Zero-Exception Baseline

- [x] 7.1 Remove the M01-to-M09 catalog exception and empty the frozen baseline; acceptance: all modules retain empty exception lists and no cycles.

## 8. Architecture Enforcement

- [x] 8.1 Add catalog, exception and AST regressions; acceptance: M01-to-M09 runtime, local and TYPE_CHECKING imports and the historical exception identity fail tests.

## 9. Documentation

- [x] 9.1 Document GUI-neutral persistence and zero debt; acceptance: architecture and M01/M09 module docs describe structural injection without claiming no runtime GUI objects.

## 10. Validation

- [x] 10.1 Run focused checks, compileall, OpenSpec strict validation, full suite and global Ruff; acceptance: execution and environment limitations are recorded in `validation.md`; all runnable in-scope checks pass and Ruff retains only the established single F401.

## 11. Archive

- [x] 11.1 Archive the validated change after review; acceptance: approved by the user and archived after successful validation.
