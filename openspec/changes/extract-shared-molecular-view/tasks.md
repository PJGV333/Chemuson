## 1. Characterization

- [x] 1.1 Finalize and validate MolView characterization tests using the
  historical ChemName import.
- [x] 1.2 Add architecture regression coverage that rejects M03 exceptions and
  M03-to-M04 imports.

## 2. Extraction

- [x] 2.1 Move the sole MolView implementation and neutral error to Core and
  expose the Core public API.
- [x] 2.2 Preserve ChemName MolView and ChemNameNotSupported compatibility
  exports with exception identity.
- [x] 2.3 Migrate ChemCalc to the Core adapter without changing its behavior.

## 3. Architecture And Documentation

- [x] 3.1 Update the module catalog to remove the M03/M04 cycle and M03
  exceptions while retaining all six M01 exceptions.
- [x] 3.2 Update module documentation for Core ownership and historical
  compatibility.

## 4. Validation And Closure

- [x] 4.1 Run focused tests, architecture tests, full suite, and Ruff.
- [ ] 4.2 Validate the OpenSpec change, complete the checklist, commit the
  scoped work, and archive the validated change.
- [ ] 4.3 Archive the validated change after pre-archive verification.
