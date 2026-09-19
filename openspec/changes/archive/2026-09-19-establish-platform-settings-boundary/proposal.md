# Establish platform settings boundary

## Why

Application configuration currently enters the GUI shell as a concrete
`QSettings` object and packaged-resource resolution is owned by `utils`, mixing
platform policy with GUI assembly and shared utilities.

## What Changes

- Add M21 `platform.settings` for application preference policy and packaged
  resource resolution.
- Preserve historical preference keys, defaults and value coercion.
- Keep `utils.resources` as an import-only compatibility shim.
- Make GUI consumers depend on the platform settings contract rather than
  constructing or interpreting QSettings policy themselves.
- Register the new ownership and dependency direction without renumbering
  M00-M20.

## Non-goals

- Do not move document state, selection, undo, canvas state, update policy,
  widgets or domain-specific chemical options.
- Do not change CMSN formats, preference keys, visible behavior or public GUI
  APIs.
- Do not introduce a dependency from platform to ChemUSON GUI modules.
