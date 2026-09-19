# Audit existing domain module boundaries

## Why

The catalog contains mature domain modules whose ownership must be verified
before introducing platform and shell boundaries. The audit must distinguish
real structural debt from responsibilities that are already correctly owned.

## What Changes

- Record an explicit audit of M05, M06, M07, M14, M16, M17, M18, M19 and M20.
- Verify paths, runtime and target dependencies, forbidden dependencies, public
  and internal APIs, tests, exceptions, cycles and documentation.
- Keep cohesive modules unchanged and preserve all existing IDs.

## Non-goals

- Do not create duplicate domain modules.
- Do not renumber M00-M20.
- Do not change chemical behavior, public APIs or runtime dependencies.
