# Proposal — ChemName / IUPAC robustness

## Why

The beta-preparation campaign deferred nomenclature changes until exact structures and independent references were available. Two cases are now reproducible: a primary amide is emitted with a spurious `amino` prefix, and a ketone containing a multisubstituted phenyl group loses both ring substituents. The existing acceptance corpus encoded the first wrong answer, so green status alone was not evidence of correct naming.

## What Changes

- Classify the nitrogen of a primary carboxamide as part of the amide group, not as a separate amine substituent. For `CC(=O)N`, emit the requested systematic name `ethanamide`.
- Preserve supported substituents and their locants when a benzene ring is named as a phenyl substituent, including the documented amino/methyl case `CC(=O)Cc1cc(N)ccc1C`.
- Add a documented, reference-backed corpus, exact-name acceptance checks, and fail-closed sentinels for unsupported ring connectivity or stereochemical decoration that must not be silently omitted.
- Correct the existing ChemName acceptance and packaged-smoke expectation for the primary-amide case.

## Capabilities

### Modified Capabilities

- `chemname-iupac-naming`: primary carboxamides are not double-counted as amines; a supported substituted phenyl group retains its substituents, locants, and unambiguous grouping; unsupported cases fail closed.

## Impact

Production scope is limited to `src/chemuson/chemname/`. Test/data changes are limited to ChemName acceptance and smoke coverage. Campaign evidence is recorded under this OpenSpec and in the existing roadmap campaign entry. No new package, runtime dependency, architecture edge, chemistry outside these cases, formula/isotope/charge presentation, `.cmsn`, GUI, Clean2D, RDKit worker, release workflow, tag, or public channel is changed.
