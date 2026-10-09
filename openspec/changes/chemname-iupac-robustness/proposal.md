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

Production scope is limited to `src/chemuson/chemname/`. Test/data changes are limited to ChemName acceptance and smoke coverage. Campaign evidence is recorded under this OpenSpec and in the existing roadmap campaign entry. No new package, runtime dependency, architecture edge, formula/isotope/charge presentation, `.cmsn`, GUI, Clean2D, RDKit worker, release workflow, tag, or public channel is changed.

## Campaign expansion (authorized 2026-10-08)

The initial two-case tranche is valid but insufficient for campaign completion. Continue on the same branch and from the existing commits. The expanded objective is to measure a broad, reference-backed set of functional-group and aromatic families; correct reusable priority, parent-chain, suffix/prefix, locant, and group-accounting defects; and test that unsupported information is never silently omitted. The independent survey currently resolves 78 unique PubChem structures (77 new beyond the prior two-structure baseline after accounting for the shared acetamide structure).

This is a staged, additive amendment to the active campaign, not a replacement OpenSpec. Stages 2–5 are tracked in `tasks.md`; exact references, previous names, and baseline classifications are in `corpus.md`. A PubChem name is not assumed to be a PIN, and unadjudicated variants remain `reference_pending` instead of being promoted to correctness.

The new work stays within `src/chemuson/chemname/`, its tests/fixtures, and campaign-specific documentation. The two simple aryl-ketone cases `CC(=O)c1ccccc1` and `CC(=O)c1ccc(O)cc1`, formerly outside the first tranche, are now explicit reference-backed stage-3 candidates. This does not authorize general aromatic ketones, fused/nested systems, or broader unrelated modules.
