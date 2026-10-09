# Design — narrowly scoped nomenclature corrections

## Context and evidence

The current engine detects the primary amide suffix for `CC(=O)N`, but also detects the same amide nitrogen as a standalone amine and renders `1-aminoethanamide`. In a second case, `CC(=O)Cc1cc(N)ccc1C`, it names the ring simply `phenyl`, returning `1-phenylpropan-2-one` and dropping the amino and methyl groups. The exact inputs and reference evidence are in `corpus.md`.

The IUPAC Blue Book, *P-66.1.1.1.1.1*, specifies substitutive amide formation by adding the suffix `amide` to the parent hydride (eliding final `e` before `a`); on the user's requested systematic output this gives `ethanamide`. The same section, *P-66.1.1.1.2.1*, identifies retained `acetamide` as the PIN. Both facts are recorded explicitly: this campaign implements the requested systematic display, not a claim that `ethanamide` is the Blue Book PIN.

## Goals / Non-Goals

**Goals:**
- Prevent an amide N from also generating an amine prefix.
- Name the documented substituted phenyl group with fixed ipso locant 1, minimum supported ring locants, all supported substituents, and parentheses around the complex substituent.
- Preserve input/output connectivity and functional-group identity; reject unsupported stereochemical, charged, isotopic, multiply attached, cross-linked, or otherwise unrecognized aryl decorations rather than flattening them to `phenyl`.
- Keep the acceptance corpus exact and auditable, with independent references and per-case status/timing evidence.

**Non-Goals:**
- General IUPAC coverage, direct acyl-substituted benzenes such as acetophenone, nested/fused/polycyclic aryl substituent expansion, or new nomenclature for unreferenced molecules.
- Changing stereo algorithms, isotope/formal-charge presentation, or any name outside the minimal regression corpus.
- Changes to ChemIO/RDKit worker isolation, GUI freshness, Clean2D, persistence, package identity, or release/publication.

## Decisions

1. **Amide N is not a second amine.** Reuse the existing carbonyl/amide recognition. When a nitrogen attached to the parent chain carbon is already the N component of a recognized amide occurrence, do not append an independent amine occurrence/prefix. Keep the amide as the principal suffix under the existing functional-group priority.
2. **Name the aryl substituent from its attachment atom.** For a benzene used as a substituent, fix the attachment carbon as phenyl locant 1; consider both ring directions from that atom, reuse the existing supported ring-substituent recognizers, and select the lowest locant sequence with the existing alphabetical tie-break. Render decorated forms as parenthesized substituents (for example `(5-amino-2-methylphenyl)`); retain bare `phenyl` for an unsubstituted ring.
3. **Fail closed before rendering unsupported information.** Validate that the named ring has one single-bond connection to the parent. Reject unsupported cross-links, nested ring decorations not covered by this tranche, formal charge/isotope/radical state not represented by the current recognizers, and any unrendered stereo annotation on the aryl group or its non-parent branches. Existing supported substituent detectors may name ordinary achiral neutral amino/methyl/halogen/hydroxy/nitro/alkoxy groups; an unrecognized group raises `ChemNameNotSupported` and returns `N/D` under the default safe API.
4. **Use exact output tests, not fuzzy performance gates.** Add the exact structures to the existing acceptance corpus and source/package smoke. Treat `pass`/`fail`/`error` counts and exact names as correctness gates; retain harness build/name/total durations as diagnostic metrics without machine-dependent millisecond thresholds.
5. **No new dependency or package.** Keep implementation inside existing `chemname` modules and reuse the current graph, ring, locant, and rendering abstractions. The explicit two-case corpus and negative sentinels bound the claim.

## Risks / Trade-offs

- The parenthetical substituent needs correct attachment-fixed locants. Tests use the independent exact name and an alternate atom-order representation of the same molecular graph.
- Shared ring recognizers may support more cases than this tranche can preserve fully. The helper must reject rather than silently fall back to `phenyl` when metadata or connectivity is not represented.
- PubChem's `IUPACName` field for acetamide returns the retained PIN `acetamide`; the campaign's exact `ethanamide` expectation is intentionally a requested systematic form supported by the Blue Book's substitutive suffix rule, not a PIN assertion.
- `N/D` for cases outside the tested support boundary is intentional and safer than a lossy name.

## Migration / rollback

No data migration or public interface change is needed. If a focused test exposes ambiguity or a lost annotation, revert the corresponding small commit rather than widening the scope. No release or distribution action is part of this campaign.
