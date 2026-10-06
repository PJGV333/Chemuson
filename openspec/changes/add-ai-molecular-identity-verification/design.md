# Design

## Architecture

Use M10's asynchronous request flow to invoke an application-side verifier after M23 returns a ChemIO-validated result. The verifier reuses M16 `resolve_name_to_structure` through an injectable resolver and M01 isolated ChemIO/RDKit worker APIs. M23 remains provider-neutral and does not resolve names; M02/Clean2D remains independent. Update the module catalog only if imports introduce an edge not already declared.

## Recognition and comparison

Recognize only full requests shaped like `Dibuja [el/la] X`, `Draw [the] X`, `Generate [the] X`, or `Genera [el/la] X`. Reject generative descriptions (e.g. “a molecule with three fused rings”). If no exact name-like subject is identified, report `not_applicable`. Compare isolated stereo-sensitive InChI strings, not SMILES text or 2D coordinates.

A trusted resolver result plus successful canonical identity for both structures yields `verified` or `mismatch`. No reference yields `unverified`; resolver/canonicalization errors yield `reference_error`. The verifier never infers success from the model's natural-language response.

## UI policy

Show ChemIO validity independently from identity status. A mismatch is not insertable by default: the user must choose an explicit “Insert anyway” action and accept a second confirmation. For `unverified`/`reference_error`, state clearly that identity was not verified and allow normal review. The Clean2D report includes requested name, identity status, and reference identifier without calling an identity mismatch a successful named-molecule evaluation.

## Tests

Use fake resolvers entirely offline; isolated ChemIO canonical identity can be exercised with small deterministic graphs. Include different text SMILES for equivalent identity, distinct structures, resolver missing/failure, open-ended descriptions, ignored model text, guarded mismatch insertion, explicit override, and the named cholesterol regression label.
