# Molecular identity verification for AI proposals

## Why

A live request for cholesterol returned a syntactically and chemically interpretable SMILES that was not cholesterol. ChemIO validation proves structural parseability, not semantic identity.

## Changes

Add conservative recognition of simple named-molecule requests, reuse the existing Name→Structure resolver, compare proposed/reference structures through isolated canonical InChI identity, and present distinct `not_applicable`, `unverified`, `verified`, `mismatch`, and `reference_error` states. A mismatch requires an explicit confirmed override before insertion. Include the `cholesterol-semantic-mismatch-01` offline regression case and expose identity provenance in the Clean2D evaluator when available.

## Boundaries

No second LLM, raw RDKit import, drawing comparison, Clean2D dependency, or semantic claims based only on ChemIO parsing. Resolver/network failure is not success; open-ended generation remains unverified/not applicable.
