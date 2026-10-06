# Improve AI structured-response recovery

## Why

Some OpenAI-compatible local endpoints reject `response_format=json_object`, and otherwise successful model responses can contain malformed JSON, Markdown, or incorrectly escaped SMILES. The current service fails closed, which is correct, but has no explicit capability fallback or bounded formatting-only retry.

## Scope

- Represent prompt-only, native JSON-object, and unknown structured-output capability.
- On the exact `response_format_not_supported` signal, make one request without `response_format`.
- On strict-decoder `invalid_json`, make at most one bounded format-repair request and decode it with the same strict schema before ChemIO.
- Record non-secret, non-persistent capability/repair diagnostics for generation and transformation.

## Out of scope

No heuristic JSON extraction/cleanup, schema or chemistry relaxation, identity-verification bypass, Clean2D changes, model training, scraping, agents, or changes to model chemical knowledge.

## Likely impact

M23 provider/config/response models/service, existing Molecular Assistant dialog diagnostics, fake-transport/provider tests, and the M23 OpenSpec delta. No new dependencies or module-boundary changes are intended.
