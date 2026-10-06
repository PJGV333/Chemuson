# Harden Molecular Assistant locality and credential transport

## Why

Identity verification currently allows network access by default, so a local-model request can implicitly query PubChem. The OpenAI-compatible provider also permits a non-TLS remote endpoint while attaching a bearer key. Both boundaries need secure defaults without changing generation, transformation, ChemIO, or Clean2D behavior.

## Scope

- Make identity verification offline-only by default, with explicit offline/external/disabled policy, a small UI control, and non-secret persisted preferences.
- Require HTTPS for API keys sent to non-loopback hosts; retain HTTP support for loopback local servers and for remote endpoints without keys.
- Add offline contract tests and record a future local chemistry specialist-model research campaign.

## Out of scope

No model training/fine-tuning, datasets, new provider protocol, changes to M23 generation/transformation or ChemIO validation, Clean2D changes, manual model evaluation, or claims that the overall Molecular Assistant project is complete.

## Likely impact

`name2structure.identity`, M23's OpenAI-compatible config, platform settings, the existing Molecular Assistant dialog/controller wiring, focused tests, OpenSpec, and `docs/history/CAMPAIGNS.md`. No dependency or module-boundary change is intended.
