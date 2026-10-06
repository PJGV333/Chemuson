# Design

## Decisions

1. `verify_molecular_identity` defaults to `allow_network=False`. It still uses the existing offline Name→Structure sources and isolated InChI canonicalization. Its injected resolver receives the same explicit `allow_network` value. An `enabled=False` option returns `not_applicable` before resolving or canonicalizing. The GUI preference defaults to enabled + offline-only and offers explicit external-reference opt-in.
2. Persist only the non-secret identity policy through `platform.settings`; never persist the transient provider API key. When offline lookup has no reference, report `unverified`, not mismatch or verified, and tell the user external sources were not queried.
3. Validate transport security in `OpenAICompatibleConfig`, below the GUI. A non-empty key is valid over HTTPS or HTTP to an unequivocal loopback host only (`localhost`, IPv4 127/8, IPv6 ::1). Never resolve DNS to classify a host. Remote HTTP without a key remains allowed. Errors are generic and do not include credentials or endpoint contents.
4. Keep existing HTTP local llama.cpp/LM Studio configuration, OpenAI-compatible protocol, JSON decoder, ChemIO validation, M23 generation/transformation and M02 Clean2D unchanged. Add tests using fakes only.
5. Record the small-model research direction as a future hypothesis in `docs/history/CAMPAIGNS.md`; do not select a model, download data, build datasets, or fine-tune.

## Affected files

Expected production changes are limited to identity policy, provider config validation, platform preference functions, existing GUI policy controls/wiring, and contract tests. The module catalog should remain unchanged because no new import edge is introduced.

## Risks and checks

- Local model users must retain HTTP loopback with no key.
- Domain names resolving to private/LAN addresses are still remote for this policy and must use HTTPS when keyed.
- Explicitly disabled verification must avoid both resolver and canonicalization calls.
- Offline absence of a reference is uncertainty, never a semantic mismatch.
- Verify M23 remains separate from identity/name resolution and M02 remains free of AI dependencies.

Use fake resolvers/transports, focused tests, architecture tests, strict OpenSpec validation, changed-file Ruff rules, and `git diff --check`. Do not run the historical monolithic suite under the requested test budget.
