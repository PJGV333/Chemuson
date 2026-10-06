## Baseline

- [x] Verify the clean branch/HEAD equals the expected origin branch HEAD before edits.
- [x] Capture compileall, collect-only, focused assistant/UI/transform/identity tests, architecture, and scoped Ruff baseline.
- [x] Respect the hard 10-minute command limit; do not run monolithic pytest.

## Structured-output capability and fallback

- [x] Add prompt-only/native/unknown capability representation with legacy config compatibility.
- [x] Recognize only exact `response_format_not_supported`; allow at most one request without `response_format`.
- [x] Strengthen the provider-neutral system prompt, including correct JSON escaping for backslash SMILES.
- [x] Return bounded output-capability diagnostics without model text or secrets.

## Strict format repair

- [x] Retry once only for strict-decoder `invalid_json`, with a byte/prompt-bounded format-only request containing the original content as quoted untrusted data.
- [x] Decode the retry through the unchanged exact schema and isolated ChemIO; preserve downstream identity verification.
- [x] Prove shared generate/transform recovery, no regex/fence cleanup, no third request, and no secret leakage with offline fakes.
- [x] Show small, non-persistent UI diagnostics without reasoning/model raw content.

## Validation and live check

- [x] Run focused tests, architecture, compileall, strict OpenSpec, changed-file Ruff, and diff checks under limits.
- [x] Check the already-running loopback `/models` endpoint; run bounded ethanol/caffeine/tetrandrine probes only because Qwen was available; no server was started and no external Internet was used.
- [ ] Commit and normally push only this branch; no merge/rebase/force-push/squash.
