## Baseline

- [x] Confirm clean branch/HEAD and matching `origin/ai/molecular-assistant-foundation` before edits.
- [x] Capture compileall, collection, focused test, architecture, and scoped Ruff baseline; do not run monolithic pytest (historical runtime exceeds the task budget).

## Identity locality

- [x] Default identity resolution to offline-only and propagate the explicit network flag to injected/default resolvers.
- [x] Support disabled verification without resolver/canonicalization work; expose a small UI opt-in for external references.
- [x] Persist only non-secret identity preferences and explain offline-unverified results.
- [x] Test default/offline/external/disabled behavior, cached offline references, settings persistence, absence of accidental requests, and identity/M23/Clean2D boundaries.

## Credential transport

- [x] Reject non-empty API keys over remote HTTP in `OpenAICompatibleConfig`; preserve HTTPS, loopback HTTP, and unkeyed HTTP.
- [x] Test remote/LAN rejection, allowed endpoints, config repr, and QSettings secret exclusion.

## Documentation and verification

- [x] Record the future Local Chemistry Specialist Model as research only; include project continuation and prudent SIGSEGV wording.
- [x] Run bounded focused tests, compileall, architecture, strict OpenSpec validation, changed-file Ruff rules, and diff check.
- [x] Commit the work in one or two clear commits; do not merge, rebase, force-push, or squash.
- [ ] Push normally to `origin/ai/molecular-assistant-foundation` (attempted once; blocked because GitHub HTTPS credentials are unavailable in this environment).
