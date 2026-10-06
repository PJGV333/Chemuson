# Provider runtime controls and UX

## Why

Manual local-model sessions showed that the existing fixed 60-second provider timeout can reject useful complex-structure requests, the llama.cpp default endpoint (8080) differs from a user's server port, and the dialog offers no persistent profile preferences or elapsed-time feedback.

## Changes

Add configurable per-profile timeout and provider-neutral output token limits, persist endpoint/model/runtime settings but never API keys, show elapsed time while requests run, and translate stable internal error codes into Spanish user messages. Keep profiles editable and retain the 8080 global llama.cpp default.

## Boundaries

Use the existing OpenAI-compatible Chat Completions adapter and QSettings infrastructure. No vendor-specific request fields, external dependencies, real cloud calls, or secret persistence.
