## ADDED Requirements

### Requirement: Identity verification network access is explicit and offline by default

Molecular identity verification SHALL default to offline-only reference resolution. It SHALL pass an explicit network policy to the Name→Structure resolver and SHALL use available local references and isolated canonicalization without external requests in offline-only mode.

#### Scenario: Local model request does not implicitly access the network
- **GIVEN** a successful proposal for a named molecule and no explicit external-reference permission
- **WHEN** identity verification runs
- **THEN** the resolver SHALL receive `allow_network=false`
- **AND** no PubChem or other external request SHALL be made
- **AND** a missing local reference SHALL yield `unverified`, not `verified` or `mismatch`
- **AND** the user-facing result SHALL state that external sources were not consulted.

#### Scenario: An available offline reference can verify identity
- **GIVEN** offline-only policy and a matching local reference
- **WHEN** isolated canonical InChI comparison completes
- **THEN** identity SHALL be `verified` only when both canonical identities match.

#### Scenario: External reference access is explicitly enabled
- **GIVEN** the user enabled external identity lookup
- **WHEN** identity verification invokes Name→Structure
- **THEN** the resolver SHALL receive `allow_network=true`.

#### Scenario: Identity verification is disabled
- **GIVEN** identity verification is disabled
- **WHEN** an assistant proposal is reviewed
- **THEN** the result SHALL be `not_applicable`
- **AND** no resolver or canonicalization work SHALL occur.

### Requirement: Provider credentials require a protected transport

`OpenAICompatibleConfig` SHALL reject a non-empty API key for an HTTP endpoint unless the endpoint host is unambiguously loopback (`localhost`, IPv4 127/8, or IPv6 `::1`). The policy SHALL NOT treat LAN/private addresses as trusted and SHALL NOT resolve hostnames to determine locality. HTTPS credentials and HTTP endpoints without credentials SHALL remain valid.

#### Scenario: Remote HTTPS credential transport is allowed
- **GIVEN** an HTTPS endpoint and non-empty API key
- **WHEN** provider configuration is validated
- **THEN** the configuration SHALL be accepted.

#### Scenario: Loopback HTTP is allowed
- **GIVEN** HTTP to `localhost`, IPv4 127/8, or IPv6 `::1`, with or without a key
- **WHEN** provider configuration is validated
- **THEN** the configuration SHALL be accepted.

#### Scenario: Remote plaintext credential transport is rejected
- **GIVEN** an HTTP endpoint on a remote hostname, private/LAN IP, or otherwise unrecognized host and a non-empty key
- **WHEN** provider configuration is validated
- **THEN** it SHALL fail with a generic error that does not reveal the key.

#### Scenario: Unkeyed HTTP behavior is preserved
- **GIVEN** an HTTP endpoint outside loopback and no API key
- **WHEN** provider configuration is validated
- **THEN** existing configuration behavior SHALL be preserved.

#### Scenario: Credentials are not represented or persisted
- **GIVEN** a provider config containing a key or provider preferences saved through QSettings
- **WHEN** the config is represented or preferences are persisted
- **THEN** the key SHALL not appear in the config representation or settings.
