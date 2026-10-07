## ADDED Requirements

### Requirement: Reference resolution is application-controlled and offline by default
The Molecular Assistant orchestration SHALL distinguish model proposals from chemical reference structures and SHALL reuse `extract_requested_molecule_name()` plus the existing `resolve_name_to_structure()` service. It MUST NOT provide model browsing, arbitrary URL selection, HTTP tools, PubChem credentials, or reference payloads to the model. External access SHALL default off; offline mode MAY use static entries and the existing local PubChem cache. PubChem PUG REST requests SHALL use the supported `SMILES`, `ConnectivitySMILES`, and `IUPACName` properties, preferring stereo-capable `SMILES` and falling back to connectivity only when necessary; legacy SMILES property names MAY remain parser fallbacks.

#### Scenario: Explicit named request with external access disabled
- **GIVEN** the prompt explicitly names one molecule and reference mode is selected
- **WHEN** the orchestration resolves it with external permission off
- **THEN** the resolver receives only the extracted name and `allow_network=false`
- **AND** static and cached references remain available without a network request.

#### Scenario: External reference permission is enabled
- **GIVEN** the user explicitly enables PubChem reference access
- **WHEN** an eligible named request is resolved
- **THEN** ChemUSON's existing resolver receives only the extracted chemical name with `allow_network=true`
- **AND** no prompt, document content, private SMILES, model reasoning, or API key is sent to PubChem.

#### Scenario: Open-ended generation does not resolve a reference
- **GIVEN** a generative description or a whole-molecule transformation rather than an explicit named-molecule request
- **WHEN** AI or AI+reference mode runs
- **THEN** the reference resolver is not called.

#### Scenario: Model has no browsing or tools
- **GIVEN** a Qwen/OpenAI-compatible generation request
- **WHEN** ChemUSON prepares the provider request
- **THEN** it uses only the explicitly configured completion endpoint and ordinary prompt/response contract
- **AND** it defines no browsing, browser, search, tool, or model-selected URL capability.

### Requirement: Verified name aliases preserve query provenance
Name→Structure MAY normalize only explicitly verified language aliases through a small deterministic alias table. The submitted name SHALL remain available as the original query, and the canonical lookup spelling SHALL be recorded separately; the resolver SHALL retain the connector source/cache provenance. Aliases SHALL NOT be inferred by an LLM or fetched from the network.

#### Scenario: Spanish tetrandrina and cholesterol names resolve through verified English aliases
- **GIVEN** PUG REST does not resolve `tetrandrina` or `colesterol`, while `tetrandrine` and `cholesterol` return usable current-property responses
- **WHEN** the user submits an explicit named-molecule request in Spanish
- **THEN** Name→Structure queries the verified canonical spelling and retains both the original query and canonical query with PubChem provenance.

### Requirement: Reference structures pass ChemIO and retain stable provenance
A reference SHALL NOT become an insertable candidate until its SMILES has passed isolated ChemIO validation. Application results and previews SHALL carry a closed structure-origin value (`ai`, `reference`, `ai_verified_by_reference`, or `ai_mismatch_reference`) and stable reference metadata including source, resolved name, and cache status when available. Provider-supplied strings SHALL NOT determine origin logic.

#### Scenario: Invalid reference is returned
- **GIVEN** a resolver result with invalid/missing SMILES or a ChemIO validation failure
- **WHEN** reference orchestration completes
- **THEN** no reference graph is available for insertion
- **AND** no raw connector exception is shown.

#### Scenario: Reference metadata is shown
- **GIVEN** a PubChem or local reference is available
- **WHEN** a preview is presented
- **THEN** it labels the structure as a reference and shows a controlled source label, resolved name, and cache indicator where applicable
- **AND** it never labels that structure as model-generated.

### Requirement: AI and reference candidates are reconciled without substitution
The orchestration SHALL compare AI and reference structures with the existing isolated, stereochemistry-aware InChI identity path, never textual SMILES equality. It SHALL keep the AI graph unchanged and require an explicit candidate choice on mismatch.

#### Scenario: AI and reference have equivalent identity
- **GIVEN** valid AI and reference SMILES that may differ textually
- **WHEN** their isolated canonical identities match
- **THEN** the unchanged AI proposal is offered as `ai_verified_by_reference` and identity is verified.

#### Scenario: AI and reference differ
- **GIVEN** valid ChemIO-accepted AI and reference structures with different isolated InChI identities
- **WHEN** the preview is shown
- **THEN** both SMILES are available for comparison
- **AND** the reference is the recommended/default candidate
- **AND** using the AI proposal requires a separate explicit override and confirmation.

#### Scenario: AI succeeds without a usable reference
- **GIVEN** a valid AI structure and a missing/unavailable/offline-only reference
- **WHEN** the result is shown
- **THEN** the AI proposal remains available with identity clearly unverified and origin `ai`.

### Requirement: Reference fallback recovers controlled AI failures
For an explicit named-molecule request in AI+reference mode, a provider failure, timeout, malformed/empty response, exhausted generation, or unusable ChemIO proposal SHALL still allow the application resolver to run. A valid ChemIO reference SHALL be reviewable as a reference-origin candidate while the AI failure remains diagnostically available.

#### Scenario: AI fails but reference exists
- **GIVEN** an AI failure and a valid resolved reference
- **WHEN** orchestration completes
- **THEN** the user sees that AI produced no usable structure, the reference is available, and AI status/reason remain bounded diagnostics
- **AND** the reference must be explicitly inserted or the flow cancelled.

#### Scenario: Both AI and reference fail
- **GIVEN** no usable AI result and no usable reference
- **WHEN** orchestration completes
- **THEN** it reports a controlled failure and creates no insertable graph.

### Requirement: Resolution modes have distinct provider semantics
The UI SHALL offer AI-only, AI+reference, and reference-only methods, defaulting to AI+reference. AI-only SHALL call only the configured model. Reference-only SHALL use Name→Structure without constructing or calling a model provider and SHALL require an explicit molecule name.

#### Scenario: Reference-only request
- **GIVEN** reference-only mode and an explicit name request
- **WHEN** the user submits
- **THEN** no model request is made and a validated reference may be previewed.

#### Scenario: Reference-only request has no molecular name
- **GIVEN** a generative/open-ended prompt in reference-only mode
- **WHEN** the user submits
- **THEN** the application reports a controlled name-required outcome without calling AI or the resolver.

### Requirement: Reference and AI insertions share normal undoable insertion
Every accepted AI or reference candidate SHALL use the existing normal canvas insertion command and one undo step. Reference insertion SHALL NOT invoke Clean2D or a parallel mutation route. Provenance SHALL remain available in preview and insertion feedback without changing the `.cmsn` format.

#### Scenario: Reference insertion undo and redo
- **GIVEN** the user accepts a validated reference candidate
- **WHEN** the canvas inserts, undoes, and redoes it
- **THEN** insertion uses the existing undoable path, Undo removes it, and Redo restores it.

### Requirement: Whole-molecule transformation stays AI-only
Reference resolution SHALL NOT run for selected-molecule transformations, which retain the current AI-only request, validation, and undo contracts.

#### Scenario: Transform request
- **GIVEN** an active whole-molecule transformation
- **WHEN** the assistant worker runs
- **THEN** it does not call a name/reference resolver regardless of prompt text.
