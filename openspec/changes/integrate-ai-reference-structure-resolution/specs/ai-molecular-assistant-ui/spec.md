## ADDED Requirements

### Requirement: Resolution method is visible and unambiguous
The Molecular Assistant window SHALL offer AI-only, AI+reference, and reference-only modes, default to AI+reference, and make the selected structure's origin explicit. Reference-only mode SHALL not require or call a model endpoint.

#### Scenario: Method selection controls execution
- **GIVEN** the user selects one of the three resolution methods
- **WHEN** the request starts
- **THEN** provider and Name→Structure calls follow that method's contract and the preview names the actual origin.

### Requirement: Mismatch preview offers explicit candidate choices
When AI and reference identities differ, the preview SHALL display both SMILES, recommend the reference, expose explicit `Use reference`, `Use AI proposal anyway`, and `Cancel` choices, and never silently replace the AI graph. AI override retains the existing confirmation guard.

#### Scenario: User selects the reference on mismatch
- **GIVEN** a mismatch preview
- **WHEN** the user chooses `Use reference`
- **THEN** only the validated reference graph is sent through normal insertion.

#### Scenario: User explicitly overrides with AI
- **GIVEN** a mismatch preview
- **WHEN** the user chooses `Use AI proposal anyway` and confirms
- **THEN** only the original AI graph is inserted and the mismatch provenance remains visible.

### Requirement: Reference fallback has a useful non-failure preview
When AI cannot provide a usable structure but a reference is valid, the dialog SHALL explain that AI failed, distinguish PubChem/reference provenance, show the resolved name and reference SMILES, and allow reference insertion or cancellation.

#### Scenario: Generation is exhausted but tetrandrine reference exists
- **GIVEN** an empty `finish_reason=length` AI result and a valid tetrandrine reference
- **WHEN** the preview appears
- **THEN** it reports exhausted generation and offers a reference-origin insertion, not a total generation failure.

### Requirement: Existing lifecycle and transform contracts are retained
Resolution SHALL use the existing Molecular Assistant worker/QThread ownership and lifecycle cleanup. Transform mode SHALL remain AI-only, and every approved candidate SHALL use the established normal undoable insertion path.

#### Scenario: Close during AI plus reference work
- **GIVEN** the modeless assistant is closed while its existing worker is processing either stage
- **WHEN** Qt destroys the dialog or the window shuts down
- **THEN** stable job-ID cleanup suppresses late results without dereferencing deleted wrappers or introducing a new worker thread.
