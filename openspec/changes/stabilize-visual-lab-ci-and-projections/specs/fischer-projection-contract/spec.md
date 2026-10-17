# Fischer projection contract

## ADDED Requirements

### Requirement: Fischer crossings encode projection conventions as metadata
Fischer geometry SHALL expose indexed vertical backbone and horizontal substituent atoms. Its typed descriptor SHALL explicitly state that horizontal bonds are toward the observer and vertical bonds recede. Graphical bonds SHALL remain ordinary drawing bonds unless an independently declared chemical stereo encoding is introduced under a separate approved contract.

#### Scenario: Preserve left/right arrangement
- **GIVEN** a Fischer center with declared left/right groups
- **WHEN** its geometry and semantics are adapted for display
- **THEN** left/right positions and horizontal/vertical viewing conventions remain available
- **AND** no wedge or E/Z metadata is invented.
