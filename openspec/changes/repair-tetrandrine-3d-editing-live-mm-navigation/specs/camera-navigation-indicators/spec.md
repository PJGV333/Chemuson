## ADDED Requirements

### Requirement: Viewport navigation remains legible and distinct from atom editing
ChemUSON SHALL display a persistent X/Y/Z indicator projected from the active camera basis in the lower-left viewport. Empty-space camera orbit drags SHALL show temporary directional guidance; camera motion and atom motion SHALL have distinct visual/status cues. Overlays SHALL not alter molecular picking, coordinates, or document chemistry.

#### Scenario: Camera orientation changes
- **GIVEN** any valid camera orientation
- **WHEN** the user orbits the camera
- **THEN** the X/Y/Z indicator rotates consistently with the camera basis and remains in the lower-left corner

#### Scenario: Orbit in empty space
- **GIVEN** the pointer begins a camera orbit in empty viewport space
- **WHEN** the user drags
- **THEN** temporary orbit guidance and a camera-motion cue appear
- **AND** the cue disappears when the drag ends

#### Scenario: Move an atom
- **GIVEN** 3D edit mode with the move tool and a picked atom
- **WHEN** the user drags the atom
- **THEN** the atom geometry preview is shown without camera-orbit guidance
- **AND** the camera axes and camera pose remain unchanged
