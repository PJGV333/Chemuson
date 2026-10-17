## ADDED Requirements

### Requirement: DrawTool previews and commits geometry in 3D
The build tool SHALL show a world-space atom sphere and bond while dragging from an atom, SHALL snap the preview to an existing atom when hit, and SHALL reuse that atom on release. A release in empty space SHALL create an atom and its source bond in one undoable action. Escape SHALL cancel without changing the graph or undo history. All commits SHALL retain existing valence, element-change, explicit-H, derived-H, undo, and redo safeguards.

#### Scenario: Snap to an existing atom
- **GIVEN** two existing atoms and the build tool
- **WHEN** the user drags from one atom toward the other and releases within its pick region
- **THEN** the preview endpoint snaps to the target's 3D coordinates
- **AND** one bond is committed without adding a duplicate atom

#### Scenario: Extend into empty space
- **GIVEN** an atom selected as the source anchor
- **WHEN** the user drags from it into empty space and releases
- **THEN** a preview sphere and bond follow a camera-plane world position
- **AND** exactly one new atom plus one bond are committed as a single undo-stack entry
- **AND** one undo removes both and redo restores both

#### Scenario: Contextual double-click
- **GIVEN** the build element is set to a valid element
- **WHEN** the user double-clicks an atom
- **THEN** only that atom's element is replaced through the existing validated undoable route
- **WHEN** the user selects one atom and double-clicks empty space
- **THEN** exactly one atom is appended with a single bond to that anchor
- **AND** Qt's first click does not create an extra mutation

#### Scenario: Cancel preview
- **GIVEN** an active 3D build drag
- **WHEN** Escape is pressed before release
- **THEN** the preview disappears and the molecule, coordinates, and undo index remain unchanged
