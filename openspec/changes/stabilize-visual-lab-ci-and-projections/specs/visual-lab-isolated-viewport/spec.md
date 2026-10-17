# Visual Lab isolated viewport

## ADDED Requirements

### Requirement: Initial centering uses the molecular scene, not the paper
The Molecular Depiction Visual Lab SHALL render both panels visibly after Qt's deferred initial-centering callbacks. Lab-local view policy SHALL center on the molecular scene bounds and SHALL NOT monkey-patch the global canvas or alter normal editor paper centering.

#### Scenario: First show displays molecular content
- **GIVEN** Current and Working are constructed for a fixture
- **WHEN** the window is shown and deferred Qt events are processed
- **THEN** actual atom and bond graphics/pixels are visible in each viewport
- **AND** their displayed geometry and scale match.

### Requirement: Case and view changes preserve synchronized visibility
After case selection, zoom, rotation, parent resize and deferred events, Current and Working SHALL remain centered on the same molecular scene with equal transforms/visible geometry and visible atoms/bonds.

#### Scenario: Resize a lab window
- **GIVEN** a visible lab with matching Current/Working panes
- **WHEN** the window is resized or the view is transformed
- **THEN** both panes retain the same visible molecular content and neither recenters on the normal paper.
