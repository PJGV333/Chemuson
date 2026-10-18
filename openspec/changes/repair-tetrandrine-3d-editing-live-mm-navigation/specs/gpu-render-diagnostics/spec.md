## ADDED Requirements

### Requirement: Viewport reports completed GPU work truthfully
The 3D viewport SHALL distinguish scene construction, CPU mesh generation, GPU buffer upload, and draw submission. A normal return from `render()` or a truthy upload result SHALL NOT be presented as proof that a draw was submitted or that pixels are visible. `render_requested` SHALL NOT be used as a success state.

#### Scenario: Renderer exits before drawing
- **GIVEN** a missing or mismatched context, a zero-sized viewport, an empty mesh, or another early-return condition
- **WHEN** the viewport paints
- **THEN** diagnostics identify the failed or skipped stage and reason
- **AND** the status does not claim a successful draw

#### Scenario: Non-empty scene is uploaded and drawn
- **GIVEN** a non-empty finite scene and a valid current OpenGL context
- **WHEN** upload and draw complete without GL errors
- **THEN** diagnostics report mesh vertices, uploaded byte count, GPU buffer size/state, VAO and shader state, GL context identity/version/vendor/renderer, draw-call count, and GL errors
- **AND** the render state means only that drawing was submitted, not that visible pixels were proven

#### Scenario: Empty or invalid mesh
- **GIVEN** a scene has atoms but mesh generation produces no drawable vertices
- **WHEN** the viewport paints
- **THEN** it reports an explicit empty-mesh diagnostic instead of a successful upload or draw

### Requirement: Camera and projection diagnostics explain molecular visibility
For a rendered scene, the viewport SHALL expose the active view and projection matrices, viewport dimensions, scene center and bounds, camera target and distance, near/far clipping planes, and the projected screen bounds and count of visible atom centers. Diagnostics SHALL make clipping, off-screen projection, and invalid viewport dimensions distinguishable from chemistry and mesh-generation failures.

#### Scenario: Scene is outside the camera frustum
- **GIVEN** a valid scene whose atom centers cannot be projected inside the active clipping range
- **WHEN** the viewport paints
- **THEN** diagnostics report zero projected atoms and the projection/clipping data
- **AND** do not label the scene as visibly rendered

### Requirement: Valid molecular geometry is visible in the real viewport
A finite non-empty molecular scene that passes mesh generation, GPU upload, and draw submission SHALL produce molecular pixels in the 3D viewport, separate from painter overlays. The scene SHALL remain centered after fit-to-view and remain visible after camera orbit and dock reparenting/context recreation.

#### Scenario: Synthetic and corrected tetrandrina scenes
- **GIVEN** a minimal synthetic molecule and the locally corrected tetrandrina graph with a valid finite conformer
- **WHEN** each scene is published through the viewport pipeline
- **THEN** molecular pixels are present in a real Qt/OpenGL framebuffer readback
- **AND** tetrandrina can be fit and rotated without changing its chemical connectivity
- **AND** a non-converged minimization result with finite coordinates remains displayable

#### Scenario: Dock context lifecycle
- **GIVEN** a visible molecular viewport with a rendered scene
- **WHEN** its dock is floated, re-docked, or resized
- **THEN** the owning context and GPU resources are safely recreated or retained
- **AND** the same scene is uploaded and rendered again without a false success status
