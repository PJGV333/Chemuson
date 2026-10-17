## ADDED Requirements

### Requirement: Live MM frames are paced visually without falsifying scientific data
ChemUSON SHALL preserve each optimizer frame and its measured energy as an immutable snapshot while interpolating only the displayed coordinates between compatible snapshots. Interpolated samples SHALL be marked visual-only, SHALL carry no measured energy or gradient, and SHALL never be persisted as optimizer results. The viewport SHALL show force field, measured E, sample-to-sample ΔE, steps/max steps, and fixed-position constraints with their units. Editing or document/job replacement SHALL clear obsolete playback and cancel stale work immediately.

#### Scenario: Frames arrive in a burst
- **GIVEN** compatible RDKit snapshots arrive faster than the display refresh rate
- **WHEN** the viewer consumes them
- **THEN** a timed visual transition shows progressive geometry rather than coalescing directly to the last pose
- **AND** HUD E and ΔE remain the latest exact optimizer sample, identified as such
- **AND** the interpolated coordinates expose no claimed energy or gradient

#### Scenario: Incompatible snapshots
- **GIVEN** adjacent snapshots differ in atom identity, force field, or generated scene topology
- **WHEN** playback validates the transition
- **THEN** it does not interpolate across that boundary and reports/resets the visual timeline safely

#### Scenario: Edit during optimization
- **GIVEN** queued playback frames for the active document
- **WHEN** a chemical or geometry edit, tab switch, or newer job invalidates the snapshot
- **THEN** pending frames are discarded before the next visual tick and the worker is cancelled
- **AND** no stale geometry is displayed or persisted

#### Scenario: Measured energy display
- **GIVEN** RDKit MMFF94/UFF frame samples and fixed atom IDs
- **WHEN** the HUD updates
- **THEN** E and ΔE use the frame's declared energy unit, steps are dimensionless counts, and fixed constraints are labeled as atom positions
- **AND** no gradient is shown unless a distinct real gradient calculation supplies it
