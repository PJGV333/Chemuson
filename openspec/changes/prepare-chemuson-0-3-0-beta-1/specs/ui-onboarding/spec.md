# UI-ONBOARDING-001 — first-run onboarding geometry

## Purpose

Keep the existing three-step welcome tour aligned with the final main-window layout on Windows and Linux, including startup, resize, move, and display-scale changes.

## ADDED Requirements

### Requirement: The first-run overlay follows finalized window geometry
The onboarding overlay MUST be created only after the main window is shown and its initial child layout has settled. It MUST cover the main-window client area, map each spotlight to its live target in overlay-local coordinates, keep the card wholly inside the overlay, and preserve the geometry of the rail, canvas, and side panel. It MUST recalculate after window/target resize, move, layout, or screen/DPI changes.

#### Scenario: First launch
- **WHEN** the main window is first shown and onboarding has not been completed
- **THEN** the existing three steps appear in order: tool rail, molecular canvas, side panel
- **AND** each spotlight aligns with its corresponding visible target and the card is fully visible.

#### Scenario: Resize or display-scale change
- **WHEN** the window is resized, moved, or receives a screen/DPI change
- **THEN** the overlay continues to cover the client area and the spotlight tracks the target's final geometry
- **AND** no main component is moved or resized by the overlay.

### Requirement: Onboarding preference semantics remain stable
The setting key `ui/onboarding/completed` MUST retain its existing QSettings meaning: completing all three steps or closing with “No volver a mostrar” stores true; closing without that option leaves it unset/false and allows the tour to reappear. Closing the tour MUST hide and release the overlay without leaving a modal blocker.

#### Scenario: User closes early
- **WHEN** the user closes the tour with or without “No volver a mostrar”
- **THEN** the existing persistence behavior is preserved and the overlay no longer intercepts the main window.

## Acceptance

The geometry correction is source/test verified only. `UI-ONBOARDING-001` remains pending the owner's manual retest of a new Windows portable preview; Linux should be checked as a cross-platform consistency retest. Automated Qt scale-factor tests do not substitute for a real Windows display/DPI check.
