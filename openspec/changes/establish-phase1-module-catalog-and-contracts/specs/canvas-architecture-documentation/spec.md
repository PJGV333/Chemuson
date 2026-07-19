# Spec: Canvas Architecture Documentation

## ADDED Requirements

### Requirement: Canvas Mixin Hierarchy Documented

The documentation SHALL specify the complete mixin inheritance chain of `ChemusonCanvas`:
- `ChemusonCanvas` extends `QGraphicsView` via 5 top-level mixins: `CanvasInputMixin`, `CanvasSelectionMixin`, `CanvasTextMixin`, `CanvasRenderMixin`, `CanvasStructureMixin`.
- `CanvasInputMixin` aggregates 6 sub-mixins: `CanvasSelectionInputMixin`, `CanvasKeyboardMixin`, `CanvasContextMenuMixin`, `CanvasToolsAnnotationsMixin`, `CanvasToolsBondingMixin`, `CanvasToolsRingsChainsMixin`.
- `CanvasToolsBondingMixin` aggregates 5 sub-mixins: bond state, model operations, drag, hit testing, and geometry.

#### Scenario: Developer traces event flow
- **GIVEN** a mouse press event occurs on the canvas
- **WHEN** the developer reads the mixin hierarchy in documentation
- **THEN** they can trace the event through `CanvasSelectionInputMixin` → `CanvasInputMixin` → `ChemusonCanvas`

### Requirement: File Inventory with Responsibilities

The documentation SHALL list all ~20 files in `src/chemuson/gui/canvas/` with a one-line responsibility description per file, organized by category: core view, input, bonding, rendering, structure, text, tools, data.

#### Scenario: Developer finds relevant file
- **GIVEN** the developer needs to understand bond hit testing
- **WHEN** they consult the file inventory
- **THEN** they find `canvas_bond_hit_testing.py` with its responsibility

### Requirement: Shared State Documented

The documentation SHALL identify key shared instance variables across mixins: `model` (MolGraph), `scene` (QGraphicsScene), `state`, `drawing_style`, drag/selection/rotation tracking variables, and any mutable state that crosses mixin boundaries.

#### Scenario: Developer avoids state conflict
- **GIVEN** two mixins both modify `self._drag_target`
- **WHEN** the developer reads the shared state section
- **THEN** they understand the interaction before modifying either mixin

### Requirement: Event Flow Documented

The documentation SHALL describe:
- Mouse event handling chain (press, move, release, double-click) through input mixins.
- Keyboard event handling through `CanvasKeyboardMixin`.
- Context menu invocation through `CanvasContextMenuMixin`.
- Tool switching via `set_current_tool` and its effect on event routing.

#### Scenario: Developer fixes click behavior
- **GIVEN** a click is not registering on bonds
- **WHEN** the developer reads the event flow documentation
- **THEN** they can trace the click from QGraphicsView through the mixin chain

### Requirement: Selection and Hit Testing Documented

The documentation SHALL cover:
- Selection model: how items are selected, deselected, and multi-selected.
- Hit testing: how `CanvasSelectionInputMixin` and bonding mixins determine which item was clicked.
- Fragment pivot rotation: how selected fragments are rotated around a pivot point.
- Bond-specific hit testing from `CanvasBondHitTestingMixin`.

#### Scenario: Developer adds new hit target
- **GIVEN** a new item type needs click detection
- **WHEN** the developer reads the hit testing section
- **THEN** they know which mixin to extend and what methods to implement

### Requirement: Undo/Redo Integration Documented

The documentation SHALL describe how the canvas integrates with `QUndoStack`:
- Commands from `gui/commands/` are pushed to the stack.
- Which canvas actions produce undoable commands vs. direct mutations.
- The relationship between canvas state changes and command execution.

#### Scenario: Developer adds undoable action
- **GIVEN** a new canvas tool needs undo support
- **WHEN** the developer reads the undo/redo section
- **THEN** they know to create a command class in `gui/commands/` and push it to the stack

### Requirement: Model-Scene Synchronization Documented

The documentation SHALL explain how `CanvasStructureMixin` maintains consistency between `MolGraph` (chemical model) and `QGraphicsScene` (visual representation), including:
- How structure changes propagate to visual items.
- How visual edits propagate back to the model.
- Import/export of structures from the canvas.

#### Scenario: Developer fixes sync bug
- **GIVEN** an atom deletion in the model doesn't remove the visual item
- **WHEN** the developer reads the synchronization section
- **THEN** they can locate the sync method in `CanvasStructureMixin`

### Requirement: Clipboard Operations Documented

The documentation SHALL cover copy/paste of:
- Chemical structures (selection → clipboard → new selection).
- Images and annotations.
- The serialization format used for clipboard data.

#### Scenario: Developer fixes paste behavior
- **GIVEN** pasted structures appear at wrong coordinates
- **WHEN** the developer reads the clipboard section
- **THEN** they understand the coordinate transformation involved

### Requirement: Text and Editing Documented

The documentation SHALL describe `CanvasTextMixin`:
- Rich text editing cycle.
- Annotation creation and formatting.
- Focus management between text editing and canvas interaction.

#### Scenario: Developer fixes text focus issue
- **GIVEN** text editing steals focus unexpectedly
- **WHEN** the developer reads the text editing section
- **THEN** they understand the focus cycle and can fix it

### Requirement: Render and Export Documented

The documentation SHALL describe `CanvasRenderMixin`:
- PNG/SVG/PDF export capabilities.
- Paper bounds and margins.
- Overlay rendering (grid, rulers).
- Opacity handling.

#### Scenario: Developer adds new export format
- **GIVEN** a new export format is needed
- **WHEN** the developer reads the render section
- **THEN** they know which methods to extend in `CanvasRenderMixin`

### Requirement: Serialization Documented

The documentation SHALL cover serialization paths available from the canvas:
- Molfile export (`molgraph_to_molfile`).
- SMILES export (`molgraph_to_smiles`).
- CMSN persistence format.
- The relationship between canvas state and serialized output.

#### Scenario: Developer traces serialization
- **GIVEN** exported molfile has incorrect bond orders
- **WHEN** the developer reads the serialization section
- **THEN** they can trace from canvas → MolGraph → molfile

### Requirement: Internal Dependencies Mapped

The documentation SHALL list all imports from outside `canvas/`:
- From `core`: `MolGraph`, `ChemState`, `Bond`, `BondStyle`.
- From `chemio`: `molgraph_to_molfile`, `molgraph_to_smiles`, `rdkit_io` functions.
- From `chemname`: `iupac_name`, `MolView`, `find_rings_simple`.
- From `gui`: `items.*`, `commands.*`, `dialogs.*`, `style.*`, `geom.*`, `numbering.*`, etc.
- From PyQt6: all Qt dependencies.

#### Scenario: Developer assesses extraction risk
- **GIVEN** a mixin imports heavily from `gui/items.py`
- **WHEN** the developer reads the dependencies section
- **THEN** they know extraction would require moving those imports

### Requirement: Risk Zones Identified

The documentation SHALL flag high-risk areas:
- `items.py` (5704 lines): geometric rendering of atoms, bonds, wedges.
- `canvas_selection.py` (5154 lines): selection, transformation, clipboard.
- `canvas_structure.py` (3875 lines): model-scene sync, import/export, analysis.
- `canvas_view.py` (1035 lines): canvas composition, shared state initialization.
- Files with intertwined model and view logic.

#### Scenario: Developer assesses change risk
- **GIVEN** the developer wants to modify wedge rendering
- **WHEN** they see `items.py` flagged as high risk
- **THEN** they know to create visual regression tests first

### Requirement: Future Extractions Mapped

The documentation SHALL propose:
- **Safe extractions**: pure calculation helpers (bounds, geometry math) that have no Qt state.
- **Risky extractions**: anything touching event order, signal emission, or visual geometry without snapshot tests.
- **Do not extract**: canvas composition, event routing, model-scene sync (too entangled).

#### Scenario: Planner designs Phase 3
- **GIVEN** the team wants to extract canvas helpers
- **WHEN** they read the extraction roadmap
- **THEN** they know which helpers are safe to move and which are not

### Requirement: Document Location and Scope

The canvas documentation SHALL be located at `docs/modules/M09-canvas.md`. It SHALL cover all topics above without copying thousands of lines of source code. It SHALL NOT explain every individual method; it SHALL explain architectural relationships and boundaries.

#### Scenario: Reviewer checks document quality
- **GIVEN** the completed canvas documentation
- **WHEN** the reviewer checks for code dumps
- **THEN** the document describes relationships, not method-by-method walkthroughs
