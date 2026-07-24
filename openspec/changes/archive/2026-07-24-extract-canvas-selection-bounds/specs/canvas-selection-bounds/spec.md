# Canvas Selection Bounds Specification

## Purpose

Define `src/chemuson/gui/canvas/selection_bounds.py` as the internal query module
for resolving selected atom IDs and computing visual selection bounds.

## ADDED Requirements

### Requirement: Module Provides resolve_selected_atom_ids

The module SHALL expose a function `resolve_selected_atom_ids(selected_atom_ids, selected_bond_ids, model)` that returns a `set[int]` of atom IDs.

#### Scenario: atoms only
- **GIVEN** `selected_atom_ids=[1, 2, 3]` and empty `selected_bond_ids`
- **WHEN** `resolve_selected_atom_ids` is called
- **THEN** the result is `{1, 2, 3}`

#### Scenario: via bonds
- **GIVEN** empty `selected_atom_ids` and `selected_bond_ids=[10]` where bond 10 connects atoms 4 and 5
- **WHEN** `resolve_selected_atom_ids` is called
- **THEN** the result is `{4, 5}`

#### Scenario: missing bond
- **GIVEN** `selected_atom_ids=[1]` and `selected_bond_ids=[999]` where bond 999 does not exist
- **WHEN** `resolve_selected_atom_ids` is called
- **THEN** the result is `{1}`

#### Scenario: no mutation
- **GIVEN** `selected_atom_ids=[1, 2]` and `selected_bond_ids=[1]`
- **WHEN** `resolve_selected_atom_ids` is called
- **THEN** the input lists are not mutated

#### Scenario: empty
- **GIVEN** empty inputs
- **WHEN** `resolve_selected_atom_ids` is called
- **THEN** the result is `set()`

### Requirement: Module Provides selection_bounds

The module SHALL expose a function `selection_bounds(scene, atom_items, bond_items, implicit_h_overlays, atom_ids, bond_ids, graphic_items)` that returns `QRectF | None`.

#### Scenario: empty returns None
- **GIVEN** empty atom IDs, bond IDs, and graphic items
- **WHEN** `selection_bounds` is called
- **THEN** the result is `None`

#### Scenario: invalid rects ignored
- **GIVEN** an atom with a null body QRectF and no visible label or charge_label
- **WHEN** `selection_bounds` is called
- **THEN** the result is `None`

#### Scenario: union of multiple rects
- **GIVEN** two atoms at `(0,0,10,10)` and `(20,0,10,10)`
- **WHEN** `selection_bounds` is called
- **THEN** the result has `x=0, y=0, width=30, height=10`

#### Scenario: atom visible body
- **GIVEN** an atom with pen/brush visible and body rect `(5,5,15,15)`
- **WHEN** `selection_bounds` is called
- **THEN** the result `x` is `5`

#### Scenario: atom hidden no pen no brush
- **GIVEN** an atom with `NoPen` and `NoBrush` and no visible label or charge_label
- **WHEN** `selection_bounds` is called
- **THEN** the result is `None`

#### Scenario: label visible included
- **GIVEN** an atom with visible label at `(10,0,5,5)`
- **WHEN** `selection_bounds` is called
- **THEN** the result `width` is at least `15`

#### Scenario: charge_label visible included
- **GIVEN** an atom with visible charge_label at `(0,10,5,5)`
- **WHEN** `selection_bounds` is called
- **THEN** the result `height` is at least `15`

#### Scenario: implicit hydrogen visible
- **GIVEN** an atom with visible implicit H overlays at `(20,20,5,5)` and `(25,25,5,5)`
- **WHEN** `selection_bounds` is called
- **THEN** the result includes the overlay rects

#### Scenario: overlays hidden or removed
- **GIVEN** overlays that are removed from scene or hidden
- **WHEN** `selection_bounds` is called
- **THEN** only the atom body rect is included

#### Scenario: bond item
- **GIVEN** a bond item at `(30,0,20,5)`
- **WHEN** `selection_bounds` is called
- **THEN** the result `x` is `30`

#### Scenario: graphic item
- **GIVEN** a graphic item at `(50,50,10,10)`
- **WHEN** `selection_bounds` is called
- **THEN** the result `x` is `50`

#### Scenario: graphic item removed from scene
- **GIVEN** a graphic item that is not in the scene
- **WHEN** `selection_bounds` is called
- **THEN** the result is `None`

#### Scenario: combination
- **GIVEN** an atom at `(0,0,10,10)`, a bond at `(20,0,10,10)`, and a graphic item at `(40,0,10,10)`
- **WHEN** `selection_bounds` is called
- **THEN** the result has `x=0, width=50`

#### Scenario: null rect ignored
- **GIVEN** a graphic item returning a null QRectF
- **WHEN** `selection_bounds` is called
- **THEN** the result is `None`

### Requirement: CanvasSelectionMixin Wrappers Delegate to Module

`CanvasSelectionMixin` SHALL retain the private methods `_selected_atom_ids_for_transform`, `_selected_items_bbox`, and `_targets_bbox` as thin delegates to the new module functions.

#### Scenario: _selected_atom_ids_for_transform delegates
- **GIVEN** `CanvasSelectionMixin` with `self.state.selected_atoms=[1,2]` and `self.state.selected_bonds=[]`
- **WHEN** `_selected_atom_ids_for_transform` is called
- **THEN** the result equals `resolve_selected_atom_ids(selected_atom_ids=[1,2], selected_bond_ids=[], model=self.model)`

#### Scenario: _selected_items_bbox delegates
- **GIVEN** `CanvasSelectionMixin` with selected atoms and bonds
- **WHEN** `_selected_items_bbox` is called
- **THEN** the result equals `selection_bounds` called with the atom IDs from `_selected_atom_ids_for_transform` and selected graphic items

#### Scenario: _targets_bbox delegates
- **GIVEN** `CanvasSelectionMixin` with explicit item groups
- **WHEN** `_targets_bbox` is called
- **THEN** the result equals `selection_bounds` called with concatenated item groups
