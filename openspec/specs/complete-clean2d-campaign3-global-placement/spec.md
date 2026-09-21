# complete-clean2d-campaign3-global-placement Specification

## Purpose
Provide topology-derived global block placement evidence and a safe internal candidate path for medium and large Clean2D multiblock layouts before local polishing.

## Requirements

### Requirement: Medium global placement is topology-derived

The Clean2D engine SHALL generate an internal global block-placement candidate for eligible medium multiblock structures before the complex `preserve-only` exit. The candidate MUST derive root, parent/child relationships, attachment directions, angular sectors, sibling occupancy, relative orientation, target separation, and child placement from topology and current coordinates, without fixture-name branches or force-directed global optimization.

#### Scenario: Multiblock medium structure receives hierarchical placement

- **WHEN** a selected graph contains connected rigid or semi-rigid blocks with attachment, bridge, or linker relationships
- **THEN** the candidate chooses a deterministic root, preserves each block's internal geometry, allocates separated child sectors, places children breadth-first, and records the plan and metrics

### Requirement: Preserve-only remains a safety fallback

The internal placement candidate MUST pass existing graph invariants and hard gates before it can compete in ranking. If it is rejected or fails to improve safely, the existing `preserve-only` behavior SHALL remain visible with a stable rejection reason. The change MUST NOT disable preserve-only globally.

#### Scenario: Difficult multiblock layout remains protected

- **WHEN** a topology-aware placement introduces a crossing, unsafe collision, ring degeneration, invalid stereo/selection state, non-finite coordinate, or excessive displacement
- **THEN** the placement candidate is rejected and the operation returns the existing controlled preserve-only result

### Requirement: Promotion evidence proves internal contribution

Corrective Campaign 3 evidence SHALL compare commits `043fe96` and the corrected HEAD for the required medium topology families. Each case MUST record baseline and current result state, selected source, candidate metrics, plus a separate internal Campaign 3 candidate record with source, before/after metrics, accepted/rejected status, reason, and selected/competitive status. The promotion gate MUST show improvement for a significant subset of the medium family caused by the internal candidate, while documenting unchanged difficult cases and simple-case behavior.

#### Scenario: Internal candidate improves a medium family member

- **WHEN** an eligible medium case is evaluated
- **THEN** its internal placement candidate records reduced hard geometric errors or better visual metrics than its own input geometry, passes the safety gates when accepted, and remains separately auditable even when an external backend is ultimately selected
