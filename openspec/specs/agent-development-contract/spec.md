# agent-development-contract Specification

## Purpose
Define the mandatory workflow, isolation, baseline verification, reporting,
scope-control and stop criteria that human and automated agents must follow
when modifying ChemUSON.
## Requirements
### Requirement: OpenSpec Workflow Mandatory

All agents SHALL follow the OpenSpec lifecycle for any change: create Proposal → Design → Specs → Tasks before writing code. No code changes SHALL occur without an active OpenSpec change.

#### Scenario: Agent starts work
- **GIVEN** an agent receives a task to add a new feature
- **WHEN** the agent reads `AGENTS.md`
- **THEN** the agent creates an OpenSpec change before modifying code

### Requirement: Baseline Capture Before Any Modification

Before editing any file, an agent SHALL execute and record: `git status --short`, `python -m compileall src tests tools packaging`, `pytest --collect-only -q`, `pytest -q`, `ruff check src tests tools packaging --select F401,F811,F821,E722,E741`. Results SHALL be saved for post-change comparison.

#### Scenario: Agent captures baseline
- **GIVEN** the repository is in a clean state
- **WHEN** the agent is about to start implementation
- **THEN** all 5 baseline commands are executed and results recorded

### Requirement: Post-Change Verification

After completing modifications, an agent SHALL re-execute the full baseline suite and compare results. New test failures SHALL be investigated and resolved. The agent SHALL NOT proceed if existing tests regressed without justification.

#### Scenario: Regression detected
- **GIVEN** baseline showed 842 passed, 55 skipped, 6 failed
- **WHEN** post-change run shows 840 passed, 55 skipped, 8 failed
- **THEN** the agent stops and investigates the 2 new failures

### Requirement: No Opportunistic Refactoring

An agent SHALL NOT modify code outside the explicit scope of the active OpenSpec change. This includes: bug fixes, style improvements, dead code removal, reorganization, and dependency updates.

#### Scenario: Temptation to fix unrelated bug
- **GIVEN** the agent's task is to add a catalog entry
- **WHEN** the agent notices a broken test in an unrelated module
- **THEN** the agent does not fix the test and reports the observation

### Requirement: No Hiding Failures

An agent SHALL NOT modify test snapshots, baselines, exception lists, or expected outputs to make tests pass. Failures SHALL be reported, not concealed.

#### Scenario: Snapshot mismatch
- **GIVEN** a visual test fails due to pixel difference
- **WHEN** the agent could update the baseline snapshot
- **THEN** the agent reports the failure instead of updating the snapshot

### Requirement: Agent Report for Autonomous Work

When working autonomously, an agent SHALL maintain `AGENT_REPORT.md` documenting: actions taken, files modified, rationale for decisions, deviations from plan, and final verification results.

#### Scenario: Multi-step autonomous task
- **GIVEN** the agent is executing a 10-step task
- **WHEN** the agent completes each step
- **THEN** the step is logged in `AGENT_REPORT.md` with outcome

### Requirement: Modified Files Listed

Every commit, report, or handoff SHALL include a complete list of modified files with brief description per file.

#### Scenario: Agent completes work
- **GIVEN** the agent modified 5 files
- **WHEN** the agent produces its final report
- **THEN** all 5 files are listed with descriptions

### Requirement: Branch or Worktree Isolation

An agent SHALL work in a dedicated branch or worktree. Only one writer agent SHALL operate per branch simultaneously.

#### Scenario: Concurrent agents
- **GIVEN** two agents are assigned different tasks
- **WHEN** both start work
- **THEN** each uses a separate branch

### Requirement: Stop Criteria

An agent SHALL stop and report if any of the following occur:
- Existing tests that were passing begin to fail without clear cause.
- The task requires moving or renaming production code (outside scope).
- New circular dependencies are introduced.
- The catalog reveals undocumented dependencies that block the task.
- The agent cannot resolve a contradiction in the specs.

#### Scenario: Tests regress unexpectedly
- **GIVEN** the agent's changes only affect documentation
- **WHEN** a Clean2D test that was previously stable now fails differently
- **THEN** the agent stops and reports

### Requirement: Protected Subsystems

The following subsystems require extra caution. Changes touching these areas SHALL be reviewed against the protection rules:
- **GUI (`gui/`)**: do not alter signals, event order, initialization sequence, or widget hierarchy.
- **Clean2D (`clean2d/`)**: do not modify heuristics or geometry without first resolving baseline failures.
- **ChemName (`chemname/`)**: do not change naming rules without acceptance test coverage.
- **Persistence (`chemio/persistence.py`)**: do not change serialization formats.

#### Scenario: Agent needs to touch GUI
- **GIVEN** the agent's task involves documenting GUI imports
- **WHEN** the agent considers modifying a GUI signal
- **THEN** the agent stops and verifies the change is within scope

### Requirement: AGENTS.md Root File Exists

The file `AGENTS.md` SHALL exist at the repository root and contain all rules defined in this specification.

#### Scenario: New agent onboarded
- **GIVEN** a new agent accesses the repository
- **WHEN** the agent looks for behavioral guidelines
- **THEN** `AGENTS.md` exists at root with complete rules
