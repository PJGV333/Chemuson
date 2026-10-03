# Repository Hygiene Specification

## ADDED Requirements

### Requirement: Historical campaign knowledge SHALL remain recoverable without retaining obsolete implementation artifacts

The project SHALL maintain a concise campaign history with objectives, relevant branch/head identifiers, strategies, affected areas, outcomes, what worked/failed, decisions, lessons and a safe restart point. Before deleting a historically significant branch or artifact, its useful technical context SHALL be captured there. Git history alone SHALL NOT substitute for project documentation.

#### Scenario: Campaign is consolidated
- **GIVEN** a completed, superseded or abandoned campaign
- **WHEN** its branch or redundant evidence is considered for deletion
- **THEN** the campaign’s intent, relevant SHAs, result and lessons are documented
- **AND** canonical contracts and uniquely useful evidence remain available.

### Requirement: Cleanup SHALL remove only proven-unused artifacts and preserve product contracts

Cleanup SHALL distinguish SAFE_DELETE, KEEP and NEEDS_REVIEW candidates. It SHALL check static and dynamic consumers, entry points, plugins, Qt actions/signals, APIs, assets, packaging, tests and documentation. It SHALL NOT change chemistry, Clean2D, canvas/editor behavior, CMSN, templates, exports or public APIs, and SHALL NOT delete current regression tests to reduce size.

#### Scenario: Candidate has no confirmed consumer
- **GIVEN** a file or symbol with no textual references
- **WHEN** it is evaluated for removal
- **THEN** it is deleted only after structural consumer checks prove it is not dynamically loaded or contractually/publicly used
- **AND** ambiguous candidates are kept and reported.

### Requirement: Branch retirement SHALL be ancestry-verified and gated by documentation and validation

A remote branch SHALL be deleted only after its exact HEAD is documented, its tip is an ancestor of `origin/main`, the maintenance branch is pushed, and the cleanup gates pass. `main`, `gh-pages`, active Clean2D branches and branches with unreviewed unique commits SHALL be retained. Git history SHALL NOT be rewritten in this campaign.

#### Scenario: Merged branch is retired
- **GIVEN** the recorded branch ref and its current remote tip
- **WHEN** `git merge-base --is-ancestor origin/<branch> origin/main` succeeds and required gates pass
- **THEN** the branch may be deleted and the action/absorbing main SHA is recorded
- **AND** all unique, active, unknown and protected branches remain.

### Requirement: Repository cleanup SHALL demonstrate functional equivalence

Cleanup SHALL compare full-test failure identity with the pre-clean baseline, run architecture/UI/package-relevant checks and OpenSpec validation, and report size changes. New functional failures or missing packaged assets SHALL block branch retirement. No Git history rewrite tool or force push SHALL be used.

#### Scenario: Validation and history-size decision
- **GIVEN** the post-clean tree
- **WHEN** validation and object accounting complete
- **THEN** compileall, architecture, targeted UI, OpenSpec, diff and Qt smoke results are recorded
- **AND** any full-suite failure is compared with baseline
- **AND** potential filter-repo savings are only proposed, never applied.
