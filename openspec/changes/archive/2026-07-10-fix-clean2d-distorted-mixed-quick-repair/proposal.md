# Proposal: Repair distorted mixed layouts safely in Quick Clean2D

## Problem

`run_clean_2d(..., mode="quick")` can leave a distorted mixed structure unchanged when its initially selected repair is rejected at application time. The remaining layout can retain structurally unsafe individual bond lengths.

## Scope

Require each structural bond in Quick and Publication repair candidates for a `needs_rebuild` baseline to remain between 70% and 145% of that bond's order-specific desired length. Rank only candidates that preserve stereo and satisfy existing safety invariants. Tighten the `block_layout` exception to the same requirements.

## Non-goals

Do not change direct `PROPOSE` behavior. It must continue rejecting non-good baselines, including `needs_polish`.
