# Proposal: Repair simple fused aromatic layouts

## Why

The fused-aromatic preserve policy correctly prevents global redraws, but leaves a collapsed simple naphthalene layout without a safe repair candidate.

## What Changes

Add a narrowly detected fused-aromatic template for an unsubstituted, neutral, stereo-free two-six-membered-ring system. It produces a deterministic two-hexagon depiction while retaining the global preserve policy for all other fused systems.
