# Design

The detector requires eight atoms, eight structural bonds, one aromatic six-cycle, two opposite ring anchors, and two allowed terminal external atoms. It excludes charge, stereo, wedge/hash, non-valence bonds, extra rings, partial targets, and all ambiguous topology.

The geometry aligns a regular aromatic hexagon to the input then places each terminal atom radially outward using its bond-order target. Candidate ranking retains all existing safety gates. `visual_score` is lower-is-better.
