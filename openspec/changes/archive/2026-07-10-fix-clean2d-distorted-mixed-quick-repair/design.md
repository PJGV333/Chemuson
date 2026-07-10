# Design: Per-bond repair safety gate

`Clean2DQualityReport` will expose the minimum and maximum ratios of actual structural-bond length to its order-specific desired length, plus the IDs below and above the allowed range. The pure safety gate will enforce the range only for Quick and Publication repair candidates when the baseline needs rebuilding.

Candidate ranking will reject stereo-changing candidates before selection. E/Z metadata remains an explicit bond invariant; endpoint orientation is not treated as a chiral-center signature. `block_layout` will no longer bypass safety or the `needs_rebuild` quality rejection. Its existing regression checks remain, and it must also pass the common safety, stereo, quality, crossing, collision, and ring-degeneracy gates.

The controller remains an application adapter; no new chemical ranking is added there. `PROPOSE` retains its early rejection for every baseline that is not `good`.
