# Design

The pure detector accepts only a whole ten-atom, eleven-aromatic-bond graph whose two six-membered cycles share one aromatic edge and have no external bonds, charge, or stereo. Each ring path around the shared edge maps to one side of two regular hexagons with aromatic target length. The complete template is aligned to the input before safety evaluation.

The candidate is evaluated in the complex-preserve route before the existing scaffold, unwrap, and conservative repair fallbacks. It must pass invariants, stereo, individual bond range, crossings, collision, ring-degeneracy, and final-quality checks. No complex-policy flags are relaxed.
