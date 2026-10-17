# Depiction Visual Lab — reproduction report

- ChemUSON version: `0.3.0-beta.1`
- Renderer source SHA at export: `0f501cc70388b429d34c7bf7132ad227849a8aad`
- Git working tree at export: `clean`
- Platform: `Linux-7.2.10-1-cachyos-x86_64-with-glibc2.44`
- Qt / PyQt: `6.11.0` / `6.11.0`
- Case: `C3` — Haworth sugar template (`template-depiction`)
- Editorial calibration eligibility: `excluded` — outside the approved AB3 editorial calibration cohort (A1, A2, A3, B1, B2, B3)
- Approved selection cohort: `A1, A2, A3, B1, B2, B3`. B4 remains a geometric fixture and is excluded from editorial preset selection/style calibration.
- Theme: `light`; panel: `1200×820`; antialiasing: on
- Profile: `working-profile.json`; schema `chemuson.depiction-style` v1; based on `current-renderer`
- SVG: exported through the same QGraphicsScene renderer

## Working profile

```json
{
  "based_on": "current-renderer",
  "name": "Working",
  "parameters": {
    "bold_double_bond_offset_ratio": null,
    "bond_stroke_ratio": 0.065,
    "cap_style": "round",
    "double_bond_inner_trim_ratio": 0.175,
    "double_bond_offset_ratio": 0.105,
    "hashed_wedge_count": 6,
    "hashed_wedge_max_width_ratio": 0.25,
    "hashed_wedge_min_width_ratio": 0.05,
    "hashed_wedge_stroke_ratio": 0.055,
    "hashed_wedge_terminal_shift_ratio": null,
    "join_style": "round",
    "label_clearance_ratio": null,
    "solid_wedge_tip_width_ratio": null,
    "solid_wedge_width_ratio": 0.225,
    "triple_bond_offset_ratio": null
  },
  "schema": "chemuson.depiction-style",
  "version": 1
}
```

## Changed values (Current → Working)

- No differences; Working matches Current.

## Artifact SHA256

- `current-C3-light.png` — SHA256 `fc51eadf9a80851403730f58236ebb5a2bfc7d0fbb22efe76259e7288f4da522`
- `working-C3-light.png` — SHA256 `fc51eadf9a80851403730f58236ebb5a2bfc7d0fbb22efe76259e7288f4da522`
- `comparison-C3-light.png` — SHA256 `d84e8751a69f3ea27f84b2044731f9500a4c97f360b01dbb37b008d1a9baa275`
- `current-C3-light.svg` — SHA256 `2ef1b92c10065c72dfe29c48e0323331e81c000c522c1fba5d9fc55b0acf3d85`
- `working-C3-light.svg` — SHA256 `2ef1b92c10065c72dfe29c48e0323331e81c000c522c1fba5d9fc55b0acf3d85`
- `working-profile.json` — SHA256 `f1312b81c9b94e9d7bbfe0744a393c47c8f3d2bc1855bb8ad02bee8f9872ae25`
