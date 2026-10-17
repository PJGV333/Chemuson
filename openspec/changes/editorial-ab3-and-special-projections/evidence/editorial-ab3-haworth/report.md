# Depiction Visual Lab — reproduction report

- ChemUSON version: `0.3.0-beta.1`
- Renderer source SHA at export: `27f39ce669bf794db1bad96c0a404c27885bb0ff`
- Git working tree at export: `clean`
- Platform: `Linux-7.2.10-1-cachyos-x86_64-with-glibc2.44`
- Qt / PyQt: `6.11.0` / `6.11.0`
- Case: `C3` — Haworth sugar template (`template-depiction`)
- Editorial calibration eligibility: `excluded` — outside the approved AB3 editorial calibration cohort (A1, A2, A3, B1, B2, B3)
- Approved selection cohort: `A1, A2, A3, B1, B2, B3`. B4 remains a geometric fixture and is excluded from editorial preset selection/style calibration.
- Theme: `light`; panel: `1200×820`; antialiasing: on
- Profile: `working-profile.json`; schema `chemuson.depiction-style` v1; based on `owner-approved valid cases A1, A2, A3, B1, B2, B3; A4 and B4 excluded; experimental only; not a production default`
- SVG: exported through the same QGraphicsScene renderer

## Working profile

```json
{
  "based_on": "owner-approved valid cases A1, A2, A3, B1, B2, B3; A4 and B4 excluded; experimental only; not a production default",
  "name": "ChemUSON Editorial AB3",
  "parameters": {
    "bold_double_bond_offset_ratio": null,
    "bond_stroke_ratio": 0.075,
    "cap_style": "round",
    "double_bond_inner_trim_ratio": 0.175,
    "double_bond_offset_ratio": 0.15,
    "hashed_wedge_count": 9,
    "hashed_wedge_max_width_ratio": 0.37,
    "hashed_wedge_min_width_ratio": 0.05,
    "hashed_wedge_stroke_ratio": 0.055,
    "hashed_wedge_terminal_shift_ratio": null,
    "join_style": "round",
    "label_clearance_ratio": 0.02,
    "solid_wedge_tip_width_ratio": null,
    "solid_wedge_width_ratio": 0.45,
    "triple_bond_offset_ratio": null
  },
  "schema": "chemuson.depiction-style",
  "version": 1
}
```

## Changed values (Current → Working)

- `bond_stroke_ratio`: `0.065` → `0.075`
- `double_bond_offset_ratio`: `0.105` → `0.15`
- `solid_wedge_width_ratio`: `0.225` → `0.45`
- `hashed_wedge_count`: `6` → `9`
- `hashed_wedge_max_width_ratio`: `0.25` → `0.37`
- `label_clearance_ratio`: `None` → `0.02`

## Artifact SHA256

- `current-C3-light.png` — SHA256 `ba9782f1641713a9496d4aa42b75a9120f57da920704c866f77be348962fe0ad`
- `working-C3-light.png` — SHA256 `35fe1cad4dfd823cbfede0a6814fa58af97c419334aac27f35ee95546e356005`
- `comparison-C3-light.png` — SHA256 `7ecf1012f47daf5dad5c8386819440fa58e063f2125109cef5328725570df776`
- `current-C3-light.svg` — SHA256 `338bca832d7b136b67a388ddfb322076873eceb20d472a7df0b61bba0b5eef3e`
- `working-C3-light.svg` — SHA256 `c3cad6a3a162f006cbd597273cba2aa87e203014d4065333f93298cfed8ecbd9`
- `working-profile.json` — SHA256 `f47643da6e651d5bd322e21926efa61a04fceae515169905a2052cde6153d39d`
