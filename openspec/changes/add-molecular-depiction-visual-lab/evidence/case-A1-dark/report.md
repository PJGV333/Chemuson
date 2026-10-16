# Depiction Visual Lab — reproduction report

- ChemUSON version: `0.3.0-beta.1`
- Renderer source SHA at export: `d4cda5f4f20e63d4eb8594e1a65d54d8b95d8aaa`
- Git working tree at export: `modified`
- Platform: `Linux-7.2.9-2-cachyos-x86_64-with-glibc2.44`
- Qt / PyQt: `6.11.0` / `6.11.0`
- Case: `A1` — Solid wedge · OH label (`molecule`)
- Theme: `dark`; panel: `1200×820`; antialiasing: on
- Profile: `working-profile.json`; schema `chemuson.depiction-style` v1; based on `current-renderer`
- SVG: exported through the same QGraphicsScene renderer

## Working profile

```json
{
  "based_on": "current-renderer",
  "name": "Unapproved experiment — evidence only",
  "parameters": {
    "bold_double_bond_offset_ratio": 0.22,
    "bond_stroke_ratio": 0.07,
    "cap_style": "round",
    "double_bond_inner_trim_ratio": 0.23,
    "double_bond_offset_ratio": 0.14,
    "hashed_wedge_count": 7,
    "hashed_wedge_max_width_ratio": 0.3,
    "hashed_wedge_min_width_ratio": 0.05,
    "hashed_wedge_stroke_ratio": 0.065,
    "hashed_wedge_terminal_shift_ratio": 0.35,
    "join_style": "round",
    "label_clearance_ratio": 0.055,
    "solid_wedge_tip_width_ratio": 0.14,
    "solid_wedge_width_ratio": 0.31,
    "triple_bond_offset_ratio": 0.22
  },
  "schema": "chemuson.depiction-style",
  "version": 1
}
```

## Changed values (Current → Working)

- `bond_stroke_ratio`: `0.065` → `0.07`
- `double_bond_offset_ratio`: `0.105` → `0.14`
- `double_bond_inner_trim_ratio`: `0.175` → `0.23`
- `triple_bond_offset_ratio`: `None` → `0.22`
- `bold_double_bond_offset_ratio`: `None` → `0.22`
- `solid_wedge_width_ratio`: `0.225` → `0.31`
- `solid_wedge_tip_width_ratio`: `None` → `0.14`
- `hashed_wedge_count`: `6` → `7`
- `hashed_wedge_max_width_ratio`: `0.25` → `0.3`
- `hashed_wedge_stroke_ratio`: `0.055` → `0.065`
- `hashed_wedge_terminal_shift_ratio`: `None` → `0.35`
- `label_clearance_ratio`: `None` → `0.055`

## Artifact SHA256

- `current-A1-dark.png` — SHA256 `d415e845e39391870bde714a49bb681d383a5ab5fc3d3a6b7cd0a1781a19d4c4`
- `working-A1-dark.png` — SHA256 `d1f48783364062ef29c88fa34e647a845d59d2644a5591bb77cc3514d130f52b`
- `comparison-A1-dark.png` — SHA256 `ed3b0a8038f2ce8ef301581ac0a0167b1e3bb140365958fd073f2396236475e8`
- `current-A1-dark.svg` — SHA256 `e6ba844ca5ac4e184c0a99ef65b63a861fabe9f1d432b49bdf5539c3c9d0ad84`
- `working-A1-dark.svg` — SHA256 `062fa08b08f6f6b55c51f36468f7ea25bc6a9af328b9a8cb8c2c2d7bbaf03d35`
- `working-profile.json` — SHA256 `ecf75b87ee50571b66ed9a83d5d7f061e99772ea511b44dcb7c064c013b2097e`
