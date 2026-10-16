# Depiction Visual Lab — gallery evidence

- ChemUSON version: `0.3.0-beta.1`
- Renderer source SHA at export: `d4cda5f4f20e63d4eb8594e1a65d54d8b95d8aaa`
- Git working tree at export: `modified`
- Platform: `Linux-7.2.9-2-cachyos-x86_64-with-glibc2.44`; Qt/PyQt `6.11.0` / `6.11.0`
- Theme: `light`; cases: 12; contact sheet: `contact-sheet-light.png`
- Profile: `working-profile.json`; based on `current-renderer`

## Cases

- `A1` — Solid wedge · OH label (`molecule`): Tetrahedral center with an OH-bearing wedge terminus.
- `A2` — Solid wedge · F label (`molecule`): Tetrahedral center with a compact halogen label.
- `A3` — Hashed wedge · OH label (`molecule`): Tetrahedral center with a hashed OH-bearing substituent.
- `A4` — Tetrahedral center · wedge + hash (`molecule`): Four distinct substituents; one solid and one hashed stereobond.
- `B1` — Hashed wedge + double bond (`molecule`): Tetrahedral stereocenter adjacent to a C=C sigma/pi depiction.
- `B2` — Solid wedge + double bond (`molecule`): Tetrahedral stereocenter adjacent to a C=C sigma/pi depiction.
- `B3` — Carbonyl · oxygen label (`molecule`): Acetone carbonyl with the actual O label and neighboring sigma bonds.
- `B4` — Three-way wedge / bold / triple (`geometric-fixture`): Abstract geometry only: the combined bond orders/styles are not a chemically valid molecule.
- `C1` — Cyclopentane · bold lower edges (`molecule`): Five-membered saturated ring with two visually bold edges.
- `C2` — Cyclohexane · bold lower edges (`molecule`): Six-membered saturated ring with two visually bold edges.
- `C3` — Haworth sugar template (`template-depiction`): Deterministic existing Haworth template; no additional stereo metadata is invented.
- `C4` — Aromatic ring · chiral side chain (`molecule`): Alternating aromatic ring with a tetrahedral fluorinated OH-bearing side chain.

## Modified parameters

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

## SHA256

- `contact-sheet-light.png` — `3d8f9a817ab561203f41a2fed068b090d8bf3fff83943f304cc9bef67ee25dc7`
- `working-profile.json` — `ecf75b87ee50571b66ed9a83d5d7f061e99772ea511b44dcb7c064c013b2097e`
