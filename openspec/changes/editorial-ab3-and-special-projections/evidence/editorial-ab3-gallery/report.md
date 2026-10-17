# Depiction Visual Lab — gallery evidence

- ChemUSON version: `0.3.0-beta.1`
- Renderer source SHA at export: `27f39ce669bf794db1bad96c0a404c27885bb0ff`
- Git working tree at export: `clean`
- Platform: `Linux-7.2.10-1-cachyos-x86_64-with-glibc2.44`; Qt/PyQt `6.11.0` / `6.11.0`
- Theme: `light`; cases: 12; contact sheet: `contact-sheet-light.png`
- Profile: `working-profile.json`; based on `owner-approved valid cases A1, A2, A3, B1, B2, B3; A4 and B4 excluded; experimental only; not a production default`
- Editorial selection cohort: `A1, A2, A3, B1, B2, B3`. B4 is explicitly excluded as a geometric fixture / non-editorial stress test.

## Cases

- `A1` — Solid wedge · OH label (`molecule`): Tetrahedral center with an OH-bearing wedge terminus.
- `A2` — Solid wedge · F label (`molecule`): Tetrahedral center with a compact halogen label.
- `A3` — Hashed wedge · OH label (`molecule`): Tetrahedral center with a hashed OH-bearing substituent.
- `A4` — Tetrahedral center · wedge + hash (`molecule`): Four distinct substituents; one solid and one hashed stereobond.
- `B1` — Hashed wedge + double bond (`molecule`): Tetrahedral stereocenter adjacent to a C=C sigma/pi depiction.
- `B2` — Solid wedge + double bond (`molecule`): Tetrahedral stereocenter adjacent to a C=C sigma/pi depiction.
- `B3` — Carbonyl · oxygen label (`molecule`): Acetone carbonyl with the actual O label and neighboring sigma bonds.
- `B4` — Three-way wedge / bold / triple (`geometric-fixture`): Geometric fixture / non-editorial: invalid combined bond-order/style stress case; excluded from preset selection and style calibration.
- `C1` — Cyclopentane · bold lower edges (`molecule`): Five-membered saturated ring with two visually bold edges.
- `C2` — Cyclohexane · bold lower edges (`molecule`): Six-membered saturated ring with two visually bold edges.
- `C3` — Haworth sugar template (`template-depiction`): Deterministic Haworth projection fixture; diagnostic only and not a source of AB3 calibration values.
- `C4` — Aromatic ring · chiral side chain (`molecule`): Alternating aromatic ring with a tetrahedral fluorinated OH-bearing side chain.

## Editorial calibration eligibility

- Eligible selection cases: `A1`, `A2`, `A3`, `B1`, `B2`, `B3`
- A4 and B4 are excluded from the approved AB3 cohort. B4: `geometric fixture`; `non-editorial`; `excluded from preset selection / style calibration`.
- C-series projection/ring cases are diagnostic only and outside the approved editorial calibration cohort (A1, A2, A3, B1, B2, B3).

## Modified parameters

- `bond_stroke_ratio`: `0.065` → `0.075`
- `double_bond_offset_ratio`: `0.105` → `0.15`
- `solid_wedge_width_ratio`: `0.225` → `0.45`
- `hashed_wedge_count`: `6` → `9`
- `hashed_wedge_max_width_ratio`: `0.25` → `0.37`
- `label_clearance_ratio`: `None` → `0.02`

## SHA256

- `contact-sheet-light.png` — `9112cd481155b187a293abfbc2f65a04ad727cbb67c7b7806b019e123ca519d2`
- `working-profile.json` — `f47643da6e651d5bd322e21926efa61a04fceae515169905a2052cde6153d39d`
