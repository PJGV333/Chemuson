# Depiction Visual Lab — reproduction report

- ChemUSON version: `0.3.0-beta.1`
- Renderer source SHA at export: `70efde4bca68adc6ffc55a10bd18727fb6fa1cec`
- Git working tree at export: `clean`
- Platform: `Linux-7.2.10-1-cachyos-x86_64-with-glibc2.44`
- Qt / PyQt: `6.11.0` / `6.11.0`
- Case: `C3` — Haworth sugar template (`template-depiction`)
- Theme: `light`; panel: `1200×820`; antialiasing: on
- Profile: `working-profile.json`; schema `chemuson.depiction-style` v1; based on `current-renderer`
- SVG: exported through the same QGraphicsScene renderer

## Working profile

```json
{
  "based_on": "current-renderer",
  "name": "Current",
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

- `current-C3-light.png` — SHA256 `58b59de56bc1e77d5230d8f942cd9e5fb3559ecc739162fe18b1a3a7c84a0c48`
- `working-C3-light.png` — SHA256 `58b59de56bc1e77d5230d8f942cd9e5fb3559ecc739162fe18b1a3a7c84a0c48`
- `comparison-C3-light.png` — SHA256 `ddfe1cd9b569fbc86ed3194a03a5e089e81a39d5e3862d46b2fefe483c4f91a8`
- `current-C3-light.svg` — SHA256 `e0a5189ea4e54a01d79022b612856a01db9341cdadab7a3f4b480cc31ae18851`
- `working-C3-light.svg` — SHA256 `e0a5189ea4e54a01d79022b612856a01db9341cdadab7a3f4b480cc31ae18851`
- `working-profile.json` — SHA256 `9590a0ff48b1fe64e135ec9cd6c27476c9473c9e421f4c8345128a4f9e60c17e`
