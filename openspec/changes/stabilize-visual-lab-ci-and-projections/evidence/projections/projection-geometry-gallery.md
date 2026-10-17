# Typed projection geometry gallery

- Source SHA: `0f501cc70388b429d34c7bf7132ad227849a8aad`
- Renderer: production `ChemusonCanvas`/`BondItem`/`AtomItem`, Current style.
- Canvas panels: 640×400 px; light theme; same scene renderer for all examples.
- Stereo safety: every graphical MolGraph bond is `BondStereo.NONE`; projection semantics are shown in the JSON metadata below and remain attached only to the in-memory `ProjectionInstance`.
- This evidence is a reproducible geometry/rendering comparison, **not owner aesthetic acceptance**.

## Examples

```json
[
  {
    "bond_stereo_values": [
      "none"
    ],
    "family": "haworth",
    "semantics": {
      "anomer": "beta",
      "anomeric_position": 1,
      "configuration": "D",
      "front_edge_indices": [
        1,
        2,
        3
      ],
      "ring_atom_indices": [
        0,
        1,
        2,
        3,
        4,
        5
      ],
      "ring_size": 6,
      "substituents": [
        {
          "group_atom_index": 6,
          "label": "OH",
          "orientation": "up",
          "ring_atom_index": 1,
          "ring_position": 1
        },
        {
          "group_atom_index": 7,
          "label": "OH",
          "orientation": "down",
          "ring_atom_index": 2,
          "ring_position": 2
        },
        {
          "group_atom_index": 8,
          "label": "OH",
          "orientation": "up",
          "ring_atom_index": 3,
          "ring_position": 3
        },
        {
          "group_atom_index": 9,
          "label": "OH",
          "orientation": "down",
          "ring_atom_index": 4,
          "ring_position": 4
        },
        {
          "group_atom_index": 10,
          "label": "CH2OH",
          "orientation": "up",
          "ring_atom_index": 5,
          "ring_position": 5
        }
      ]
    },
    "title": "Haworth · β-D-glucopyranose"
  },
  {
    "bond_stereo_values": [
      "none"
    ],
    "family": "haworth",
    "semantics": {
      "anomer": null,
      "anomeric_position": null,
      "configuration": null,
      "front_edge_indices": [
        1,
        2,
        3
      ],
      "ring_atom_indices": [
        0,
        1,
        2,
        3,
        4
      ],
      "ring_size": 5,
      "substituents": [
        {
          "group_atom_index": 5,
          "label": "OH",
          "orientation": "up",
          "ring_atom_index": 1,
          "ring_position": 1
        },
        {
          "group_atom_index": 6,
          "label": "OH",
          "orientation": "down",
          "ring_atom_index": 2,
          "ring_position": 2
        },
        {
          "group_atom_index": 7,
          "label": "OH",
          "orientation": "up",
          "ring_atom_index": 3,
          "ring_position": 3
        },
        {
          "group_atom_index": 8,
          "label": "OH",
          "orientation": "down",
          "ring_atom_index": 4,
          "ring_position": 4
        }
      ]
    },
    "title": "Haworth · furanose"
  },
  {
    "bond_stereo_values": [
      "none"
    ],
    "family": "fischer",
    "semantics": {
      "centers": [
        {
          "absolute_configuration": null,
          "center_atom_index": 1,
          "left_atom_index": 7,
          "lower_atom_index": 2,
          "right_atom_index": 6,
          "upper_atom_index": 0
        },
        {
          "absolute_configuration": null,
          "center_atom_index": 2,
          "left_atom_index": 8,
          "lower_atom_index": 3,
          "right_atom_index": 9,
          "upper_atom_index": 1
        },
        {
          "absolute_configuration": null,
          "center_atom_index": 3,
          "left_atom_index": 11,
          "lower_atom_index": 4,
          "right_atom_index": 10,
          "upper_atom_index": 2
        },
        {
          "absolute_configuration": null,
          "center_atom_index": 4,
          "left_atom_index": 12,
          "lower_atom_index": 5,
          "right_atom_index": 13,
          "upper_atom_index": 3
        }
      ],
      "horizontal_toward_observer": true,
      "vertical_away_from_observer": true
    },
    "title": "Fischer · horizontal toward / vertical away"
  },
  {
    "bond_stereo_values": [
      "none"
    ],
    "family": "chair",
    "semantics": {
      "conformation": "chair_a",
      "ring_atom_indices": [
        0,
        1,
        2,
        3,
        4,
        5
      ],
      "substituents": [
        {
          "group_atom_index": 6,
          "label": "OH",
          "orientation": "up",
          "placement": "axial",
          "ring_atom_index": 0,
          "ring_position": 0
        }
      ]
    },
    "title": "Chair A · OH axial up"
  },
  {
    "bond_stereo_values": [
      "none"
    ],
    "family": "chair",
    "semantics": {
      "conformation": "chair_b",
      "ring_atom_indices": [
        0,
        1,
        2,
        3,
        4,
        5
      ],
      "substituents": [
        {
          "group_atom_index": 6,
          "label": "OH",
          "orientation": "up",
          "placement": "equatorial",
          "ring_atom_index": 0,
          "ring_position": 0
        }
      ]
    },
    "title": "Chair B · same OH equatorial up"
  }
]
```
