from __future__ import annotations

import hashlib
import json
from pathlib import Path

from chemuson.gui.diagram_presets import (
    list_atomic_presets,
    list_ligand_field_presets,
    list_molecular_orbital_presets,
)


def main() -> None:
    presets = [
        *list_atomic_presets(),
        *list_molecular_orbital_presets(),
        *list_ligand_field_presets(),
    ]
    fingerprints = {}
    for preset in presets:
        diagram = preset.build().to_json_dict()
        canonical = json.dumps(
            diagram,
            sort_keys=True,
            separators=(",", ":"),
            ensure_ascii=False,
        ).encode("utf-8")
        fingerprints[preset.name] = hashlib.sha256(canonical).hexdigest()
    output = Path(__file__).with_name("semantic-preset-fingerprints.json")
    output.write_text(
        json.dumps(fingerprints, ensure_ascii=False, indent=2) + "\n",
        encoding="utf-8",
    )


if __name__ == "__main__":
    main()
