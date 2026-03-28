"""Genera una hoja PNG con todos los orbitales usando los presets actuales."""

from __future__ import annotations

import argparse
import os
from pathlib import Path
import sys

os.environ.setdefault("QT_QPA_PLATFORM", "offscreen")

ROOT = Path(__file__).resolve().parents[1]
SRC = ROOT / "src"
if str(SRC) not in sys.path:
    sys.path.insert(0, str(SRC))

from PyQt6.QtWidgets import QApplication

from chemuson.gui.orbitals import (
    ORBITAL_FAMILY_ORDER,
    ORBITAL_PALETTE_MODEL,
    ORBITAL_PRESET_CONFIG_PATH,
    build_orbital_renderer,
    load_orbital_presets_payload,
)


DEFAULT_OUTPUT = ROOT / "tests" / "data" / "orbitals" / "palette_preview.png"
DEFAULT_TRIPTYCH_DIR = ROOT / "tests" / "data" / "orbitals" / "family_triptychs"


def _build_arg_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--config",
        type=Path,
        default=ORBITAL_PRESET_CONFIG_PATH,
        help=f"Archivo de presets a usar (default: {ORBITAL_PRESET_CONFIG_PATH})",
    )
    parser.add_argument(
        "--output",
        type=Path,
        default=DEFAULT_OUTPUT,
        help=f"PNG de salida para la hoja completa (default: {DEFAULT_OUTPUT})",
    )
    parser.add_argument(
        "--triptych-dir",
        type=Path,
        default=DEFAULT_TRIPTYCH_DIR,
        help=f"Directorio opcional para previews por familia (default: {DEFAULT_TRIPTYCH_DIR})",
    )
    parser.add_argument(
        "--skip-triptychs",
        action="store_true",
        help="No generar previews individuales outline/shaded/solid por familia.",
    )
    return parser


def main(argv: list[str] | None = None) -> int:
    args = _build_arg_parser().parse_args(argv)
    app = QApplication.instance() or QApplication([])
    del app

    payload = load_orbital_presets_payload(args.config, include_docs=True)
    renderer = build_orbital_renderer(payload)

    args.output.parent.mkdir(parents=True, exist_ok=True)
    palette = renderer.render_palette_image(ORBITAL_PALETTE_MODEL)
    if not palette.save(str(args.output)):
        raise RuntimeError(f"No se pudo guardar {args.output}")

    print(f"palette_preview={args.output}")

    if not args.skip_triptychs:
        args.triptych_dir.mkdir(parents=True, exist_ok=True)
        for family in ORBITAL_FAMILY_ORDER:
            if family not in payload:
                continue
            path = args.triptych_dir / f"{family}.png"
            image = renderer.render_style_triptych(family, panel_size=180, gap=18)
            if not image.save(str(path)):
                raise RuntimeError(f"No se pudo guardar {path}")
            print(f"triptych[{family}]={path}")

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
