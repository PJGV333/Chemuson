from __future__ import annotations

import argparse
from dataclasses import asdict
import json
import os
from pathlib import Path
import subprocess
import sys

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "src"))
sys.path.insert(0, str(ROOT))
os.environ.setdefault("QT_QPA_PLATFORM", "offscreen")

from PyQt6.QtGui import QColor, QFont, QImage, QPainter
from PyQt6.QtWidgets import QApplication

from chemuson.core.model import BondStereo
from chemuson.gui.projections import (
    ChairParameters,
    ChairSubstituent,
    FischerParameters,
    HaworthParameters,
    HaworthSubstituent,
    build_chair_geometry,
    build_fischer_geometry,
    build_haworth_geometry,
)
from tools.visual.depiction_lab_app.renderer import ProductionCanvasAdapter


def _source_sha() -> str:
    result = subprocess.run(
        ["git", "rev-parse", "HEAD"], cwd=ROOT, capture_output=True, text=True, check=True
    )
    return result.stdout.strip()


def main() -> int:
    parser = argparse.ArgumentParser(description="Render typed projection examples through ChemusonCanvas")
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    output = args.output.resolve()
    output.mkdir(parents=True, exist_ok=True)
    app = QApplication.instance() or QApplication([])

    chair_parameters = ChairParameters(
        substituents=(ChairSubstituent(0, "OH", "up", "axial"),)
    )
    examples = (
        (
            "Haworth · β-D-glucopyranose",
            build_haworth_geometry(HaworthParameters()),
        ),
        (
            "Haworth · furanose",
            build_haworth_geometry(
                HaworthParameters(
                    ring_size=5,
                    sugar_family="generic",
                    substituents=(
                        HaworthSubstituent(1, "OH", "up"),
                        HaworthSubstituent(2, "OH", "down"),
                        HaworthSubstituent(3, "OH", "up"),
                        HaworthSubstituent(4, "OH", "down"),
                    ),
                )
            ),
        ),
        (
            "Fischer · horizontal toward / vertical away",
            build_fischer_geometry(
                FischerParameters(substituent_sides=("right", "left", "right", "left"))
            ),
        ),
        ("Chair A · OH axial up", build_chair_geometry(chair_parameters)),
        ("Chair B · same OH equatorial up", build_chair_geometry(chair_parameters.flipped())),
    )

    columns = 2
    card_width, card_height = 680, 450
    rows = (len(examples) + columns - 1) // columns
    sheet = QImage(
        columns * card_width,
        rows * card_height,
        QImage.Format.Format_ARGB32_Premultiplied,
    )
    sheet.fill(QColor("#E2E8F0"))
    painter = QPainter(sheet)
    painter.setRenderHint(QPainter.RenderHint.Antialiasing, True)
    painter.setRenderHint(QPainter.RenderHint.TextAntialiasing, True)
    metadata = []
    for index, (title, geometry) in enumerate(examples):
        instance = geometry.to_projection_instance()
        if any(bond.stereo != BondStereo.NONE for bond in instance.render_graph.bonds.values()):
            raise AssertionError(f"projection rendering invented bond stereo: {title}")
        adapter = ProductionCanvasAdapter.create(instance.render_graph, theme="light")
        try:
            panel = adapter.render_image(width=640, height=400)
        finally:
            adapter.close()
            app.processEvents()
        row, column = divmod(index, columns)
        x, y = column * card_width, row * card_height
        painter.fillRect(x + 8, y + 8, card_width - 16, card_height - 16, QColor("#F8FAFC"))
        painter.setPen(QColor("#0F172A"))
        painter.setFont(QFont("DejaVu Sans", 12, QFont.Weight.DemiBold))
        painter.drawText(x + 18, y + 30, title)
        painter.setFont(QFont("DejaVu Sans", 8))
        painter.setPen(QColor("#475569"))
        painter.drawText(
            x + 18,
            y + 48,
            f"{len(instance.render_graph.atoms)} atoms · {len(instance.render_graph.bonds)} bonds · rendering-only MolGraph; semantics retained in ProjectionInstance",
        )
        painter.drawImage(x + 18, y + 56, panel)
        metadata.append(
            {
                "title": title,
                "family": geometry.family,
                "semantics": asdict(geometry.semantics) if geometry.semantics else None,
                "bond_stereo_values": sorted({bond.stereo.value for bond in instance.render_graph.bonds.values()}),
            }
        )
    painter.end()

    image_path = output / "projection-geometry-gallery.png"
    if not sheet.save(str(image_path), "PNG"):
        raise OSError(f"could not save {image_path}")
    report = (
        "# Typed projection geometry gallery\n\n"
        f"- Source SHA: `{_source_sha()}`\n"
        "- Renderer: production `ChemusonCanvas`/`BondItem`/`AtomItem`, Current style.\n"
        "- Canvas panels: 640×400 px; light theme; same scene renderer for all examples.\n"
        "- Stereo safety: every graphical MolGraph bond is `BondStereo.NONE`; projection semantics are shown in the JSON metadata below and remain attached only to the in-memory `ProjectionInstance`.\n"
        "- This evidence is a reproducible geometry/rendering comparison, **not owner aesthetic acceptance**.\n\n"
        "## Examples\n\n```json\n"
        + json.dumps(metadata, ensure_ascii=False, indent=2, sort_keys=True)
        + "\n```\n"
    )
    (output / "projection-geometry-gallery.md").write_text(report, encoding="utf-8")
    print(json.dumps({"image": str(image_path), "report": str(output / "projection-geometry-gallery.md")}))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
