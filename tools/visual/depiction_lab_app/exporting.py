"""Deterministic PNG/SVG comparison and profile evidence export."""
from __future__ import annotations

import hashlib
import json
import platform
from pathlib import Path
import subprocess
from typing import Any

from PyQt6.QtCore import QRect, QRectF, QSize, QT_VERSION_STR, PYQT_VERSION_STR, Qt
from PyQt6.QtGui import QColor, QFont, QImage, QPainter
from PyQt6.QtWidgets import QApplication

try:
    from PyQt6.QtSvg import QSvgGenerator
except ImportError:  # pragma: no cover - depends on the installed Qt build
    QSvgGenerator = None  # type: ignore[assignment,misc]

from chemuson import __version__

from .parameters import DepictionStyleParameters
from .profiles import encode_profile, profile_payload
from .renderer import ProductionCanvasAdapter
from .scenes import CASES, DepictionCase, case_graph

PANEL_SIZE = (1200, 820)
REPOSITORY_ROOT = Path(__file__).resolve().parents[3]


def source_sha() -> str:
    try:
        completed = subprocess.run(
            ["git", "rev-parse", "HEAD"],
            check=True,
            capture_output=True,
            text=True,
            timeout=2,
            cwd=REPOSITORY_ROOT,
        )
    except (OSError, subprocess.SubprocessError):
        return "unavailable"
    return completed.stdout.strip()


def source_worktree_state() -> str:
    try:
        completed = subprocess.run(
            ["git", "status", "--porcelain", "--untracked-files=all"],
            check=True,
            capture_output=True,
            text=True,
            timeout=2,
            cwd=REPOSITORY_ROOT,
        )
    except (OSError, subprocess.SubprocessError):
        return "unavailable"
    return "clean" if not completed.stdout.strip() else "modified"


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def _save_image(image: QImage, path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    if not image.save(str(path), "PNG"):
        raise OSError(f"could not write PNG image: {path}")


def _render_svg(adapter: ProductionCanvasAdapter, path: Path, width: int, height: int) -> bool:
    if QSvgGenerator is None:
        return False
    path.parent.mkdir(parents=True, exist_ok=True)
    generator = QSvgGenerator()
    generator.setFileName(str(path))
    generator.setSize(QSize(width, height))
    generator.setViewBox(QRect(0, 0, width, height))
    generator.setTitle("ChemUSON Molecular Depiction Visual Lab")
    generator.setDescription("Production ChemusonCanvas renderer; profile and source SHA are in report.md")
    painter = QPainter(generator)
    if not painter.isActive():
        return False
    painter.setRenderHint(QPainter.RenderHint.Antialiasing, True)
    painter.setRenderHint(QPainter.RenderHint.TextAntialiasing, True)
    adapter.canvas.scene.render(
        painter,
        QRectF(0.0, 0.0, float(width), float(height)),
        adapter.scene_rect,
        Qt.AspectRatioMode.KeepAspectRatio,
    )
    painter.end()
    return path.is_file() and path.stat().st_size > 0


def render_comparison_sheet(current: QImage, working: QImage, title: str) -> QImage:
    width, height = current.width(), current.height()
    sheet = QImage(width * 2, height + 42, QImage.Format.Format_ARGB32_Premultiplied)
    sheet.fill(QColor("#E2E8F0"))
    painter = QPainter(sheet)
    painter.setRenderHint(QPainter.RenderHint.TextAntialiasing, True)
    painter.setFont(QFont("DejaVu Sans", 12, QFont.Weight.DemiBold))
    painter.setPen(QColor("#0F172A"))
    painter.drawText(16, 28, f"CURRENT · {title}")
    painter.drawText(width + 16, 28, f"WORKING · {title}")
    painter.drawImage(0, 42, current)
    painter.drawImage(width, 42, working)
    painter.end()
    return sheet


def export_case_bundle(
    output_dir: Path,
    case_id: str,
    current: DepictionStyleParameters,
    working: DepictionStyleParameters,
    *,
    theme: str = "light",
    profile_name: str = "Working",
) -> dict[str, Any]:
    """Export a matched Current/Working case and provenance bundle."""
    app = QApplication.instance()
    if app is None:
        raise RuntimeError("QApplication must exist before exporting depiction evidence")
    output = Path(output_dir)
    output.mkdir(parents=True, exist_ok=True)
    case = next((item for item in CASES if item.case_id == case_id), None)
    if case is None:
        raise KeyError(f"unknown depiction case ID: {case_id}")
    graph = case_graph(case_id)
    adapters: list[ProductionCanvasAdapter] = []
    generated: list[Path] = []
    try:
        current_adapter = ProductionCanvasAdapter.create(
            graph, parameters=current, theme=theme
        )
        adapters.append(current_adapter)
        working_adapter = ProductionCanvasAdapter.create(
            graph, parameters=working, theme=theme
        )
        adapters.append(working_adapter)
        current_png = output / f"current-{case_id}-{theme}.png"
        working_png = output / f"working-{case_id}-{theme}.png"
        comparison_png = output / f"comparison-{case_id}-{theme}.png"
        _save_image(current_adapter.render_image(width=PANEL_SIZE[0], height=PANEL_SIZE[1]), current_png)
        _save_image(working_adapter.render_image(width=PANEL_SIZE[0], height=PANEL_SIZE[1]), working_png)
        comparison = render_comparison_sheet(
            QImage(str(current_png)), QImage(str(working_png)), f"{case_id} · {case.title}"
        )
        _save_image(comparison, comparison_png)
        generated.extend((current_png, working_png, comparison_png))

        svg_status = "not available (PyQt6.QtSvg is not installed)"
        if QSvgGenerator is not None:
            current_svg = output / f"current-{case_id}-{theme}.svg"
            working_svg = output / f"working-{case_id}-{theme}.svg"
            if _render_svg(current_adapter, current_svg, *PANEL_SIZE) and _render_svg(
                working_adapter, working_svg, *PANEL_SIZE
            ):
                generated.extend((current_svg, working_svg))
                svg_status = "exported through the same QGraphicsScene renderer"
            else:
                svg_status = "Qt SVG generator failed; PNG export remains available"

        profile_path = output / "working-profile.json"
        profile_path.write_text(
            encode_profile(working, name=profile_name), encoding="utf-8"
        )
        generated.append(profile_path)
        report_path = output / "report.md"
        report = _build_report(
            case,
            current,
            working,
            theme=theme,
            profile_name=profile_name,
            artifacts=generated,
            svg_status=svg_status,
            profile_path=profile_path,
        )
        report_path.write_text(report, encoding="utf-8")
        generated.append(report_path)
        return {
            "case_id": case_id,
            "artifacts": [str(path) for path in generated],
            "report": str(report_path),
            "svg": svg_status,
        }
    finally:
        for adapter in adapters:
            adapter.close()
        app.processEvents()


def _build_report(
    case: DepictionCase,
    current: DepictionStyleParameters,
    working: DepictionStyleParameters,
    *,
    theme: str,
    profile_name: str,
    artifacts: list[Path],
    svg_status: str,
    profile_path: Path,
) -> str:
    payload = profile_payload(working, name=profile_name)
    hashes = "\n".join(
        f"- `{path.name}` — SHA256 `{_sha256(path)}`" for path in artifacts if path.is_file()
    )
    return (
        "# Depiction Visual Lab — reproduction report\n\n"
        f"- ChemUSON version: `{__version__}`\n"
        f"- Renderer source SHA at export: `{source_sha()}`\n"
        f"- Git working tree at export: `{source_worktree_state()}`\n"
        f"- Platform: `{platform.platform()}`\n"
        f"- Qt / PyQt: `{QT_VERSION_STR}` / `{PYQT_VERSION_STR}`\n"
        f"- Case: `{case.case_id}` — {case.title} (`{case.kind}`)\n"
        f"- Theme: `{theme}`; panel: `{PANEL_SIZE[0]}×{PANEL_SIZE[1]}`; antialiasing: on\n"
        f"- Profile: `{profile_path.name}`; schema `{payload['schema']}` v{payload['version']}; based on `{payload['based_on']}`\n"
        f"- SVG: {svg_status}\n\n"
        "## Working profile\n\n```json\n"
        f"{json.dumps(payload, ensure_ascii=False, indent=2, sort_keys=True)}\n``\n\n"
        "## Changed values (Current → Working)\n\n"
        + (
            "\n".join(
                f"- `{key}`: `{value['current']}` → `{value['working']}`"
                for key, value in current.changed_values(working).items()
            )
            or "- No differences; Working matches Current."
        )
        + "\n\n## Artifact SHA256\n\n"
        + hashes
        + "\n"
    )


def export_gallery(
    output_dir: Path,
    current: DepictionStyleParameters,
    working: DepictionStyleParameters,
    *,
    theme: str = "light",
    profile_name: str = "Working",
) -> dict[str, Any]:
    """Generate a deterministic side-by-side contact sheet for every case."""
    app = QApplication.instance()
    if app is None:
        raise RuntimeError("QApplication must exist before exporting depiction evidence")
    output = Path(output_dir)
    output.mkdir(parents=True, exist_ok=True)
    card_width, card_height = 1040, 360
    columns = 3
    rows = (len(CASES) + columns - 1) // columns
    sheet = QImage(
        card_width * columns,
        (card_height + 36) * rows + 60,
        QImage.Format.Format_ARGB32_Premultiplied,
    )
    sheet.fill(QColor("#E2E8F0"))
    painter = QPainter(sheet)
    painter.setRenderHint(QPainter.RenderHint.Antialiasing, True)
    painter.setRenderHint(QPainter.RenderHint.TextAntialiasing, True)
    painter.setFont(QFont("DejaVu Sans", 18, QFont.Weight.DemiBold))
    painter.setPen(QColor("#0F172A"))
    painter.drawText(18, 32, f"ChemUSON Molecular Depiction Visual Lab · {theme.upper()}")
    try:
        for index, case in enumerate(CASES):
            graph = case_graph(case.case_id)
            left = ProductionCanvasAdapter.create(graph, parameters=current, theme=theme)
            right = ProductionCanvasAdapter.create(graph, parameters=working, theme=theme)
            try:
                left_image = left.render_image(width=500, height=300)
                right_image = right.render_image(width=500, height=300)
            finally:
                left.close()
                right.close()
                app.processEvents()
            row, column = divmod(index, columns)
            x = column * card_width
            y = 60 + row * (card_height + 36)
            painter.fillRect(x, y, card_width - 12, card_height, QColor("#F8FAFC"))
            painter.setPen(QColor("#0F172A"))
            painter.setFont(QFont("DejaVu Sans", 11, QFont.Weight.DemiBold))
            painter.drawText(x + 10, y + 20, f"{case.case_id} · {case.title} · {case.kind}")
            painter.setFont(QFont("DejaVu Sans", 8))
            painter.setPen(QColor("#475569"))
            painter.drawText(x + 10, y + 34, case.description)
            painter.drawText(x + 10, y + 52, "CURRENT")
            painter.drawText(x + 520, y + 52, "WORKING")
            painter.drawImage(x + 10, y + 58, left_image.scaled(500, 290))
            painter.drawImage(x + 520, y + 58, right_image.scaled(500, 290))
    finally:
        painter.end()

    contact_path = output / f"contact-sheet-{theme}.png"
    _save_image(sheet, contact_path)
    profile_path = output / "working-profile.json"
    profile_path.write_text(encode_profile(working, name=profile_name), encoding="utf-8")
    summary_path = output / "report.md"
    changes = current.changed_values(working)
    artifact_hash = _sha256(contact_path)
    profile_hash = _sha256(profile_path)
    summary_path.write_text(
        "# Depiction Visual Lab — gallery evidence\n\n"
        f"- ChemUSON version: `{__version__}`\n"
        f"- Renderer source SHA at export: `{source_sha()}`\n"
        f"- Git working tree at export: `{source_worktree_state()}`\n"
        f"- Platform: `{platform.platform()}`; Qt/PyQt `{QT_VERSION_STR}` / `{PYQT_VERSION_STR}`\n"
        f"- Theme: `{theme}`; cases: {len(CASES)}; contact sheet: `{contact_path.name}`\n"
        f"- Profile: `{profile_path.name}`; based on `current-renderer`\n\n"
        "## Cases\n\n"
        + "\n".join(f"- `{case.case_id}` — {case.title} (`{case.kind}`): {case.description}" for case in CASES)
        + "\n\n## Modified parameters\n\n"
        + (
            "\n".join(f"- `{key}`: `{value['current']}` → `{value['working']}`" for key, value in changes.items())
            or "- Working matches Current."
        )
        + f"\n\n## SHA256\n\n- `{contact_path.name}` — `{artifact_hash}`\n"
        + f"- `{profile_path.name}` — `{profile_hash}`\n",
        encoding="utf-8",
    )
    return {
        "contact_sheet": str(contact_path),
        "profile": str(profile_path),
        "report": str(summary_path),
        "sha256": artifact_hash,
        "case_count": len(CASES),
    }
