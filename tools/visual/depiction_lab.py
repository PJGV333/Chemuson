#!/usr/bin/env python3
"""Launch ChemUSON Molecular Depiction Visual Lab or export reproducible evidence."""
from __future__ import annotations

import argparse
import json
import os
from pathlib import Path
import sys

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "src"))
sys.path.insert(0, str(ROOT))


def _parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description="Calibrate the production ChemUSON molecular depiction renderer."
    )
    parser.add_argument("--smoke", action="store_true", help="construct, show, and close the Qt window without entering the event loop")
    parser.add_argument("--export-gallery", type=Path, metavar="DIR", help="export a complete side-by-side PNG contact sheet and profile report")
    parser.add_argument("--export-case", metavar="CASE_ID", help="export Current/Working PNG and SVG for one case")
    parser.add_argument("--output", type=Path, default=Path("depiction-lab-evidence"), help="output directory for batch exports")
    parser.add_argument("--profile", type=Path, help="load a validated Working JSON profile")
    parser.add_argument("--theme", choices=("light", "dark"), default="light")
    parser.add_argument("--profile-name", default="Working")
    return parser


def main(argv: list[str] | None = None) -> int:
    args = _parser().parse_args(argv)
    batch = bool(args.export_gallery or args.export_case or args.smoke)
    if batch:
        os.environ.setdefault("QT_QPA_PLATFORM", "offscreen")

    from PyQt6.QtCore import QCoreApplication
    from PyQt6.QtWidgets import QApplication

    QCoreApplication.setOrganizationName("ChemUSON")
    QCoreApplication.setApplicationName("DepictionVisualLab")
    app = QApplication.instance() or QApplication([sys.argv[0]])

    from tools.visual.depiction_lab_app.parameters import DepictionStyleParameters
    from tools.visual.depiction_lab_app.profiles import ProfileError, load_profile

    profile = DepictionStyleParameters.current()
    profile_name = args.profile_name
    if args.profile is not None:
        try:
            profile_name, profile = load_profile(args.profile)
        except (OSError, ProfileError) as exc:
            print(f"Profile error: {exc}", file=sys.stderr)
            return 2

    if args.export_gallery:
        from tools.visual.depiction_lab_app.exporting import export_gallery

        result = export_gallery(
            args.export_gallery,
            DepictionStyleParameters.current(),
            profile,
            theme=args.theme,
            profile_name=profile_name,
        )
        print(json.dumps(result, ensure_ascii=False, indent=2))
        return 0

    if args.export_case:
        from tools.visual.depiction_lab_app.exporting import export_case_bundle

        result = export_case_bundle(
            args.output,
            args.export_case,
            DepictionStyleParameters.current(),
            profile,
            theme=args.theme,
            profile_name=profile_name,
        )
        print(json.dumps(result, ensure_ascii=False, indent=2))
        return 0

    from tools.visual.depiction_lab_app.ui import DepictionLabWindow

    window = DepictionLabWindow()
    if args.profile is not None:
        window.working_name = profile_name
        window.profile_name.setText(profile_name)
        window.working_parameters = profile
        if window.working_panel.adapter is not None:
            window.working_panel.adapter.apply_parameters(profile)
        window._refresh_parameter_controls()
        window._update_change_summary()
    window.show()
    app.processEvents()
    if args.smoke:
        window.close()
        app.processEvents()
        return 0
    return int(app.exec())


if __name__ == "__main__":
    raise SystemExit(main())
