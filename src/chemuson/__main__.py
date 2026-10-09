"""CLI entry point para ejecutar Chemuson como módulo."""

from __future__ import annotations

import argparse
import io
import json
import os
import sys
from pathlib import Path

from chemuson.version import get_app_version


def _build_parser() -> argparse.ArgumentParser:
    """Construye parser de CLI de Chemuson."""
    parser = argparse.ArgumentParser(prog="chemuson")
    parser.add_argument(
        "--version",
        action="store_true",
        help="Muestra la versión de Chemuson y termina.",
    )
    parser.add_argument("--icon-smoke-test", action="store_true", help=argparse.SUPPRESS)
    parser.add_argument("--rdkit-packaged-smoke-test", action="store_true", help=argparse.SUPPRESS)
    parser.add_argument("--rdkit-smoke-report", type=Path, help=argparse.SUPPRESS)
    parser.add_argument("--chemname-packaged-smoke-test", action="store_true", help=argparse.SUPPRESS)
    parser.add_argument("--chemname-smoke-report", type=Path, help=argparse.SUPPRESS)
    return parser


def _run_internal_rdkit_worker(arguments: list[str]) -> int:
    """Handle private JSON-file IPC before the frozen GUI bootstrap."""
    if (
        not getattr(sys, "frozen", False)
        or os.environ.get("CHEMUSON_INTERNAL_RDKIT_WORKER") != "1"
        or len(arguments) != 2
    ):
        return 64

    request_path, response_path = (Path(value) for value in arguments)
    try:
        from chemuson.chemio._rdkit_worker import main as worker_main

        with request_path.open("r", encoding="utf-8") as request_file:
            request_text = request_file.read()
        with response_path.open("w", encoding="utf-8") as response_file:
            result = worker_main(io.StringIO(request_text), response_file)
            response_file.flush()
        return int(result or 0)
    except Exception as exc:
        try:
            response_path.write_text(
                json.dumps({"ok": False, "error": "worker_bootstrap_failed", "detail": str(exc)}),
                encoding="utf-8",
            )
        except OSError:
            return 70
        return 0


def main() -> int:
    """Punto de entrada principal de CLI."""
    if len(sys.argv) > 1 and sys.argv[1] == "--chemuson-internal-rdkit-worker":
        return _run_internal_rdkit_worker(sys.argv[2:])

    parser = _build_parser()
    args = parser.parse_args()
    if args.icon_smoke_test:
        if os.environ.get("CHEMUSON_ICON_SMOKE_TEST") != "1":
            parser.error("--icon-smoke-test is reserved for the packaging smoke workflow.")
        from chemuson.gui.theme.icon_smoke import run_icon_smoke_test

        try:
            print(json.dumps(run_icon_smoke_test(), sort_keys=True))
        except Exception as exc:
            parser.exit(1, f"Packaged icon smoke failed: {exc}\n")
        return 0
    if args.rdkit_packaged_smoke_test:
        if (
            not getattr(sys, "frozen", False)
            or os.environ.get("CHEMUSON_RDKIT_PACKAGED_SMOKE") != "1"
        ):
            parser.error("--rdkit-packaged-smoke-test is reserved for the packaging smoke workflow.")
        if args.rdkit_smoke_report is None:
            parser.error("--rdkit-smoke-report is required for packaged RDKit validation.")
        try:
            from chemuson.chemio.rdkit_packaged_smoke import run_smoke

            report = run_smoke()
        except Exception as exc:
            report = {"ok": False, "error": "packaged_smoke_failed", "detail": str(exc)}
        try:
            args.rdkit_smoke_report.write_text(
                json.dumps(report, sort_keys=True), encoding="utf-8"
            )
        except OSError:
            return 1
        return 0 if report.get("ok") is True else 1
    if args.chemname_packaged_smoke_test:
        if (
            not getattr(sys, "frozen", False)
            or os.environ.get("CHEMUSON_CHEMNAME_PACKAGED_SMOKE") != "1"
        ):
            parser.error("--chemname-packaged-smoke-test is reserved for the packaging smoke workflow.")
        if args.chemname_smoke_report is None:
            parser.error("--chemname-smoke-report is required for packaged ChemName validation.")
        try:
            from chemuson.chemname.packaged_smoke import run_smoke

            report = run_smoke()
        except Exception as exc:
            report = {
                "ok": False,
                "error": "packaged_chemname_smoke_failed",
                "exception_type": type(exc).__name__,
                "detail": str(exc),
            }
        try:
            args.chemname_smoke_report.write_text(
                json.dumps(report, sort_keys=True), encoding="utf-8"
            )
        except OSError:
            return 1
        return 0 if report.get("ok") is True else 1
    if args.version:
        print(get_app_version())
        return 0

    from chemuson.app.bootstrap import run_app

    run_app()
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
