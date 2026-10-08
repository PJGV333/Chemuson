"""CLI entry point para ejecutar Chemuson como módulo."""

from __future__ import annotations

import argparse

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
    return parser


def main() -> None:
    """Punto de entrada principal de CLI."""
    parser = _build_parser()
    args = parser.parse_args()
    if args.icon_smoke_test:
        import os

        if os.environ.get("CHEMUSON_ICON_SMOKE_TEST") != "1":
            parser.error("--icon-smoke-test is reserved for the packaging smoke workflow.")
        from chemuson.gui.theme.icon_smoke import run_icon_smoke_test

        try:
            import json

            print(json.dumps(run_icon_smoke_test(), sort_keys=True))
        except Exception as exc:
            parser.exit(1, f"Packaged icon smoke failed: {exc}\n")
        return
    if args.version:
        print(get_app_version())
        return

    from chemuson.app.bootstrap import run_app

    run_app()


if __name__ == "__main__":
    main()
