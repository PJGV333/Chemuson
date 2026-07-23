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
    return parser


def main() -> None:
    """Punto de entrada principal de CLI."""
    parser = _build_parser()
    args = parser.parse_args()
    if args.version:
        print(get_app_version())
        return

    from chemuson.app.bootstrap import run_app

    run_app()


if __name__ == "__main__":
    main()
