"""Valida el layout esperado de artifacts Flatpak locales y de GitHub Pages."""

from __future__ import annotations

import argparse
from pathlib import Path


def _require_file(path: Path) -> None:
    if not path.is_file():
        raise FileNotFoundError(f"Missing file: {path}")
    if path.stat().st_size == 0:
        raise ValueError(f"Empty file: {path}")


def _require_dir_with_files(path: Path) -> None:
    if not path.is_dir():
        raise FileNotFoundError(f"Missing directory: {path}")
    if not any(child.is_file() for child in path.rglob("*")):
        raise ValueError(f"Directory has no files: {path}")


def validate_build_output(*, root: Path, basename: str, channel: str) -> None:
    _require_file(root / f"{basename}-{channel}.flatpakref")
    _require_file(root / f"{basename}-{channel}.flatpakrepo")

    repo_dir = root / "repo"
    _require_file(repo_dir / "config")
    _require_file(repo_dir / "summary")
    _require_dir_with_files(repo_dir / "objects")
    _require_dir_with_files(repo_dir / "refs")


def validate_publish_payload(*, root: Path, basename: str, channel: str) -> None:
    _require_file(root / "icon.svg")
    channel_dir = root / channel

    _require_file(channel_dir / f"{basename}-{channel}.flatpakref")
    _require_file(channel_dir / f"{basename}-{channel}.flatpakrepo")

    repo_dir = channel_dir / "repo"
    _require_file(repo_dir / "config")
    _require_file(repo_dir / "summary")
    _require_dir_with_files(repo_dir / "objects")
    _require_dir_with_files(repo_dir / "refs")


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description="Valida artifacts Flatpak para build local y payload de GitHub Pages."
    )
    parser.add_argument(
        "--mode",
        required=True,
        choices=("build-output", "publish-payload"),
        help="Tipo de layout esperado.",
    )
    parser.add_argument("--root", required=True, help="Directorio raiz a validar.")
    parser.add_argument("--channel", required=True, help="Canal Flatpak a validar.")
    parser.add_argument("--basename", default="Chemuson", help="Prefijo de artifacts.")
    return parser


def main() -> None:
    args = build_parser().parse_args()
    root = Path(args.root).resolve()

    if args.mode == "build-output":
        validate_build_output(root=root, basename=args.basename, channel=args.channel)
    else:
        validate_publish_payload(root=root, basename=args.basename, channel=args.channel)

    print(f"Validated {args.mode} at {root}")


if __name__ == "__main__":
    main()
