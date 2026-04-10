"""Genera el index.html de GitHub Pages para los canales Flatpak publicados."""

from __future__ import annotations

import argparse
from pathlib import Path

CHANNELS = (
    ("stable", "Stable"),
    ("beta", "Beta"),
)


def collect_channels(*, root: Path, basename: str) -> list[dict[str, str]]:
    channels: list[dict[str, str]] = []
    for channel, label in CHANNELS:
        ref_rel = Path("flatpak") / channel / f"{basename}-{channel}.flatpakref"
        repo_rel = Path("flatpak") / channel / f"{basename}-{channel}.flatpakrepo"
        summary_rel = Path("flatpak") / channel / "repo" / "summary"
        if not (root / ref_rel).is_file():
            continue
        if not (root / repo_rel).is_file():
            continue
        if not (root / summary_rel).is_file():
            continue
        channels.append(
            {
                "label": label,
                "ref_href": f"./{ref_rel.as_posix()}",
                "repo_href": f"./{repo_rel.as_posix()}",
            }
        )
    return channels


def build_index_html(channels: list[dict[str, str]]) -> str:
    if channels:
        items = "\n".join(
            (
                f'      <li><a href="{entry["ref_href"]}">{entry["label"]} install (.flatpakref)</a> '
                f'- <a href="{entry["repo_href"]}">{entry["label"]} remote (.flatpakrepo)</a></li>'
            )
            for entry in channels
        )
        intro = "Install from the published Flatpak channels."
    else:
        items = "      <li>No Flatpak channels published yet.</li>"
        intro = "Published channels will appear here after the first successful release."

    return f"""<!doctype html>
<html lang="en">
<head>
  <meta charset="utf-8">
  <meta name="viewport" content="width=device-width, initial-scale=1">
  <title>Chemuson Flatpak</title>
</head>
<body>
  <h1>Chemuson Flatpak</h1>
  <p>{intro}</p>
  <ul>
{items}
  </ul>
</body>
</html>
"""


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description="Genera index.html para el remoto Flatpak publicado en GitHub Pages."
    )
    parser.add_argument("--root", required=True, help="Directorio raiz de gh-pages.")
    parser.add_argument("--basename", default="Chemuson", help="Prefijo de artifacts.")
    parser.add_argument(
        "--output",
        default="index.html",
        help="Ruta relativa a root donde se escribira el index.",
    )
    return parser


def main() -> None:
    args = build_parser().parse_args()
    root = Path(args.root).resolve()
    channels = collect_channels(root=root, basename=args.basename)
    output_path = root / args.output
    output_path.parent.mkdir(parents=True, exist_ok=True)
    output_path.write_text(build_index_html(channels), encoding="utf-8")
    print(f"Wrote {output_path}")


if __name__ == "__main__":
    main()
