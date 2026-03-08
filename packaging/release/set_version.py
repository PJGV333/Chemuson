"""Sincroniza la versión canónica de Chemuson en `_version.py`."""

from __future__ import annotations

import argparse
import re
import sys
from pathlib import Path


REPO_ROOT = Path(__file__).resolve().parents[2]
SRC_ROOT = REPO_ROOT / "src"
if str(SRC_ROOT) not in sys.path:
    sys.path.insert(0, str(SRC_ROOT))

from chemuson.update.semver import parse_semver


VERSION_ASSIGN_RE = re.compile(r'^(?P<prefix>\s*__version__\s*=\s*")(?P<value>[^"]+)(?P<suffix>".*)$')


def validate_version(version: str) -> str:
    """Normaliza y valida una versión SemVer soportada por el updater."""
    normalized = parse_semver(str(version or "").strip()).normalized()
    return normalized


def read_version(version_file: Path) -> str:
    """Lee la versión actual desde el archivo indicado."""
    for line in version_file.read_text(encoding="utf-8").splitlines():
        match = VERSION_ASSIGN_RE.match(line)
        if match:
            return str(match.group("value"))
    raise ValueError(f"No se encontró __version__ en {version_file}")


def write_version(version_file: Path, version: str) -> str:
    """Reemplaza `__version__` por la nueva versión validada."""
    normalized = validate_version(version)
    lines = version_file.read_text(encoding="utf-8").splitlines()
    updated_lines: list[str] = []
    replaced = False
    for line in lines:
        match = VERSION_ASSIGN_RE.match(line)
        if match:
            updated_lines.append(
                f'{match.group("prefix")}{normalized}{match.group("suffix")}'
            )
            replaced = True
        else:
            updated_lines.append(line)
    if not replaced:
        raise ValueError(f"No se encontró __version__ en {version_file}")
    version_file.write_text("\n".join(updated_lines) + "\n", encoding="utf-8")
    return normalized


def build_parser() -> argparse.ArgumentParser:
    """Construye el parser CLI."""
    parser = argparse.ArgumentParser(
        description="Sincroniza la versión canónica de Chemuson en src/chemuson/_version.py."
    )
    parser.add_argument("--version", required=True, help="Versión SemVer destino.")
    parser.add_argument(
        "--file",
        default=str(REPO_ROOT / "src" / "chemuson" / "_version.py"),
        help="Ruta del archivo de versión.",
    )
    return parser


def main() -> None:
    """Punto de entrada CLI."""
    parser = build_parser()
    args = parser.parse_args()
    version_file = Path(args.file).resolve()
    normalized = write_version(version_file, args.version)
    print(f"Updated {version_file} to {normalized}")


if __name__ == "__main__":
    main()
