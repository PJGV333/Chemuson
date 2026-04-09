"""Sincroniza versión y metadatos públicos de release de Chemuson."""

from __future__ import annotations

import argparse
import re
import sys
import xml.etree.ElementTree as ET
from datetime import date, datetime, timezone
from pathlib import Path


REPO_ROOT = Path(__file__).resolve().parents[2]
SRC_ROOT = REPO_ROOT / "src"
if str(SRC_ROOT) not in sys.path:
    sys.path.insert(0, str(SRC_ROOT))

from chemuson.update.semver import parse_semver


VERSION_ASSIGN_RE = re.compile(r'^(?P<prefix>\s*__version__\s*=\s*")(?P<value>[^"]+)(?P<suffix>".*)$')
DEFAULT_PROJECT_LICENSE = "AGPL-3.0-only"
DEFAULT_METAINFO_PATH = (
    REPO_ROOT / "packaging" / "flatpak" / "io.github.PJGV333.Chemuson.metainfo.xml"
)


def validate_version(version: str) -> str:
    """Normaliza y valida una versión SemVer soportada por el updater."""
    normalized = parse_semver(str(version or "").strip()).normalized()
    return normalized


def default_release_date() -> str:
    """Devuelve la fecha UTC actual en formato ISO."""
    return datetime.now(timezone.utc).date().isoformat()


def validate_release_date(value: str) -> str:
    """Valida una fecha de release en formato ISO `YYYY-MM-DD`."""
    normalized = str(value or "").strip()
    try:
        date.fromisoformat(normalized)
    except ValueError as exc:
        raise ValueError(f"Fecha de release inválida: {value!r}") from exc
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


def update_flatpak_metainfo(
    metainfo_file: Path,
    version: str,
    release_date: str | None = None,
    project_license: str = DEFAULT_PROJECT_LICENSE,
) -> tuple[str, str]:
    """Sincroniza versión, fecha y licencia del metainfo de Flatpak/AppStream."""
    normalized = validate_version(version)
    normalized_date = validate_release_date(release_date or default_release_date())
    tree = ET.parse(metainfo_file)
    root = tree.getroot()

    license_element = root.find("project_license")
    if license_element is None:
        license_element = ET.SubElement(root, "project_license")
    license_element.text = str(project_license or DEFAULT_PROJECT_LICENSE).strip()

    releases_element = root.find("releases")
    if releases_element is None:
        releases_element = ET.SubElement(root, "releases")

    release_element = None
    for candidate in releases_element.findall("release"):
        if candidate.attrib.get("version") == normalized:
            release_element = candidate
            break
    if release_element is None:
        release_element = ET.Element("release")

    release_element.set("version", normalized)
    release_element.set("date", normalized_date)

    remaining_children = [
        child for child in list(releases_element) if child is not release_element
    ]
    releases_element[:] = [release_element, *remaining_children]

    ET.indent(tree, space="  ")
    tree.write(metainfo_file, encoding="utf-8", xml_declaration=True)
    return normalized, normalized_date


def build_parser() -> argparse.ArgumentParser:
    """Construye el parser CLI."""
    parser = argparse.ArgumentParser(
        description=(
            "Sincroniza versión canónica y metadatos Flatpak/AppStream de Chemuson."
        )
    )
    parser.add_argument("--version", required=True, help="Versión SemVer destino.")
    parser.add_argument(
        "--file",
        default=str(REPO_ROOT / "src" / "chemuson" / "_version.py"),
        help="Ruta del archivo de versión.",
    )
    parser.add_argument(
        "--metainfo-file",
        default=str(DEFAULT_METAINFO_PATH),
        help="Ruta del archivo metainfo XML a sincronizar.",
    )
    parser.add_argument(
        "--release-date",
        default=default_release_date(),
        help="Fecha ISO del release para AppStream (default: hoy UTC).",
    )
    parser.add_argument(
        "--project-license",
        default=DEFAULT_PROJECT_LICENSE,
        help="Licencia SPDX publicada en AppStream/Flatpak.",
    )
    return parser


def main() -> None:
    """Punto de entrada CLI."""
    parser = build_parser()
    args = parser.parse_args()
    version_file = Path(args.file).resolve()
    metainfo_file = Path(args.metainfo_file).resolve()
    normalized = write_version(version_file, args.version)
    _, normalized_date = update_flatpak_metainfo(
        metainfo_file=metainfo_file,
        version=normalized,
        release_date=args.release_date,
        project_license=args.project_license,
    )
    print(f"Updated {version_file} to {normalized}")
    print(
        f"Updated {metainfo_file} to version {normalized}, "
        f"date {normalized_date}, license {args.project_license}"
    )


if __name__ == "__main__":
    main()
