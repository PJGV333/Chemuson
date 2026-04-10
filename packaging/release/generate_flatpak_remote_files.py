"""Genera archivos `.flatpakrepo` y `.flatpakref` para el remoto oficial."""

from __future__ import annotations

import argparse
import base64
from pathlib import Path


def normalize_repo_url(url: str) -> str:
    """Normaliza URLs de repo OSTree con slash final."""
    normalized = str(url or "").strip()
    if not normalized:
        raise ValueError("La URL del repo Flatpak no puede estar vacía.")
    return normalized.rstrip("/") + "/"


def _single_line(value: str) -> str:
    return " ".join(str(value or "").split())


def _required_url(url: str) -> str:
    normalized = str(url or "").strip()
    if not normalized:
        raise ValueError("La URL requerida no puede estar vacía.")
    return normalized


def encode_gpg_key(gpg_key_file: str) -> str:
    """Codifica una clave pública GPG al formato esperado por Flatpak."""
    key_path = Path(gpg_key_file).resolve()
    payload = key_path.read_bytes()
    if not payload:
        raise ValueError(f"El archivo de clave GPG está vacío: {key_path}")
    return base64.b64encode(payload).decode("ascii")


def build_flatpak_repo_config(
    *,
    title: str,
    repo_url: str,
    homepage: str = "",
    comment: str = "",
    description: str = "",
    icon_url: str = "",
    default_branch: str = "",
    gpg_key: str = "",
) -> str:
    """Construye el contenido de un archivo `.flatpakrepo`."""
    lines = [
        "[Flatpak Repo]",
        f"Title={_single_line(title)}",
        f"Url={normalize_repo_url(repo_url)}",
    ]
    if homepage:
        lines.append(f"Homepage={homepage.strip()}")
    if comment:
        lines.append(f"Comment={_single_line(comment)}")
    if description:
        lines.append(f"Description={_single_line(description)}")
    if icon_url:
        lines.append(f"Icon={icon_url.strip()}")
    if default_branch:
        lines.append(f"DefaultBranch={default_branch.strip()}")
    if gpg_key:
        lines.append(f"GPGKey={gpg_key.strip()}")
    return "\n".join(lines) + "\n"


def build_flatpak_ref_config(
    *,
    title: str,
    app_id: str,
    branch: str,
    repo_url: str,
    runtime_repo: str,
    suggest_remote_name: str = "",
    gpg_key: str = "",
) -> str:
    """Construye el contenido de un archivo `.flatpakref`."""
    lines = [
        "[Flatpak Ref]",
        f"Title={_single_line(title)}",
        f"Name={app_id.strip()}",
        f"Branch={branch.strip()}",
        "IsRuntime=false",
        f"Url={normalize_repo_url(repo_url)}",
        f"RuntimeRepo={_required_url(runtime_repo)}",
    ]
    if suggest_remote_name:
        lines.append(f"SuggestRemoteName={suggest_remote_name.strip()}")
    if gpg_key:
        lines.append(f"GPGKey={gpg_key.strip()}")
    return "\n".join(lines) + "\n"


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description="Genera archivos Flatpak para instalar Chemuson desde un remoto oficial."
    )
    parser.add_argument("--repo-url", required=True, help="URL pública del repo OSTree.")
    parser.add_argument("--app-id", required=True, help="App ID Flatpak.")
    parser.add_argument("--branch", required=True, help="Branch/canal Flatpak.")
    parser.add_argument("--out-dir", required=True, help="Directorio destino.")
    parser.add_argument("--basename", default="Chemuson", help="Prefijo de archivos.")
    parser.add_argument("--repo-title", default="Chemuson", help="Título del remoto.")
    parser.add_argument(
        "--ref-title",
        default="Chemuson",
        help="Título visible del archivo .flatpakref.",
    )
    parser.add_argument("--homepage", default="", help="Homepage opcional.")
    parser.add_argument("--comment", default="", help="Comentario corto opcional.")
    parser.add_argument("--description", default="", help="Descripción opcional.")
    parser.add_argument("--icon-url", default="", help="URL pública del icono.")
    parser.add_argument(
        "--default-branch",
        default="",
        help="Branch por defecto del remoto.",
    )
    parser.add_argument(
        "--runtime-repo",
        default="https://dl.flathub.org/repo/flathub.flatpakrepo",
        help="Repositorio Flatpak del runtime.",
    )
    parser.add_argument(
        "--suggest-remote-name",
        default="",
        help="Nombre sugerido del remote al instalar desde .flatpakref.",
    )
    parser.add_argument(
        "--gpg-key-file",
        default="",
        help="Ruta opcional de clave pública GPG para incrustar en .flatpakrepo y .flatpakref.",
    )
    return parser


def main() -> None:
    parser = build_parser()
    args = parser.parse_args()
    out_dir = Path(args.out_dir).resolve()
    out_dir.mkdir(parents=True, exist_ok=True)

    gpg_key = encode_gpg_key(args.gpg_key_file) if args.gpg_key_file else ""
    repo_text = build_flatpak_repo_config(
        title=args.repo_title,
        repo_url=args.repo_url,
        homepage=args.homepage,
        comment=args.comment,
        description=args.description,
        icon_url=args.icon_url,
        default_branch=args.default_branch,
        gpg_key=gpg_key,
    )
    ref_text = build_flatpak_ref_config(
        title=args.ref_title,
        app_id=args.app_id,
        branch=args.branch,
        repo_url=args.repo_url,
        runtime_repo=args.runtime_repo,
        suggest_remote_name=args.suggest_remote_name,
        gpg_key=gpg_key,
    )

    repo_path = out_dir / f"{args.basename}-{args.branch}.flatpakrepo"
    ref_path = out_dir / f"{args.basename}-{args.branch}.flatpakref"
    repo_path.write_text(repo_text, encoding="utf-8")
    ref_path.write_text(ref_text, encoding="utf-8")
    print(f"Wrote {repo_path}")
    print(f"Wrote {ref_path}")


if __name__ == "__main__":
    main()
