"""Valida el layout esperado de artifacts Flatpak locales y de GitHub Pages."""

from __future__ import annotations

import argparse
import base64
import binascii
import configparser
from pathlib import Path
import subprocess
import tempfile


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


def _load_ini(path: Path) -> configparser.ConfigParser:
    parser = configparser.ConfigParser(interpolation=None)
    try:
        parser.read_string(path.read_text(encoding="utf-8"))
    except (configparser.Error, UnicodeDecodeError) as exc:
        raise ValueError(f"Invalid Flatpak config file: {path}") from exc
    return parser


def _read_gpg_key(repo_config: Path, ref_config: Path) -> bytes:
    repo_parser = _load_ini(repo_config)
    ref_parser = _load_ini(ref_config)
    repo_key = repo_parser.get("Flatpak Repo", "GPGKey", fallback="").strip()
    ref_key = ref_parser.get("Flatpak Ref", "GPGKey", fallback="").strip()

    if repo_key and ref_key and repo_key != ref_key:
        raise ValueError(
            "Flatpak repo/ref artifacts embed different GPG keys; refusing to validate."
        )

    encoded_key = repo_key or ref_key
    if not encoded_key:
        return b""

    try:
        return base64.b64decode(encoded_key, validate=True)
    except binascii.Error as exc:
        raise ValueError("Embedded Flatpak GPG key is not valid base64.") from exc


def _run_command(args: list[str]) -> subprocess.CompletedProcess[str]:
    try:
        return subprocess.run(args, check=True, capture_output=True, text=True)
    except FileNotFoundError as exc:
        raise RuntimeError(f"Required command not found: {args[0]}") from exc
    except subprocess.CalledProcessError as exc:
        details = (exc.stderr or exc.stdout or str(exc)).strip()
        raise RuntimeError(f"Command failed: {' '.join(args)}\n{details}") from exc


def _commitmeta_path(repo_dir: Path, commit: str) -> Path:
    checksum = commit.strip()
    if len(checksum) < 3:
        raise ValueError(f"Invalid OSTree commit checksum: {checksum!r}")
    return repo_dir / "objects" / checksum[:2] / f"{checksum[2:]}.commitmeta"


def _verify_signed_app_ref(
    *,
    repo_dir: Path,
    repo_config: Path,
    ref_config: Path,
    channel: str,
    app_id: str,
    arch: str,
) -> None:
    gpg_key = _read_gpg_key(repo_config, ref_config)
    if not gpg_key:
        return

    _require_file(repo_dir / "summary.sig")

    app_ref = f"app/{app_id}/{arch}/{channel}"
    commit = _run_command(
        ["ostree", "rev-parse", f"--repo={repo_dir}", app_ref]
    ).stdout.strip()
    if not commit:
        raise ValueError(f"Unable to resolve OSTree ref: {app_ref}")

    _require_file(_commitmeta_path(repo_dir, commit))

    with tempfile.TemporaryDirectory(prefix="chemuson-flatpak-verify-") as temp_dir:
        temp_root = Path(temp_dir)
        key_path = temp_root / "remote.gpg"
        verify_repo = temp_root / "verify-repo"
        key_path.write_bytes(gpg_key)

        _run_command(["ostree", "init", f"--repo={verify_repo}", "--mode=bare-user-only"])
        _run_command(
            [
                "ostree",
                "remote",
                "add",
                f"--repo={verify_repo}",
                "--force",
                f"--gpg-import={key_path}",
                "chemuson-verify",
                repo_dir.resolve().as_uri(),
            ]
        )
        _run_command(["ostree", "pull", f"--repo={verify_repo}", "chemuson-verify", app_ref])

        verified_commit = _run_command(
            ["ostree", "rev-parse", f"--repo={verify_repo}", app_ref]
        ).stdout.strip()
        if verified_commit != commit:
            raise ValueError(
                f"Verified commit mismatch for {app_ref}: expected {commit}, got {verified_commit}"
            )


def validate_build_output(
    *,
    root: Path,
    basename: str,
    channel: str,
    app_id: str = "io.github.PJGV333.Chemuson",
    arch: str = "x86_64",
) -> None:
    ref_config = root / f"{basename}-{channel}.flatpakref"
    repo_config = root / f"{basename}-{channel}.flatpakrepo"
    _require_file(ref_config)
    _require_file(repo_config)

    repo_dir = root / "repo"
    _require_file(repo_dir / "config")
    _require_file(repo_dir / "summary")
    _require_dir_with_files(repo_dir / "objects")
    _require_dir_with_files(repo_dir / "refs")
    _verify_signed_app_ref(
        repo_dir=repo_dir,
        repo_config=repo_config,
        ref_config=ref_config,
        channel=channel,
        app_id=app_id,
        arch=arch,
    )


def validate_publish_payload(
    *,
    root: Path,
    basename: str,
    channel: str,
    app_id: str = "io.github.PJGV333.Chemuson",
    arch: str = "x86_64",
) -> None:
    _require_file(root / "icon.svg")
    channel_dir = root / channel

    ref_config = channel_dir / f"{basename}-{channel}.flatpakref"
    repo_config = channel_dir / f"{basename}-{channel}.flatpakrepo"
    _require_file(ref_config)
    _require_file(repo_config)

    repo_dir = channel_dir / "repo"
    _require_file(repo_dir / "config")
    _require_file(repo_dir / "summary")
    _require_dir_with_files(repo_dir / "objects")
    _require_dir_with_files(repo_dir / "refs")
    _verify_signed_app_ref(
        repo_dir=repo_dir,
        repo_config=repo_config,
        ref_config=ref_config,
        channel=channel,
        app_id=app_id,
        arch=arch,
    )


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
    parser.add_argument(
        "--app-id",
        default="io.github.PJGV333.Chemuson",
        help="App ID Flatpak a verificar.",
    )
    parser.add_argument(
        "--arch",
        default="x86_64",
        help="Arquitectura Flatpak a verificar.",
    )
    return parser


def main() -> None:
    args = build_parser().parse_args()
    root = Path(args.root).resolve()

    if args.mode == "build-output":
        validate_build_output(
            root=root,
            basename=args.basename,
            channel=args.channel,
            app_id=args.app_id,
            arch=args.arch,
        )
    else:
        validate_publish_payload(
            root=root,
            basename=args.basename,
            channel=args.channel,
            app_id=args.app_id,
            arch=args.arch,
        )

    print(f"Validated {args.mode} at {root}")


if __name__ == "__main__":
    main()
