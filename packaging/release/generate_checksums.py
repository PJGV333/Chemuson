"""Genera sidecars `.sha256` y archivo agregado `checksums.txt`."""

from __future__ import annotations

import argparse
import hashlib
import os
from pathlib import Path


def sha256_of(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as fh:
        while True:
            chunk = fh.read(8192)
            if not chunk:
                break
            digest.update(chunk)
    return digest.hexdigest()


def iter_artifacts(root: Path):
    for path in sorted(root.iterdir()):
        if not path.is_file():
            continue
        if path.name.endswith(".sha256") or path.name.endswith(".sig"):
            continue
        if path.name in {"checksums.txt", "stable.json", "beta.json"}:
            continue
        yield path


def main() -> None:
    parser = argparse.ArgumentParser(description="Genera checksums SHA-256 de artifacts.")
    parser.add_argument("--dir", required=True, help="Directorio de artifacts.")
    args = parser.parse_args()
    root = Path(args.dir).resolve()
    root.mkdir(parents=True, exist_ok=True)
    lines: list[str] = []
    for artifact in iter_artifacts(root):
        checksum = sha256_of(artifact)
        sidecar = artifact.with_name(f"{artifact.name}.sha256")
        sidecar.write_text(f"{checksum}  {artifact.name}\n", encoding="utf-8")
        lines.append(f"{checksum}  {artifact.name}")
    (root / "checksums.txt").write_text("\n".join(lines) + ("\n" if lines else ""), encoding="utf-8")
    print(f"Generated {len(lines)} checksums in {root}")


if __name__ == "__main__":
    main()

