"""Firma artifacts usando HMAC-SHA256 (MVP)."""

from __future__ import annotations

import argparse
import hashlib
import hmac
import os
from pathlib import Path


def sign_file(path: Path, key: str) -> str:
    mac = hmac.new(key.encode("utf-8"), digestmod=hashlib.sha256)
    with path.open("rb") as fh:
        while True:
            chunk = fh.read(8192)
            if not chunk:
                break
            mac.update(chunk)
    return mac.hexdigest()


def iter_signable(root: Path):
    for path in sorted(root.iterdir()):
        if not path.is_file():
            continue
        if path.name.endswith(".sig"):
            continue
        if path.name.endswith(".sha256"):
            continue
        yield path


def main() -> None:
    parser = argparse.ArgumentParser(description="Firma artifacts con HMAC-SHA256.")
    parser.add_argument("--dir", required=True, help="Directorio con artifacts.")
    parser.add_argument("--key", default="", help="Clave HMAC (si no, usa env CHEMUSON_SIGN_KEY).")
    args = parser.parse_args()

    key = args.key or os.getenv("CHEMUSON_SIGN_KEY", "")
    if not key:
        raise SystemExit("Missing HMAC key. Provide --key or CHEMUSON_SIGN_KEY.")

    root = Path(args.dir).resolve()
    count = 0
    for artifact in iter_signable(root):
        signature = sign_file(artifact, key)
        sig_path = artifact.with_name(f"{artifact.name}.sig")
        sig_path.write_text(signature + "\n", encoding="utf-8")
        count += 1
    print(f"Signed {count} files in {root}")


if __name__ == "__main__":
    main()

