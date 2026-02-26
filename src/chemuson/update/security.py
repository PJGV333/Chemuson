"""Verificación de integridad y firma para updates."""

from __future__ import annotations

import base64
import hashlib
import hmac
import os
import time
from urllib.request import Request, urlopen


def _validate_https_url(url: str) -> str:
    raw = str(url or "").strip()
    if not raw.lower().startswith("https://"):
        raise ValueError("Solo se permiten URLs HTTPS para descarga de updates.")
    return raw


def sha256_file(path: str) -> str:
    """Calcula SHA-256 hex de un archivo."""
    digest = hashlib.sha256()
    with open(path, "rb") as fh:
        while True:
            chunk = fh.read(8192)
            if not chunk:
                break
            digest.update(chunk)
    return digest.hexdigest()


def sha256_bytes(payload: bytes) -> str:
    """Calcula SHA-256 hex de bytes."""
    return hashlib.sha256(payload).hexdigest()


def download_binary(
    url: str,
    target_path: str,
    timeout: float = 20.0,
    retries: int = 2,
    backoff_sec: float = 0.35,
    max_bytes: int = 1024 * 1024 * 1024,
) -> str:
    """Descarga binario a disco y devuelve ruta."""
    os.makedirs(os.path.dirname(os.path.abspath(target_path)), exist_ok=True)
    download_url = _validate_https_url(url)
    attempts = max(1, int(retries) + 1)
    last_error: Exception | None = None
    for attempt in range(attempts):
        try:
            request = Request(url=download_url, headers={"User-Agent": "chemuson-updater"})
            with urlopen(request, timeout=timeout) as response, open(target_path, "wb") as out:
                total = 0
                while True:
                    chunk = response.read(1024 * 256)
                    if not chunk:
                        break
                    total += len(chunk)
                    if total > int(max_bytes):
                        raise ValueError("Respuesta excede tamaño máximo permitido.")
                    out.write(chunk)
            return target_path
        except Exception as exc:
            last_error = exc
            if attempt >= attempts - 1:
                raise
            time.sleep(max(0.0, float(backoff_sec)) * (2 ** attempt))
    if last_error is not None:
        raise last_error
    return target_path


def download_text(
    url: str,
    timeout: float = 20.0,
    retries: int = 2,
    backoff_sec: float = 0.35,
    max_bytes: int = 4 * 1024 * 1024,
) -> str:
    """Descarga texto UTF-8."""
    download_url = _validate_https_url(url)
    attempts = max(1, int(retries) + 1)
    last_error: Exception | None = None
    for attempt in range(attempts):
        try:
            request = Request(url=download_url, headers={"User-Agent": "chemuson-updater"})
            with urlopen(request, timeout=timeout) as response:
                data = response.read(int(max_bytes) + 1)
            if len(data) > int(max_bytes):
                raise ValueError("Respuesta de texto excede tamaño máximo permitido.")
            return data.decode("utf-8", errors="replace")
        except Exception as exc:
            last_error = exc
            if attempt >= attempts - 1:
                raise
            time.sleep(max(0.0, float(backoff_sec)) * (2 ** attempt))
    if last_error is not None:
        raise last_error
    return ""


def parse_sha256_text(payload: str) -> str:
    """Extrae hash SHA-256 de archivos tipo `sha256sum`."""
    for line in str(payload or "").splitlines():
        stripped = line.strip()
        if not stripped:
            continue
        token = stripped.split()[0]
        if len(token) == 64 and all(ch in "0123456789abcdefABCDEF" for ch in token):
            return token.lower()
    return ""


def verify_file_hash(path: str, expected_sha256: str) -> bool:
    """Valida hash SHA-256 esperado."""
    expected = str(expected_sha256 or "").strip().lower()
    if len(expected) != 64:
        return False
    return hmac.compare_digest(sha256_file(path), expected)


class SignatureVerifier:
    """Verificador de firma (MVP: HMAC-SHA256, extensible a Ed25519)."""

    def __init__(self, key: str = "") -> None:
        self._key = str(key or "")

    def verify_payload(
        self,
        payload: bytes,
        signature: str,
        algorithm: str = "hmac-sha256",
    ) -> bool:
        """Verifica firma sobre bytes."""
        algo = str(algorithm or "").strip().lower()
        if algo == "none":
            return False
        if algo == "hmac-sha256":
            if not self._key:
                return False
            mac = hmac.new(self._key.encode("utf-8"), payload, hashlib.sha256).hexdigest()
            candidate = str(signature or "").strip().lower()
            return hmac.compare_digest(mac, candidate)
        if algo == "ed25519":
            try:
                from cryptography.hazmat.primitives.asymmetric.ed25519 import (
                    Ed25519PublicKey,
                )
            except Exception:
                return False
            if not self._key:
                return False
            try:
                key_bytes = base64.b64decode(self._key)
                sig_bytes = base64.b64decode(str(signature or "").strip())
                pub = Ed25519PublicKey.from_public_bytes(key_bytes)
                pub.verify(sig_bytes, payload)
                return True
            except Exception:
                return False
        return False

    def is_configured(self, algorithm: str = "hmac-sha256") -> bool:
        """Indica si hay material criptográfico para verificar `algorithm`."""
        algo = str(algorithm or "").strip().lower()
        if algo in {"hmac-sha256", "ed25519"}:
            return bool(self._key)
        return False

    def verify_file(
        self,
        file_path: str,
        signature: str,
        algorithm: str = "hmac-sha256",
    ) -> bool:
        """Verifica firma de un archivo completo."""
        with open(file_path, "rb") as fh:
            payload = fh.read()
        return self.verify_payload(payload, signature, algorithm=algorithm)
