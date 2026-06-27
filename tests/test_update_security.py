"""Pruebas de endurecimiento de descargas del updater."""

from io import BytesIO


import chemuson.update.security as security


class _FakeResponse:
    def __init__(self, payload: bytes, status: int = 200):
        self._buffer = BytesIO(payload)
        self.status = status

    def read(self, size: int = -1) -> bytes:
        return self._buffer.read(size)

    def __enter__(self):
        return self

    def __exit__(self, exc_type, exc, tb):
        _ = exc_type, exc, tb
        return False


def test_download_binary_rejects_non_https(tmp_path) -> None:
    out = tmp_path / "out.bin"
    try:
        security.download_binary("http://example.com/file.bin", str(out))
        assert False, "Debe rechazar URL no HTTPS"
    except ValueError:
        pass


def test_download_text_retries_and_succeeds(monkeypatch) -> None:
    attempts = {"n": 0}

    def _fake_urlopen(_request, timeout=0, context=None):
        _ = timeout, context
        attempts["n"] += 1
        if attempts["n"] < 2:
            raise OSError("transient")
        return _FakeResponse(b"ok")

    monkeypatch.setattr(security, "urlopen", _fake_urlopen)
    text = security.download_text("https://example.com/file.txt", retries=2, backoff_sec=0.0)
    assert text == "ok"
    assert attempts["n"] == 2


def test_download_text_passes_ssl_context(monkeypatch) -> None:
    sentinel = object()
    seen: dict[str, object] = {}

    def _fake_create_ssl_context():
        return sentinel

    def _fake_urlopen(_request, timeout=0, context=None):
        _ = timeout
        seen["context"] = context
        return _FakeResponse(b"ok")

    monkeypatch.setattr(security, "_create_ssl_context", _fake_create_ssl_context)
    monkeypatch.setattr(security, "urlopen", _fake_urlopen)

    assert security.download_text("https://example.com/file.txt", retries=0) == "ok"
    assert seen["context"] is sentinel
