"""Pruebas del proveedor GitHub Releases para update."""

import json
from datetime import datetime, timedelta, timezone


from chemuson.update.provider import GitHubReleasesProvider
from chemuson.update.types import UpdateChannel


def _sample_releases_payload() -> list[dict]:
    return [
        {
            "tag_name": "v1.2.0-beta.2",
            "name": "Chemuson 1.2.0-beta.2",
            "prerelease": True,
            "draft": False,
            "published_at": "2026-02-20T12:00:00Z",
            "html_url": "https://example.com/release/beta",
            "assets": [
                {
                    "name": "Chemuson-v1.2.0-beta.2-windows-x86_64-portable.exe",
                    "browser_download_url": "https://example.com/Chemuson-v1.2.0-beta.2-windows-x86_64-portable.exe",
                    "size": 123,
                    "content_type": "application/octet-stream",
                },
                {
                    "name": "Chemuson-v1.2.0-beta.2-windows-x86_64-portable.exe.sha256",
                    "browser_download_url": "https://example.com/Chemuson-v1.2.0-beta.2-windows-x86_64-portable.exe.sha256",
                    "size": 64,
                    "content_type": "text/plain",
                },
            ],
        },
        {
            "tag_name": "v1.1.4",
            "name": "Chemuson 1.1.4",
            "prerelease": False,
            "draft": False,
            "published_at": "2026-02-10T12:00:00Z",
            "html_url": "https://example.com/release/stable",
            "assets": [
                {
                    "name": "Chemuson-v1.1.4-linux-x86_64.AppImage",
                    "browser_download_url": "https://example.com/Chemuson-v1.1.4-linux-x86_64.AppImage",
                    "size": 456,
                    "content_type": "application/octet-stream",
                },
                {
                    "name": "Chemuson-v1.1.4-linux-x86_64.AppImage.sha256",
                    "browser_download_url": "https://example.com/Chemuson-v1.1.4-linux-x86_64.AppImage.sha256",
                    "size": 64,
                    "content_type": "text/plain",
                },
                {
                    "name": "Chemuson-v1.1.4-linux-x86_64.AppImage.sig",
                    "browser_download_url": "https://example.com/Chemuson-v1.1.4-linux-x86_64.AppImage.sig",
                    "size": 64,
                    "content_type": "text/plain",
                },
            ],
        },
    ]


def test_latest_release_selection_by_channel(monkeypatch) -> None:
    provider = GitHubReleasesProvider("PJGV333", "Chemuson")
    monkeypatch.setattr(provider, "fetch_releases_payload", _sample_releases_payload)

    stable = provider.latest_for_channel(UpdateChannel.STABLE)
    beta = provider.latest_for_channel(UpdateChannel.BETA)

    assert stable is not None
    assert stable.version == "1.1.4"
    assert beta is not None
    assert beta.version == "1.2.0-beta.2"


def test_asset_parser_picks_platform_and_links_sidecars(monkeypatch) -> None:
    provider = GitHubReleasesProvider("PJGV333", "Chemuson")
    monkeypatch.setattr(provider, "fetch_releases_payload", _sample_releases_payload)
    release = provider.latest_for_channel(UpdateChannel.STABLE)
    assert release is not None

    asset = provider.find_platform_asset(release, platform_tag="linux-x86_64")
    assert asset is not None
    assert asset.name.endswith(".AppImage")
    assert asset.sha256_url.endswith(".sha256")
    assert asset.signature_url.endswith(".sig")


def test_asset_parser_prefers_flavor_when_requested() -> None:
    provider = GitHubReleasesProvider("PJGV333", "Chemuson")
    release = provider.parse_release(
        {
            "tag_name": "v1.3.0",
            "name": "Chemuson 1.3.0",
            "prerelease": False,
            "draft": False,
            "published_at": "2026-02-22T12:00:00Z",
            "html_url": "https://example.com/release/1.3.0",
            "assets": [
                {
                    "name": "Chemuson-v1.3.0-windows-x86_64-portable.exe",
                    "browser_download_url": "https://example.com/portable.exe",
                },
                {
                    "name": "Chemuson-v1.3.0-windows-x86_64-setup.exe",
                    "browser_download_url": "https://example.com/setup.exe",
                },
            ],
        }
    )
    assert release is not None

    installer = provider.find_platform_asset(
        release,
        platform_tag="windows-x86_64",
        preferred_flavor="installer",
    )
    portable = provider.find_platform_asset(
        release,
        platform_tag="windows-x86_64",
        preferred_flavor="portable",
    )
    assert installer is not None
    assert portable is not None
    assert installer.name.endswith("-setup.exe")
    assert portable.name.endswith("-portable.exe")


def test_provider_falls_back_to_cache_when_github_unavailable(monkeypatch, tmp_path) -> None:
    provider = GitHubReleasesProvider(
        "PJGV333",
        "Chemuson",
        cache_dir=str(tmp_path),
        cache_max_age_hours=24,
    )
    payload = _sample_releases_payload()
    calls = {"count": 0}

    def _mock_read_json(_url):
        calls["count"] += 1
        if calls["count"] == 1:
            return payload
        raise OSError("offline")

    monkeypatch.setattr(provider, "_read_json", _mock_read_json)

    first = provider.fetch_releases_payload()
    assert first
    assert provider.last_fetch_source == "remote"

    second = provider.fetch_releases_payload()
    assert second
    assert provider.last_fetch_source == "cache"


def test_provider_rejects_stale_cache(monkeypatch, tmp_path) -> None:
    provider = GitHubReleasesProvider(
        "PJGV333",
        "Chemuson",
        cache_dir=str(tmp_path),
        cache_max_age_hours=1,
    )
    cache_path = provider._cache_path()
    old_stamp = (datetime.now(timezone.utc) - timedelta(hours=48)).isoformat()
    with open(cache_path, "w", encoding="utf-8") as fh:
        fh.write(json.dumps({"fetched_at": old_stamp, "payload": _sample_releases_payload()}))

    def _raise(_url):
        raise OSError("offline")

    monkeypatch.setattr(provider, "_read_json", _raise)

    try:
        provider.fetch_releases_payload()
        assert False, "Debe fallar sin usar caché vencida"
    except OSError:
        pass
    assert provider.last_fetch_source == "error"


def test_provider_can_disable_cache_fallback(monkeypatch, tmp_path) -> None:
    provider = GitHubReleasesProvider(
        "PJGV333",
        "Chemuson",
        cache_dir=str(tmp_path),
        cache_max_age_hours=24,
        allow_cached_fallback=False,
    )
    provider._write_cache(_sample_releases_payload())

    def _raise(_url):
        raise OSError("offline")

    monkeypatch.setattr(provider, "_read_json", _raise)

    try:
        provider.fetch_releases_payload()
        assert False, "No debe usar caché si el fallback está deshabilitado"
    except OSError:
        pass
    assert provider.last_fetch_source == "error"


def test_provider_rejects_insecure_api_base() -> None:
    try:
        GitHubReleasesProvider("PJGV333", "Chemuson", api_base="http://api.github.com")
        assert False, "Debe rechazar api_base sin HTTPS"
    except ValueError:
        pass


def test_provider_skips_non_https_assets() -> None:
    provider = GitHubReleasesProvider("PJGV333", "Chemuson")
    release = provider.parse_release(
        {
            "tag_name": "v1.9.0",
            "name": "Chemuson 1.9.0",
            "prerelease": False,
            "draft": False,
            "published_at": "2026-02-22T12:00:00Z",
            "html_url": "https://example.com/release/1.9.0",
            "assets": [
                {
                    "name": "Chemuson-v1.9.0-windows-x86_64-portable.exe",
                    "browser_download_url": "http://example.com/portable.exe",
                },
                {
                    "name": "Chemuson-v1.9.0-windows-x86_64-setup.exe",
                    "browser_download_url": "https://example.com/setup.exe",
                },
            ],
        }
    )
    assert release is not None
    names = [asset.name for asset in release.assets]
    assert "Chemuson-v1.9.0-windows-x86_64-portable.exe" not in names
    assert "Chemuson-v1.9.0-windows-x86_64-setup.exe" in names
