"""Pruebas del núcleo de actualización (hash/firma/rollback)."""

import hashlib
import hmac
import shutil


from chemuson.update.core import AutoUpdateCore, VerificationResult
from chemuson.update.types import (
    DownloadedUpdate,
    ReleaseAsset,
    ReleaseInfo,
    UpdateCandidate,
    UpdateChannel,
    UpdateMode,
    UpdateSettings,
)
from chemuson.update.security import SignatureVerifier


class _FakeProvider:
    def __init__(self, release: ReleaseInfo):
        self._release = release

    def latest_for_channel(self, _channel):
        return self._release

    def find_platform_asset(self, release, platform_tag=None, preferred_flavor=None):
        _ = platform_tag
        _ = preferred_flavor
        return release.assets[0] if release.assets else None


def _make_release(version: str = "1.4.0") -> ReleaseInfo:
    asset = ReleaseAsset(
        name=f"Chemuson-v{version}-windows-x86_64-portable.exe",
        url="https://example.invalid/fake.exe",
    )
    return ReleaseInfo(
        tag_name=f"v{version}",
        version=version,
        prerelease=False,
        published_at="2026-02-26T00:00:00Z",
        html_url="https://example.invalid/release",
        assets=[asset],
    )


def test_check_for_updates_returns_candidate() -> None:
    release = _make_release("1.4.0")
    provider = _FakeProvider(release)
    settings = UpdateSettings(
        enabled=True,
        channel=UpdateChannel.STABLE,
        mode=UpdateMode.NOTIFY,
        check_interval_hours=24,
    )
    core = AutoUpdateCore(provider, settings)
    result = core.check_for_updates(current_version="1.3.0", platform_tag="windows-x86_64")

    assert result.available is True
    assert result.latest_version == "1.4.0"
    assert result.candidate is not None
    assert result.candidate.artifact.name.endswith(".exe")


def test_check_for_updates_blocks_downgrade_replay() -> None:
    release = _make_release("1.4.0")
    provider = _FakeProvider(release)
    settings = UpdateSettings(
        enabled=True,
        channel=UpdateChannel.STABLE,
        mode=UpdateMode.NOTIFY,
        check_interval_hours=24,
        highest_seen_version="1.5.0",
    )
    core = AutoUpdateCore(provider, settings)
    result = core.check_for_updates(current_version="1.3.0", platform_tag="windows-x86_64")
    assert result.available is False
    assert "downgrade" in result.reason.lower() or "replay" in result.reason.lower()


def test_verify_download_with_sha256_and_hmac_signature(tmp_path) -> None:
    artifact_path = tmp_path / "update.bin"
    artifact_bytes = b"chemuson-update"
    artifact_path.write_bytes(artifact_bytes)

    sha_value = hashlib.sha256(artifact_bytes).hexdigest()
    checksum_path = tmp_path / "update.bin.sha256"
    checksum_path.write_text(f"{sha_value}  update.bin\n", encoding="utf-8")

    key = "test-signing-key"
    signature = hmac.new(key.encode("utf-8"), artifact_bytes, hashlib.sha256).hexdigest()
    signature_path = tmp_path / "update.bin.sig"
    signature_path.write_text(signature, encoding="utf-8")

    release = _make_release("1.5.0")
    candidate = UpdateCandidate(release=release, artifact=release.assets[0])
    downloaded = DownloadedUpdate(
        candidate=candidate,
        artifact_path=str(artifact_path),
        checksum_path=str(checksum_path),
        signature_path=str(signature_path),
    )
    settings = UpdateSettings(enabled=True)
    core = AutoUpdateCore(
        provider=_FakeProvider(release),
        settings=settings,
        verifier=SignatureVerifier(key),
        signature_algorithm="hmac-sha256",
    )

    check = core.verify_download(downloaded)
    assert check.ok is True


def test_verify_download_skips_signature_when_no_verifier_key(tmp_path) -> None:
    artifact_path = tmp_path / "update.bin"
    artifact_bytes = b"chemuson-update"
    artifact_path.write_bytes(artifact_bytes)

    sha_value = hashlib.sha256(artifact_bytes).hexdigest()
    checksum_path = tmp_path / "update.bin.sha256"
    checksum_path.write_text(f"{sha_value}  update.bin\n", encoding="utf-8")

    # Firma inválida intencional: el MVP permite omitir firma si no hay clave.
    signature_path = tmp_path / "update.bin.sig"
    signature_path.write_text("invalid-signature", encoding="utf-8")

    release = _make_release("1.5.1")
    candidate = UpdateCandidate(release=release, artifact=release.assets[0])
    downloaded = DownloadedUpdate(
        candidate=candidate,
        artifact_path=str(artifact_path),
        checksum_path=str(checksum_path),
        signature_path=str(signature_path),
    )
    settings = UpdateSettings(enabled=True)
    core = AutoUpdateCore(
        provider=_FakeProvider(release),
        settings=settings,
        verifier=SignatureVerifier(""),
        signature_algorithm="hmac-sha256",
    )

    check = core.verify_download(downloaded)
    assert check.ok is True


def test_verify_download_requires_checksum_minimum(tmp_path) -> None:
    artifact_path = tmp_path / "update.bin"
    artifact_path.write_bytes(b"payload")
    release = _make_release("1.5.2")
    candidate = UpdateCandidate(release=release, artifact=release.assets[0])
    downloaded = DownloadedUpdate(candidate=candidate, artifact_path=str(artifact_path))
    core = AutoUpdateCore(
        provider=_FakeProvider(release),
        settings=UpdateSettings(enabled=True, require_sha256=True),
    )
    check = core.verify_download(downloaded)
    assert check.ok is False
    assert "sha-256 requerido" in check.reason.lower()


def test_verify_download_requires_signature_when_enabled(tmp_path) -> None:
    artifact_path = tmp_path / "update.bin"
    artifact_bytes = b"chemuson-update"
    artifact_path.write_bytes(artifact_bytes)
    sha_value = hashlib.sha256(artifact_bytes).hexdigest()
    checksum_path = tmp_path / "update.bin.sha256"
    checksum_path.write_text(f"{sha_value}  update.bin\n", encoding="utf-8")

    release = _make_release("1.5.3")
    candidate = UpdateCandidate(release=release, artifact=release.assets[0])
    downloaded = DownloadedUpdate(
        candidate=candidate,
        artifact_path=str(artifact_path),
        checksum_path=str(checksum_path),
        signature_path="",
    )
    core = AutoUpdateCore(
        provider=_FakeProvider(release),
        settings=UpdateSettings(enabled=True, require_sha256=True, require_signature=True),
        verifier=SignatureVerifier(""),
        signature_algorithm="hmac-sha256",
    )
    check = core.verify_download(downloaded)
    assert check.ok is False
    assert "firma requerida" in check.reason.lower()


def test_apply_update_rolls_back_when_replace_fails(tmp_path, monkeypatch) -> None:
    target_path = tmp_path / "Chemuson.exe"
    artifact_path = tmp_path / "download.exe"
    target_path.write_bytes(b"old-version")
    artifact_path.write_bytes(b"new-version")

    release = _make_release("2.0.0")
    candidate = UpdateCandidate(release=release, artifact=release.assets[0])
    downloaded = DownloadedUpdate(candidate=candidate, artifact_path=str(artifact_path))

    class _FakeRollbackManager:
        def __init__(self, rollback_dir: str):
            self.rollback_dir = rollback_dir

        def stage_backup(self, target: str) -> str:
            backup = str(tmp_path / "Chemuson.exe.bak")
            shutil.copy2(target, backup)
            return backup

        def restore(self, backup: str, target: str) -> None:
            shutil.copy2(backup, target)

    def _raise_replace(_src, _dst):
        raise OSError("replace failed")

    monkeypatch.setattr("chemuson.update.core.RollbackManager", _FakeRollbackManager)
    monkeypatch.setattr("chemuson.update.core.os.replace", _raise_replace)

    core = AutoUpdateCore(_FakeProvider(release), UpdateSettings(enabled=True))
    result = core.apply_update(downloaded, str(target_path), str(tmp_path / "rollback"))

    assert isinstance(result, VerificationResult)
    assert result.ok is False
    assert target_path.read_bytes() == b"old-version"


def test_apply_update_rolls_back_when_validation_fails(tmp_path) -> None:
    target_path = tmp_path / "Chemuson.exe"
    artifact_path = tmp_path / "download.exe"
    rollback_dir = tmp_path / "rollback"
    target_path.write_bytes(b"old-version")
    artifact_path.write_bytes(b"new-version")

    release = _make_release("2.1.0")
    candidate = UpdateCandidate(release=release, artifact=release.assets[0])
    downloaded = DownloadedUpdate(candidate=candidate, artifact_path=str(artifact_path))

    core = AutoUpdateCore(_FakeProvider(release), UpdateSettings(enabled=True))
    result = core.apply_update(
        downloaded,
        str(target_path),
        str(rollback_dir),
        validate_target=lambda _path: False,
    )

    assert result.ok is False
    assert target_path.read_bytes() == b"old-version"


def test_apply_update_cleans_backup_on_success(tmp_path) -> None:
    target_path = tmp_path / "Chemuson.exe"
    artifact_path = tmp_path / "download.exe"
    rollback_dir = tmp_path / "rollback"
    target_path.write_bytes(b"old-version")
    artifact_path.write_bytes(b"new-version")

    release = _make_release("2.2.0")
    candidate = UpdateCandidate(release=release, artifact=release.assets[0])
    downloaded = DownloadedUpdate(candidate=candidate, artifact_path=str(artifact_path))

    core = AutoUpdateCore(_FakeProvider(release), UpdateSettings(enabled=True))
    result = core.apply_update(downloaded, str(target_path), str(rollback_dir))

    assert result.ok is True
    assert target_path.read_bytes() == b"new-version"
    backups = list(rollback_dir.glob("*.bak"))
    assert backups == []
