"""Núcleo de auto-update (MVP) para Chemuson."""

from __future__ import annotations

import os
from dataclasses import dataclass
from typing import Callable, Optional

from chemuson.update.policy import can_offer_update, is_downgrade_suspected
from chemuson.update.semver import compare_versions
from chemuson.update.provider import GitHubReleasesProvider, detect_platform_tag
from chemuson.update.rollback import RollbackManager
from chemuson.update.security import (
    SignatureVerifier,
    download_binary,
    parse_sha256_text,
    verify_file_hash,
)
from chemuson.update.types import (
    DownloadedUpdate,
    UpdateCandidate,
    UpdateCheckResult,
    UpdateSettings,
)


@dataclass(slots=True)
class VerificationResult:
    """Resultado de verificación de descarga."""

    ok: bool
    reason: str = ""


class AutoUpdateCore:
    """Orquesta chequeo, descarga, verificación y aplicación con rollback."""

    def __init__(
        self,
        provider: GitHubReleasesProvider,
        settings: UpdateSettings,
        verifier: Optional[SignatureVerifier] = None,
        signature_algorithm: str = "hmac-sha256",
    ) -> None:
        self.provider = provider
        self.settings = settings
        self.verifier = verifier or SignatureVerifier("")
        self.signature_algorithm = signature_algorithm

    def check_for_updates(
        self,
        current_version: str,
        platform_tag: str | None = None,
        preferred_asset_flavor: str | None = None,
    ) -> UpdateCheckResult:
        """Busca actualización disponible para versión/plataforma."""
        channel = self.settings.channel.value
        release = self.provider.latest_for_channel(channel)
        if release is None:
            return UpdateCheckResult(
                available=False,
                current_version=current_version,
                channel=channel,
                reason="No se encontraron releases candidatas.",
            )
        latest_version = release.version
        if is_downgrade_suspected(latest_version, self.settings.highest_seen_version):
            return UpdateCheckResult(
                available=False,
                current_version=current_version,
                latest_version=latest_version,
                channel=channel,
                reason="Se detectó posible downgrade/replay del feed remoto.",
            )
        previous_highest = str(self.settings.highest_seen_version or "").strip()
        if not previous_highest:
            self.settings.highest_seen_version = latest_version
        else:
            try:
                if compare_versions(latest_version, previous_highest) > 0:
                    self.settings.highest_seen_version = latest_version
            except Exception:
                pass
        if not can_offer_update(
            channel=channel,
            current_version=current_version,
            candidate_version=latest_version,
        ):
            return UpdateCheckResult(
                available=False,
                current_version=current_version,
                latest_version=latest_version,
                channel=channel,
                reason="La versión remota no es elegible para el canal o no es más nueva.",
            )
        asset = self.provider.find_platform_asset(
            release,
            platform_tag=platform_tag,
            preferred_flavor=preferred_asset_flavor,
        )
        if asset is None:
            return UpdateCheckResult(
                available=False,
                current_version=current_version,
                latest_version=latest_version,
                channel=channel,
                reason=f"No hay asset compatible con {platform_tag or detect_platform_tag()}.",
            )
        candidate = UpdateCandidate(release=release, artifact=asset)
        return UpdateCheckResult(
            available=True,
            current_version=current_version,
            latest_version=latest_version,
            channel=channel,
            reason="",
            candidate=candidate,
        )

    def download_candidate(self, candidate: UpdateCandidate, download_dir: str) -> DownloadedUpdate:
        """Descarga asset y sidecars de verificación."""
        os.makedirs(download_dir, exist_ok=True)
        artifact_path = os.path.join(download_dir, candidate.artifact.name)
        download_binary(candidate.artifact.url, artifact_path)
        checksum_path = ""
        signature_path = ""

        if candidate.artifact.sha256_url:
            checksum_path = artifact_path + ".sha256"
            download_binary(candidate.artifact.sha256_url, checksum_path)
        if candidate.artifact.signature_url:
            signature_path = artifact_path + ".sig"
            download_binary(candidate.artifact.signature_url, signature_path)
        return DownloadedUpdate(
            candidate=candidate,
            artifact_path=artifact_path,
            checksum_path=checksum_path,
            signature_path=signature_path,
        )

    def verify_download(self, downloaded: DownloadedUpdate) -> VerificationResult:
        """Valida hash y firma del update descargado."""
        require_sha256 = bool(self.settings.require_sha256)
        if not downloaded.checksum_path and require_sha256:
            return VerificationResult(False, "SHA-256 requerido y no se encontró sidecar de checksum.")
        if downloaded.checksum_path:
            try:
                with open(downloaded.checksum_path, "r", encoding="utf-8") as fh:
                    checksum_text = fh.read()
            except Exception as exc:
                return VerificationResult(False, f"No se pudo leer checksum: {exc}")
            expected_sha256 = parse_sha256_text(checksum_text)
            if not expected_sha256:
                return VerificationResult(False, "Archivo de checksum inválido.")
            if not verify_file_hash(downloaded.artifact_path, expected_sha256):
                return VerificationResult(False, "Hash SHA-256 no coincide.")

        require_signature = bool(self.settings.require_signature)
        if not downloaded.signature_path and require_signature:
            return VerificationResult(False, "Firma requerida y no se encontró sidecar de firma.")
        if downloaded.signature_path:
            # En MVP la firma es opcional en cliente cuando no hay clave pública
            # configurada. Si hay clave, la firma se valida estrictamente.
            if self.verifier.is_configured(self.signature_algorithm):
                try:
                    with open(downloaded.signature_path, "r", encoding="utf-8") as fh:
                        signature_value = fh.read().strip()
                except Exception as exc:
                    return VerificationResult(False, f"No se pudo leer firma: {exc}")
                if not self.verifier.verify_file(
                    downloaded.artifact_path,
                    signature_value,
                    algorithm=self.signature_algorithm,
                ):
                    return VerificationResult(False, "Firma inválida.")
            elif require_signature:
                return VerificationResult(
                    False,
                    "Firma requerida, pero no hay verificador criptográfico configurado.",
                )
        return VerificationResult(True, "")

    def apply_update(
        self,
        downloaded: DownloadedUpdate,
        target_path: str,
        rollback_dir: str,
        validate_target: Callable[[str], bool] | None = None,
        cleanup_backup_on_success: bool = True,
    ) -> VerificationResult:
        """Aplica update sobre `target_path` con rollback básico."""
        manager = RollbackManager(rollback_dir)
        backup_path = ""
        target_abs = os.path.abspath(target_path)
        artifact_abs = os.path.abspath(downloaded.artifact_path)
        try:
            if os.path.exists(target_abs):
                backup_path = manager.stage_backup(target_abs)
            os.replace(artifact_abs, target_abs)
            if validate_target is not None:
                if not bool(validate_target(target_abs)):
                    raise RuntimeError("Validación posterior a update falló.")
            if backup_path and cleanup_backup_on_success:
                try:
                    manager.cleanup(backup_path)
                except Exception:
                    pass
            return VerificationResult(True, "")
        except Exception as exc:
            if backup_path:
                try:
                    manager.restore(backup_path, target_abs)
                except Exception:
                    pass
            return VerificationResult(False, f"Error aplicando update: {exc}")
