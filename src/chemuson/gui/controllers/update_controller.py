from __future__ import annotations

import os
from dataclasses import dataclass, field
from datetime import datetime
from typing import Callable

from PyQt6.QtCore import QSettings
from PyQt6.QtWidgets import QMessageBox, QWidget

from chemuson.update import (
    AutoUpdateCore,
    GitHubReleasesProvider,
    PortableUpdateContext,
    UpdateChannel,
    UpdateMode,
    UpdateSettings,
    UpdateTelemetryLogger,
    choose_windows_asset_flavor,
    detect_portable_update_context,
    detect_windows_install_context,
    is_portable_target_writable,
    is_windows_installer_asset,
    launch_portable_update_script,
    launch_inno_installer,
    mark_checked,
    should_check_now,
    write_portable_update_script,
)


FLATPAK_APP_ID = "io.github.PJGV333.Chemuson"


def _empty_portable_context() -> PortableUpdateContext:
    return PortableUpdateContext(
        is_portable=False,
        is_windows=False,
        is_appimage=False,
        target_path="",
        executable_path="",
        display_name="",
    )


def channel_display_name(channel: str) -> str:
    """Devuelve un nombre legible para el canal de updates."""
    value = str(channel or "").strip().lower()
    if value == "stable":
        return "estable"
    if value == "beta":
        return "beta"
    return value or "desconocido"


def update_source_display_name(source: str) -> str:
    """Devuelve un nombre legible para el origen del feed de updates."""
    value = str(source or "").strip().lower()
    if value == "cache":
        return "caché local"
    if value == "remote":
        return "GitHub"
    if value == "error":
        return "error"
    return ""


def format_no_update_message(
    current_version: str,
    channel: str,
    reason: str = "",
    latest_version: str = "",
    source: str = "",
) -> str:
    """Construye mensaje explicativo cuando no hay update elegible."""
    lines = [
        "No hay una version publicada mas nueva para tu canal actual.",
        "",
        f"Version instalada: {str(current_version or '').strip() or 'desconocida'}",
        f"Canal: {channel_display_name(channel)}",
    ]
    latest = str(latest_version or "").strip()
    if latest:
        lines.append(f"Ultima version consultada: {latest}")
    source_name = update_source_display_name(source)
    if source_name:
        lines.append(f"Origen de datos: {source_name}")
    lines.extend(
        [
            "",
            "Este verificador compara releases publicadas; no distribuye commits sueltos.",
        ]
    )
    if str(source or "").strip().lower() == "cache":
        lines.extend(
            [
                "",
                "Aviso: el resultado proviene de la caché local. Si acabas de publicar una release, "
                "vuelve a intentarlo cuando GitHub esté accesible.",
            ]
        )
    detail = str(reason or "").strip()
    if detail:
        lines.extend(["", f"Detalle: {detail}"])
    return "\n".join(lines)


def is_running_in_flatpak() -> bool:
    """Detecta si Chemuson corre dentro de un sandbox Flatpak."""
    flatpak_id = str(os.getenv("FLATPAK_ID", "") or "").strip()
    if flatpak_id:
        return True
    return os.path.exists("/.flatpak-info")


def format_update_disabled_message(flatpak: bool = False, app_id: str = FLATPAK_APP_ID) -> str:
    """Construye mensaje cuando el chequeo interno de updates está deshabilitado."""
    if flatpak:
        return (
            "Esta edicion Flatpak no usa auto-update real dentro de Chemuson.\n"
            "Se actualiza con Flatpak, no desde la propia app.\n\n"
            f"Usa:\nflatpak update {app_id}\n\n"
            "Si instalaste desde un bundle local sin un remote configurado, "
            "instala el bundle mas reciente manualmente."
        )
    return "El chequeo de actualizaciones está deshabilitado por entorno."


@dataclass(slots=True)
class UpdateControllerContext:
    """Dependencias mínimas de UI para el flujo de updates."""

    parent: QWidget
    show_status: Callable[[str, int], None]
    close_window: Callable[[], bool]


@dataclass(slots=True)
class UpdatePendingState:
    """Estado transitorio de updates preparados para aplicar al salir."""

    windows_installer_path: str = ""
    windows_installer_version: str = ""
    windows_download: object | None = None
    portable_target_path: str = ""
    portable_version: str = ""
    portable_download: object | None = None
    portable_relaunch: bool = False
    portable_context: PortableUpdateContext = field(default_factory=_empty_portable_context)


class UpdateController:
    """Gestiona chequeo, cola y aplicación diferida de actualizaciones."""

    def __init__(self, settings: QSettings, app_version: str) -> None:
        self._settings_store = settings
        self._app_version = app_version
        self._settings = self._load_preferences()
        self._windows_install_context = detect_windows_install_context()
        self._portable_update_context = detect_portable_update_context(
            windows_context=self._windows_install_context
        )
        self._pending = UpdatePendingState()
        self._telemetry = UpdateTelemetryLogger()

    @property
    def settings(self) -> UpdateSettings:
        return self._settings

    @staticmethod
    def setting_bool(value, default: bool) -> bool:
        """Normaliza valores de QSettings a booleano."""
        if value is None:
            return bool(default)
        if isinstance(value, bool):
            return value
        if isinstance(value, str):
            lowered = value.strip().lower()
            if lowered in {"1", "true", "yes", "on", "si", "sí"}:
                return True
            if lowered in {"0", "false", "no", "off"}:
                return False
        try:
            return bool(int(value))
        except Exception:
            return bool(value)

    def settings_payload(self) -> dict:
        """Devuelve preferencias de update para precargar en diálogo."""
        return {
            "enabled": bool(self._settings.enabled),
            "channel": str(self._settings.channel.value),
            "mode": str(self._settings.mode.value),
            "check_interval_hours": int(max(1, self._settings.check_interval_hours)),
        }

    def save_preferences(self) -> None:
        """Persiste preferencias de actualización en QSettings."""
        self._settings_store.setValue("update/enabled", bool(self._settings.enabled))
        self._settings_store.setValue("update/channel", str(self._settings.channel.value))
        self._settings_store.setValue("update/mode", str(self._settings.mode.value))
        self._settings_store.setValue(
            "update/check_interval_hours",
            int(max(1, self._settings.check_interval_hours)),
        )
        self._settings_store.setValue("update/last_check_iso", str(self._settings.last_check_iso or ""))
        self._settings_store.setValue(
            "update/highest_seen_version",
            str(self._settings.highest_seen_version or ""),
        )
        self._settings_store.setValue("update/require_sha256", bool(self._settings.require_sha256))
        self._settings_store.setValue(
            "update/require_signature",
            bool(self._settings.require_signature),
        )

    def apply_preferences(self, prefs: dict) -> None:
        """Actualiza settings de update desde el diálogo de preferencias."""
        update_enabled = self.setting_bool(
            prefs.get("update_enabled", self._settings.enabled),
            self._settings.enabled,
        )
        update_channel_raw = str(
            prefs.get("update_channel", self._settings.channel.value) or self._settings.channel.value
        ).strip().lower()
        update_mode_raw = str(
            prefs.get("update_mode", self._settings.mode.value) or self._settings.mode.value
        ).strip().lower()
        try:
            update_interval = int(
                prefs.get("update_check_interval_hours", self._settings.check_interval_hours)
            )
        except Exception:
            update_interval = self._settings.check_interval_hours
        try:
            channel = UpdateChannel(update_channel_raw)
        except Exception:
            channel = UpdateChannel.STABLE
        try:
            mode = UpdateMode(update_mode_raw)
        except Exception:
            mode = UpdateMode.NOTIFY
        self._settings.enabled = bool(update_enabled)
        self._settings.channel = channel
        self._settings.mode = mode
        self._settings.check_interval_hours = max(1, int(update_interval))
        self.save_preferences()

    def check_for_updates(
        self,
        context: UpdateControllerContext,
        *,
        force: bool,
        interactive: bool,
    ) -> None:
        """Ejecuta chequeo de updates aplicando política y canal."""
        self._log_event(
            "check_start",
            force=bool(force),
            interactive=bool(interactive),
            channel=self._settings.channel.value,
            mode=self._settings.mode.value,
        )
        if os.getenv("CHEMUSON_DISABLE_UPDATE_CHECK", "").strip().lower() in {"1", "true", "yes"}:
            flatpak_runtime = is_running_in_flatpak()
            self._log_event(
                "check_skipped",
                reason_code="disabled_flatpak" if flatpak_runtime else "disabled_env",
                channel=self._settings.channel.value,
                mode=self._settings.mode.value,
            )
            if interactive:
                QMessageBox.information(
                    context.parent,
                    "Actualizaciones",
                    format_update_disabled_message(flatpak=flatpak_runtime),
                )
            return
        if not force and not should_check_now(self._settings):
            self._log_event(
                "check_skipped",
                reason_code="interval_not_elapsed",
                channel=self._settings.channel.value,
                mode=self._settings.mode.value,
            )
            return
        mark_checked(self._settings)
        self.save_preferences()
        provider = None
        try:
            manual_check = bool(force and interactive)
            provider = GitHubReleasesProvider(
                "PJGV333",
                "Chemuson",
                timeout=8.0 if manual_check else 4.0,
                allow_cached_fallback=not manual_check,
            )
            updater = AutoUpdateCore(provider, self._settings)
            result = updater.check_for_updates(
                self._app_version,
                preferred_asset_flavor=self._preferred_update_asset_flavor(),
            )
        except Exception as exc:
            self._log_event(
                "check_error",
                error_type=exc.__class__.__name__,
                channel=self._settings.channel.value,
                mode=self._settings.mode.value,
            )
            if self._settings.mode == UpdateMode.NOTIFY or interactive:
                context.show_status(f"No se pudo comprobar actualización: {exc}", 10000)
                QMessageBox.warning(
                    context.parent,
                    "Actualizaciones",
                    f"No se pudo comprobar actualización:\n{exc}",
                )
            return
        source = ""
        if provider is not None:
            source = str(getattr(provider, "last_fetch_source", "") or "")
        if source == "cache":
            context.show_status(
                "GitHub no disponible: se usó caché local de actualizaciones.",
                12000,
            )
        self._log_event(
            "check_result",
            available=bool(getattr(result, "available", False)),
            current_version=str(getattr(result, "current_version", "") or ""),
            latest_version=str(getattr(result, "latest_version", "") or ""),
            source=source or "unknown",
            channel=self._settings.channel.value,
            mode=self._settings.mode.value,
        )
        self.save_preferences()
        self._handle_update_check_result(
            context,
            result,
            interactive=interactive,
            source=source,
        )

    def apply_pending_windows_update_on_exit(self, context: UpdateControllerContext) -> bool:
        """Ejecuta instalador pendiente antes de salir cuando corresponde."""
        installer_path = str(self._pending.windows_installer_path or "").strip()
        if not installer_path:
            return True
        clear_pending = False
        try:
            downloaded = self._pending.windows_download
            if downloaded is None:
                raise RuntimeError("No hay metadata de verificación del instalador pendiente.")
            provider = GitHubReleasesProvider("PJGV333", "Chemuson", timeout=6.0)
            updater = AutoUpdateCore(provider, self._settings)
            verify_result = updater.verify_download(downloaded)
            if not verify_result.ok:
                raise RuntimeError(
                    f"No se ejecutará artefacto no verificado: {verify_result.reason}"
                )
            launch_inno_installer(installer_path, silent=True)
            self._log_event(
                "apply_queued",
                latest_version=self._pending.windows_installer_version or "",
                context="windows_installer",
                result="launched",
            )
            clear_pending = True
            return True
        except Exception as exc:
            self._log_event(
                "apply_failed",
                latest_version=self._pending.windows_installer_version or "",
                context="windows_installer",
                error_type=exc.__class__.__name__,
            )
            reply = QMessageBox.question(
                context.parent,
                "Actualización pendiente",
                (
                    "No se pudo lanzar el instalador de actualización pendiente.\n"
                    f"Error: {exc}\n\n"
                    "¿Deseas salir de todas formas sin aplicar la actualización?"
                ),
                QMessageBox.StandardButton.Yes | QMessageBox.StandardButton.No,
                QMessageBox.StandardButton.No,
            )
            if reply == QMessageBox.StandardButton.Yes:
                clear_pending = True
                return True
            return False
        finally:
            if clear_pending:
                self._pending.windows_download = None
                self._pending.windows_installer_path = ""
                self._pending.windows_installer_version = ""

    def apply_pending_portable_update_on_exit(self, context: UpdateControllerContext) -> bool:
        """Lanza helper para reemplazar AppImage/binario portable al cerrar."""
        target_path = str(self._pending.portable_target_path or "").strip()
        if not target_path:
            return True
        clear_pending = False
        try:
            downloaded = self._pending.portable_download
            if downloaded is None:
                raise RuntimeError("No hay metadata de verificación del binario pendiente.")
            provider = GitHubReleasesProvider("PJGV333", "Chemuson", timeout=6.0)
            updater = AutoUpdateCore(provider, self._settings)
            verify_result = updater.verify_download(downloaded)
            if not verify_result.ok:
                raise RuntimeError(
                    f"No se reemplazará binario no verificado: {verify_result.reason}"
                )
            version = str(self._pending.portable_version or "latest").strip() or "latest"
            script_dir = os.path.dirname(str(getattr(downloaded, "artifact_path", "") or ""))
            rollback_dir = os.path.join(os.path.expanduser("~"), ".chemuson", "rollback")
            os.makedirs(rollback_dir, exist_ok=True)
            backup_name = (
                f"{os.path.basename(target_path)}."
                f"{datetime.utcnow().strftime('%Y%m%d_%H%M%S')}.bak"
            )
            backup_path = os.path.join(rollback_dir, backup_name)
            log_path = os.path.join(script_dir, "portable-update.log")
            script_name = (
                "apply-portable-update.cmd"
                if self._pending.portable_context.is_windows
                else "apply-portable-update.sh"
            )
            script_path = os.path.join(script_dir, script_name)
            write_portable_update_script(
                script_path,
                source_path=str(getattr(downloaded, "artifact_path", "") or ""),
                target_path=target_path,
                backup_path=backup_path,
                log_path=log_path,
                pid=os.getpid(),
                relaunch=bool(self._pending.portable_relaunch),
                platform_name="windows" if self._pending.portable_context.is_windows else "linux",
            )
            launch_portable_update_script(
                script_path,
                platform_name="windows" if self._pending.portable_context.is_windows else "linux",
            )
            self._log_event(
                "apply_queued",
                latest_version=version,
                context="portable_binary",
                result="helper_launched",
                target_kind="appimage" if self._pending.portable_context.is_appimage else "portable",
            )
            clear_pending = True
            return True
        except Exception as exc:
            self._log_event(
                "apply_failed",
                latest_version=self._pending.portable_version or "",
                context="portable_binary",
                error_type=exc.__class__.__name__,
            )
            reply = QMessageBox.question(
                context.parent,
                "Actualización pendiente",
                (
                    "No se pudo preparar el reemplazo del ejecutable pendiente.\n"
                    f"Error: {exc}\n\n"
                    "¿Deseas salir de todas formas sin aplicar la actualización?"
                ),
                QMessageBox.StandardButton.Yes | QMessageBox.StandardButton.No,
                QMessageBox.StandardButton.No,
            )
            if reply == QMessageBox.StandardButton.Yes:
                clear_pending = True
                return True
            return False
        finally:
            if clear_pending:
                self._pending = UpdatePendingState()

    def _load_preferences(self) -> UpdateSettings:
        enabled = self.setting_bool(self._settings_store.value("update/enabled", False), False)
        raw_channel = str(self._settings_store.value("update/channel", "stable") or "stable").strip().lower()
        raw_mode = str(self._settings_store.value("update/mode", "notify") or "notify").strip().lower()
        try:
            interval_hours = int(self._settings_store.value("update/check_interval_hours", 24) or 24)
        except Exception:
            interval_hours = 24
        require_sha256 = self.setting_bool(self._settings_store.value("update/require_sha256", True), True)
        require_signature = self.setting_bool(
            self._settings_store.value("update/require_signature", False),
            False,
        )
        highest_seen_version = str(
            self._settings_store.value("update/highest_seen_version", "") or ""
        ).strip()
        try:
            channel = UpdateChannel(raw_channel)
        except Exception:
            channel = UpdateChannel.STABLE
        try:
            mode = UpdateMode(raw_mode)
        except Exception:
            mode = UpdateMode.NOTIFY
        last_check_iso = str(self._settings_store.value("update/last_check_iso", "") or "").strip()
        return UpdateSettings(
            enabled=bool(enabled),
            channel=channel,
            mode=mode,
            check_interval_hours=max(1, interval_hours),
            last_check_iso=last_check_iso,
            highest_seen_version=highest_seen_version,
            require_sha256=bool(require_sha256),
            require_signature=bool(require_signature),
        )

    def _preferred_update_asset_flavor(self) -> str | None:
        self._windows_install_context = detect_windows_install_context()
        self._portable_update_context = detect_portable_update_context(
            windows_context=self._windows_install_context
        )
        flavor = choose_windows_asset_flavor(self._windows_install_context)
        return flavor or None

    def _current_portable_update_context(self) -> PortableUpdateContext:
        self._windows_install_context = detect_windows_install_context()
        self._portable_update_context = detect_portable_update_context(
            windows_context=self._windows_install_context
        )
        return self._portable_update_context

    def _log_event(self, event: str, **fields) -> None:
        try:
            self._telemetry.log_event(event, **fields)
        except Exception:
            return

    def _download_update_candidate(self, candidate):
        version = str(getattr(candidate.release, "version", "latest") or "latest")
        self._log_event(
            "download_start",
            latest_version=version,
            channel=self._settings.channel.value,
            mode=self._settings.mode.value,
        )
        provider = GitHubReleasesProvider("PJGV333", "Chemuson", timeout=12.0)
        updater = AutoUpdateCore(provider, self._settings)
        download_dir = os.path.join(os.path.expanduser("~"), ".chemuson", "updates", version)
        try:
            downloaded = updater.download_candidate(candidate, download_dir)
            verification = updater.verify_download(downloaded)
            if not verification.ok:
                reason = verification.reason or "Falló verificación del paquete descargado."
                self._log_event(
                    "download_failed",
                    latest_version=version,
                    reason_code="verification_failed",
                )
                raise RuntimeError(reason)
            self._log_event("download_ok", latest_version=version)
            return downloaded
        except Exception as exc:
            self._log_event(
                "download_failed",
                latest_version=version,
                error_type=exc.__class__.__name__,
            )
            raise

    def _queue_windows_installer_update(
        self,
        context: UpdateControllerContext,
        candidate,
        show_errors: bool = True,
    ) -> bool:
        candidate_version = str(getattr(candidate.release, "version", "") or "")
        if (
            self._pending.windows_installer_path
            and self._pending.windows_installer_version == candidate_version
            and os.path.exists(self._pending.windows_installer_path)
            and self._pending.windows_download is not None
        ):
            self._log_event(
                "queue_reused",
                latest_version=candidate_version,
                context="windows_installer",
            )
            return True
        try:
            downloaded = self._download_update_candidate(candidate)
        except Exception as exc:
            self._log_event(
                "queue_failed",
                latest_version=candidate_version,
                context="windows_installer",
                error_type=exc.__class__.__name__,
            )
            if show_errors:
                QMessageBox.warning(
                    context.parent,
                    "Actualización",
                    f"No se pudo preparar el instalador de actualización:\n{exc}",
                )
            else:
                context.show_status(
                    f"No se pudo preparar la instalación silenciosa de actualización: {exc}",
                    12000,
                )
            return False
        self._pending.windows_download = downloaded
        self._pending.windows_installer_path = str(getattr(downloaded, "artifact_path", "") or "")
        self._pending.windows_installer_version = candidate_version
        self._log_event(
            "queue_ok",
            latest_version=candidate_version,
            context="windows_installer",
        )
        return True

    def _queue_portable_binary_update(
        self,
        context: UpdateControllerContext,
        candidate,
        portable_context: PortableUpdateContext,
        show_errors: bool = True,
    ) -> bool:
        candidate_version = str(getattr(candidate.release, "version", "") or "")
        target_path = str(getattr(portable_context, "target_path", "") or "").strip()
        if (
            self._pending.portable_target_path
            and self._pending.portable_target_path == target_path
            and self._pending.portable_version == candidate_version
            and self._pending.portable_download is not None
            and os.path.exists(str(getattr(self._pending.portable_download, "artifact_path", "") or ""))
        ):
            self._log_event(
                "queue_reused",
                latest_version=candidate_version,
                context="portable_binary",
            )
            return True
        try:
            if not portable_context.can_self_update:
                raise RuntimeError("No se pudo determinar la ruta del ejecutable actual.")
            if not is_portable_target_writable(target_path):
                raise PermissionError(
                    "Chemuson no tiene permisos para reemplazar el ejecutable actual."
                )
            downloaded = self._download_update_candidate(candidate)
        except Exception as exc:
            self._log_event(
                "queue_failed",
                latest_version=candidate_version,
                context="portable_binary",
                error_type=exc.__class__.__name__,
            )
            if show_errors:
                QMessageBox.warning(
                    context.parent,
                    "Actualización",
                    f"No se pudo preparar la actualización portable:\n{exc}",
                )
            else:
                context.show_status(
                    f"No se pudo preparar el auto-update real: {exc}",
                    12000,
                )
            return False
        self._pending.portable_download = downloaded
        self._pending.portable_target_path = target_path
        self._pending.portable_version = candidate_version
        self._pending.portable_relaunch = False
        self._pending.portable_context = portable_context
        self._log_event(
            "queue_ok",
            latest_version=candidate_version,
            context="portable_binary",
            target_kind="appimage" if portable_context.is_appimage else "portable",
        )
        return True

    def _offer_windows_installer_update(
        self,
        context: UpdateControllerContext,
        result,
        interactive: bool,
    ) -> None:
        candidate = getattr(result, "candidate", None)
        if candidate is None:
            return
        version = str(getattr(result, "latest_version", "") or "")
        if self._settings.mode == UpdateMode.SILENT and not interactive:
            if self._queue_windows_installer_update(context, candidate, show_errors=False):
                context.show_status(
                    f"Instalación silenciosa {version} lista para aplicarse al cerrar Chemuson.",
                    20000,
                )
            return

        reply = QMessageBox.question(
            context.parent,
            "Actualización disponible",
            (
                f"Hay una nueva versión disponible ({version}).\n\n"
                "Esta edición no usa auto-update real.\n"
                "Chemuson descargará el instalador oficial y lo ejecutará en silencio al cerrar.\n\n"
                "¿Quieres prepararlo ahora?"
            ),
            QMessageBox.StandardButton.Yes | QMessageBox.StandardButton.No,
            QMessageBox.StandardButton.Yes,
        )
        if reply != QMessageBox.StandardButton.Yes:
            return
        if not self._queue_windows_installer_update(context, candidate):
            return
        context.show_status(
            f"Instalador {version} preparado. Se ejecutará en silencio al cerrar Chemuson.",
            20000,
        )
        if interactive:
            close_now = QMessageBox.question(
                context.parent,
                "Aplicar actualización",
                "El instalador de actualización está listo.\n"
                "¿Deseas cerrar Chemuson ahora para ejecutar la instalación silenciosa?",
                QMessageBox.StandardButton.Yes | QMessageBox.StandardButton.No,
                QMessageBox.StandardButton.No,
            )
            if close_now == QMessageBox.StandardButton.Yes:
                context.close_window()

    def _offer_portable_binary_update(
        self,
        context: UpdateControllerContext,
        result,
        interactive: bool,
        portable_context: PortableUpdateContext,
    ) -> None:
        candidate = getattr(result, "candidate", None)
        if candidate is None:
            return
        version = str(getattr(result, "latest_version", "") or "")
        target_label = "AppImage" if portable_context.is_appimage else "ejecutable portable"
        if self._settings.mode == UpdateMode.SILENT and not interactive:
            if self._queue_portable_binary_update(context, candidate, portable_context, show_errors=False):
                context.show_status(
                    f"Auto-update real {version} listo: Chemuson reemplazará el {target_label} al cerrar.",
                    20000,
                )
            return

        reply = QMessageBox.question(
            context.parent,
            "Actualización disponible",
            (
                f"Hay una nueva versión disponible ({version}).\n\n"
                "Esto es auto-update real.\n"
                f"Chemuson descargará la actualización y reemplazará el {target_label} actual al cerrar.\n\n"
                "¿Quieres prepararla ahora?"
            ),
            QMessageBox.StandardButton.Yes | QMessageBox.StandardButton.No,
            QMessageBox.StandardButton.Yes,
        )
        if reply != QMessageBox.StandardButton.Yes:
            return
        if not self._queue_portable_binary_update(context, candidate, portable_context):
            return
        context.show_status(
            f"Auto-update real {version} preparado. Se aplicará al cerrar Chemuson.",
            20000,
        )
        if interactive:
            close_now = QMessageBox.question(
                context.parent,
                "Aplicar actualización",
                "El auto-update real está listo.\n"
                "¿Deseas cerrar Chemuson ahora para completar el reemplazo?",
                QMessageBox.StandardButton.Yes | QMessageBox.StandardButton.No,
                QMessageBox.StandardButton.No,
            )
            if close_now == QMessageBox.StandardButton.Yes:
                self._pending.portable_relaunch = True
                if not context.close_window():
                    self._pending.portable_relaunch = False

    def _handle_update_check_result(
        self,
        context: UpdateControllerContext,
        result,
        *,
        interactive: bool = False,
        source: str = "",
    ) -> None:
        if not getattr(result, "available", False):
            if interactive:
                QMessageBox.information(
                    context.parent,
                    "Actualizaciones",
                    format_no_update_message(
                        str(getattr(result, "current_version", "") or self._app_version),
                        str(getattr(result, "channel", "") or self._settings.channel.value),
                        str(getattr(result, "reason", "") or ""),
                        str(getattr(result, "latest_version", "") or ""),
                        source,
                    ),
                )
            return
        candidate = getattr(result, "candidate", None)
        artifact_name = ""
        if candidate is not None and getattr(candidate, "artifact", None) is not None:
            artifact_name = str(candidate.artifact.name)

        if (
            self._windows_install_context.is_windows
            and self._windows_install_context.installed
            and artifact_name
            and is_windows_installer_asset(artifact_name)
        ):
            self._offer_windows_installer_update(context, result, interactive=interactive)
            return

        portable_context = self._current_portable_update_context()
        if candidate is not None and portable_context.can_self_update:
            self._offer_portable_binary_update(
                context,
                result,
                interactive=interactive,
                portable_context=portable_context,
            )
            return

        message = (
            f"Actualización disponible: {result.latest_version}"
            + (f" ({artifact_name})" if artifact_name else "")
        )
        context.show_status(message, 15000)
        if self._settings.mode == UpdateMode.NOTIFY or interactive:
            release_url = ""
            if candidate is not None and getattr(candidate, "release", None) is not None:
                release_url = str(candidate.release.html_url or "")
            source_line = f"\n\nOrigen de datos: {'caché local' if source == 'cache' else 'GitHub'}"
            extra = f"\n\nRelease: {release_url}" if release_url else ""
            QMessageBox.information(
                context.parent,
                "Actualización disponible",
                f"Hay una nueva versión disponible ({result.latest_version}).{extra}{source_line}",
            )
