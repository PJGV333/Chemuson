"""Rollback básico para reemplazos de binarios."""

from __future__ import annotations

import os
import shutil
from datetime import datetime, timezone


class RollbackManager:
    """Gestiona backups de actualización y restauración de emergencia."""

    def __init__(self, rollback_dir: str) -> None:
        self.rollback_dir = os.path.abspath(rollback_dir)
        os.makedirs(self.rollback_dir, exist_ok=True)

    def stage_backup(self, target_path: str) -> str:
        """Crea backup del objetivo y retorna su ruta."""
        src = os.path.abspath(target_path)
        if not os.path.exists(src):
            raise FileNotFoundError(src)
        stamp = datetime.now(timezone.utc).strftime("%Y%m%d_%H%M%S")
        base = os.path.basename(src)
        backup_path = os.path.join(self.rollback_dir, f"{base}.{stamp}.bak")
        shutil.copy2(src, backup_path)
        return backup_path

    def restore(self, backup_path: str, target_path: str) -> None:
        """Restaura backup sobre el objetivo."""
        src = os.path.abspath(backup_path)
        dst = os.path.abspath(target_path)
        if not os.path.exists(src):
            raise FileNotFoundError(src)
        os.replace(src, dst)

    def cleanup(self, backup_path: str) -> None:
        """Elimina backup luego de validar update."""
        path = os.path.abspath(backup_path)
        if os.path.exists(path):
            os.remove(path)

