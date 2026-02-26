"""Utilidades de versión de la aplicación."""

from __future__ import annotations

from importlib.metadata import PackageNotFoundError, version as package_version

from chemuson._version import __version__


def get_app_version() -> str:
    """Devuelve la versión de la app desde una fuente única."""
    return __version__


def get_installed_version() -> str:
    """Devuelve la versión instalada del paquete o fallback local."""
    try:
        return package_version("chemuson")
    except PackageNotFoundError:
        return __version__

