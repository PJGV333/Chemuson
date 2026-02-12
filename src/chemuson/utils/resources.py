"""Helpers para resolver recursos empaquetados."""

from importlib.resources import as_file, files


def open_resource_path(*parts: str, package: str = "chemuson"):
    """Retorna un context manager con path real a un recurso empacado."""
    target = files(package).joinpath(*parts)
    return as_file(target)
