"""Helpers for resolving packaged resources."""

from importlib.resources import as_file, files


def open_resource_path(*parts: str, package: str = "chemuson"):
    """Return a context manager yielding a filesystem path to a resource."""
    target = files(package).joinpath(*parts)
    return as_file(target)
