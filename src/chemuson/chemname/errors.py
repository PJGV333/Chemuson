"""Excepciones específicas del motor de nomenclatura."""


class ChemNameNotSupported(Exception):
    """Se lanza cuando la molécula está fuera del alcance soportado."""


class ChemNameInternalError(Exception):
    """Se lanza ante errores internos inesperados del motor."""
