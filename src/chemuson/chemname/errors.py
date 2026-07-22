"""Excepciones específicas del motor de nomenclatura."""

from chemuson.core.molview import MolecularViewNotSupported


ChemNameNotSupported = MolecularViewNotSupported


class ChemNameInternalError(Exception):
    """Se lanza ante errores internos inesperados del motor."""
