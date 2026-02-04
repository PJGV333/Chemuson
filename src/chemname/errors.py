class ChemNameNotSupported(Exception):
    """Raised when the naming engine does not support the input."""


class ChemNameInternalError(Exception):
    """Raised for unexpected internal errors in the naming engine."""
