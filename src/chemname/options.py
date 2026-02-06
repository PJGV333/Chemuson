"""Opciones de configuración para el motor de nombres IUPAC-lite."""

from dataclasses import dataclass


@dataclass
class NameOptions:
    """Opciones de control del algoritmo de nomenclatura."""

    # Si hay empate en locantes, ordenar por nombre del sustituyente.
    prefer_alphabetical_tiebreak: bool = True
    # Si falla el soporte, devolver "N/D" en lugar de lanzar excepción.
    return_nd_on_fail: bool = True
    # Esquema de numeración para sistemas fusionados (p. ej., pireno).
    fused_numbering_scheme: str = "iupac2004"
