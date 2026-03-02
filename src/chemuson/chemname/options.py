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
    # Modo estricto: ante ambigüedad, preferir excepción en lugar de degradar.
    strict: bool = False
    # Habilita nomenclatura de coordinación (experimental).
    allow_coordination: bool = True
    # Habilita heteroaromáticos exóticos (P/Si/B) en rutas experimentales.
    enable_exotic_hetero: bool = True
    # Habilita plantillas especiales (carbohidratos/esteroides).
    enable_special_templates: bool = True
    # Habilita estereoquímica avanzada (R_a/S_a, M/P, endo/exo, si/re).
    enable_advanced_stereo: bool = True
    # Master switch para funcionalidades experimentales.
    enable_experimental: bool = True
    # Ejecuta extracción estereo RDKit en subproceso para aislar fallos nativos.
    rdkit_isolated: bool = True
