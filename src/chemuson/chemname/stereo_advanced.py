"""Stubs para estereoquímica avanzada (fase 6, best-effort pendiente)."""

from __future__ import annotations

from .molview import MolView


def extract_advanced_stereo(view: MolView) -> list[str]:
    """Hook para R_a/S_a, M/P y si/re cuando exista información confiable.

    TODO(PR28/Phase6):
    - Integrar asignación axial/helicoidal desde RDKit (modo aislado).
    - Requerir 3D o información suficiente para si/re en carbonilos.
    """
    _ = view
    return []
