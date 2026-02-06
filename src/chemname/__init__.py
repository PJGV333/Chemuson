"""Motor de nomenclatura IUPAC-lite para Chemuson.

Soporta (v0) un subconjunto de hidrocarburos y aromáticos sencillos:
- Cadenas acíclicas (alcanos), halógenos y sustituyentes alquilo lineales.
- Una sola insaturación C=C o C#C en la cadena principal.
- Un único grupo funcional simple (-ol o -amine).
- Cicloalcanos simples (C3–C8) con sustituyentes sencillos.
- Benceno con sustituyentes simples (halógenos, alquilo, hidroxi, nitro).
- Heteroaromáticos: furano, tiofeno, pirrol, piridina, diazinas, imidazol,
  oxazol, tiazol.
- Naftaleno con derivados mono/di-sustituidos (sustituyentes simples).
- Sustituyentes alquilo ramificados: isopropilo, sec-butilo, tert-butilo,
  neopentilo.
- Sustituyentes de anillo: fenilo, bencilo, ciclohexilo.
- Aromáticos fusionados: antraceno y fenantreno (mono-sustituidos).
- Aromáticos fusionados: pireno (mono/di-sustituidos con numeración por plantilla).
- Heteroaromáticos fusionados: quinolina, isoquinolina, indol (mono-sustituidos).

Fuera de alcance se devuelve "N/D" si está habilitada la opción.
"""

from .engine import iupac_name
from .options import NameOptions

__all__ = ["iupac_name", "NameOptions"]
