# Auditoría del composition root

Fecha: 2026-09-19

## Decisión

M19 (`bootstrap`) está **consolidated / no structural change required**. `__main__.py`
se limita al parser y al defer de `run_app`; `app/bootstrap.py` crea QApplication,
instala resiliencia, configura la ventana, muestra la aplicación y entrega el
control al event loop.

## Ownership

M19 conserva `src/chemuson/__main__.py` y `src/chemuson/app/`. La composición
conoce los adaptadores GUI y M22 por diseño; no se mueve lógica de dominio ni se
crea una segunda raíz de composición.

## M24

M24 queda reservado. No se crea un módulo vacío ni se fragmenta un composition
root ya pequeño y explícito.
