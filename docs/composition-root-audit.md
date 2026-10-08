# Auditoría del composition root

Fecha: 2026-09-19

## Decisión

M19 (`bootstrap`) está **consolidated / no structural change required**. `__main__.py`
contiene el parser, un despacho interno congelado para el worker aislado de RDKit
(M01, antes del bootstrap GUI) y el defer de `run_app`; `app/bootstrap.py` crea
QApplication, instala resiliencia, configura la ventana, muestra la aplicación y
entrega el control al event loop. El despacho RDKit no importa Qt ni RDKit en el
proceso padre y no crea una segunda raíz de composición.

## Ownership

M19 conserva `src/chemuson/__main__.py` y `src/chemuson/app/`. La composición
conoce los adaptadores GUI, M22 y el worker aislado de M01 por diseño; no se mueve
lógica de dominio ni se crea una segunda raíz de composición. M19 registra M01
como dependencia de packaging/worker en `architecture/modules.yml`.

## M24

M24 queda reservado. No se crea un módulo vacío ni se fragmenta un composition
root ya pequeño y explícito.
