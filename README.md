# Chemuson

Chemuson es un editor molecular 2D libre y open source, desarrollado como **proyecto en marcha del Grupo de Química Supramolecular de la Universidad de Sonora**.

## Estado del proyecto

> **Proyecto en construcción.**
>
> Chemuson está en desarrollo activo. Las funciones, formatos y resultados pueden cambiar mientras el proyecto madura.

## ¿Qué hace Chemuson?

Chemuson permite dibujar estructuras químicas 2D, analizarlas y exportarlas en formatos comunes para trabajo académico y de investigación.

Capacidades actuales (resumen):
- Editor gráfico 2D con herramientas para átomos, enlaces, anillos, cadenas y anotaciones.
- Soporte de estilos de enlace (simple, doble, triple, aromático, cuña, punteado, coordinativo, etc.).
- Limpieza geométrica (`Limpiar 2D`) y opciones visuales (hidrógenos, carbonos, aromáticos como círculo, rejilla/reglas).
- Numeración visual no destructiva (átomos, estructuras o ambos), configurable para exportación.
- Importación y exportación química: `SMILES`, `Molfile` y archivos propios `.cmsn`.
- Autosave automático con backups rotativos por documento y recuperación de sesiones al iniciar.
- Exportación gráfica: `PNG`, `SVG`, `PDF`.
- Herramientas de análisis integradas: nombre, fórmula, masa exacta, peso molecular, `m/z` y análisis elemental.
- Motor de nomenclatura **IUPAC-lite** para un subconjunto de estructuras orgánicas.
- Biblioteca de plantillas (incluye base predefinida + categorías de usuario importables/exportables).

## Autosave y recuperación

- Chemuson guarda autosaves automáticamente cada ~2 minutos y también después de unos segundos de inactividad cuando hay cambios.
- Al iniciar, si existen autosaves pendientes, aparece un diálogo para **recuperar** o **descartar** cada sesión.
- Al guardar manualmente, se limpian autosaves obsoletos del documento y se conserva un respaldo rotativo reciente.

Rutas locales de trabajo:
- Autosaves: `~/.chemuson/autosave/`
- Autosaves archivados tras recuperar/descartar: `~/.chemuson/autosave/old/`
- Logs de crash: `~/.chemuson/crash_logs/crash_YYYYMMDD_HHMMSS.txt`

## Fortalezas

- Arquitectura modular (núcleo químico, GUI, nomenclatura, IO y persistencia separados).
- Integración con RDKit y rutas de respaldo internas cuando RDKit no cubre ciertos casos.
- Cobertura de pruebas amplia para una base en evolución (49 pruebas en `tests/`).
- Orientado a uso docente y de laboratorio con enfoque práctico de dibujo y análisis rápido.

## Debilidades y límites actuales

- El motor de nombres es **IUPAC-lite**: no cubre toda la nomenclatura IUPAC ni todos los sistemas cíclicos/fusionados.
- Algunas rutas de cálculo/exportación usan subconjuntos de elementos o aproximaciones, según el caso.
- Al ser un proyecto en construcción, puede haber cambios de comportamiento entre versiones.
- No sustituye herramientas de validación química regulatoria o flujos de producción validados.

## Alcance de IUPAC-lite (actual)

Incluye soporte para:
- Cadenas lineales y algunas insaturaciones.
- Sustituyentes comunes (halógenos, alquilos simples, algunos grupos funcionales).
- Cicloalcanos simples, benceno y un subconjunto de heteroaromáticos.
- Algunos sistemas aromáticos fusionados (por ejemplo: naftaleno, antraceno, fenantreno, pireno; y casos concretos heterofusionados).

Limitaciones relevantes:
- Cobertura parcial en anillos fusionados complejos.
- Cobertura parcial en heterociclos fuera del subconjunto soportado.
- Casos con múltiples grupos funcionales y escenarios avanzados aún no completamente soportados.

## Stack técnico

- **GUI**: PyQt6
- **Química**: RDKit
- **Lenguaje**: Python 3.10+

## Instalación

```bash
python3 -m venv chemuson
./chemuson/bin/pip install -r requirements.txt
```

Opcional (modo editable):

```bash
./chemuson/bin/pip install -e .
```

## Ejecución

Desde el entorno virtual:

```bash
./chemuson/bin/python src/main.py
```

Si instalaste en modo editable:

```bash
./chemuson/bin/chemuson
```
