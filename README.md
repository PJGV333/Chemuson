# ChemUSON

**Editor molecular libre para la docencia y la investigación química.**<br>
Desarrollado en el **Grupo de Química Supramolecular de la Universidad de Sonora, México**.

**Idioma:** Español · [Read this page in English](README.en.md)

ChemUSON es una aplicación de código abierto para **dibujar, visualizar, analizar y comunicar estructuras químicas**. Integra un editor 2D, visualización molecular 3D y herramientas de apoyo para química computacional, nomenclatura y exportación de figuras.

> **Versión pública más reciente: [ChemUSON 0.3.0-beta.3](https://github.com/PJGV333/Chemuson/releases/tag/v0.3.0-beta.3)** (25 de octubre de 2026). Es una **versión preliminar**: puede contener errores o funciones incompletas. La última versión etiquetada como estable es [v0.2.5](https://github.com/PJGV333/Chemuson/releases/tag/v0.2.5), anterior a esta beta.
>
> Para obtener **el código exacto de la beta publicada**, consulta [su etiqueta](https://github.com/PJGV333/Chemuson/tree/v0.3.0-beta.3). La rama `main` puede estar en una etapa diferente de desarrollo.

## Descargar ChemUSON

| Sistema | Modalidad | Descarga oficial |
| --- | --- | --- |
| Windows (64 bits) | Instalador | [Descargar `.exe`](https://github.com/PJGV333/Chemuson/releases/download/v0.3.0-beta.3/Chemuson-v0.3.0-beta.3-windows-x86_64-setup.exe) |
| Windows (64 bits) | Portable, sin instalación | [Descargar `.exe`](https://github.com/PJGV333/Chemuson/releases/download/v0.3.0-beta.3/Chemuson-v0.3.0-beta.3-windows-x86_64-portable.exe) |
| Linux (64 bits) | AppImage Type 2 | [Descargar `.AppImage`](https://github.com/PJGV333/Chemuson/releases/download/v0.3.0-beta.3/Chemuson-v0.3.0-beta.3-linux-x86_64.AppImage) |
| Linux (64 bits) | Flatpak | [Descargar `.flatpak`](https://github.com/PJGV333/Chemuson/releases/download/v0.3.0-beta.3/Chemuson-v0.3.0-beta.3-linux-x86_64.flatpak) |

Para Linux, también está disponible el [canal beta de Flatpak](https://pjgv333.github.io/Chemuson/flatpak/beta/Chemuson-beta.flatpakref), que permite recibir actualizaciones mediante el remoto publicado. No se requiere una cuenta de GitHub para descargar los archivos de una versión pública.

[**Ver todos los archivos de la versión, sus notas y sumas SHA-256**](https://github.com/PJGV333/Chemuson/releases/tag/v0.3.0-beta.3) · [**Consultar las sumas SHA-256**](https://github.com/PJGV333/Chemuson/releases/download/v0.3.0-beta.3/checksums.txt)

### Instalación rápida

**Windows:** ejecuta el instalador descargado o abre el ejecutable portable.

**Linux, AppImage:** concede permiso de ejecución y ábrelo desde el entorno gráfico o la terminal:

```bash
chmod +x Chemuson-v0.3.0-beta.3-linux-x86_64.AppImage
./Chemuson-v0.3.0-beta.3-linux-x86_64.AppImage
```

**Linux, Flatpak:** instala desde el remoto beta oficial:

```bash
flatpak install --user https://pjgv333.github.io/Chemuson/flatpak/beta/Chemuson-beta.flatpakref
flatpak run io.github.PJGV333.Chemuson
```

Consulta el [manual de usuario](docs/MANUAL_USUARIO.md) para instrucciones detalladas, uso de herramientas y solución de problemas.

## ¿Qué puedes hacer con ChemUSON?

- **Edición molecular 2D:** dibujar átomos, enlaces, anillos, cadenas, anotaciones, flechas y representaciones estereoquímicas.
- **Mejora de diagramas:** utilizar **Limpiar 2D** para reorganizar estructuras y preparar figuras para clases, informes y publicaciones.
- **Visualización molecular 3D:** generar conformeros, explorar geometría molecular y utilizar funciones de mecánica molecular basadas en RDKit.
- **Análisis químico:** consultar fórmula, masa molecular, masa exacta, análisis elemental y descriptores disponibles.
- **Nomenclatura y conversión:** utilizar herramientas de nomenclatura de alcance parcial (**IUPAC-lite**) y flujos de nombre a estructura.
- **Diagramación científica:** trabajar con orbitales, diagramas de energía y otras anotaciones de laboratorio.
- **Archivos y exportación:** guardar documentos editables `.cmsn`, intercambiar estructuras compatibles mediante SMILES/Molfile y exportar figuras a PNG, SVG o PDF.
- **Recuperación de trabajo:** disponer de guardado automático y mecanismos de recuperación de sesiones.

Las funciones avanzadas pueden depender del formato de entrada, el tipo de estructura, la plataforma o las dependencias disponibles.

## ChemUSON en imágenes

Las siguientes son **capturas de referencia de la interfaz**, conservadas en la documentación de modernización visual. Pueden diferir en algunos detalles de la versión beta publicada.

| Tema claro | Tema oscuro |
| --- | --- |
| ![Interfaz de ChemUSON en tema claro](docs/ui-modernization/after/main-light-1440x900.png) | ![Interfaz de ChemUSON en tema oscuro](docs/ui-modernization/after/main-dark-1440x900.png) |

[Ver más capturas de la interfaz](docs/ui-modernization/after/README.md).

## Primeros pasos

1. Abre ChemUSON y crea un documento con **Archivo → Nuevo**.
2. Selecciona herramientas de enlace, átomo o anillo en la barra lateral izquierda.
3. Dibuja una molécula en el lienzo; usa **Seleccionar** para manipular sus componentes.
4. Cuando lo necesites, aplica **Limpiar 2D** (`Ctrl+K`).
5. Abre el visor 3D para explorar la geometría de una estructura compatible.
6. Guarda el documento (`Ctrl+S`) en formato `.cmsn` o exporta una figura.

**Guías:** [Manual de usuario](docs/MANUAL_USUARIO.md) · [Nomenclatura ChemName](docs/chemname.md) · [Química computacional 3D](docs/geometry3d_compchem.md)

## Alcance y limitaciones

ChemUSON está en desarrollo activo. Su motor IUPAC-lite **no cubre toda la nomenclatura IUPAC**; algunas estructuras complejas, coordinaciones y asignaciones estereoquímicas pueden quedar fuera de su alcance. Las predicciones espectrales disponibles son orientativas y **no sustituyen mediciones experimentales**. Conviene verificar independientemente los resultados antes de utilizarlos en decisiones científicas relevantes.

Esta beta no constituye una validación regulatoria ni una garantía de ausencia de errores. Si encuentras un problema reproducible, puedes [abrir una incidencia](https://github.com/PJGV333/Chemuson/issues) con la versión, el sistema operativo y los pasos para reproducirlo, evitando compartir datos personales o confidenciales.

## Tecnologías utilizadas

- **Interfaz gráfica:** Python y PyQt6.
- **Herramientas químicas:** RDKit ejecutado en procesos aislados, Open Babel opcional y módulos internos de validación y geometría.
- **Distribución:** instalador y ejecutable portable para Windows; AppImage Type 2 y Flatpak para Linux.

Los detalles de arquitectura, diagnóstico de calidad de Clean2D, contratos de pruebas y preparación de paquetes están en la [documentación técnica](docs/architecture.md), la [guía de desarrollo](docs/distribution-plan.md) y las [instrucciones de distribución](docs/release/PREVIEW_BUILDS.md).

### Para desarrolladores

Con Python 3.10 o posterior, en un entorno virtual:

```bash
python -m venv .venv
source .venv/bin/activate
python -m pip install -e .
python -m pip install -r requirements-dev.txt
```

En Windows activa `.venv\Scripts\activate`; en fish usa `source .venv/bin/activate.fish`. El código fuente que corresponde a un paquete concreto se identifica por la etiqueta de su versión.

## Origen del proyecto y autores

ChemUSON nace del trabajo del **Grupo de Química Supramolecular de la Universidad de Sonora (México)**, con el propósito de contribuir a la enseñanza, investigación y comunicación de la química mediante herramientas libres. Su desarrollo busca hacer más accesible la representación y el análisis molecular, sin sustituir el juicio científico.

**Autores:** Pedro Jancarlo Gomez Vega*, José Octavio Juárez Sánchez, Alondra Karolina Montaño Preciado, Ramón Moreno Corral, Felipe Medrano Valenzuela, David Octavio Corona Martínez y Karen L. Ochoa Lara*.

\* Autores de correspondencia.

**Licencia:** [AGPL-3.0](LICENSE). Se agradecen informes de errores, propuestas y contribuciones a través de [Issues](https://github.com/PJGV333/Chemuson/issues) y [Pull Requests](https://github.com/PJGV333/Chemuson/pulls).

---

**ChemUSON — Universidad de Sonora, México.** · [English version](README.en.md)
