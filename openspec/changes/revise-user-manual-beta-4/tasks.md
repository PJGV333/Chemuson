# Tareas — revisión del manual de usuario contra 0.3.0-beta.4

- [x] 1.1 Registrar `git status`, rama, HEAD y la equivalencia `v0.3.0-beta.4..main` en `baseline.md`.
- [x] 1.2 Volcar menús, acciones con atajo, paletas del rail y menú contextual de la aplicación ejecutada sin pantalla.
- [x] 1.3 Levantar las discrepancias A1–A9, duplicaciones B1–B6, hallazgos C1–C10 y omisiones D (no corregir código).
- [x] 2.1 Generar capturas auténticas de `0.3.0-beta.4` (ventana, menús, paletas, lienzo, contextual, paneles) y versionar el script `docs/manual-figures/capturar_figuras.py`.
- [x] 2.2 Redactar `docs/manual-figures/README.md` con procedencia, entorno, regeneración e índice de archivos.
- [x] 3.1 Reescribir el manual con 22 capítulos, índice y anclas internas actualizadas.
- [x] 3.2 Corregir §7/Estructura: retirar la matriz de uniones, enlazar la figura del menú real y desarrollar limpieza 2D.
- [x] 3.3 Presentar `Ctrl+K` como método recomendado y documentar `Limpiar 2D` como entrada heredada equivalente (C1).
- [x] 3.4 Añadir el enlace público al manual y declarar que **no** existe «Ayuda → Manual de usuario» en la interfaz (propuesta futura).
- [x] 3.5 Eliminar los falsos accesos directos de elementos y documentar la paleta de átomos y la tabla periódica.
- [x] 3.6 Separar giro 2D, ramas, fragmento con pivote, perspectiva *trackball* y modelo 3D (capítulos 6 y 15).
- [x] 3.7 Ilustrar y corregir flechas, corchetes, símbolos, placas, diagramas de energía y orbitales.
- [x] 3.8 Sustituir las matrices comparativas por operaciones reales de enlace y fragmento.
- [x] 3.9 Absorber las secciones duplicadas de paneles, barras flotantes, TLC/gel y plantillas.
- [x] 3.10 Reescribir el menú contextual por tipo de objeto, con capturas de cada estado.
- [x] 3.11 Unificar plantillas en un solo capítulo con su biblioteca de fábrica real.
- [x] 3.12 Auditoría final: atajos inexistentes, menús inexistentes, afirmaciones de química sin soporte, enlaces y referencias cruzadas.
- [x] 4.1 Validar rutas de imagen, anclas, correlación imagen↔leyenda y consecutividad de figuras.
- [x] 4.2 `git diff --check` y comprobación de que los README siguen enlazando al manual principal.
- [x] 4.3 `openspec validate revise-user-manual-beta-4 --strict`.
- [x] 4.4 Pasada cruzada de las 55 rutas de menú citadas contra los volcados de la interfaz.
- [x] 4.5 Barrido inverso de los 229 elementos de la interfaz para detectar funciones sin documentar; las cuatro lagunas halladas se documentan.
- [x] 4.6 Regeneración y verificación del inventario de figuras (37/37, cero desajustes).
## 5. Entrega

- [x] 5.1 Publicar la rama `docs/manual-usuario-beta4-campana` y abrir la PR hacia `main` (**#13**), sin combinarla.
- [x] 5.2 Entregar el informe de campaña en español con los once apartados solicitados, incluidos los inventarios de ilustraciones y de funciones no verificables.
- [ ] 5.3 Revisión editorial del propietario en la PR: el cambio no se combina sin aprobación expresa.
- [ ] 5.4 Campaña posterior, fuera de este alcance: aprobación del Markdown y generación del PDF definitivo.
