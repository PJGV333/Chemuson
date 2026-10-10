# Matriz de aceptación manual — ChemUSON 0.3.0-beta.1

> **P1 — Missing UI icons in packaged Windows and Linux builds — FAILED — Blocks beta acceptance.** El propietario observó iconos esenciales ausentes en los paquetes portable Windows y Linux del Actions run [#37826597134](https://github.com/PJGV333/Chemuson/actions/runs/37826597134), en tema claro y oscuro. Los casos históricos `UI-FAIL-WIN-01` y `UI-FAIL-LINUX-01` quedan `FAILED — P1 blocks beta acceptance`. Los nuevos smokes de iconos Windows/AppImage pasan en [Build Preview #38006731375](https://github.com/PJGV333/Chemuson/actions/runs/38006731375), SHA `57d481d21fee2db6ecacf6a862839caa63a1c56d`; no reemplazan la verificación visual manual, que sigue pendiente.
>
> **P1 — RDKit isolated backend unavailable in packaged executable — FAILED — Blocks beta acceptance.** El propietario informó que Windows portable calcula fórmula, masa y espectros estimados, pero no muestra descriptores RDKit. Los nuevos ejecutables congelados Windows/AppImage pasan worker, extensiones nativas, descriptores, SMILES y 3D en [Build Preview #38006731375](https://github.com/PJGV333/Chemuson/actions/runs/38006731375), SHA `57d481d21fee2db6ecacf6a862839caa63a1c56d`; el retest manual de Windows/Linux sigue pendiente.
>
> **P1 — ChemName MOL templates omitted from frozen packages — FAILED — Blocks beta acceptance.** Los CArchive reales de Windows portable y AppImage del run `37855524295` carecían de los nueve `.mol` necesarios; no se observó un defecto de reglas de nomenclatura. Build Preview [#38006731375](https://github.com/PJGV333/Chemuson/actions/runs/38006731375), SHA `57d481d21fee2db6ecacf6a862839caa63a1c56d`, ahora valida 9/9 templates y 6/6 controles ChemName en Windows portable/AppImage/Flatpak, comparados con Python. Los casos históricos permanecen `FAILED`; el retest manual del propietario sigue pendiente.

> Las filas históricas ya iniciadas conservan su estado. Los casos nuevos empiezan en `NOT TESTED`; la matriz no acredita aceptación hasta completar pruebas con paquetes instalados y SHA identificables.

## Registro de ejecución

- Artefacto: `Windows portable / Windows setup / Linux AppImage Type 2 / Linux Flatpak`
- Versión y SHA completo: `NOT RECORDED`
- Fuente: `Actions run URL + artifact name`
- Sistema/versión/arquitectura: `NOT RECORDED`
- Tester/fecha: `NOT RECORDED`
- Configuración (escala, tema, red, proveedor IA): `NOT RECORDED`
- Evidencia/directorio de logs: `NOT RECORDED`

Resultado permitido por caso: `PASS`, `FAIL`, `BLOCKED` o `NOT TESTED`. Registrar pasos observados, no sólo la conclusión. Si el mismo caso se prueba en varias plataformas, crear una copia de la tabla por ejecución.

## A. Inicio, salida y ciclo de vida

| ID | Pasos y datos | Resultado esperado | Severidad | Resultado / evidencia |
|---|---|---|---|---|
| START-01 | Instalar/extraer el paquete en un perfil limpio y abrir ChemUSON desde su entrypoint. | Ventana principal visible, versión `0.3.0-beta.1`, sin diálogo fatal. | P0 | NOT TESTED — |
| START-02 | Iniciar desde ruta con espacios y caracteres no ASCII. | Arranque y acceso a recursos sin error de ruta. | P2 | NOT TESTED — |
| START-03 | Abrir dos documentos y cambiar repetidamente entre pestañas. | El documento activo y sus propiedades corresponden a la pestaña elegida. | P1 | NOT TESTED — |
| START-04 | Abrir/cerrar el panel de análisis y el de Assistant tres veces. | No hay wrapper/diálogo obsoleto, bloqueo ni crash. | P0 | NOT TESTED — |
| START-05 | Cerrar la ventana principal con cálculo/descriptores aún activos. | Cierre controlado; no hay QThread activo destruido ni proceso colgado. | P0 | NOT TESTED — |
| START-06 | Salir desde menú y repetir mediante el botón cerrar de ventana. | Ambos caminos liberan la app sin crash ni pérdida inesperada. | P0 | NOT TESTED — |
| START-07 | Abrir dos instancias según el modo de instalación disponible. | La segunda instancia falla de forma controlada o se comporta según el contrato del paquete. | P2 | NOT TESTED — |
| START-08 | Ejecutar `chemuson --version` si la instalación ofrece CLI. | La versión coincide con el manifest y el binario instalado. | P2 | NOT TESTED — |

## B. Edición, dibujo y acciones

| ID | Pasos y datos | Resultado esperado | Severidad | Resultado / evidencia |
|---|---|---|---|---|
| DRAW-01 | Dibujar dos carbonos y un enlace simple. | Se crea la identidad/conectividad exacta de dos átomos y un enlace, sin enlaces extra o perdidos. | P0 | NOT TESTED — |
| DRAW-02 | Cambiar el enlace a doble y triple. | Orden químico y representación corresponden a cada elección, sin cambio silencioso de identidad. | P0 | NOT TESTED — |
| DRAW-03 | Dibujar enlace aromático y revisar las etiquetas/valencias. | Aromaticidad y grafo son correctos, sin enlaces duplicados accidentales. | P0 | NOT TESTED — |
| DRAW-04 | Crear un anillo de seis miembros y un anillo fusionado. | Los anillos quedan conectados con el número esperado de átomos/enlaces y son editables. | P0 | NOT TESTED — |
| DRAW-05 | Añadir O, N, S y halógeno desde selector de elemento; comprobar valencias ordinarias. | Elemento, valencia y etiqueta persisten tras cambiar herramienta; no se altera el grafo en silencio. | P0 | NOT TESTED — |
| DRAW-06 | Seleccionar átomo, enlace y región; mover cada selección. | Sólo el objeto seleccionado cambia de posición. | P1 | NOT TESTED — |
| DRAW-07 | Copiar/pegar una estructura pequeña. | La copia conserva topología y propiedades químicas relevantes. | P1 | NOT TESTED — |
| DRAW-08 | Ejecutar Undo y Redo después de dibujar y mover una molécula. | Cada acción restaura exactamente el grafo, propiedades y estado anterior. | P0 | NOT TESTED — |
| DRAW-09 | Eliminar un átomo terminal y después restaurarlo con Undo. | El grafo y enlace se restauran exactamente, sin pérdida o elementos fantasma. | P0 | NOT TESTED — |
| DRAW-10 | Cambiar zoom, pan y rejilla; volver al zoom inicial. | La vista responde y no cambia la química del documento. | P2 | NOT TESTED — |
| DRAW-11 | Dibujar cuña sólida y discontinua en centro estereogénico. | La notación indicada se muestra y conserva en operaciones soportadas. | P1 | NOT TESTED — |
| DRAW-12 | Crear texto/anotación, seleccionarlo y borrarlo. | Sólo cambia la anotación; la estructura permanece. | P2 | NOT TESTED — |

## C. Persistencia, importación y estructura química

| ID | Pasos y datos | Resultado esperado | Severidad | Resultado / evidencia |
|---|---|---|---|---|
| DATA-01 | Dibujar etanol `CCO`, guardar como `.cmsn`, cerrar y reabrir. | Identidad y SMILES canónico isomérico `CCO` iguales antes/después; se conservan fórmula, 3 átomos, 2 enlaces, etiquetas, carga/isótopos y coordenadas relevantes. | P0 | NOT TESTED — |
| DATA-02 | Guardar documento modificado con un nombre nuevo usando “Guardar como”. | El original no se sobrescribe y el nuevo documento abre completo. | P1 | NOT TESTED — |
| DATA-03 | Guardar, editar y usar Undo hasta el estado guardado. | El indicador dirty/asterisco refleja el estado limpio. | P2 | NOT TESTED — |
| DATA-04 | Abrir `.cmsn` válido de una versión anterior disponible. | Se carga sin perder estructuras soportadas; avisos son claros. | P1 | NOT TESTED — |
| DATA-05 | Intentar abrir un `.cmsn` truncado/corrupto de copia de prueba. | Error controlado; no modifica ni elimina el archivo original. | P1 | NOT TESTED — |
| DATA-06 | Importar SMILES `CCO`. | Identidad de etanol (`CCO`), fórmula, conectividad y elementos coinciden con la referencia; una cadena válida por sí sola no basta. | P0 | NOT TESTED — |
| DATA-07 | Importar SMILES `c1ccccc1`. | Identidad de benceno, seis carbonos, fórmula/conectividad y aromaticidad esperadas. | P0 | NOT TESTED — |
| DATA-08 | Importar un Molfile del catión metilamonio isotópico `[13CH3][NH3+]` (dos elementos, carga +1 e isótopo 13C). | Identidad/fórmula `CH6N+`, conectividad C–N, carga e isótopo 13C coinciden con la referencia independiente; no se altera silenciosamente. | P0 | NOT TESTED — |
| DATA-09 | Importar SMILES inválido `C1CC`. | Rechazo controlado con explicación; documento activo intacto. | P1 | NOT TESTED — |
| DATA-10 | Exportar `CCO` a Molfile y reimportar en documento nuevo; repetir con casos de `DATA-11` en las rutas SMILES/MOL/SDF soportadas. | Conservar identidad molecular: grafo, fórmula, enlaces, carga/isótopos y estereoquímica representable iguales antes/después, comprobados con referencia independiente; no basta que el texto se pueda volver a parsear. | P0 | NOT TESTED — |
| DATA-11 | Importar los enantiómeros `C[C@H](O)F` y `C[C@@H](O)F`, el aminoácido `N[C@@H](C)C(=O)O` y los alquenos `F/C=C/F` y `F/C=C\F`; hacer round-trip SMILES↔MOL/SDF en las rutas soportadas. | Se conserva la identidad tetraédrica/E/Z de cada entrada y los enantiómeros siguen siendo distintos. Comparar con referencia estereoquímica independiente; no comparar sólo el texto `@`/`@@`. Si un formato no puede representarla, debe informar/fallar explícitamente, nunca borrarla en silencio. | P0 | NOT TESTED — |
| DATA-12 | Guardar estructura con carga formal y aromaticidad; volver a abrir. | Propiedades químicas no cambian silenciosamente. | P0 | NOT TESTED — |

## D. Exportación gráfica y archivos

| ID | Pasos y datos | Resultado esperado | Severidad | Resultado / evidencia |
|---|---|---|---|---|
| EXP-01 | Exportar etanol a PNG. | Archivo existe, no vacío y muestra la estructura. | P2 | NOT TESTED — |
| EXP-02 | Exportar anillo aromático a SVG. | XML abre y el gráfico contiene estructura completa. | P2 | NOT TESTED — |
| EXP-03 | Exportar a PDF y abrir en visor externo. | Página legible y estructura completa. | P2 | NOT TESTED — |
| EXP-04 | Exportar varias páginas/estructuras si la UI lo permite. | Orden y contenido de página corresponden al documento. | P2 | NOT TESTED — |
| EXP-05 | Exportar con numeración visual activada y desactivada. | Numeración afecta la salida gráfica, no el grafo químico. | P2 | NOT TESTED — |
| EXP-06 | Exportar a un directorio sin permiso o con destino inválido. | Error informado; documento fuente permanece intacto. | P1 | NOT TESTED — |

## E. ChemName, Clean2D y propiedades

| ID | Pasos y datos | Resultado esperado | Severidad | Resultado / evidencia |
|---|---|---|---|---|
| CHEMNAME-WIN-OLD-01 | Inspección del portable Windows del run `37855524295`, SHA-256 `f2648e29f3557c301fe8b2a094bd4ca5cd1bbb0b20c69090419d2e51fd6897cd`: revisar CArchive PyInstaller. | Los nueve `.mol` de ChemName deben estar en el paquete. La ausencia confirma el P1 del portable anterior, sin atribuirlo a reglas químicas. | P1 | **FAILED — P1 blocks beta acceptance** — CArchive TOC contiene 0 `.mol` de `chemuson/chemname/templates`. |
| CHEMNAME-LINUX-OLD-01 | Inspección del ejecutable extraído del AppImage preview `37855524295`, SHA-256 AppImage `e2d76a6af54367e8ce1dddf1b5b8d3c8a2ace0d1bc00e71cce9262bf0e1e2e43`. | Los nueve `.mol` deben existir dentro del CArchive y resolverse por el ejecutable. | P1 | **FAILED — P1 blocks beta acceptance** — CArchive TOC contiene 0 templates `.mol`. |
| CHEMNAME-FLATPAK-OLD-01 | Bundle preview `37855524295`, SHA-256 `e319cd2cb73d7eea7e979cb6447ad99b6c2e659e4789340e4c832503ea208700`; importarlo sólo en un OSTree temporal y enumerar package data. | Los nueve `.mol` aparecen bajo `site-packages/chemuson/chemname/templates`; el smoke de nombres requiere nueva build. | P1 | **NOT TESTED —** 9/9 recursos están en el bundle; no se ejecutó el smoke de ChemName de la app instalada. |
| CHEMNAME-WIN-AUTO-01 | Preview [#38006731375](https://github.com/PJGV333/Chemuson/actions/runs/38006731375), Windows portable [artifact 11651882043](https://github.com/PJGV333/Chemuson/actions/runs/38006731375/artifacts/11651882043), `Chemuson-v0.3.0-beta.1-preview-57d481d2-windows-x86_64-portable.exe`, SHA-256 `dde7c8f94fef2fd344756c8aa25eac53bd998c2ff8a007021f4e05f014e00bc6`. `checksums.sha256` y provenance verificados; provenance SHA fuente `57d481d21fee2db6ecacf6a862839caa63a1c56d`, `publication=false`. | Ejecutable congelado resuelve/nombra 9/9 templates desde `_MEIPASS` y 6/6 moléculas igual a Python; no importa GUI. | P1 | **PASS — build smoke únicamente**; inspección manual no hecha. |
| CHEMNAME-WIN-SETUP-AUTO-01 | Preview [#38006731375](https://github.com/PJGV333/Chemuson/actions/runs/38006731375), Windows setup [artifact 11651667473](https://github.com/PJGV333/Chemuson/actions/runs/38006731375/artifacts/11651667473), `Chemuson-v0.3.0-beta.1-preview-57d481d2-windows-x86_64-setup.exe`, SHA-256 `95604f16ab7600d41e997c60d4ed485cf86782246a15faec8720dd6383d8eeed`. `checksums.sha256` y provenance verificados; source SHA `57d481d21fee2db6ecacf6a862839caa63a1c56d`, `publication=false`. | Inno genera el instalador desde el executable validado y checksums/provenance coinciden; no se instala ni prueba upgrade/uninstall aquí. | P1 | **PASS — artifact/build only**; instalación y retest manual siguen NOT TESTED. |
| CHEMNAME-APPIMAGE-AUTO-01 | Preview [#38006731375](https://github.com/PJGV333/Chemuson/actions/runs/38006731375), AppImage Type 2 [artifact 11650889195](https://github.com/PJGV333/Chemuson/actions/runs/38006731375/artifacts/11650889195), `Chemuson-v0.3.0-beta.1-preview-57d481d2-linux-x86_64.AppImage`, SHA-256 `ec26ce53bfd6202d72e310e160c0bbc5a32c091deafcf6cc5b6fdeb066826e9c`. Header Type 2 `AI\x02`; provenance source SHA `57d481d21fee2db6ecacf6a862839caa63a1c56d`, `publication=false`; checksum/provenance verificados. | Ejecutable extraído verifica 9/9 templates y 6/6 nombres frente a Python; valida también RDKit y un headless launch acotado. | P1 | **PASS — build smoke únicamente**; aceptación GUI manual no hecha. |
| CHEMNAME-FLATPAK-AUTO-01 | Preview [#38006731375](https://github.com/PJGV333/Chemuson/actions/runs/38006731375), Flatpak [artifact 11651388730](https://github.com/PJGV333/Chemuson/actions/runs/38006731375/artifacts/11651388730), `Chemuson-v0.3.0-beta.1-preview-57d481d2-linux-x86_64.flatpak`, SHA-256 `f4562e2a5cb6676042043eb8c0efc3d4f6006f7013ca58067da2924ddfeeea84`. `checksums.sha256` y provenance verificados; source SHA `57d481d21fee2db6ecacf6a862839caa63a1c56d`, `publication=false`. | Smoke dentro del build bajo `/app/lib/python3.13/site-packages`: 9/9 templates y 6/6 nombres esperados. | P1 | **PASS — build smoke únicamente**; aceptación GUI manual no hecha. |
| CHEMNAME-RETEST-01 | En los cuatro paquetes del preview [#38006731375](https://github.com/PJGV333/Chemuson/actions/runs/38006731375), SHA fuente `57d481d21fee2db6ecacf6a862839caa63a1c56d`, crear/abrir etanol `CCO`, benceno, acetamida, etano y ciclohexano; revisar nombre en barra de estado y anotación; probar elemento no soportado. Artifact IDs y checksums están en los cuatro casos `*-AUTO-01` anteriores. | Los nombres coinciden con Python: `ethan-1-ol`, `benzene`, `1-aminoethanamide`, `ethane`, `cyclohexane`; el caso no soportado muestra `N/D` sin excepción ni cambio químico. | P1 | NOT TESTED — owner retest required; automation is not manual acceptance. |
| CHEMNAME-UPDATE-01 | Sobre cada paquete exacto del Build Preview [#38006731375](https://github.com/PJGV333/Chemuson/actions/runs/38006731375), source SHA `57d481d21fee2db6ecacf6a862839caa63a1c56d`, dibujar etano, añadir O terminal para producir etanol, consultar barra/anotación y hacer Undo/Redo. Usar los artifact IDs/checksums de `*-AUTO-01` y registrar el SHA descargado. | El nombre sigue la estructura inmediatamente: `ethane` → `ethan-1-ol` → `ethane` → `ethan-1-ol`; el texto de una anotación nueva corresponde al grafo actual y ninguna actualización tardía restaura el nombre anterior. Elemento no soportado sigue `N/D`. | P1 | NOT TESTED — owner manual retest required; package smoke passes do not replace acceptance. |
| CHEMNAME-ROBUSTNESS-01 | Campaña diferida: reproducir moléculas de nombres discutidos con SMILES/mol y opciones exactas; comparar contra nombre preferido IUPAC actual y referencia normativa/versionada. Investigar por separado de la actualización UI. | La fuente, estructura, nombre esperado, salida ChemName y discrepancia están registrados antes de cambiar reglas; salida no verificada no se presenta como aceptación normativa. | P2 — review futura; no defecto confirmado | BLOCKED — no se suministró un caso propietario ni se completó la validación normativa. Una sonda exploratoria de `C=C` produjo `eth-1-ene`; requiere adjudicación antes de clasificar como defecto. |
| CLEAN-01 | Ejecutar Clean2D en etanol. | Geometría mejora y número/conectividad de átomos no cambia. | P1 | NOT TESTED — |
| CLEAN-02 | Ejecutar Clean2D en benceno. | Anillo queda legible; aromaticidad y conectividad permanecen. | P1 | NOT TESTED — |
| CLEAN-03 | Ejecutar Clean2D en `N[C@@H](C)C(=O)O`. | No se altera centro estéreo/topología sin aviso. | P0 | NOT TESTED — |
| CLEAN-04 | Aplicar Clean2D y después Undo/Redo. | La operación se revierte/restaura como una acción coherente. | P1 | NOT TESTED — |
| CLEAN-05 | Ejecutar Clean2D en estructura con sustituyentes cercanos. | No desaparecen enlaces ni átomos; resultado revisable. | P1 | NOT TESTED — |
| CLEAN-06 | Abrir reporte/calidad tras una estructura válida y otra problemática. | Diagnóstico corresponde a la molécula y no altera su geometría. | P2 | NOT TESTED — |
| RDKIT-WIN-01 | **Histórico — etanol `CCO` en Windows portable previo al fix:** abrir Propiedades y consultar Descriptores RDKit. | logP ≈ -0.0014, TPSA 20.23, HBD 1 y HBA 1; fórmula/masa/espectros estimados siguen disponibles; worker usa el runtime empaquetado. | P1 | FAILED — owner report: “RDKit no disponible; resultado parcial”; retest bloqueado hasta preview corregido |
| RDKIT-WIN-02 | **Preview exacto para retest:** [Windows portable, artifact 11651882043](https://github.com/PJGV333/Chemuson/actions/runs/38006731375/artifacts/11651882043), SHA-256 `dde7c8f94fef2fd344756c8aa25eac53bd998c2ff8a007021f4e05f014e00bc6`, provenance source SHA `57d481d21fee2db6ecacf6a862839caa63a1c56d`, `publication=false`. En etanol `CCO`, abrir Propiedades y consultar Descriptores RDKit. | Los descriptores conocidos aparecen; no requiere Python/RDKit del sistema y la aplicación conserva los resultados parciales si el worker falla. | P1 | NOT TESTED — frozen worker/build smoke passed; owner GUI retest pending. |
| RDKIT-LINUX-01 | **Preview exacto para retest:** [AppImage, artifact 11650889195](https://github.com/PJGV333/Chemuson/actions/runs/38006731375/artifacts/11650889195), SHA-256 `ec26ce53bfd6202d72e310e160c0bbc5a32c091deafcf6cc5b6fdeb066826e9c`, provenance source SHA `57d481d21fee2db6ecacf6a862839caa63a1c56d`, `publication=false`. En etanol `CCO`, abrir Propiedades y consultar Descriptores RDKit. | Los mismos descriptores que Windows; no requiere Python/RDKit del sistema. | P1 | NOT TESTED — frozen worker/build smoke passed; owner GUI retest pending. |
| RDKIT-3D-SMILES-01 | Ejecutar worker/función 3D con etanol en los previews exactos de Windows portable y AppImage Type 2, artifacts 11651882043 y 11650889195 del run 38006731375 (SHA de fuente `57d481d21fee2db6ecacf6a862839caa63a1c56d`). | Worker aislado produce SMILES canónico `CCO` y coordenadas 3D finitas; proceso padre no carga RDKit; timeout y errores no bloquean la GUI. | P1 | NOT TESTED — build smoke passed; owner manual retest pending. |

## F. Molecular Assistant (experimental)

| ID | Pasos y datos | Resultado esperado | Severidad | Resultado / evidencia |
|---|---|---|---|---|
| AI-01 | Abrir Assistant sin configurar proveedor. | Estado offline/configuración explicado; no hace lookup implícito. | P1 | NOT TESTED — |
| AI-02 | Consultar con proveedor autorizado usando solicitud de molécula pequeña. | Se indica fuente/proveedor; resultado revisable antes de aplicar. | P1 | NOT TESTED — |
| AI-03 | Introducir clave temporal y cerrar Assistant. | La clave no aparece en preferencias persistentes ni logs. | P0 | NOT TESTED — |
| AI-04 | Solicitar estructura con proveedor desconectado. | Error controlado, sin bloqueo de UI ni afirmación de éxito. | P1 | NOT TESTED — |
| AI-05 | Solicitar respuesta con SMILES inválido, si el harness permite fixture. | La propuesta inválida se rechaza antes de modificar documento. | P1 | NOT TESTED — |
| AI-06 | Revisar propuesta válida pero distinta de la identidad solicitada. | Validez sintáctica y verificación de identidad se muestran como estados separados. | P0 | NOT TESTED — |
| AI-07 | Previsualizar y cancelar una propuesta completa. | Cancelar no modifica la molécula activa ni UndoStack. | P1 | NOT TESTED — |
| AI-08 | Aplicar propuesta a copia de prueba, luego Undo/Redo. | Sustitución es exactamente reversible; origen y preview son identificables. | P1 | NOT TESTED — |
| AI-09 | Cerrar ventana mientras petición asíncrona está activa. | Cierre controlado; sin crash, diálogo huérfano ni tarea Qt liberada prematuramente. | P0 | NOT TESTED — |
| AI-10 | Repetir petición con proveedor que devuelve timeout/límite agotado. | Se informa el fallo y no se fabrica una estructura de reemplazo. | P1 | NOT TESTED — |

## G. Interfaz, preferencias y accesibilidad

| ID | Pasos y datos | Resultado esperado | Severidad | Resultado / evidencia |
|---|---|---|---|---|
| UI-01 | Cambiar tema claro/oscuro y reiniciar. | Preferencia persiste y controles siguen legibles. | P2 | NOT TESTED — |
| UI-02 | Usar pantalla 1280×720 y luego una ventana estrecha. | Barra, paneles y canvas siguen accesibles, sin controles críticos fuera de pantalla. | P2 | NOT TESTED — |
| UI-03 | Usar escala del sistema 150–200%. | Texto/íconos no se solapan críticamente; canvas conserva operación. | P2 | NOT TESTED — |
| UI-04 | Navegar menú/herramientas con teclado y revisar atajos visibles. | Acciones disponibles por teclado y foco perceptible donde corresponda. | P2 | NOT TESTED — |
| UI-05 | Cambiar intervalos/canal del updater sin guardar credenciales. | Preferencias permitidas persisten; secretos no se guardan. | P1 | NOT TESTED — |
| UI-06 | Forzar excepción de una acción y continuar usando la app. | Mensaje controlado; estado del documento coherente. | P1 | NOT TESTED — |
| UI-07 | **Retest Windows portable** del preview [#38006731375](https://github.com/PJGV333/Chemuson/actions/runs/38006731375), [artifact 11651882043](https://github.com/PJGV333/Chemuson/actions/runs/38006731375/artifacts/11651882043), SHA-256 `dde7c8f94fef2fd344756c8aa25eac53bd998c2ff8a007021f4e05f014e00bc6`. Abrir en temas claro/oscuro; comprobar puntero/selección, enlace simple, anillo aromático, buscar, deshacer, rehacer, documento nuevo y limpieza. | Cada control muestra un icono visible y legible en ambos temas; adjuntar captura y SHA del portable. | P1 | NOT TESTED — owner retest required; automated frozen icon smoke passed. |
| UI-08 | **Retest Linux AppImage Type 2** del preview [#38006731375](https://github.com/PJGV333/Chemuson/actions/runs/38006731375), [artifact 11650889195](https://github.com/PJGV333/Chemuson/actions/runs/38006731375/artifacts/11650889195), SHA-256 `ec26ce53bfd6202d72e310e160c0bbc5a32c091deafcf6cc5b6fdeb066826e9c`. Abrir/extraer en temas claro/oscuro; comprobar puntero/selección, enlace simple, anillo aromático, buscar, deshacer, rehacer, documento nuevo y limpieza. | Cada control muestra un icono visible y legible en ambos temas; adjuntar captura, SHA y confirmar tipo Type 2. | P1 | NOT TESTED — owner retest required; automated frozen icon smoke passed. |
| UI-ONBOARDING-001 | **Retest del onboarding en Windows portable y Linux.** En primera ejecución, comprobar rail, lienzo y panel lateral en 980×600, 1440×900 y 1600×900; probar escalas 100%, 125%, 150% y 200%, mover/redimensionar la ventana, avanzar/retroceder/cerrar y reiniciar con/sin «No volver a mostrar». | La máscara cubre la ventana cliente, cada agujero coincide con su objetivo, la tarjeta queda visible, el rail no se desplaza y las preferencias/cierre funcionan igual en Windows y Linux. | P2 | NOT TESTED — previous Windows report; owner retest required |
| BRANDING-001 | **Retest de marca en Windows portable/setup y Linux.** Revisar título y barra superior, Acerca de/Ayuda, mensajes visibles, nombre mostrado del instalador y desinstalador, launcher y AppStream. | Todas las superficies presentan `ChemUSON`; el título es `ChemUSON 0.3.0-beta.1 — Editor Molecular Libre` y Acerca de muestra la descripción aprobada; instalación/actualización/desinstalación detectan la identidad existente y permanecen operativas. | P2 | NOT TESTED — owner retest required |

## H. Paquetes, actualización y separación de canales

| ID | Pasos y datos | Resultado esperado | Severidad | Resultado / evidencia |
|---|---|---|---|---|
| DIST-01 | Verificar SHA-256 del Windows portable descargado contra `checksums.sha256`. | Hash coincide con el paquete del mismo artifact group/provenance. | P1 | NOT TESTED — |
| DIST-02 | Ejecutar Windows portable en VM/perfil aislado. | Arranca con versión del manifest; no instala ni actualiza canal público. | P1 | NOT TESTED — |
| DIST-03 | Instalar Windows setup en VM, iniciar, cerrar y desinstalar. | Setup muestra versión preparada y se instala/desinstala sin pérdida inesperada. | P1 | NOT TESTED — |
| DIST-04 | Comparar SHA/versión entre portable y setup. | Ambos declaran misma versión y source SHA del run; nombres dicen preview. | P1 | NOT TESTED — |
| DIST-05 | Ejecutar el Linux portable/AppImage Type 2 corregido en distro/VM soportada. | Lanza con dependencias esperadas; verificar firma `AI\x02` y extracción sin FUSE antes de la prueba visual. | P1 | NOT TESTED — |
| DIST-06 | Inspeccionar preview portable: `.updateinfo`, `.update.json`, `.zsync`. | Los tres sidecars de updater público están ausentes. | P1 | NOT TESTED — |
| DIST-07 | Instalar bundle Flatpak preview desde archivo local, sin añadir remoto ChemUSON. | Bundle se instala y abre en rama/nombre de preview; no crea remoto beta/stable. | P1 | NOT TESTED — |
| DIST-08 | Confirmar que instalar preview no cambia manifiestos beta/stable ni el updater instalado. | Canales públicos y versiones disponibles siguen idénticos al iniciar la prueba. | P0 | NOT TESTED — |
| DIST-09 | Abrir los cuatro manifests y comprobar versión, rama, SHA y `publication=false`. | Los cuatro grupos apuntan a la misma ejecución/commit y checksums verificables. | P1 | NOT TESTED — |
| DIST-10 | Probar actualización de una instalación estable con preview presente en otra VM. | Instalación estable no ofrece ni instala el preview. | P0 | NOT TESTED — |

### Incidencia P1 observada en el primer preview (estado histórico, no retest)

| ID | Artifact/run observado | Pasos y resultado | Severidad | Estado / evidencia |
|---|---|---|---|---|
| UI-FAIL-WIN-01 | Windows portable, Actions run `37826597134` | En tema claro y oscuro, botones de barra lateral y controles superiores aparecen sin numerosos iconos; C/T/esfera y otros símbolos sí aparecen. Ventana y bienvenida abren, versión `0.3.0-beta.1`. | P1 | **FAILED — P1 blocks beta acceptance** — propietario |
| UI-FAIL-LINUX-01 | Linux portable del run `37826597134` (binario PyInstaller renombrado, no Type 2) | En tema claro y oscuro faltan numerosos iconos de barras de dibujo, herramientas y controles superiores. | P1 | **FAILED — P1 blocks beta acceptance** — propietario |

Estos defectos confirman la baseline del candidato, no el estado de los paquetes corregidos. UI-07/UI-08 siguen `NOT TESTED` hasta el retest manual propietario.

## Correcciones UI/identidad pendientes de retest

- **UI-ONBOARDING-001:** el propietario observó la tarjeta desalineada y el spotlight/rail incorrectos en Windows portable. La fuente ahora difiere el inicio hasta el primer layout visible y recalcula geometría; el resultado manual del candidato anterior permanece **FAILED**, y el paquete corregido queda **NOT TESTED — owner retest required**.
- **BRANDING-001:** las cadenas y metadatos visibles se normalizaron a `ChemUSON`; aún no se ha revisado el paquete resultante en Windows/Linux. Resultado manual: **NOT TESTED — owner retest required**.
- **P1 — RDKit isolated backend unavailable in packaged executable:** el fallo portable Windows informado queda **FAILED**; `RDKIT-LINUX-01` y `RDKIT-3D-SMILES-01` permanecen **NOT TESTED**. Los tests del ejecutable congelado son gate de build, no sustituyen la aceptación manual del propietario.
- **P1 — ChemName templates absent from frozen executables:** Windows portable y AppImage anteriores quedan **FAILED** por cero `.mol`; el Flatpak anterior tenía 9/9 archivos, pero no ejecutó su smoke. Los tres formatos del Build Preview `37866732971` pasan sus gates automatizados (9/9 templates, 6/6 nombres); `CHEMNAME-RETEST-01` sigue **NOT TESTED — owner retest required** para los artifacts exactos.

UI-ONBOARDING-001, BRANDING-001 y RDKit sólo pueden actualizarse con la versión y SHA exactos de los nuevos artifacts, entorno real y evidencia del propietario. Los tests estáticos, smoke del ejecutable congelado y Qt scale-factors simulados no sustituyen la comprobación manual.

## Severidad, promoción y cierre

- **P0:** pérdida/corrupción química o `.cmsn`, exposición de credenciales, canal/update público alterado por preview, no arranque en plataformas soportadas, identidad química crítica incorrecta. Bloquea beta/stable hasta resolver.
- **P1:** flujo esencial roto, crash reproducible de app empaquetada, instalación/actualización incorrecta, datos no reversibles o limitación de privacidad. En general bloquea stable; cualquier excepción beta requiere decisión explícita del propietario y riesgo acotado.
- **Regla específica de esta campaña:** los P1 confirmados (iconos ausentes, fallo RDKit empaquetado y templates ChemName ausentes) bloquean la aceptación y publicación beta hasta que el propietario reteste los paquetes corregidos en Windows y Linux, registre resultados por formato y decida cualquier incidencia. Los smokes automatizados no satisfacen el gate manual; estado actual `BETA PUBLICATION: BLOCKED — AWAITING MANUAL RETEST`.
- **P2:** defecto no crítico con workaround claro. Puede aceptarse sólo con decisión del propietario y seguimiento.
- **P3:** cosmético/menor. Registrar; aceptación requiere decisión explícita para stable.

Para declarar un caso finalizado, reemplazar `NOT TESTED`, añadir SHA del artifact, entorno, tester, resultado observado y enlace a evidencia. Los fallos históricos del test harness Qt no son equivalentes a un crash reproducido de la aplicación empaquetada, y tampoco se consideran cerrados por no reproducirse en una prueba aislada. Stable requiere aprobación del propietario y cero P0/P1 atribuibles al candidato.
