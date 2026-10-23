# Propuesta — endurecer dependencias y permisos Flatpak de beta.1

## Por qué

El manifiesto Flatpak de beta.1 construía e incluía Pillow 12.0.0, afectado por múltiples advisories OSV (incluidos CVE-2026-42311 y CVE-2026-54058). ChemUSON no importa PIL: las anotaciones usan `QImage`/`QPixmap` de Qt y el módulo nativo `rdMolDraw2D` de RDKit funciona con PIL bloqueado. Pillow no es una dependencia directa necesaria; sin embargo, la metadata PyPI de RDKit declara Pillow transitivamente para los entornos PyInstaller. La inspección de Build Preview #39899351520 encontró cinco extensiones nativas `PIL/*` en los CArchive Windows y AppImage; su entorno de build resolvió Pillow 13.0.0 (sin advisories en la consulta OSV del 2026-10-23). Esta revisión elimina la dependencia directa y el módulo Flatpak, excluye PIL de los paquetes congelados y añade un gate que valida el CArchive. La beta queda sin recursos Pillow distribuibles aunque el entorno de build siga resolviendo la dependencia transitiva; no se reescribe el requisito histórico de conservarla.

El manifiesto también concede `--filesystem=home`, aunque los únicos accesos directos al home identificados son rutas propias de ChemUSON (`~/.chemuson` y el fallback `~/.config/Chemuson`); la interacción de archivos del usuario usa selectores Qt, que deben depender del portal Flatpak y requiere retest manual. `--share=network` sí permite las funciones elegibles de Assistant remoto y consulta opcional PubChem; se conserva y se hace explícito el riesgo residual.

## Qué cambia

- Eliminar Pillow de las dependencias directas y el módulo/wheel Flatpak; excluir los recursos `PIL` transitivos de Windows/AppImage PyInstaller y hacer que la inspección del CArchive falle si reaparecen. No se reemplaza por otro paquete de imágenes; la versión transitiva que resuelva RDKit en build se documenta y audita por separado.
- Limitar permisos de escritura/lectura directa a los directorios de datos/configuración propios de ChemUSON, en lugar de todo `home`.
- Mantener permiso de red por las funciones remotas existentes y documentar que el tráfico a proveedores configurados/consultas externas puede salir del dispositivo.
- Añadir contratos estáticos focalizados y capturar baseline antes de cambios.
- Crear una rama aislada desde `08a90c199f526639453a8ae8b82da957a1458d08`, producir un único Build Preview con los cuatro paquetes y verificar identidad, provenance y hashes del mismo SHA. No ejecutar pytest completo ni crear tag/Release/publicar Flatpak.

## Fuera de alcance

No cambiar persistencia `.cmsn`, workflows o permisos de publicación, lógica química, rendering Qt, proveedores, portales, updater ni arquitectura de módulos. No se amplía la aceptación manual: queda pendiente sobre los artefactos exactos nuevos.

## Aceptación

1. Pillow no aparece como dependencia directa ni módulo Flatpak; los CArchive Windows/AppImage no contienen entradas top-level `PIL`/`Pillow`, y el árbol Flatpak no contiene `site-packages/PIL` ni `Pillow-*.dist-info`. La dependencia transitiva de RDKit resuelta durante build se registra y no tiene advisory OSV para el SHA construido. Las rutas de imagen/`rdMolDraw2D` relevantes pasan sin que PIL sea importable.
2. El manifiesto no concede `--filesystem=home`; conserva únicamente las rutas de datos/configuración justificadas y `--share=network`. La build Flatpak valida el manifiesto y empaqueta correctamente.
3. Tests focalizados de dependencia, Flatpak, imágenes, ChemName, RDKit y 2D→3D pasan; compileall, OpenSpec y controles estáticos pasan. La suite pytest completa no se ejecuta.
4. El nuevo Build Preview contiene Windows portable, Windows installer, AppImage Type 2 y Flatpak con el SHA completo de esta rama, `publication=false`, checksums y provenance concordantes.
5. Los casos manuales P0/P1 siguen `NOT TESTED` hasta que el propietario pruebe esos archivos exactos. No se crea tag, Release ni despliegue público; dictamen permanece NO-GO mientras la aceptación o los controles externos estén pendientes.
