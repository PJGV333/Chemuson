# RFC: Estrategia de Distribución Híbrida y Auto-Update para Chemuson

Estado: Propuesto  
Autor: Equipo Chemuson  
Fecha: 2026-02-26

## 1. Resumen ejecutivo

Se propone una estrategia híbrida de distribución para Chemuson que mantenga los binarios portables actuales (`.AppImage` y `.exe`) y agregue instaladores oficiales:

- Windows: instalador oficial + actualizador con rollback.
- Linux: Flatpak como canal instalable recomendado + mantenimiento de AppImage portable.

El sistema de actualización se centraliza en GitHub Releases con manifiestos por canal (`stable` y `beta`), verificación criptográfica y política configurable (silenciosa o con aviso).

## 2. Objetivos y no objetivos

### Objetivos

- Mantener experiencia portable sin instalación.
- Ofrecer instalación oficial en Windows y Linux.
- Implementar auto-update robusto con:
  - versionado semántico,
  - canales `stable`/`beta`,
  - verificación de integridad/autenticidad,
  - rollback en fallo.
- Permitir política de actualización configurable por usuario/organización.
- Firmar artefactos de release.

### No objetivos (en esta RFC)

- Cambiar la arquitectura funcional de Chemuson.
- Implementar telemetría obligatoria.
- Definir empaquetado nativo de todas las distribuciones Linux (`deb/rpm`) en primera etapa.

## 3. Contexto actual

- Ya existen binarios portables para Windows (`.exe`) y Linux (`.AppImage`).
- No existe aún un flujo formal de instalación oficial con actualización automática y firma end-to-end.

## 4. Opciones evaluadas

## 4.1 Instalador y update en Windows

### Opción A: Inno Setup + actualizador propio (helper externo)

Pros:
- Control total del flujo.
- Compatible con GitHub Releases sin infraestructura adicional.
- Fácil convivencia con binario portable.

Contras:
- Mayor esfuerzo de implementación/mantenimiento.
- Requiere manejo explícito de privilegios y bloqueo de archivos.

### Opción B: NSIS + plugin/update framework

Pros:
- Muy usado en ecosistema Windows.
- Flexibilidad alta.

Contras:
- Similar complejidad operativa que Inno.
- Mayor riesgo de fragmentación en scripts.

### Opción C: MSIX + App Installer / winget

Pros:
- Modelo moderno de instalación y update.
- Integración nativa con firma y políticas empresariales.

Contras:
- Curva operativa superior (certificados, packaging, restricciones sandbox).
- Menor alineación con el patrón actual de binarios portables.

## 4.2 Linux instalable

### Opción A: Flatpak (preferido)

Pros:
- Estándar de facto para distribución desktop transversal.
- Actualización gestionada por runtime/remote.
- Buen aislamiento y reproducibilidad.

Contras:
- Requiere manifest y pipeline específico.
- Integración host depende de permisos/portals.

### Opción B: Snap

Pros:
- Actualizaciones automáticas integradas.

Contras:
- Dependencia de infraestructura Snap Store.
- Menor preferencia en parte de la comunidad Linux desktop científica.

### Opción C: Solo AppImage

Pros:
- Muy simple para usuario final.

Contras:
- No cubre bien flujos instalables administrados.
- Auto-update más frágil y heterogéneo.

## 4.3 Motor de actualización

### Opción A: Librería third-party de auto-update

Pros:
- Menor tiempo inicial.

Contras:
- Riesgo de dependencia no mantenida.
- Menor control sobre canales/firma/rollback.

### Opción B: Cliente de update propio + manifiesto firmado

Pros:
- Ajuste exacto a requerimientos (`stable/beta`, rollback, políticas).
- Control de seguridad y trazabilidad.

Contras:
- Coste inicial mayor.

## 5. Decisión recomendada

1. Mantener binarios portables existentes (`.AppImage`, `.exe`).
2. Windows: adoptar **Inno Setup** (instalador) + **updater helper propio**.
3. Linux: adoptar **Flatpak** como canal instalable principal y mantener AppImage.
4. Auto-update: **cliente propio** consumiendo manifiesto firmado desde GitHub Releases.
5. Firma: esquema de doble capa:
   - firma del manifiesto,
   - verificación de hash de cada artefacto.

## 6. Arquitectura de actualización

## 6.1 Componentes

- `UpdateClient` (en app): chequeo de versión, canal, política, descarga y orquestación.
- `UpdaterHelper` (proceso externo): aplica reemplazo/instalación y rollback.
- `Release Manifest` por canal (`stable.json`, `beta.json`).
- `Signature Verifier`: valida firma del manifiesto y hash del artefacto.
- `Rollback Manager`: conserva versión previa y revierte si el arranque post-update falla.

## 6.2 Flujo de actualización

1. App inicia y lee política (`silent`/`notify`) + canal (`stable`/`beta`).
2. `UpdateClient` descarga manifiesto del canal desde GitHub Releases.
3. Verifica firma del manifiesto con clave pública embebida.
4. Compara versión local vs remota con SemVer (incluyendo prereleases).
5. Si hay update elegible:
   - `notify`: solicita confirmación.
   - `silent`: descarga en segundo plano y aplica según ventana de mantenimiento.
6. Descarga artefacto, valida hash SHA-256 contra manifiesto.
7. Lanza `UpdaterHelper` para aplicar update atómico.
8. Al reiniciar, se ejecuta health-check mínimo de arranque.
9. Si falla, `Rollback Manager` restaura versión previa.

## 6.3 SemVer y canales

Reglas propuestas:

- `stable` acepta versiones sin sufijo prerelease (`X.Y.Z`).
- `beta` acepta prerelease (`X.Y.Z-beta.N`) y finales.
- Usuario en `stable` no recibe prereleases.
- Usuario en `beta` sí puede recibir prereleases más nuevas.

## 6.4 Manifiesto de actualización (propuesto)

```json
{
  "channel": "stable",
  "latest": "1.4.0",
  "published_at": "2026-02-26T18:00:00Z",
  "min_supported": "1.2.0",
  "artifacts": {
    "windows-x86_64-installer": {
      "version": "1.4.0",
      "url": "https://github.com/PJGV333/Chemuson/releases/download/v1.4.0/Chemuson-1.4.0-setup.exe",
      "sha256": "..."
    },
    "windows-x86_64-portable": {
      "version": "1.4.0",
      "url": "https://github.com/PJGV333/Chemuson/releases/download/v1.4.0/Chemuson-1.4.0-portable.exe",
      "sha256": "..."
    },
    "linux-x86_64-appimage": {
      "version": "1.4.0",
      "url": "https://github.com/PJGV333/Chemuson/releases/download/v1.4.0/Chemuson-1.4.0-x86_64.AppImage",
      "sha256": "..."
    }
  },
  "signature": {
    "algorithm": "ed25519",
    "value": "base64(...)"
  }
}
```

Nota: Flatpak se actualiza por su propio mecanismo (`flatpak update`), pero puede usar el mismo canal de release para etiquetado y comunicación.

## 6.5 Política de actualización configurable

Parámetros sugeridos:

- `update.channel`: `stable` | `beta`
- `update.mode`: `silent` | `notify`
- `update.check_on_startup`: `true`/`false`
- `update.check_interval_hours`: entero (ej. 24)
- `update.allow_prerelease`: bool (sincronizado con canal, opcional avanzado)
- `update.allow_rollback`: bool (default `true`)

Reglas:

- En `silent`, se permite descarga automática y aplicación diferida.
- En `notify`, siempre se solicita acción del usuario.
- En entornos administrados, valores pueden bloquearse por archivo/política de sistema.

## 7. Matriz de compatibilidad objetivo

| Plataforma | Formato | Estado objetivo | Auto-update | Firma |
|---|---|---|---|---|
| Windows x86_64 | Portable `.exe` | Mantener | Sí (cliente + helper) | Sí |
| Windows x86_64 | Instalador `.exe` | Nuevo oficial | Sí (cliente + helper) | Sí |
| Linux x86_64 | `.AppImage` | Mantener | Sí (cliente) | Sí |
| Linux x86_64 | Flatpak | Nuevo oficial preferido | Gestionado por Flatpak | Sí (artefacto/repo) |

## 8. Firma de artefactos

Estrategia recomendada:

- Firmar manifiesto por canal con clave offline (Ed25519) y distribuir clave pública en la app.
- Publicar `sha256` de cada artefacto dentro del manifiesto firmado.
- En Windows, además firmar binarios/instalador con certificado de code signing para reducir fricción con SmartScreen.

Modelo de confianza:

- Si falla firma de manifiesto: abortar update.
- Si falla hash de artefacto: abortar update.
- Sin fallback inseguro por defecto.

## 9. Riesgos y mitigaciones

- Riesgo: bloqueo de ejecutable durante update en Windows.  
  Mitigación: aplicar update desde proceso externo (`UpdaterHelper`) con swap atómico y reinicio.

- Riesgo: corrupción de descarga o MITM.  
  Mitigación: TLS + firma de manifiesto + hash por artefacto.

- Riesgo: actualización defectuosa deja app inutilizable.  
  Mitigación: backup de versión previa + health-check de arranque + rollback automático.

- Riesgo: desalineación entre AppImage updater y Flatpak updater.  
  Mitigación: separar claramente rutas de update por formato y documentarlas.

- Riesgo: compromiso de clave de firma.  
  Mitigación: rotación de claves, revocación y lista de claves confiables versionada.

- Riesgo: rate limiting/API changes en GitHub.  
  Mitigación: consumir URL de manifiesto estática como asset; cache local con backoff.

## 10. Plan incremental por fases

### Fase 0: Base de release y versionado

- Normalizar versión SemVer única en app y pipeline.
- Definir canales `stable`/`beta`.
- Entregables: convención de tags/releases, checklist de publicación.

### Fase 1: Seguridad y metadatos de update

- Diseñar/emitir manifiestos por canal.
- Integrar firma de manifiesto y hashes en pipeline.
- Entregables: `stable.json`, `beta.json`, firmas verificables.

### Fase 2: Update client (modo notify)

- Implementar chequeo, comparación SemVer y descarga verificada.
- Aplicación manual asistida (sin silencioso aún).
- Entregables: actualización funcional con confirmación de usuario.

### Fase 3: Instalador Windows + helper + rollback

- Publicar instalador oficial Inno Setup.
- Implementar `UpdaterHelper` con reemplazo atómico y rollback.
- Entregables: update automático Windows robusto.

### Fase 4: Flatpak oficial

- Crear manifiesto Flatpak y pipeline de build/publicación.
- Documentar convivencia Flatpak/AppImage.
- Entregables: canal Linux instalable oficial.

### Fase 5: Políticas avanzadas y modo silent

- Activar política `silent` configurable.
- Añadir controles de administración local (archivo/política).
- Entregables: estrategia de actualización empresarial/doméstica configurable.

## 11. Criterios de aceptación

- Binarios portables siguen disponibles en cada release.
- Instalador Windows oficial publicado y actualizado automáticamente.
- Flatpak publicado como canal Linux recomendado.
- Update rechaza artefactos sin firma/hash válidos.
- Rollback comprobado en al menos un escenario de fallo inducido.

## 12. Preguntas abiertas

- ¿El instalador Windows será per-user por defecto o machine-wide?
- ¿La distribución Flatpak inicial irá a Flathub o a remote propio?
- ¿Se exigirá firma Authenticode desde la primera release con updater?
- ¿Se habilitará telemetría opcional de éxito/fallo de update?
