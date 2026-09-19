# M21 — Configuración y recursos de plataforma

## Responsabilidad

M21 (`platform.settings`) posee la política de configuración persistente de la
aplicación y la resolución de recursos empaquetados. No importa módulos
ChemUSON, `chemuson.gui` ni `QtWidgets`.

## API

- `application_settings`
- `setting_bool`
- `NamingPreferences` y `NumberingPreferences`
- `load_naming_preferences` / `save_naming_preferences`
- `load_numbering_preferences` / `save_numbering_preferences`
- `open_resource_path` en el módulo interno `platform.resources`

## Compatibilidad

`chemuson.utils.resources` conserva un shim de importación para consumidores
históricos. Las claves `naming/*` y `numbering/*`, sus defaults y su coerción se
mantienen sin cambios.

## Límites

M21 no posee estado documental, selección, undo, actualización, widgets ni
composición de la ventana. La GUI puede depender de M21; M21 no depende de la
GUI.
