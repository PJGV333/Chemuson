window.CHEMUSON_MODULES = [
 {
  "id": "M00",
  "name": "core",
  "title": "Modelo molecular fundamental",
  "responsibility": "Modelo molecular independiente de UI, estado químico, validación básica, análisis elemental y grafo químico multicapa.",
  "status": "stable",
  "risk": "low",
  "deps": [],
  "public_api": [
   "Atom",
   "Bond",
   "BondStyle",
   "BondStereo",
   "ChemState",
   "MolGraph",
   "MolView",
   "MolecularViewNotSupported",
   "ValidationCorrectionAction",
   "ValidationIssue",
   "BlockEdge",
   "BlockEdgeKind"
  ]
 },
 {
  "id": "M01",
  "name": "chemio",
  "title": "Import/export químico y persistencia",
  "responsibility": "Import/export químico y persistencia: SMILES, Molfile, CML, CMSN y acceso seguro a RDKit.",
  "status": "evolving",
  "risk": "medium",
  "deps": [
   "M00"
  ],
  "public_api": []
 },
 {
  "id": "M02",
  "name": "clean2d",
  "title": "Limpieza/depiction 2D",
  "responsibility": "Limpieza/depiction 2D desacoplada de la GUI: candidatos, ranking, invariantes, seguridad geométrica, reparación local y políticas de moléculas complejas.",
  "status": "evolving",
  "risk": "high",
  "deps": [
   "M00",
   "M01"
  ],
  "public_api": [
   "Clean2DParameters",
   "optimize_clean2d_positions",
   "Clean2DQualityReport",
   "evaluate_clean2d_layout",
   "is_clean2d_candidate_safe",
   "has_cycles",
   "count_new_bond_crossings",
   "min_nonbonded_distance",
   "ring_degeneracy_score",
   "max_atom_displacement",
   "bond_length_stats",
   "length_only_polish"
  ]
 },
 {
  "id": "M03",
  "name": "chemcalc",
  "title": "Cálculos químicos auxiliares",
  "responsibility": "Fórmula, masa molecular, valencias típicas e hidrógenos implícitos auxiliares.",
  "status": "stable",
  "risk": "low",
  "deps": [
   "M00"
  ],
  "public_api": [
   "molecular_formula",
   "format_formula",
   "molecular_weight",
   "implicit_h_count",
   "TYPICAL_VALENCE"
  ]
 },
 {
  "id": "M04",
  "name": "chemname",
  "title": "Nomenclatura IUPAC-lite",
  "responsibility": "IUPAC-lite y nomenclatura semisistemática: selección de cadenas/anillos, plantillas, grupos funcionales, estereoquímica y coordinación.",
  "status": "evolving",
  "risk": "medium",
  "deps": [
   "M00",
   "M03",
   "M01",
   "M21"
  ],
  "public_api": [
   "iupac_name",
   "NameOptions"
  ]
 },
 {
  "id": "M05",
  "name": "geometry3d",
  "title": "Servicios 3D",
  "responsibility": "Modelos y servicios 3D, generación/proyección de conformeros, cache y export XYZ.",
  "status": "stable",
  "risk": "medium",
  "deps": [
   "M00",
   "M01"
  ],
  "public_api": [
   "Conformer3DResult",
   "CoordinateSet3D",
   "cache_key_for_3d",
   "DepthCue",
   "ForceField",
   "OptimizationFrame",
   "OptimizationResult",
   "OptimizationSettings",
   "SceneAtom3D",
   "SceneBond3D",
   "SceneMolecule3D",
   "ProjectedAtom3D"
  ]
 },
 {
  "id": "M06",
  "name": "compchem",
  "title": "Exportación química computacional",
  "responsibility": "Exportadores de química computacional y especificaciones de entrada para Gaussian/ORCA/NWChem.",
  "status": "stable",
  "risk": "low",
  "deps": [
   "M00",
   "M05"
  ],
  "public_api": [
   "export_gaussian_input",
   "export_nwchem_input",
   "export_orca_input"
  ]
 },
 {
  "id": "M07",
  "name": "spectroscopy",
  "title": "Predicción espectral",
  "responsibility": "Predicción MVP de NMR/masa y registro de predictores espectrales.",
  "status": "evolving",
  "risk": "medium",
  "deps": [
   "M00",
   "M01"
  ],
  "public_api": [
   "CarbonNmrPeak",
   "MassPeak",
   "ProtonNmrPeak",
   "SpectralPrediction",
   "SpectrumPredictor",
   "predict_spectra",
   "register_predictor"
  ]
 },
 {
  "id": "M08",
  "name": "gui",
  "title": "Orquestación de interfaz PyQt6",
  "responsibility": "UI PyQt6: ventana principal, canvas/editor 2D, acciones, docks, controllers, comandos undo/redo, rendering y herramientas visuales.",
  "status": "evolving",
  "risk": "high",
  "deps": [
   "M00",
   "M01",
   "M05",
   "M06",
   "M07",
   "M09",
   "M10",
   "M11",
   "M12",
   "M13",
   "M14",
   "M16",
   "M17",
   "M18",
   "M21",
   "M22"
  ],
  "public_api": [
   "ChemusonWindow"
  ]
 },
 {
  "id": "M09",
  "name": "gui.canvas",
  "title": "Canvas de edición molecular",
  "responsibility": "Canvas de edición molecular interactivo basado en QGraphicsView con arquitectura de mixins.",
  "status": "evolving",
  "risk": "high",
  "deps": [
   "M00",
   "M01",
   "M04",
   "M05",
   "M08",
   "M11",
   "M12",
   "M13",
   "M20"
  ],
  "public_api": [
   "AROMATIC_CIRCLE_ATOMS_ROLE",
   "BRANCH_ROTATION_NOOP_TOLERANCE_DEG",
   "BRANCH_ROTATION_STEP_DEG",
   "ChemusonCanvas",
   "FRAGMENT_ROTATION_STEP_DEG",
   "molgraph_to_molfile",
   "molgraph_to_smiles",
   "QInputDialog"
  ]
 },
 {
  "id": "M10",
  "name": "gui.controllers",
  "title": "Controllers de la GUI",
  "responsibility": "Controllers de la GUI: orquestación de canvas, documentos, exportación, actualización y validación.",
  "status": "evolving",
  "risk": "medium",
  "deps": [
   "M00",
   "M01",
   "M02",
   "M05",
   "M08",
   "M09",
   "M11",
   "M14",
   "M21",
   "M22"
  ],
  "public_api": [
   "Clean2DController",
   "CompChem3DController",
   "CompChem3DWorker",
   "CompChemJobSpec",
   "DocumentController",
   "DocumentDiscardContext",
   "DocumentTabsContext",
   "ExportController",
   "FileController",
   "FileWorkflowContext",
   "RecentFilesContext",
   "RecoveryController"
  ]
 },
 {
  "id": "M11",
  "name": "gui.commands",
  "title": "Comandos undo/redo",
  "responsibility": "Comandos undo/redo para el canvas: operaciones atómicas sobre átomos, enlaces, texto y diagramas.",
  "status": "evolving",
  "risk": "medium",
  "deps": [
   "M00",
   "M08"
  ],
  "public_api": [
   "AddAtomCommand",
   "ChangeAtomCommand",
   "ChangeChargeCommand",
   "ChangeNoImplicitCommand",
   "ChangeAtomLabelScaleCommand",
   "SetCoordinationCenterCommand",
   "ChangeCoordinationSphereStyleCommand",
   "AddBondCommand",
   "ChangeBondCommand",
   "ChangeBondLengthCommand",
   "ChangeBondStrokeCommand",
   "ChangeBondColorCommand"
  ]
 },
 {
  "id": "M12",
  "name": "gui.dialogs",
  "title": "Diálogos de la GUI",
  "responsibility": "Diálogos de la GUI: preferencias, estilos, selección, inserción de placas/gels, guía rápida.",
  "status": "evolving",
  "risk": "low",
  "deps": [
   "M00",
   "M08"
  ],
  "public_api": []
 },
 {
  "id": "M13",
  "name": "gui.items",
  "title": "Items gráficos (átomos, enlaces, etc.)",
  "responsibility": "Items gráficos (átomos, enlaces, etc.) para la escena de Qt: subclases de QGraphicsItem.",
  "status": "evolving",
  "risk": "medium",
  "deps": [
   "M00",
   "M08"
  ],
  "public_api": []
 },
 {
  "id": "M14",
  "name": "update",
  "title": "Subsistema de auto-actualización",
  "responsibility": "Subsistema de auto-actualización: política, proveedor, seguridad, rollback, portable y telemetría.",
  "status": "stable",
  "risk": "low",
  "deps": [],
  "public_api": [
   "AutoUpdateCore",
   "VerificationResult",
   "GitHubReleasesProvider",
   "detect_platform_tag",
   "RollbackManager",
   "SignatureVerifier",
   "UpdateTelemetryLogger",
   "parse_semver",
   "compare_versions",
   "is_newer_version",
   "is_prerelease",
   "channel_accepts_version"
  ]
 },
 {
  "id": "M15",
  "name": "utils",
  "title": "Utilidades compartidas",
  "responsibility": "Shims históricos de compatibilidad para servicios canónicos de plataforma y resiliencia.",
  "status": "stable",
  "risk": "low",
  "deps": [
   "M21",
   "M22"
  ],
  "public_api": []
 },
 {
  "id": "M16",
  "name": "name2structure",
  "title": "Resolución nombre a estructura",
  "responsibility": "Resolver nombres químicos a estructuras con conectores estáticos/PubChem y fallback seguro.",
  "status": "stable",
  "risk": "medium",
  "deps": [
   "M00",
   "M01"
  ],
  "public_api": [
   "NameToStructureResult",
   "PubChemNameConnector",
   "StaticNameConnector",
   "resolve_name_to_structure"
  ]
 },
 {
  "id": "M17",
  "name": "markush",
  "title": "Estructuras Markush y polímeros",
  "responsibility": "Estructuras Markush y polímeros: R-groups, grupos repetitivos, resumen y sanitización.",
  "status": "stable",
  "risk": "low",
  "deps": [
   "M00"
  ],
  "public_api": [
   "MarkushSummary",
   "PolymerRepeat",
   "RGroupAtom",
   "sanitize_r_group_substituents",
   "set_r_group_substituents",
   "summarize_markush"
  ]
 },
 {
  "id": "M18",
  "name": "version",
  "title": "Gestión de versión",
  "responsibility": "Gestión de versión: fuente única de __version__, metadata de packaging.",
  "status": "stable",
  "risk": "low",
  "deps": [],
  "public_api": [
   "__version__",
   "get_app_version"
  ]
 },
 {
  "id": "M19",
  "name": "bootstrap",
  "title": "Arranque y composición de la aplicación",
  "responsibility": "Punto de entrada CLI, parsing inicial y arranque diferido de la GUI.",
  "status": "stable",
  "risk": "medium",
  "deps": [
   "M18",
   "M08",
   "M22"
  ],
  "public_api": [
   "main"
  ]
 },
 {
  "id": "M20",
  "name": "gui.editor2d.selection",
  "title": "Selección del editor 2D",
  "responsibility": "Geometría, hit testing, overlays y política de clipboard consultivas de la selección del editor 2D.",
  "status": "evolving",
  "risk": "low",
  "deps": [],
  "public_api": []
 },
 {
  "id": "M21",
  "name": "platform.settings",
  "title": "Configuración y recursos de plataforma",
  "responsibility": "Configuración persistente de aplicación y resolución de recursos empaquetados sin dependencia de GUI.",
  "status": "evolving",
  "risk": "medium",
  "deps": [],
  "public_api": [
   "NamingPreferences",
   "NumberingPreferences",
   "UI_THEME_CHOICES",
   "UiPreferences",
   "application_settings",
   "load_naming_preferences",
   "load_numbering_preferences",
   "load_ui_preferences",
   "save_naming_preferences",
   "save_numbering_preferences",
   "save_ui_preferences",
   "setting_bool"
  ]
 },
 {
  "id": "M22",
  "name": "resilience",
  "title": "Resiliencia y recuperación de runtime",
  "responsibility": "Persistencia de recuperación, autosave rotativo, logging de crashes y aislamiento de fallos de runtime.",
  "status": "evolving",
  "risk": "medium",
  "deps": [],
  "public_api": [
   "AutosaveManager",
   "archive_autosave",
   "install",
   "list_autosave_entries",
   "read_autosave_metadata",
   "write_crash_log"
  ]
 }
];
