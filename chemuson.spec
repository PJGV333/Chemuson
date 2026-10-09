"""PyInstaller spec para Chemuson con recolección automática de Qt/RDKit."""

from pathlib import Path

from PyInstaller.utils.hooks import collect_all, collect_dynamic_libs

PROJECT_ROOT = Path(SPECPATH).resolve()
STATIC_ICON_DIR = PROJECT_ROOT / "src" / "chemuson" / "gui" / "theme" / "icons"
STATIC_ICON_FILES = sorted(STATIC_ICON_DIR.glob("i-*.svg"))
if len(STATIC_ICON_FILES) != 69 or any(path.stat().st_size == 0 for path in STATIC_ICON_FILES):
    raise SystemExit(
        f"Expected 69 nonempty ChemUSON static SVG icons in {STATIC_ICON_DIR}; "
        f"found {len(STATIC_ICON_FILES)}."
    )
if not (STATIC_ICON_DIR / "LICENSE.txt").is_file():
    raise SystemExit(f"ChemUSON static icon license is missing: {STATIC_ICON_DIR}")

datas_c, binaries_c, hidden_c = collect_all("chemuson")
# These entries are enumerated explicitly below, even when collect_all can see
# an installed/editable copy of the source package.
datas_c = [
    (source, destination)
    for source, destination in datas_c
    if not destination.replace(chr(92), "/").startswith(
        ("chemuson/gui/theme/icons/", "chemuson/chemname/templates/")
    )
]
datas_icons = [
    (str(path), "chemuson/gui/theme/icons") for path in STATIC_ICON_FILES
]
datas_icons.append((str(STATIC_ICON_DIR / "LICENSE.txt"), "chemuson/gui/theme/icons"))
CHEMNAME_TEMPLATE_DIR = PROJECT_ROOT / "src" / "chemuson" / "chemname" / "templates"
EXPECTED_CHEMNAME_TEMPLATES = {
    "fused/pyrene_cas.mol",
    "fused/pyrene_iupac2004.mol",
    "simple/benzene.mol",
    "special/alpha_d_glucopyranose.mol",
    "special/androstane_core.mol",
    "special/beta_d_fructofuranose.mol",
    "special/beta_d_glucopyranose.mol",
    "special/cholestane_core.mol",
    "special/d_ribose.mol",
}
CHEMNAME_TEMPLATE_FILES = sorted(CHEMNAME_TEMPLATE_DIR.glob("*/*.mol"))
CHEMNAME_TEMPLATE_INVENTORY = {
    path.relative_to(CHEMNAME_TEMPLATE_DIR).as_posix() for path in CHEMNAME_TEMPLATE_FILES
}
if CHEMNAME_TEMPLATE_INVENTORY != EXPECTED_CHEMNAME_TEMPLATES or any(
    path.stat().st_size == 0 for path in CHEMNAME_TEMPLATE_FILES
):
    raise SystemExit(
        "ChemName template inventory is incomplete or unexpected: "
        f"{sorted(CHEMNAME_TEMPLATE_INVENTORY)}"
    )
datas_chemname_templates = [
    (str(path), f"chemuson/chemname/templates/{path.parent.name}")
    for path in CHEMNAME_TEMPLATE_FILES
]
datas_qt, binaries_qt, hidden_qt = collect_all("PyQt6")
binaries_rdkit = collect_dynamic_libs("rdkit")

datas = datas_c + datas_chemname_templates + datas_icons + datas_qt
binaries = binaries_c + binaries_qt + binaries_rdkit
hiddenimports = sorted(
    set(
        hidden_c
        + hidden_qt
        + [
            "PyQt6.QtSvg",
            "PyQt6.QtPrintSupport",
            "chemuson.chemio._rdkit_worker",
            "chemuson.chemio.rdkit_packaged_smoke",
            "chemuson.chemname.packaged_smoke",
            "rdkit",
            "rdkit.Chem.AllChem",
            "rdkit.Chem.rdchem",
            "rdkit.Chem.rdDistGeom",
            "rdkit.Chem.rdMolDescriptors",
            "rdkit.Chem.rdForceFieldHelpers",
        ]
    )
)

entry_script = str(Path("src") / "chemuson" / "__main__.py")

a = Analysis(
    [entry_script],
    pathex=["src"],
    binaries=binaries,
    datas=datas,
    hiddenimports=hiddenimports,
    runtime_hooks=["packaging/pyinstaller/rthook_qt.py"],
    noarchive=False,
)

pyz = PYZ(a.pure)

exe = EXE(
    pyz,
    a.scripts,
    a.binaries,
    a.datas,
    [],
    name="Chemuson",
    console=False,
    debug=False,
    strip=False,
    upx=False,
    disable_windowed_traceback=False,
)
