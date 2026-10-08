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
    if not destination.replace(chr(92), "/").startswith("chemuson/gui/theme/icons/")
]
datas_icons = [
    (str(path), "chemuson/gui/theme/icons") for path in STATIC_ICON_FILES
]
datas_icons.append((str(STATIC_ICON_DIR / "LICENSE.txt"), "chemuson/gui/theme/icons"))
datas_qt, binaries_qt, hidden_qt = collect_all("PyQt6")
binaries_rdkit = collect_dynamic_libs("rdkit")

datas = datas_c + datas_icons + datas_qt
binaries = binaries_c + binaries_qt + binaries_rdkit
hiddenimports = sorted(
    set(hidden_c + hidden_qt + ["PyQt6.QtSvg", "PyQt6.QtPrintSupport", "rdkit"])
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
