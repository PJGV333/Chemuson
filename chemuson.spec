"""PyInstaller spec para Chemuson."""

from pathlib import Path

from PyInstaller.utils.hooks import collect_all, collect_dynamic_libs, collect_submodules

datas, binaries, hiddenimports = collect_all("chemuson")

hiddenimports += [
    "PyQt6.QtSvg",
    "PyQt6.QtPrintSupport",
]
hiddenimports += collect_submodules("PyQt6")
binaries += collect_dynamic_libs("rdkit")

entry_script = str(Path("src") / "chemuson" / "__main__.py")

a = Analysis(
    [entry_script],
    pathex=["src"],
    binaries=binaries,
    datas=datas,
    hiddenimports=hiddenimports,
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
)
