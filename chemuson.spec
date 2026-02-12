"""PyInstaller spec para Chemuson con recolección automática de Qt/RDKit."""

from pathlib import Path

from PyInstaller.utils.hooks import collect_all, collect_dynamic_libs

datas_c, binaries_c, hidden_c = collect_all("chemuson")
datas_qt, binaries_qt, hidden_qt = collect_all("PyQt6")
binaries_rdkit = collect_dynamic_libs("rdkit")

datas = datas_c + datas_qt
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
