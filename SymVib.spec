# PyInstaller 設定: GUI 版 (SymVib.exe) をフォルダ形式で作る。
#   ビルド:  .\build_exe.ps1   (または  pyinstaller SymVib.spec --noconfirm)
#   出力:    dist\SymVib\SymVib.exe  (フォルダごと zip で配布)
from PyInstaller.utils.hooks import collect_data_files, collect_dynamic_libs

datas = collect_data_files("pymatgen") + collect_data_files("spglib") + collect_data_files("spgrep")
binaries = collect_dynamic_libs("spglib")

a = Analysis(
    ["vibration_gui.py"],
    pathex=["."],
    binaries=binaries,
    datas=datas,
    hiddenimports=["matplotlib.backends.backend_tkagg"],
    excludes=["pytest", "IPython", "notebook", "jupyter", "PyQt5", "PyQt6", "PySide2",
              "PySide6", "plotly", "seaborn", "sympy"],
    noarchive=False,
)
pyz = PYZ(a.pure)
exe = EXE(
    pyz,
    a.scripts,
    [],
    exclude_binaries=True,
    name="SymVib",
    console=False,          # コンソール窓を出さない
    upx=False,
)
coll = COLLECT(exe, a.binaries, a.datas, name="SymVib", upx=False)
