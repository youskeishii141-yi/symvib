"""symvib — 対称性 (群論) に基づく分子・結晶の振動モード解析。

主な入口:
    analyze_molecule(cart, numbers)        分子 (点群)
    analyze_crystal(lattice, frac, numbers) 結晶 (Γ 点の因子群)
"""
from .common import AnalysisError, AnalysisResult, Mode

__all__ = ["AnalysisError", "AnalysisResult", "Mode",
           "analyze_molecule", "load_molecule", "analyze_crystal", "load_crystal"]
__version__ = "2.2.1"


def __getattr__(name):
    # 重い依存 (pymatgen, spglib) を必要になるまで読み込まない
    if name in ("analyze_molecule", "load_molecule"):
        from . import molecule
        return getattr(molecule, name)
    if name in ("analyze_crystal", "load_crystal"):
        from . import crystal
        return getattr(crystal, name)
    raise AttributeError(name)
