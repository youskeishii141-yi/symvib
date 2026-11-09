import os
import sys
import warnings

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, ROOT)
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
for mod in ("spgrep", "spglib", "pymatgen"):
    warnings.filterwarnings("ignore", category=DeprecationWarning, module=mod)

# 既存のテストは日本語の注記・エラー文を確認するので、テスト中は日本語表示にする
from symvib.i18n import set_language  # noqa: E402

set_language("ja")
