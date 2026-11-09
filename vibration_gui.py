#!/usr/bin/env python3
"""振動モード解析 GUI の起動スクリプト。

    python vibration_gui.py [構造ファイル]

配布用の exe は build_exe.ps1 で作成する (README の「GUI 版」を参照)。
"""
from symvib.gui import main

if __name__ == "__main__":
    main()
