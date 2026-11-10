#!/usr/bin/env python3
"""SymVib GUI の起動スクリプト。

    python vibration_gui.py [構造ファイル]
    SymVib.exe --self-test <結果ファイル>    (配布版の動作確認。画面を出さずに解析まで試す)

配布用の exe は build_exe.ps1 で作成する (README の「GUI 版」を参照)。

配布版 (PyInstaller) では、起動に失敗したとき (zip の展開が途中で終わった、セキュリティソフトが
DLL を隔離した、など) に、ビルド時に記録したファイル一覧 (_internal/symvib_manifest.json) と
照合して、欠けている・壊れているファイルと対処法をダイアログで示す。
"""
import json
import os
import sys
import traceback

MANIFEST = "symvib_manifest.json"
ERROR_LOG = os.path.join(os.path.expanduser("~"), "symvib_error.log")
VC_REDIST = "https://learn.microsoft.com/cpp/windows/latest-supported-vc-redist"


def bundle_problems():
    """配布版で、ビルド時の一覧と比べて欠けている・大きさが違うファイルを返す (開発版では空)。"""
    internal = getattr(sys, "_MEIPASS", None)
    if not internal:
        return []
    root = os.path.dirname(os.path.abspath(sys.executable))
    try:
        with open(os.path.join(internal, MANIFEST), encoding="utf-8") as f:
            files = json.load(f)["files"]
    except (OSError, ValueError, KeyError):
        return [f"missing: _internal/{MANIFEST}"]
    problems = []
    for rel, size in files.items():
        try:
            actual = os.path.getsize(os.path.join(root, rel))
        except OSError:
            problems.append(f"missing: {rel}")
            continue
        if actual != size:
            problems.append(f"size differs: {rel} ({actual} bytes, expected {size})")
    return problems


def _write_log(text):
    import datetime
    try:
        with open(ERROR_LOG, "a", encoding="utf-8") as f:
            f.write(f"--- {datetime.datetime.now():%Y-%m-%d %H:%M:%S} (SymVib startup)\n{text}\n")
    except OSError:
        pass


def startup_failure_message(problems):
    """起動に失敗したときに示す文 (英語・日本語)。"""
    if problems:
        shown = problems[:8] + ([f"... and {len(problems) - 8} more"] if len(problems) > 8 else [])
        found = ("Missing or damaged program files / 欠けている・壊れているファイル:\n  "
                 + "\n  ".join(shown))
    else:
        found = ("All program files are present, so Windows or security software is probably\n"
                 "blocking them. / ファイルはそろっているので、Windows やセキュリティソフトが\n"
                 "読み込みを止めている可能性があります。")
    return f"""SymVib could not start because some of its program files could not be loaded.
SymVib を起動できませんでした (プログラムの一部のファイルを読み込めません)。

{found}

How to fix / 対処:
1. Delete this folder and extract SymVib.zip again with "Extract All". Wait until the extraction
   has finished, and use a short location such as C:\\SymVib.
   このフォルダを削除し、SymVib.zip を「すべて展開」で展開し直してください。展開が終わるまで
   待ってから起動し、C:\\SymVib のような短い場所を使ってください。
   If "Extract All" stops partway, extract with PowerShell (README: "Troubleshooting").
   「すべて展開」が途中で止まる場合は PowerShell で展開してください (README の「起動しないとき」)。
2. Open Windows Security > Virus & threat protection > Protection history and check whether
   files of SymVib were quarantined. On a PC managed by your organization, ask the administrator.
   Windows セキュリティ >「ウイルスと脅威の防止」>「保護の履歴」で SymVib のファイルが隔離されて
   いないか確認してください。組織が管理する PC では管理者に相談してください。
3. Install the Microsoft Visual C++ Redistributable (x64) and try again:
   Microsoft Visual C++ 再頒布可能パッケージ (x64) をインストールしてから再度お試しください:
   {VC_REDIST}

Details / 詳細: {ERROR_LOG}
See also the "Troubleshooting" section of README.md / README.ja.md の「起動しないとき」も参照。"""


def _show_error(text):
    if sys.platform == "win32":
        try:
            import ctypes
            ctypes.windll.user32.MessageBoxW(None, text, "SymVib", 0x10)   # MB_ICONERROR
            return
        except Exception:                                                  # noqa: BLE001
            pass
    print(text, file=sys.stderr)


def self_test(out_path):
    """画面を出さずに、配布版のファイル・依存ライブラリ・解析が動くかを確かめる (0 = 成功)。"""
    lines, ok = [], True
    problems = bundle_problems()
    if problems:
        ok = False
        lines += ["bundle check: FAILED"] + problems
    else:
        lines.append("bundle check: OK")
    try:
        import tkinter
        root = tkinter.Tk()                    # Tcl/Tk のデータ (_tcl_data, _tk_data) まで確かめる
        root.withdraw()
        root.destroy()
        lines.append("tcl/tk: OK")
        import matplotlib
        matplotlib.use("Agg")
        from symvib import gui                                             # noqa: F401
        from symvib.crystal import analyze_crystal, graphite_2h
        from symvib.molecule import analyze_molecule, water_demo
        lat, pos, num, sym = graphite_2h()
        lines.append("graphite: " + str(analyze_crystal(lat, pos, num, sym).decomposition))
        cart, z, sym = water_demo()
        lines.append("water: " + str(analyze_molecule(cart, z, sym).decomposition))
    except BaseException:                                                  # noqa: BLE001
        ok = False
        lines += ["analysis: FAILED", traceback.format_exc()]
    lines.append("RESULT: " + ("OK" if ok else "FAILED"))
    with open(out_path, "w", encoding="utf-8") as f:
        f.write("\n".join(lines) + "\n")
    return 0 if ok else 1


def _startup_failed():
    """起動 (読み込み・画面の初期化) の失敗を、欠けたファイルの診断つきで知らせて終了する。"""
    details = traceback.format_exc()
    problems = bundle_problems()
    _write_log(details + ("\n".join(problems) if problems else "(no missing files)"))
    _show_error(startup_failure_message(problems))
    sys.exit(1)


def main():
    argv = sys.argv[1:]
    if argv[:1] == ["--self-test"]:
        sys.exit(self_test(argv[1] if len(argv) > 1 else "symvib_self_test.txt"))
    frozen = getattr(sys, "frozen", False)
    try:
        from symvib.gui import main as gui_main
    except Exception:                                 # noqa: BLE001
        if not frozen:
            raise                                     # 開発環境ではそのまま例外を見せる
        _startup_failed()
    try:
        gui_main(argv)
    except Exception:                                 # noqa: BLE001  (Tcl/Tk の初期化失敗など)
        if not frozen:
            raise
        _startup_failed()


if __name__ == "__main__":
    main()
