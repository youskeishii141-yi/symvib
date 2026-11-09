#!/usr/bin/env python3
"""qe_frozen_phonon.py — 力ベース Γ 点ブロック凍結フォノン (Quantum ESPRESSO 連携) の CLI。

因子群解析の対称性適応変位基底を使い、pw.x の SCF を数回行うだけで、
各既約表現の厳密な Γ 点振動数と真の固有ベクトルを求める。
物理の説明は symvib/frozen_phonon.py の冒頭を参照。

使い方 (2 フェーズ):
  [1] 変位構造の生成:
      python qe_frozen_phonon.py generate mycrystal.cif --config qe_config.json
        → qe_fp/ に scf_*.in (pw.x 入力) と manifest.json。これらを pw.x にかける
  [2] 力の回収と解析 (全 pw.x 計算の終了後):
      python qe_frozen_phonon.py analyze --dir qe_fp
        → 振動数表 + 固有ベクトルの XSF
  自己検証 (QE 不要):
  (出力は英語が既定。--lang ja で日本語。例: python qe_frozen_phonon.py --lang ja analyze)
      python qe_frozen_phonon.py validate [mycrystal.cif]

qe_config.json の例は examples/qe_config_example.json。

注意:
  ・必ず緩和済み (vc-relax) の構造を入力すること。
  ・力の精度のため conv_thr は十分小さく (1e-10 程度)。
  ・δ は 0.01〜0.02 Å 程度。小さすぎると力のノイズ、大きすぎると非調和が混入。
"""
import argparse
import json
import sys

from symvib.common import AnalysisError, setup_console
from symvib.i18n import LANGUAGES, set_language, t
from symvib.frozen_phonon import analyze, generate, validate


def main(argv=None):
    setup_console()
    p = argparse.ArgumentParser(description="Force-based Γ-point block frozen phonon (QE)")
    p.add_argument("--lang", choices=sorted(LANGUAGES), default="en",
                   help="output language: en (default) or ja")
    # --lang はサブコマンドの前でも後でも受け付ける
    lang = argparse.ArgumentParser(add_help=False)
    lang.add_argument("--lang", choices=sorted(LANGUAGES), default=argparse.SUPPRESS,
                      help="output language: en (default) or ja")
    sub = p.add_subparsers(dest="cmd")
    g = sub.add_parser("generate", parents=[lang],
                       help="generate pw.x inputs for displaced structures")
    g.add_argument("cif", nargs="?", help="CIF file (graphite 2H demo if omitted)")
    g.add_argument("--config", help="QE settings JSON (manifest only if omitted)")
    g.add_argument("--delta", type=float, default=0.01, help="displacement δ in Å (default 0.01)")
    g.add_argument("--labels", nargs="*", help="irrep labels to include (e.g. A1g Eg)")
    g.add_argument("--dir", default="qe_fp", help="output folder (default qe_fp)")
    g.add_argument("--include", default="raman", choices=["raman", "optical_all"])
    g.add_argument("--symprec", type=float, default=1e-3)
    a = sub.add_parser("analyze", parents=[lang], help="analyze frequencies from pw.x outputs")
    a.add_argument("--dir", default="qe_fp")
    v = sub.add_parser("validate", parents=[lang],
                       help="self-check of the method (no QE needed)")
    v.add_argument("cif", nargs="?")
    args = p.parse_args(argv)
    set_language(args.lang)
    try:
        if args.cmd == "generate":
            cfg = None
            if args.config:
                with open(args.config, encoding="utf-8-sig") as f:   # BOM 付きも読む
                    cfg = json.load(f)
            generate(cif=args.cif, config=cfg, delta=args.delta, labels=args.labels,
                     outdir=args.dir, include=args.include, symprec=args.symprec)
        elif args.cmd == "analyze":
            analyze(outdir=args.dir)
        elif args.cmd == "validate":
            ok = validate(cif=args.cif)
            sys.exit(0 if ok else 1)
        else:
            p.print_help()
    except (AnalysisError, OSError, ValueError) as e:      # ValueError: 設定 JSON の誤りなど
        sys.exit(t("cli.error", e=e))


if __name__ == "__main__":
    main()
