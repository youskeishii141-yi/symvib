#!/usr/bin/env python3
"""factor_group_analysis.py — 結晶の因子群解析 (Γ 点) の CLI。

CIF → Γ 点振動モードの既約分解 + Raman/IR 活性 + 偏光ラマン選択則 + XSF/VESTA 出力。

使い方:
    python factor_group_analysis.py                  # 引数なし → 黒鉛 2H のデモ
    python factor_group_analysis.py mycrystal.cif
    python factor_group_analysis.py mycrystal.cif --outdir out --symprec 1e-3 --lang ja

出力は英語が既定 (--lang ja で日本語)。
解析本体は symvib/crystal.py。GUI は vibration_gui.py。
"""
import argparse
import os
import sys

from symvib.common import AnalysisError, setup_console
from symvib.crystal import (analyze_crystal, graphite_2h, load_crystal,
                                        write_mode_files)
from symvib.i18n import LANGUAGES, set_language, t

load_from_cif = load_crystal          # 旧名 (互換用)


def run(lattice, positions, numbers, species=None, outdir="output_modes", symprec=1e-3,
        notes=None, formats=("xsf", "vesta")):
    """解析してレポートを表示し、光学モードの XSF/VESTA を書き出す。AnalysisResult を返す。"""
    result = analyze_crystal(lattice, positions, numbers, species, symprec=symprec, notes=notes)
    print(result.report)
    paths = write_mode_files(result, outdir, formats=formats)
    print()
    print(t("cli.output", d=outdir, n=len(paths)))
    print(t("cli.components"))
    print(t("cli.vesta_hint"))
    return result


def main(argv=None):
    setup_console()
    p = argparse.ArgumentParser(description="Factor group analysis (Γ-point lattice vibrations)")
    p.add_argument("cif", nargs="?", help="CIF file (graphite 2H demo if omitted)")
    p.add_argument("--outdir", help="output folder for XSF/VESTA files (default: modes_<name>)")
    p.add_argument("--symprec", type=float, default=1e-3,
                   help="symmetry tolerance for the space group in Å (default 1e-3)")
    p.add_argument("--format", choices=["both", "xsf", "vesta"], default="both",
                   help="mode files to write: both (default), xsf or vesta")
    p.add_argument("--lang", choices=sorted(LANGUAGES), default="en",
                   help="report language: en (default) or ja")
    args = p.parse_args(argv)
    set_language(args.lang)
    try:
        notes = []
        if args.cif:
            print(f"CIF: {args.cif}")
            print()
            structure = load_crystal(args.cif, notes)
            outdir = args.outdir or f"modes_{os.path.splitext(os.path.basename(args.cif))[0]}"
        else:
            print(t("cli.cry_demo"))
            print()
            structure = graphite_2h()
            outdir = args.outdir or "modes_graphite"
        formats = ("xsf", "vesta") if args.format == "both" else (args.format,)
        run(*structure, outdir=outdir, symprec=args.symprec, notes=notes, formats=formats)
    except (AnalysisError, OSError) as e:
        sys.exit(t("cli.error", e=e))


if __name__ == "__main__":
    main()
