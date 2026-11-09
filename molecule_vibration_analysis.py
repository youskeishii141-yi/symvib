#!/usr/bin/env python3
"""molecule_vibration_analysis.py — 孤立分子の振動モード解析 (DFT 不要) の CLI。

使い方:
    python molecule_vibration_analysis.py                 # 引数なし → 水のデモ
    python molecule_vibration_analysis.py mymolecule.xyz  # xyz / CIF
    python molecule_vibration_analysis.py mol.xyz --outdir out --tol 0.05 --lang ja

出力は英語が既定 (--lang ja で日本語)。
解析本体は symvib/molecule.py。GUI は vibration_gui.py。
"""
import argparse
import os
import sys

from symvib.common import AnalysisError, setup_console
from symvib.i18n import LANGUAGES, set_language, t
from symvib.molecule import (analyze_molecule, load_molecule, water_demo,
                                         write_mode_files)


def run(cart, Z, sym=None, outdir="modes_molecule", tol=0.1, formats=("xsf", "vesta"),
        notes=None):
    """解析してレポートを表示し、XSF/VESTA を書き出す。AnalysisResult を返す。"""
    result = analyze_molecule(cart, Z, sym, tol=tol, notes=notes)
    print(result.report)
    if result.modes:
        paths = write_mode_files(result, outdir, formats=formats)
        print()
        print(t("cli.output", d=outdir, n=len(paths)))
        print(t("cli.vesta_hint"))
    return result


def main(argv=None):
    setup_console()
    p = argparse.ArgumentParser(description="Molecular vibration analysis (no DFT required)")
    p.add_argument("input", nargs="?", help="xyz or CIF file (water demo if omitted)")
    p.add_argument("--outdir", help="output folder for XSF/VESTA files (default: modes_<name>)")
    p.add_argument("--tol", type=float, default=0.1,
                   help="symmetry tolerance in Å (default 0.1); smaller deviations are "
                        "symmetrized away")
    p.add_argument("--format", choices=["both", "xsf", "vesta"], default="both",
                   help="mode files to write: both (default), xsf or vesta")
    p.add_argument("--lang", choices=sorted(LANGUAGES), default="en",
                   help="report language: en (default) or ja")
    args = p.parse_args(argv)
    set_language(args.lang)
    try:
        notes = []
        if args.input:
            cart, Z, sym = load_molecule(args.input, notes)
            outdir = args.outdir or f"modes_{os.path.splitext(os.path.basename(args.input))[0]}"
        else:
            print(t("cli.mol_demo"))
            print()
            cart, Z, sym = water_demo()
            outdir = args.outdir or "modes_water"
        formats = ("xsf", "vesta") if args.format == "both" else (args.format,)
        run(cart, Z, sym, outdir=outdir, tol=args.tol, formats=formats, notes=notes)
    except (AnalysisError, OSError) as e:
        sys.exit(t("cli.error", e=e))


if __name__ == "__main__":
    main()
