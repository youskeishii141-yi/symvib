"""Write THIRD_PARTY_LICENSES.txt for the distributed exe (PyInstaller bundle).

The exe bundles the Python runtime and the packages installed in the build environment
(numpy, scipy, pymatgen, spglib, spgrep, matplotlib, Tcl/Tk, ...). Their licenses (BSD, MIT,
PSF, ...) require the copyright notices to accompany binary redistributions, so this script
collects, for every installed distribution, its name, version, declared license and the full
text of its license / notice files, plus the license of the Python runtime itself.

Build- and test-only tools that are not bundled (pytest and its plugins, pip) are skipped;
anything else is included, so the list errs on the side of over-inclusion.

Usage:  python tools/third_party_notices.py [output path]   (default: THIRD_PARTY_LICENSES.txt)
"""
from __future__ import annotations

import os
import re
import sys
from importlib import metadata

# Installed in the build environment but not part of the exe.
NOT_BUNDLED = {"pip", "pytest", "iniconfig", "pluggy", "pygments", "packaging-test",
               "exceptiongroup", "tomli", "colorama-test"}
LICENSE_FILE = re.compile(r"(^|/)(licen[cs]e|copying|notice|authors|copyright)[^/]*$", re.I)
RULE = "=" * 78


def _license_texts(dist):
    """(file name, text) of the license-like files a distribution ships."""
    out, seen = [], set()
    for f in dist.files or []:
        name = str(f)
        if not LICENSE_FILE.search(name) or name.endswith((".py", ".pyc")):
            continue
        # .dist-info/ (or .dist-info/licenses/) files, plus top-level files of the package
        if ".dist-info" not in name and name.count("/") > 1:
            continue
        try:
            with open(f.locate(), encoding="utf-8", errors="replace") as fh:
                text = fh.read()
        except OSError:
            continue
        if text and text not in seen:
            seen.add(text)
            out.append((name, text.strip()))
    return out


def _declared_license(dist):
    md = dist.metadata
    expr = md.get("License-Expression")
    if expr:
        return expr
    lic = (md.get("License") or "").strip()
    if lic and "\n" not in lic and len(lic) < 120:
        return lic
    classifiers = [c.split("::")[-1].strip() for c in md.get_all("Classifier") or []
                   if c.startswith("License ::")]
    return ", ".join(classifiers) or "(see license text below)"


def main(argv=None):
    argv = sys.argv[1:] if argv is None else argv
    path = argv[0] if argv else "THIRD_PARTY_LICENSES.txt"
    dists = {}
    for d in metadata.distributions():
        name = d.metadata["Name"]
        if name and name.lower() not in NOT_BUNDLED:
            dists.setdefault(name.lower(), d)
    parts = [
        "Third-party software bundled with SymVib.exe",
        RULE,
        "SymVib.exe bundles the Python runtime and the open-source packages listed",
        "below. Each is distributed under its own license, reproduced here as required by",
        "those licenses. They are not covered by the license of SymVib itself",
        "(see LICENSE).",
        "",
        "Packages: " + ", ".join(f"{dists[k].metadata['Name']} {dists[k].version}"
                                 for k in sorted(dists)),
        "",
    ]
    py_license = os.path.join(sys.base_prefix, "LICENSE.txt")
    parts += [RULE, f"Python {sys.version.split()[0]} (runtime; includes the licenses of the "
              "libraries bundled with the Windows build of Python, e.g. Tcl/Tk, OpenSSL, zlib)", RULE]
    if os.path.exists(py_license):
        with open(py_license, encoding="utf-8", errors="replace") as f:
            parts.append(f.read().strip())
    else:
        parts.append("Python Software Foundation License Version 2: "
                     "https://docs.python.org/3/license.html")
        print(f"warning: {py_license} not found; wrote a link to the Python license instead",
              file=sys.stderr)
    missing = []
    for key in sorted(dists):
        d = dists[key]
        parts += ["", RULE, f"{d.metadata['Name']} {d.version}",
                  f"License: {_declared_license(d)}"]
        home = d.metadata.get("Home-page") or next(
            (u.split(",", 1)[-1].strip() for u in d.metadata.get_all("Project-URL") or []), "")
        if home:
            parts.append(f"Project: {home}")
        parts.append(RULE)
        texts = _license_texts(d)
        if not texts:
            declared = (d.metadata.get("License") or "").strip()
            if "\n" in declared or len(declared) >= 120:       # full text in the metadata
                texts = [("License (package metadata)", declared)]
            else:
                missing.append(d.metadata["Name"])
                parts.append("(No license file is shipped with this package; "
                             "see the project page.)")
        for name, text in texts:
            parts += [f"--- {name} ---", text, ""]
    with open(path, "w", encoding="utf-8", newline="\n") as f:
        f.write("\n".join(parts).rstrip() + "\n")
    print(f"wrote {path}: {len(dists)} packages"
          + (f"; no license file shipped by: {', '.join(missing)}" if missing else ""))
    return 0


if __name__ == "__main__":
    sys.exit(main())
