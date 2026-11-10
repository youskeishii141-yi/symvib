"""Record every file of the built exe folder (dist/SymVib) with its size.

The list is written to <folder>/_internal/symvib_manifest.json. At startup, SymVib.exe compares the
files on disk with it, so that an incomplete extraction or a file removed by security software is
reported by name instead of as an obscure DLL error.

Usage:  python tools/bundle_manifest.py dist/SymVib
"""
import json
import os
import sys

NAME = "symvib_manifest.json"


def main(argv=None):
    argv = sys.argv[1:] if argv is None else argv
    root = argv[0] if argv else os.path.join("dist", "SymVib")
    internal = os.path.join(root, "_internal")
    if not os.path.isdir(internal):
        print(f"{internal} not found (build the exe first)")
        return 1
    files = {}
    for dirpath, _, names in os.walk(root):
        for name in names:
            path = os.path.join(dirpath, name)
            rel = os.path.relpath(path, root).replace(os.sep, "/")
            if rel != f"_internal/{NAME}":
                files[rel] = os.path.getsize(path)
    with open(os.path.join(internal, NAME), "w", encoding="utf-8") as f:
        json.dump({"files": dict(sorted(files.items()))}, f, indent=0)
    print(f"wrote {NAME}: {len(files)} files")
    return 0


if __name__ == "__main__":
    sys.exit(main())
