"""配布版の起動時チェック (vibration_gui.py) を、PyInstaller の配置をまねたフォルダで確かめる。"""
import os
import subprocess
import sys

import vibration_gui as launcher

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))


def _fake_bundle(tmp_path, monkeypatch):
    app = tmp_path / "SymVib"
    (app / "_internal" / "numpy.libs").mkdir(parents=True)
    (app / "SymVib.exe").write_bytes(b"exe")
    (app / "_internal" / "numpy.libs" / "msvcp140-abc.dll").write_bytes(b"dll" * 10)
    subprocess.run([sys.executable, os.path.join(ROOT, "tools", "bundle_manifest.py"), str(app)],
                   check=True, capture_output=True)
    monkeypatch.setattr(sys, "_MEIPASS", str(app / "_internal"), raising=False)
    monkeypatch.setattr(sys, "executable", str(app / "SymVib.exe"))
    return app


def test_complete_bundle_has_no_problems(tmp_path, monkeypatch):
    _fake_bundle(tmp_path, monkeypatch)
    assert launcher.bundle_problems() == []


def test_missing_and_truncated_files_are_named(tmp_path, monkeypatch):
    app = _fake_bundle(tmp_path, monkeypatch)
    (app / "_internal" / "numpy.libs" / "msvcp140-abc.dll").unlink()
    (app / "SymVib.exe").write_bytes(b"")
    problems = launcher.bundle_problems()
    assert "missing: _internal/numpy.libs/msvcp140-abc.dll" in problems
    assert any(p.startswith("size differs: SymVib.exe") for p in problems)
    text = launcher.startup_failure_message(problems)
    assert "msvcp140-abc.dll" in text and "すべて展開" in text and "Extract All" in text


def test_development_checkout_is_not_checked():
    assert not hasattr(sys, "_MEIPASS")
    assert launcher.bundle_problems() == []
