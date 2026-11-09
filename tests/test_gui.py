"""GUI の描画部品のうち Tk なしで確かめられる部分。"""
import numpy as np
import pytest

matplotlib = pytest.importorskip("matplotlib")


def _axes3d():
    from matplotlib.figure import Figure
    return Figure().add_subplot(projection="3d")


def test_crystal_without_table_bonds_builds_a_scene():
    """VESTA の結合表にない組だけの結晶 (AsI3 (COD 9008098) や fcc Al) では結合が 1 本もない。
    空の Line3DCollection を add_collection3d すると matplotlib が例外を出し、exe の起動確認で
    結果表示の直前に落ちていた。外部データに頼らないよう、fcc Al の単位胞で確かめる。"""
    from symvib.common import bond_pairs
    from symvib.gui import add_bond_collection
    frac = np.array([[0, 0, 0], [.5, .5, 0], [.5, 0, .5], [0, .5, .5]])
    assert bond_pairs(frac @ (4.05 * np.eye(3)), [13] * 4) == []
    ax = _axes3d()
    assert add_bond_collection(ax, np.zeros((0, 2, 3))) is None
    assert len(ax.collections) == 0


def test_bond_collection_is_added_when_bonds_exist():
    from symvib.gui import add_bond_collection
    ax = _axes3d()
    col = add_bond_collection(ax, np.array([[[0, 0, 0], [1, 0, 0]]], float))
    assert col is not None and list(ax.collections) == [col]


def test_error_log_returns_its_path(tmp_path, monkeypatch):
    from symvib import gui
    target = tmp_path / "error.log"
    monkeypatch.setattr(gui, "_ERROR_LOG", str(target))
    assert gui._write_error_log("test traceback") == str(target)
    assert "test traceback" in target.read_text(encoding="utf-8")
