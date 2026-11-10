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


def _tk_root():
    tk = pytest.importorskip("tkinter")
    try:
        return tk.Tk()
    except tk.TclError:
        pytest.skip("no display")


def test_mode_list_rows_are_taller_than_the_font_at_high_dpi():
    """高 DPI (Tk の拡大率 2) でも、モード一覧の行の高さがフォントの高さを上回る
    (固定のピクセル値の行で、日本語・英語の文字の上下が切れていた)。"""
    from tkinter import ttk
    root = _tk_root()
    try:
        root.tk.call("tk", "scaling", 2.0)
        from symvib.gui import VibrationApp
        app = VibrationApp(root)
        rowheight = int(ttk.Style(root).lookup("Modes.Treeview", "rowheight"))
        assert rowheight >= app.text_font.metrics("linespace") + 4
        assert str(app.report.cget("font")) == str(app.text_font)   # 一覧とレポートは同じフォント
        app._fit_columns()
        heading = app.tree.heading("label", "text")
        assert int(app.tree.column("label", "width")) >= app.head_font.measure(heading)
    finally:
        root.destroy()


def test_scrollbar_hides_when_the_content_fits():
    from symvib.gui import AutoScrollbar
    root = _tk_root()
    try:
        sb = AutoScrollbar(root, orient="vertical")
        sb.grid(row=0, column=0)
        sb.set(0.0, 1.0)
        assert not sb.grid_info()
        sb.set(0.0, 0.5)
        assert sb.grid_info()
    finally:
        root.destroy()


def test_font_choice_is_applied_and_remembered(tmp_path, monkeypatch):
    """フォントの変更が一覧の行の高さに反映され、次回の起動時にも使われる。"""
    from tkinter import ttk
    from symvib import gui
    monkeypatch.setattr(gui, "_SETTINGS", str(tmp_path / "settings.json"))
    root = _tk_root()
    try:
        app = gui.VibrationApp(root)
        family = gui._font_families(root)[0]
        app.apply_font(family, 20)
        rowheight = int(ttk.Style(root).lookup("Modes.Treeview", "rowheight"))
        assert rowheight >= app.text_font.metrics("linespace") + 4
        assert gui._load_settings() == {"font_family": family, "font_size": 20}
        assert gui._report_font(root, gui._load_settings()) == (family, 20)
        # 存在しないフォント・範囲外の大きさは無視して既定に戻る
        assert gui._report_font(root, {"font_family": "No Such Font", "font_size": 99}) == \
            gui._report_font(root)
    finally:
        root.destroy()


def test_japanese_fonts_are_listed_first():
    from symvib import gui
    root = _tk_root()
    try:
        families = gui._font_families(root)
        japanese = [f for f in gui.JAPANESE_FONTS if f in families]
        assert families[:len(japanese)] == japanese
        assert not any(f.startswith("@") for f in families)
    finally:
        root.destroy()
