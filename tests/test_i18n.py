"""表示言語 (英語 / 日本語) の切り替え。"""
import pytest

from structures import molecules
from symvib.common import AnalysisError
from symvib.crystal import analyze_crystal, graphite_2h
from symvib.i18n import MESSAGES, Msg, get_language, set_language, width
from symvib.molecule import analyze_molecule, water_demo


@pytest.fixture
def english():
    previous = get_language()
    set_language("en")
    yield
    set_language(previous)


def test_every_message_has_both_languages():
    for key, (en, ja) in MESSAGES.items():
        assert en and ja, key


def test_reports_switch_language_without_reanalysis(english):
    r = analyze_crystal(*graphite_2h())
    assert "Factor Group Analysis" in r.report and "Space group" in r.report
    assert not any(0x3040 <= ord(c) <= 0x9fff for c in r.report)      # 日本語を含まない
    set_language("ja")
    assert "空間群" in r.report                                          # 同じ結果を日本語で
    set_language("en")
    m = analyze_molecule(*water_demo())
    assert "Axis convention" in m.report and "Internal vibrations" in m.report


def test_notes_and_errors_follow_language(english):
    cart, Z, sym = molecules()["H2O"]
    note = analyze_molecule(cart + 0.02, Z, sym).notes[0]
    assert isinstance(note, Msg)
    assert "Axis convention" in str(note)
    set_language("ja")
    assert "軸の取り方" in str(note)
    set_language("en")
    with pytest.raises(AnalysisError, match="At least 2 atoms"):
        analyze_molecule([[0, 0, 0]], [1])


def test_table_columns_are_aligned_in_japanese():
    previous = get_language()
    set_language("ja")
    try:
        lines = analyze_crystal(*graphite_2h()).report.splitlines()
        rows = [ln for ln in lines if ln.startswith(" A2u") or ln.startswith(" E2g")]
        # 「方向」列の開始位置 (表示幅) がそろっている
        starts = {width(ln[:ln.index(word)]) for ln in rows
                  for word in ("対称性確定", "推測(網目)", "—") if word in ln}
        assert len(starts) == 1
    finally:
        set_language(previous)
