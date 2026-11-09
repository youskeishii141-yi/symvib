[English](validation_report.md) | 日本語

# 公開前検査の報告 (v2.1, 2026-11-08)

対称性に基づく振動モード解析 (`symvib/`) を、外部の基準 (Bilbao Crystallographic
Server、Crystallography Open Database、文献値、別実装 phonopy) と照合して検査した結果です。
検査に使ったスクリプトはすべてリポジトリに含めてあり、再実行できます。外部のデータ (Bilbao の指標表、
COD の CIF) は配布物に含めていません。`python tests/fetch_bilbao_tables.py` と
`python tests/fetch_cod.py` で各公式サイトから取得してください (COD は `tests/data/cod/index.json` に
記録した ID と同じ CIF を取得します)。取得していない場合、それらを使うテストはスキップされます。

| 検査 | 内容 | 実行 |
|---|---|---|
| 独立検査 (不変量) | 既約性・直交性・自由度・音響 3・重心不動・活性・ラマンテンソル成分 (射影演算子)・相関法の重複度・Bilbao 点群表のラベル | `tests/invariants.py` (`pytest tests/test_invariants.py`) |
| 2a 全空間群スイープ | spglib の 530 Hall 設定すべてに合成構造 (一般位置 2 軌道 + 特殊位置) を作って解析。標準設定に直した構造との一致、phonopy との一致 | `python tests/sweep_space_groups.py --phonopy` → [validation_space_groups.ja.md](validation_space_groups.ja.md) |
| 2b COD の実在 CIF | 230 空間群中 227 群の代表 CIF (ID は `tests/data/cod/index.json`) を解析。読み込みの注記、同じ Wyckoff 位置の合成構造との一致 | `python tests/check_cod.py` → [validation_cod.ja.md](validation_cod.ja.md) |
| 2c 文献値 | 33 物質 (7 晶系) の Raman / IR / 不活性モード数 | `pytest tests/test_literature.py` |
| 分子 | 45 点群すべての指標表の群論的検証、代表分子の文献値、任意回転・原子順入替・微小歪みへの不変性 | `pytest tests/test_tables.py tests/test_molecule.py tests/test_point_groups.py`, `python tests/sweep_point_groups.py` (45/45) |

## 1. 外部基準

- **Bilbao Crystallographic Server, Point Group Tables** (Bradley & Cracknell 1972; Altmann & Herzig 1994):
  `https://www.cryst.ehu.eus/cgi-bin/rep/programs/sam/point.py?sg=<空間群番号>`。32 点群の指標表
  (クラスの操作記号と軸方向、Mulliken ラベル、x, y, z と 2 次形式の帰属) を `tests/fetch_bilbao_tables.py`
  で取得し `tests/data/bilbao/point_tables.json` に保存して使う (配布物には含めない)。結晶の各既約表現について、標準設定の操作を
  軸方向で Bilbao のクラスに対応づけ、指標が一致する行のラベル・活性・テンソル成分と照合した。
  - Bilbao の **SAM 本体** (`nph-sam`: Wyckoff 位置ごとの既約分解) は Cloudflare Turnstile による
    人間の確認が必要で、自動取得できなかった。Wyckoff 位置ごとの重複度は、代わりに相関法
    (Bilbao の指標 × 各操作で不動な原子数) と phonopy で独立に確かめた。
- **COD** (`https://www.crystallography.net/cod/`, パブリックドメイン): 検索 API
  (`result?format=json&space_group_number=N`) で、秩序化 (化学式に小数係数なし)・元素数 ≤ 6・
  ファイル ≤ 120 kB (No. 48 のみ条件を緩めて 182 kB)・原子数が少ない順に各空間群 1 件を選んだ
  (`tests/fetch_cod.py`)。
  No. 89 (P422), 106 (P4₂bc), 168 (P6) はこの条件で見つからず「入手不可」。
- **phonopy 4.10** (`pip install phonopy`; テスト専用の追加依存): 対称化した乱数の力定数で動力学
  行列を作り、phonopy の irreps (Bilbao 由来の指標表による Mulliken ラベル) と既約分解を比較。
- 文献: 各物質の出典は `tests/test_literature.py` のコメントに記載。

## 2. 見つかった誤りと修正

| # | 重大度 | 箇所 | 内容 | 再現 | 修正 |
|---|---|---|---|---|---|
| 1 | 中 | `crystal.py` `_RAMAN_FORMS` | 4/m, -4, 4, 6/m, 6, -6 の E 系 (xz, yz) / (x²−y², xy) と 3, -3 の E のラマンテンソルの形を、4/mmm などと同じ 1 パラメータの形で表示していた。正しくは 2 パラメータ (3, -3 は 4 パラメータ) (Loudon, Adv. Phys. 13, 423 (1964))。成分の集合と活性は正しかった | `tests/test_literature.py` / `test_invariants.py::test_raman_tensor_forms_match_loudon` (P4/m の Eg など) | 表を廃止し、各モードの変位と行射影演算子 P = (1/h) Σ ⟨u\|D(g)\|u⟩ T(g) から形を厳密に計算 (`raman_tensor_forms`)。成分の集合も射影演算子の対角要素で判定 |
| 2 | 中 | `frozen_phonon.py` `generate` | 変位 δ を 3N 次元で規格化した基底ベクトルに掛けていたため、原子 1 個あたりの変位が δ/√(動く原子数) になり、大きなセルで力が数値ノイズに埋もれる | `tests/test_frozen_phonon.py::test_generate_and_analyze_with_synthetic_forces` (最大原子変位 = δ を確認) | 最大原子変位が δ になる振幅 (`displacement_amplitudes`) を使い、manifest に記録。解析側は記録した振幅で中心差分 |
| 3 | 低 | `molecule.py` `analyze_molecule` | 歪んだ入力の対称化を「操作数が変わらなくなるまで」で止めていたため、操作が厳密な群にならず残差 ~1e-7 Å が残り、縮退モードの部分空間が 1e-6 程度ずれていた | 全点群の一般分子に 0.02 Å の歪みを加えた検査 (`tests/invariants.check_molecule`) | 厳密な許容差で操作を取り直して平均する反復を収束 (1e-12 Å) まで行う |
| 4 | 低 | `cif_reader.py` `_symmetry_candidates` | 対称操作の記載がなく短い記号 ("P 21/c", "C 2/c") のとき、格子角 (unique axis) を見ずに常に b 軸設定で展開していた | `tests/test_cif_reader.py::test_monoclinic_unique_axis_follows_the_cell_angles` | 90° でない角から unique axis (α→a, β→b, γ→c) を決めて設定を絞る |
| 5 | 低 (潜在) | `crystal.py` `mulliken_labels` | 主軸を「最高次の真回転軸のうち最初のもの」にしており、D2d / S4 では C2' を主軸にとる可能性があった (spglib の操作順では実害なし) | `tests/test_invariants.py::test_mulliken_labels_do_not_depend_on_operation_order` | 回映軸 (−4, −6) の次数が高ければそれを主軸にとる |
| 6 | 中 | `crystal.py` `symmetry_adapted_blocks` | 格子定数が検出した対称性と厳密に整合しない入力 (例: AsI₃, COD 9008098, R-3, 菱面体角 51.683/51.683/51.685°) では直交座標の回転行列が直交でなくなり (8e-5)、網目モデルの力定数が対称操作と可換にならず、重複度 m>1 の「向きの推測」が同じ既約表現の別のコピーを混ぜていた (閉性の誤差 1e-5)。既約分解とラベルは正しく、縮退成分の表示だけが影響を受ける | `tests/test_cod_regression.py` (9008098) | 解析前に格子の計量を群で平均し、座標を軌道で平均して理想化 (`idealize_cell`)。変化が 1e-5 を超えるときはレポートに注記 |
| 7 | 中 | `common.py` `network_hessian` | 網目モデルが隣接セルを ±1 像しか数えておらず、ゆがんだ primitive セル (体心・面心・菱面体) ではカットオフ内の隣接原子を取りこぼして力定数が対称操作と可換でなくなっていた (COD 1010030, I-42m: 閉性の誤差 0.06)。影響は 6 と同じく m>1 の縮退成分の向きのみ | `tests/test_cod_regression.py` (1010030) | カットオフと格子面間隔から必要な像の範囲を決める |
| 8 | 低 | `molecule.py` `find_point_group_operations` | 対称操作の探索の基準原子を「同種原子の少ない順」だけで選ぶと重心のすぐそばの原子が選ばれ、0.02 Å 程度の歪みで回転の推定が大きく外れて対称性を見失う (C3 分子が C1 と判定された) | `tests/test_molecule.py::test_distorted_molecule_with_atom_near_centre_of_mass` | 外殻 (重心からの距離が最大の 40% 以上) の原子から基準を選ぶ |
| 9 | 中 (GUI) | `gui.py` `_build_scene` | VESTA の結合表にない元素の組だけの結晶 (AsI₃ など) では結合が 1 本もなく、空の Line3DCollection を `add_collection3d` すると matplotlib が例外を出して、解析結果が画面に表示されなかった (exe の起動確認で発見) | `tests/test_gui.py::test_crystal_without_table_bonds_builds_a_scene` | 結合がないときは線の集合を作らない |

誤りではないが、慣例として明記したもの:
- **-4m2 系 (No. 115–120)**: Bilbao の点群表は -42m の向き (2 回軸 ∥ [100]) で書かれている。-4m2 では
  2 回軸が [110] 方向なので、コードの B1 (2 回軸に対称: Mulliken の規則) は xy、B2 (z) は x²−y² になり、
  Bilbao の表の B1 = x²−y² と入れ替わって見える。群論的に B2 ∋ z は両設定で同じであり、慣例の違い。
  SAM での直接確認はできていない (上記)。
- **h-BN** (P6₃/mmc, B 2c, N 2d) の不活性モードは文献 (Geick 1966) で B1g と書かれるが、Bilbao の慣例
  (添字 1 = [120] 方向の 2 回軸に対称) では黒鉛と同じ B2g。指標を直接計算して確認した
  (`tests/test_literature.py::test_b2g_of_graphite_and_hbn_follow_bilbao_convention`)。
- **単斜晶**の Bilbao 点群表は z を 2 回軸にとるが、結晶の標準設定は unique axis b。コードは
  標準設定 (Ag ∋ xz, Bg ∋ xy, yz) で表す (文献の慣例と同じ)。
- 分子のラマンテンソル成分は指標表の標準の向き (z = 主軸) で表されており、入力座標の向きとは
  無関係であることをレポートに注記した。

## 3. 検査結果

### 3a. 230 空間群 (530 Hall 設定) の合成構造
`tests/sweep_space_groups.py --phonopy` (最終コードで実行、約 1 時間)。設定ごとの結果は
[validation_space_groups.ja.md](validation_space_groups.ja.md)。

- **独立検査と標準設定との一致: 530 / 530 設定、230 / 230 空間群が合格。** 各設定について、
  自由度 3N・音響 3・重心不動・各ブロックの D(g) での閉性・|χ|²/h ∈ {1, 2}・ラベルの整合 (直交性)・
  活性 (ベクトル表現・対称 2 階テンソル表現との内積)・ラマンテンソル成分 (射影演算子)・相関法に
  よる Wyckoff 位置ごとの重複度・Bilbao 点群表のラベルと関数 (x, y, z, 2 次形式) の帰属、および
  同じ構造を spglib で標準設定に直して解析した結果との一致を確かめた。
- Bilbao の表と**向きが異なる**設定 (ラベルの入れ替わりとして記録、誤りではない): -4m2 系の 6 設定
  (No. 115–120; 第 2 節参照)。
- **phonopy との比較: 451 設定で既約分解が一致、78 設定は phonopy 側がラベルを付けられず照合不可
  (P3₁2-型の "Rotation symbols not found"、A/B や 1 次元複素表現の区別がつかない群など)、
  1 設定 (No. 193 P6₃/mcm) で B1g/B1u ↔ B2g/B2u が入れ替わった。** No. 193 については、Bilbao の点群表
  (添字 1 = [100] 方向の 2 回軸に対称) と直接に指標を照合してコードのラベルが正しいことを確認して
  おり、phonopy 側の 2 回軸の向きの取り方の違いと判断した (No. 194 P6₃/mmc は phonopy が照合不可)。
- 直方晶 (P222, Pmmm 型など) の B1/B2/B3 は spglib の標準化 (a < b < c) で決まる軸の順に従う。入力の
  軸の取り方が違えば文献とラベルが入れ替わることがある (入力が標準設定でないときはレポートに
  「ラベルとテンソル成分は標準設定の軸で表す」と注記する)。

### 3b. COD の実在 CIF
`tests/check_cod.py` (最終コードで実行、約 2 時間半)。各 CIF の結果は [validation_cod.ja.md](validation_cod.ja.md)。

- **227 / 227 件 (230 空間群中 227 群) が合格、問題 0 件。** 各 CIF について、読み込み → 解析 →
  独立検査 (3a と同じ) を通し、さらに同じ空間群・同じ Wyckoff 位置に乱数の自由座標を置いた合成構造
  (元素は軌道ごとに別) を解析して既約分解の重複度が一致することを確かめた (CIF の読み込みと
  Wyckoff 位置の扱いが正しいことの検査)。原子数は 1〜1344 (単位胞)。
- 読み込み時の注記が付いたもの 73 件: 丸めた分数座標 (0.333 → 1/3 など) の読み替え 50 件、
  部分占有サイトを全占有として扱ったもの、同じ位置を共有するサイトの擬原子化。注記の内容は
  表の「読み込み時の注記」列と GUI / CLI のレポートに出る。
- 49 件では、CIF に書かれた空間群より**高い対称性**が検出された (表の「検出した空間群」が COD の
  番号と異なる行)。すべて symprec = 1e-6 Å でも同じ結果で、座標がちょうど高対称位置にある
  (理想化された構造、水素化物・金属間化合物の仮想的な秩序構造、部分占有を全占有として扱った
  無秩序構造 (例: Bi₂O₃ COD 1528292 Fm-3 → Pm-3, セル a/2) など) ことによる。宣言より**低い**
  対称性になった CIF はなかった (`tests/test_cod_regression.py` でも 135 件について点群の位数が
  下がらないことを確認)。
- 小さい 135 件 (単位胞 12 原子以下) は通常のテスト (`tests/test_cod_regression.py`) としても回る。

### 3c. 文献値 (33 物質)
`tests/test_literature.py` の全項目が一致。対象: ダイヤモンド, NaCl, CsCl, CaF₂, ZnS, SrTiO₃,
スピネル MgAl₂O₄, YAG, 黄鉄鉱, Cu₂O (立方); ルチル, アナターゼ, ジルコン, シーライト, 黄銅鉱, PbTiO₃
(正方); α-石英, コランダム, 方解石, LiNbO₃, Bi₂Se₃, イルメナイト (三方); 黒鉛, h-BN, 2H-MoS₂, ZnO
(六方); フォルステライト (Pbnm → Pnma 読み替え), KTiOPO₄ (直方); バデレアイト ZrO₂, β-Ga₂O₃,
透輝石 (単斜); カイヤナイト (三斜)。

入力の与え方の不変性: 非標準設定 (Pbnm / Pnma, R / H 軸, 原点 1 / 2, unique axis a/b/c)、
primitive / 超格子、原子順の入れ替え、座標の書き方 (1/3, 0.333, 0.3333) は
`tests/sweep_space_groups.py` (530 設定それぞれを標準設定と比較), `tests/test_crystal.py`,
`tests/test_cif_reader.py` で確認。

### 3d. 分子
- 45 点群の指標表 (実 + 複素): 正則表現上で射影演算子になること、IR/Raman の帰属と成分が
  ベクトル表現・対称 2 階テンソル表現の分解と一致することを確認 (`tests/test_tables.py`)。
- 文献値: 水, ホルムアルデヒド, NH₃, CH₄, CO₂, HCN, C₂H₂, エチレン, ベンゼン, SF₆, BF₃, XeF₄, エタン,
  アレン, フェロセン (D5h / D5d), B(OH)₃, B₁₂H₁₂²⁻, C₆₀, PCl₅, XeF₂, N₂O, アダマンタン, Ni(CO)₄,
  Mo(CO)₆, トロピリウム, IF₇, スピロペンタン, ギ酸, HNCS ほか (`tests/test_molecule.py`,
  `tests/test_point_groups.py`)。
- 全 45 点群の一般位置の分子で、任意回転・原子順入れ替え・0.02 Å の歪みを加えても点群・分解が
  変わらないこと、独立検査に合格することを確認 (検査スクリプトは本報告作成時に実行; 歪み後の
  構造が厳密に対称化されることは修正 3 で改善)。

## 4. リファクタリング
- `crystal.py`: ラマンテンソルの形の表 (`_RAMAN_FORMS`, 13 項目の lambda) を廃止し、射影演算子に
  よる計算 (`sym2_rep`, `raman_tensor_forms`, `format_tensor`) に置き換え。音響 / 光学の分離を
  `split_acoustic_optical` に切り出し。`tensor_components` を厳密な射影演算子に。重複 import の削除。
- `frozen_phonon.py`: 変位振幅を `displacement_amplitudes` に集約。
- i18n: 新しい表示文字列 (ラマンテンソルの注記、成分ごとの形、分子のテンソルの向き) を英日で追加。
- テスト: 独立検査を `tests/invariants.py` に集約し、各テストから再利用。

## 5. 残っている懸念・未検証事項
- Bilbao SAM (Wyckoff 位置ごとの分解) との直接照合は未実施 (Turnstile)。相関法・phonopy・文献値で
  代替した。-4m2 系の B1/B2 の慣例も SAM では未確認。
- `symprec = 1e-3 Å` の既定値: COD の 227 件では、宣言より高い対称性が検出された 49 件はすべて
  座標がちょうど高対称位置にあるもので (1e-6 Å でも同じ)、許容差によって対称性が上がった例も
  下がった例もなかった。実験構造の座標の精度 (1e-4 程度) に対して 1e-3 Å は妥当だが、精度の低い
  古い CIF では `--symprec 1e-2` が必要になりうる。逆に、わずかに歪んだ低対称相 (擬立方のペロブスカイト
  など) では 1e-3 Å でも高対称に丸められることがあるので、検出した空間群がレポートの先頭に出る。
- `cif_reader._coordinate`: 不確かさの付かない 3 桁の小数 (0.333, 0.667, 0.167 …) を分数 (1/3 …) に
  読み替える。精密化された自由パラメータがたまたま 0.333 と書かれていた場合は 0.0003 ほど動かす
  ことになる (symprec 以下では対称性に影響しない; 1/3 丸めが多い古い CIF で対称性が下がる方が実害が
  大きいと判断)。注記はレポートに出る。
- 網目モデル (m>1 の向きの推測) と縮退成分の選び方は物理的な予測ではなく表示規則であり、レポートに
  その旨を明示している。
- 分子のテンソル成分は点群の標準の向きで表す (入力座標系への変換は行わない)。
- GUI は数学的処理を持たないため、本検査では起動確認と結果表示の確認のみ行った。
