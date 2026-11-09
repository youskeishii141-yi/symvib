[English](README_stage1_directions.md) | 日本語

# Stage 1: 振動モードの「方向」可視化 (DFT 不要)

`factor_group_analysis.py` は、CIF を入れるだけで各振動モードについて
「どの原子がどの方向に動くか」を XSF (VESTA 用) で出力します。
**方向の確かさには 2 種類**あり、出力で明示的に区別されます。

## 「対称性確定」と「網目推測」の違い

### 対称性確定 (厳密)
同じ既約表現が **1 回だけ** 出現する (m=1) モードは、変位の方向が
**群論だけで一意に決まります**。DFT もモデルも不要で、出力がそのまま
真の固有ベクトルの方向です。ダイヤモンドの T2g、多くの高対称結晶の
ラマンモードがこれに当たります。

### 網目推測 (m>1 のときのみ; 真の向きは Stage 2 で確定)
同じ既約表現が **複数回** 出現する (m>1) モードは、対称性だけでは
個々の向きが決まりません (2m 次元空間内の線形結合の自由度)。
そこで **1/d² 減衰の中心力バネ網目モデル** で力定数行列を作り、その
既約表現の部分空間内で対角化して「推測方向」を提示します。

**重要な注意:**
- これは**幾何学的な当て推量**であり、真の固有ベクトルとは限りません。
- バネ定数 (1/d²) は実際の結合の強さ (電子状態が決める) ではありません。
- 真の向きが必要なら **Stage 2 (qe_frozen_phonon.py)** で QE の力を使い、
  ブロック対角化して厳密に確定します。
- ただし「結合が短い=強い」という幾何学は反映するので、向きの**当たり**を
  つける目的には実用上有効です (例: 黒鉛の 2 つの E2g モードを、層間
  せん断モードと層内 G バンドに正しく分離できます)。

## 出力の読み方

モード表の「方向」列:
- `対称性確定` … 群論で一意 (厳密)
- `推測(網目)` … m>1 のためバネ網目モデルによる推測
- `—`          … 音響モード (並進; 方向の概念なし)

XSF ファイルのヘッダ (`# Direction:`) にも同じ区別が記録されます:
- `symmetry-determined (exact)`
- `NETWORK-MODEL GUESS (m>1; 1/d^2 springs) - true eigenvector requires DFT`

## 使い方
```
python factor_group_analysis.py mycrystal.cif
```
→ `modes_mycrystal/` に各モードの XSF が出力されます。
VESTA で開き、Properties → Vectors → Show vectors で矢印表示。

## Stage 2 への接続
m>1 モードの真の向きと振動数が欲しくなったら:
```
python qe_frozen_phonon.py generate mycrystal.cif --config qe_config.json
# (HPC で pw.x 実行)
python qe_frozen_phonon.py analyze --dir qe_fp
```
Stage 1 の対称性適応基底が、QE の力を取り込んで真の固有ベクトルに昇華します。
