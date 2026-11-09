[English](validation_space_groups.md) | 日本語

# 230 空間群の網羅検査 (合成構造)

`tests/sweep_space_groups.py` の結果。各空間群について spglib の全 Hall 設定で合成構造 (一般位置 + 特殊位置) を解析し、独立検査 (`tests/invariants.py`: 既約性・活性・ラマンテンソル成分・相関法の重複度・Bilbao 点群表のラベル) と標準設定との一致を確かめた。

合格: 230 / 230 空間群, 530 / 530 設定

| No. | 設定 (Hall 設定数) | 既約分解 (標準設定・合成構造) | 結果 |
|---|---|---|---|
| 1 | P 1 (1) | 6A | OK; phonopy 一致 |
| 2 | P -1 (1) | 6Ag + 15Au | OK; phonopy 一致 |
| 3 | P 1 2 1 (3) | 9A + 12B | OK; phonopy 一致 |
| 4 | P 1 2_1 1 (3) | 6A + 6B | OK; phonopy 一致 |
| 5 | C 1 2 1 (9) | 9A + 12B | OK; phonopy 一致 |
| 6 | P 1 m 1 (3) | 12A' + 9A'' | OK; phonopy 一致 |
| 7 | P 1 c 1 (9) | 6A' + 6A'' | OK; phonopy 一致 |
| 8 | C 1 m 1 (9) | 12A' + 9A'' | OK; phonopy 一致 |
| 9 | C 1 c 1 (18) | 6A' + 6A'' | OK; phonopy 一致 |
| 10 | P 1 2/m 1 (3) | 6Ag + 9Au + 6Bg + 12Bu | OK; phonopy 一致 |
| 11 | P 1 2_1/m 1 (3) | 6Ag + 15Au + 6Bg + 15Bu | OK; phonopy 一致 |
| 12 | C 1 2/m 1 (9) | 6Ag + 9Au + 6Bg + 12Bu | OK; phonopy 一致 |
| 13 | P 1 2/c 1 (9) | 6Ag + 15Au + 6Bg + 15Bu | OK; phonopy 一致 |
| 14 | P 1 2_1/c 1 (9) | 6Ag + 15Au + 6Bg + 15Bu | OK; phonopy 一致 |
| 15 | C 1 2/c 1 (18) | 6Ag + 12Au + 6Bg + 12Bu | OK; phonopy 照合不可 |
| 16 | P 2 2 2 (1) | 6A + 9B1 + 9B2 + 9B3 | OK; phonopy 一致 |
| 17 | P 2 2 2_1 (3) | 9A + 12B1 + 12B2 + 9B3 | OK; phonopy 一致 |
| 18 | P 2_1 2_1 2 (3) | 9A + 9B1 + 12B2 + 12B3 | OK; phonopy 照合不可 |
| 19 | P 2_1 2_1 2_1 (1) | 6A + 6B1 + 6B2 + 6B3 | OK; phonopy 一致 |
| 20 | C 2 2 2_1 (3) | 9A + 12B1 + 12B2 + 9B3 | OK; phonopy 一致 |
| 21 | C 2 2 2 (3) | 6A + 9B1 + 9B2 + 9B3 | OK; phonopy 一致 |
| 22 | F 2 2 2 (1) | 4A + 6B1 + 7B2 + 7B3 | OK; phonopy 一致 |
| 23 | I 2 2 2 (1) | 6A + 9B1 + 9B2 + 9B3 | OK; phonopy 一致 |
| 24 | I 2_1 2_1 2_1 (1) | 7A + 8B1 + 8B2 + 7B3 | OK; phonopy 照合不可 |
| 25 | P m m 2 (3) | 9A1 + 6A2 + 9B1 + 9B2 | OK; phonopy 一致 |
| 26 | P m c 2_1 (6) | 12A1 + 9A2 + 9B1 + 12B2 | OK; phonopy 一致 |
| 27 | P c c 2 (3) | 9A1 + 9A2 + 12B1 + 12B2 | OK; phonopy 照合不可 |
| 28 | P m a 2 (6) | 9A1 + 9A2 + 12B1 + 12B2 | OK; phonopy 照合不可 |
| 29 | P c a 2_1 (6) | 6A1 + 6A2 + 6B1 + 6B2 | OK; phonopy 一致 |
| 30 | P n c 2 (6) | 9A1 + 9A2 + 12B1 + 12B2 | OK; phonopy 照合不可 |
| 31 | P m n 2_1 (6) | 12A1 + 9A2 + 9B1 + 12B2 | OK; phonopy 一致 |
| 32 | P b a 2 (3) | 9A1 + 9A2 + 12B1 + 12B2 | OK; phonopy 照合不可 |
| 33 | P n a 2_1 (6) | 6A1 + 6A2 + 6B1 + 6B2 | OK; phonopy 一致 |
| 34 | P n n 2 (3) | 9A1 + 9A2 + 12B1 + 12B2 | OK; phonopy 照合不可 |
| 35 | C m m 2 (3) | 9A1 + 6A2 + 9B1 + 9B2 | OK; phonopy 一致 |
| 36 | C m c 2_1 (6) | 12A1 + 9A2 + 9B1 + 12B2 | OK; phonopy 一致 |
| 37 | C c c 2 (3) | 9A1 + 9A2 + 12B1 + 12B2 | OK; phonopy 照合不可 |
| 38 | A m m 2 (6) | 9A1 + 6A2 + 9B1 + 9B2 | OK; phonopy 一致 |
| 39 | A e m 2 (6) | 9A1 + 9A2 + 12B1 + 12B2 | OK; phonopy 照合不可 |
| 40 | A m a 2 (6) | 8A1 + 8A2 + 10B1 + 10B2 | OK; phonopy 一致 |
| 41 | A e a 2 (6) | 9A1 + 9A2 + 12B1 + 12B2 | OK; phonopy 照合不可 |
| 42 | F m m 2 (3) | 9A1 + 6A2 + 9B1 + 9B2 | OK; phonopy 一致 |
| 43 | F d d 2 (3) | 8A1 + 8A2 + 10B1 + 10B2 | OK; phonopy 一致 |
| 44 | I m m 2 (3) | 9A1 + 6A2 + 9B1 + 9B2 | OK; phonopy 一致 |
| 45 | I b a 2 (3) | 9A1 + 9A2 + 12B1 + 12B2 | OK; phonopy 照合不可 |
| 46 | I m a 2 (6) | 8A1 + 8A2 + 10B1 + 10B2 | OK; phonopy 一致 |
| 47 | P 2/m 2/m 2/m (1) | 6Ag + 6Au + 6B1g + 9B1u + 6B2g + 9B2u + 6B3g + 9B3u | OK; phonopy 一致 |
| 48 | P 2/n 2/n 2/n (2) | 6Ag + 6Au + 9B1g + 9B1u + 9B2g + 9B2u + 9B3g + 9B3u | OK; phonopy 一致 |
| 49 | P 2/c 2/c 2/m (3) | 6Ag + 9Au + 6B1g + 9B1u + 6B2g + 12B2u + 6B3g + 12B3u | OK; phonopy 照合不可 |
| 50 | P 2/b 2/a 2/n (6) | 6Ag + 6Au + 9B1g + 9B1u + 9B2g + 9B2u + 9B3g + 9B3u | OK; phonopy 一致 |
| 51 | P 2_1/m 2/m 2/a (6) | 6Ag + 9Au + 6B1g + 12B1u + 6B2g + 9B2u + 6B3g + 12B3u | OK; phonopy 一致 |
| 52 | P 2/n 2_1/n 2/a (6) | 6Ag + 12Au + 6B1g + 12B1u + 6B2g + 12B2u + 6B3g + 12B3u | OK; phonopy 一致 |
| 53 | P 2/m 2/n 2_1/a (6) | 3Ag + 6Au + 3B1g + 9B1u + 3B2g + 9B2u + 3B3g + 6B3u | OK; phonopy 一致 |
| 54 | P 2_1/c 2/c 2/a (6) | 6Ag + 12Au + 6B1g + 12B1u + 6B2g + 12B2u + 6B3g + 12B3u | OK; phonopy 一致 |
| 55 | P 2_1/b 2_1/a 2/m (3) | 6Ag + 9Au + 6B1g + 9B1u + 6B2g + 12B2u + 6B3g + 12B3u | OK; phonopy 照合不可 |
| 56 | P 2_1/c 2_1/c 2/n (3) | 6Ag + 12Au + 6B1g + 12B1u + 6B2g + 12B2u + 6B3g + 12B3u | OK; phonopy 一致 |
| 57 | P 2/b 2_1/c 2_1/m (6) | 8Ag + 13Au + 8B1g + 13B1u + 7B2g + 14B2u + 7B3g + 14B3u | OK; phonopy 照合不可 |
| 58 | P 2_1/n 2_1/n 2/m (3) | 6Ag + 9Au + 6B1g + 9B1u + 6B2g + 12B2u + 6B3g + 12B3u | OK; phonopy 照合不可 |
| 59 | P 2_1/m 2_1/m 2/n (6) | 9Ag + 6Au + 6B1g + 9B1u + 9B2g + 9B2u + 9B3g + 9B3u | OK; phonopy 一致 |
| 60 | P 2_1/b 2/c 2_1/n (6) | 6Ag + 12Au + 6B1g + 12B1u + 6B2g + 12B2u + 6B3g + 12B3u | OK; phonopy 一致 |
| 61 | P 2_1/b 2_1/c 2_1/a (2) | 6Ag + 12Au + 6B1g + 12B1u + 6B2g + 12B2u + 6B3g + 12B3u | OK; phonopy 一致 |
| 62 | P 2_1/n 2_1/m 2_1/a (6) | 6Ag + 12Au + 6B1g + 12B1u + 6B2g + 12B2u + 6B3g + 12B3u | OK; phonopy 一致 |
| 63 | C 2/m 2/c 2_1/m (6) | 8Ag + 9Au + 7B1g + 12B1u + 7B2g + 12B2u + 8B3g + 9B3u | OK; phonopy 一致 |
| 64 | C 2/m 2/c 2_1/e (6) | 5Ag + 6Au + 4B1g + 9B1u + 4B2g + 9B2u + 5B3g + 6B3u | OK; phonopy 一致 |
| 65 | C 2/m 2/m 2/m (3) | 6Ag + 6Au + 6B1g + 9B1u + 6B2g + 9B2u + 6B3g + 9B3u | OK; phonopy 一致 |
| 66 | C 2/c 2/c 2/m (3) | 7Ag + 9Au + 7B1g + 9B1u + 8B2g + 12B2u + 8B3g + 12B3u | OK; phonopy 一致 |
| 67 | C 2/m 2/m 2/e (6) | 7Ag + 9Au + 8B1g + 12B1u + 8B2g + 12B2u + 7B3g + 9B3u | OK; phonopy 一致 |
| 68 | C 2/c 2/c 2/e (12) | 7Ag + 7Au + 10B1g + 10B1u + 10B2g + 10B2u + 9B3g + 9B3u | OK; phonopy 一致 |
| 69 | F 2/m 2/m 2/m (1) | 4Ag + 3Au + 3B1g + 6B1u + 4B2g + 6B2u + 4B3g + 6B3u | OK; phonopy 一致 |
| 70 | F 2/d 2/d 2/d (2) | 5Ag + 5Au + 7B1g + 7B1u + 8B2g + 8B2u + 7B3g + 7B3u | OK; phonopy 一致 |
| 71 | I 2/m 2/m 2/m (1) | 6Ag + 6Au + 6B1g + 9B1u + 6B2g + 9B2u + 6B3g + 9B3u | OK; phonopy 一致 |
| 72 | I 2/b 2/a 2/m (3) | 7Ag + 9Au + 7B1g + 9B1u + 8B2g + 12B2u + 8B3g + 12B3u | OK; phonopy 一致 |
| 73 | I 2/b 2/c 2/a (2) | 7Ag + 10Au + 8B1g + 11B1u + 8B2g + 11B2u + 7B3g + 10B3u | OK; phonopy 一致 |
| 74 | I 2/m 2/m 2/a (6) | 5Ag + 6Au + 4B1g + 9B1u + 4B2g + 9B2u + 5B3g + 6B3u | OK; phonopy 一致 |
| 75 | P 4 (1) | 9A + 6B + 9E | OK; phonopy 一致 |
| 76 | P 4_1 (1) | 6A + 6B + 6E | OK; phonopy 一致 |
| 77 | P 4_2 (1) | 9A + 9B + 12E | OK; phonopy 一致 |
| 78 | P 4_3 (1) | 6A + 6B + 6E | OK; phonopy 一致 |
| 79 | I 4 (1) | 9A + 7B + 10E | OK; phonopy 一致 |
| 80 | I 4_1 (1) | 9A + 9B + 12E | OK; phonopy 一致 |
| 81 | P -4 (1) | 6A + 9B + 9E | OK; phonopy 一致 |
| 82 | I -4 (1) | 7A + 9B + 10E | OK; phonopy 一致 |
| 83 | P 4/m (1) | 6Ag + 9Au + 6Bg + 6Bu + 6Eg + 9Eu | OK; phonopy 一致 |
| 84 | P 4_2/m (1) | 6Ag + 9Au + 6Bg + 9Bu + 6Eg + 12Eu | OK; phonopy 一致 |
| 85 | P 4/n (2) | 7Ag + 9Au + 8Bg + 6Bu + 9Eg + 9Eu | OK; phonopy 一致 |
| 86 | P 4_2/n (2) | 7Ag + 9Au + 9Bg + 7Bu + 10Eg + 10Eu | OK; phonopy 一致 |
| 87 | I 4/m (1) | 6Ag + 9Au + 6Bg + 7Bu + 6Eg + 10Eu | OK; phonopy 一致 |
| 88 | I 4_1/a (2) | 7Ag + 8Au + 8Bg + 7Bu + 9Eg + 9Eu | OK; phonopy 一致 |
| 89 | P 4 2 2 (1) | 6A1 + 9A2 + 6B1 + 6B2 + 15E | OK; phonopy 一致 |
| 90 | P 4 2_1 2 (1) | 7A1 + 9A2 + 8B1 + 6B2 + 18E | OK; phonopy 一致 |
| 91 | P 4_1 2 2 (1) | 9A1 + 12A2 + 9B1 + 12B2 + 21E | OK; phonopy 一致 |
| 92 | P 4_1 2_1 2 (1) | 9A1 + 12A2 + 12B1 + 9B2 + 21E | OK; phonopy 一致 |
| 93 | P 4_2 2 2 (1) | 6A1 + 9A2 + 6B1 + 9B2 + 18E | OK; phonopy 一致 |
| 94 | P 4_2 2_1 2 (1) | 7A1 + 9A2 + 9B1 + 7B2 + 20E | OK; phonopy 一致 |
| 95 | P 4_3 2 2 (1) | 9A1 + 12A2 + 9B1 + 12B2 + 21E | OK; phonopy 一致 |
| 96 | P 4_3 2_1 2 (1) | 9A1 + 12A2 + 12B1 + 9B2 + 21E | OK; phonopy 一致 |
| 97 | I 4 2 2 (1) | 6A1 + 9A2 + 6B1 + 7B2 + 16E | OK; phonopy 一致 |
| 98 | I 4_1 2 2 (1) | 8A1 + 10A2 + 10B1 + 8B2 + 21E | OK; phonopy 一致 |
| 99 | P 4 m m (1) | 9A1 + 6A2 + 6B1 + 6B2 + 15E | OK; phonopy 一致 |
| 100 | P 4 b m (1) | 9A1 + 8A2 + 6B1 + 7B2 + 18E | OK; phonopy 照合不可 |
| 101 | P 4_2 c m (1) | 9A1 + 7A2 + 7B1 + 9B2 + 20E | OK; phonopy 一致 |
| 102 | P 4_2 n m (1) | 9A1 + 7A2 + 7B1 + 9B2 + 20E | OK; phonopy 一致 |
| 103 | P 4 c c (1) | 9A1 + 9A2 + 7B1 + 7B2 + 20E | OK; phonopy 一致 |
| 104 | P 4 n c (1) | 9A1 + 9A2 + 7B1 + 7B2 + 20E | OK; phonopy 一致 |
| 105 | P 4_2 m c (1) | 9A1 + 6A2 + 9B1 + 6B2 + 18E | OK; phonopy 一致 |
| 106 | P 4_2 b c (1) | 9A1 + 9A2 + 9B1 + 9B2 + 24E | OK; phonopy 一致 |
| 107 | I 4 m m (1) | 9A1 + 6A2 + 7B1 + 6B2 + 16E | OK; phonopy 一致 |
| 108 | I 4 c m (1) | 9A1 + 8A2 + 6B1 + 7B2 + 18E | OK; phonopy 照合不可 |
| 109 | I 4_1 m d (1) | 9A1 + 6A2 + 9B1 + 6B2 + 18E | OK; phonopy 一致 |
| 110 | I 4_1 c d (1) | 7A1 + 7A2 + 7B1 + 7B2 + 16E | OK; phonopy 一致 |
| 111 | P -4 2 m (1) | 6A1 + 6A2 + 6B1 + 9B2 + 15E | OK; phonopy 一致 |
| 112 | P -4 2 c (1) | 7A1 + 7A2 + 9B1 + 9B2 + 20E | OK; phonopy 一致 |
| 113 | P -4 2_1 m (1) | 7A1 + 6A2 + 8B1 + 9B2 + 18E | OK; phonopy 照合不可 |
| 114 | P -4 2_1 c (1) | 7A1 + 7A2 + 9B1 + 9B2 + 20E | OK; phonopy 一致 |
| 115 | P -4 m 2 (1) | 6A1 + 6A2 + 6B1 + 9B2 + 15E | OK (Bilbao 点群表と向きが異なる設定: B1: tensor ['xy'] vs Bilbao ['xx', 'yy']; B2: tensor ['xx', 'yy'] vs Bilbao ['xy']); phonopy 一致 |
| 116 | P -4 c 2 (1) | 7A1 + 7A2 + 9B1 + 9B2 + 20E | OK (Bilbao 点群表と向きが異なる設定: B1: tensor ['xy'] vs Bilbao ['xx', 'yy']; B2: tensor ['xx', 'yy'] vs Bilbao ['xy']); phonopy 一致 |
| 117 | P -4 b 2 (1) | 6A1 + 7A2 + 8B1 + 9B2 + 18E | OK (Bilbao 点群表と向きが異なる設定: B1: tensor ['xy'] vs Bilbao ['xx', 'yy']; B2: tensor ['xx', 'yy'] vs Bilbao ['xy']); phonopy 照合不可 |
| 118 | P -4 n 2 (1) | 7A1 + 7A2 + 9B1 + 9B2 + 20E | OK (Bilbao 点群表と向きが異なる設定: B1: tensor ['xy'] vs Bilbao ['xx', 'yy']; B2: tensor ['xx', 'yy'] vs Bilbao ['xy']); phonopy 一致 |
| 119 | I -4 m 2 (1) | 7A1 + 6A2 + 6B1 + 9B2 + 16E | OK (Bilbao 点群表と向きが異なる設定: B1: tensor ['xy'] vs Bilbao ['xx', 'yy']; B2: tensor ['xx', 'yy'] vs Bilbao ['xy']); phonopy 一致 |
| 120 | I -4 c 2 (1) | 7A1 + 8A2 + 8B1 + 9B2 + 20E | OK (Bilbao 点群表と向きが異なる設定: B1: tensor ['xy'] vs Bilbao ['xx', 'yy']; B2: tensor ['xx', 'yy'] vs Bilbao ['xy']); phonopy 一致 |
| 121 | I -4 2 m (1) | 6A1 + 7A2 + 6B1 + 9B2 + 16E | OK; phonopy 一致 |
| 122 | I -4 2 d (1) | 7A1 + 7A2 + 8B1 + 8B2 + 18E | OK; phonopy 一致 |
| 123 | P 4/m 2/m 2/m (1) | 6A1g + 6A1u + 6A2g + 9A2u + 6B1g + 6B1u + 6B2g + 6B2u + 12Eg + 15Eu | OK; phonopy 一致 |
| 124 | P 4/m 2/c 2/c (1) | 6A1g + 9A1u + 6A2g + 9A2u + 6B1g + 7B1u + 6B2g + 7B2u + 12Eg + 20Eu | OK; phonopy 照合不可 |
| 125 | P 4/n 2/b 2/m (2) | 6A1g + 6A1u + 8A2g + 9A2u + 6B1g + 6B1u + 7B2g + 6B2u + 15Eg + 15Eu | OK; phonopy 一致 |
| 126 | P 4/n 2/n 2/c (2) | 6A1g + 6A1u + 9A2g + 9A2u + 6B1g + 6B1u + 7B2g + 7B2u + 16Eg + 16Eu | OK; phonopy 一致 |
| 127 | P 4/m 2_1/b m (1) | 6A1g + 8A1u + 6A2g + 9A2u + 6B1g + 7B1u + 6B2g + 6B2u + 12Eg + 18Eu | OK; phonopy 一致 |
| 128 | P 4/m 2_1/n c (1) | 6A1g + 9A1u + 6A2g + 9A2u + 6B1g + 7B1u + 6B2g + 7B2u + 12Eg + 20Eu | OK; phonopy 照合不可 |
| 129 | P 4/n 2_1/m m (2) | 7A1g + 6A1u + 6A2g + 9A2u + 8B1g + 6B1u + 6B2g + 6B2u + 15Eg + 15Eu | OK; phonopy 照合不可 |
| 130 | P 4/n 2_1/c c (2) | 8A1g + 9A1u + 8A2g + 9A2u + 8B1g + 7B1u + 8B2g + 7B2u + 20Eg + 20Eu | OK; phonopy 一致 |
| 131 | P 4_2/m 2/m 2/c (1) | 6A1g + 6A1u + 6A2g + 9A2u + 6B1g + 6B1u + 6B2g + 9B2u + 12Eg + 18Eu | OK; phonopy 一致 |
| 132 | P 4_2/m 2/c 2/m (1) | 6A1g + 7A1u + 6A2g + 9A2u + 6B1g + 9B1u + 6B2g + 7B2u + 12Eg + 20Eu | OK; phonopy 一致 |
| 133 | P 4_2/n 2/b 2/c (2) | 7A1g + 8A1u + 8A2g + 9A2u + 9B1g + 8B1u + 8B2g + 7B2u + 20Eg + 20Eu | OK; phonopy 一致 |
| 134 | P 4_2/n 2/n 2/m (2) | 6A1g + 6A1u + 7A2g + 9A2u + 6B1g + 6B1u + 9B2g + 7B2u + 16Eg + 16Eu | OK; phonopy 一致 |
| 135 | P 4_2/m 2_1/b 2/c (1) | 7A1g + 9A1u + 7A2g + 9A2u + 7B1g + 9B1u + 7B2g + 9B2u + 16Eg + 24Eu | OK; phonopy 照合不可 |
| 136 | P 4_2/m 2_1/n 2/m (1) | 6A1g + 7A1u + 6A2g + 9A2u + 6B1g + 9B1u + 6B2g + 7B2u + 12Eg + 20Eu | OK; phonopy 一致 |
| 137 | P 4_2/n 2_1/m 2/c (2) | 7A1g + 6A1u + 6A2g + 9A2u + 9B1g + 6B1u + 6B2g + 7B2u + 16Eg + 16Eu | OK; phonopy 照合不可 |
| 138 | P 4_2/n 2_1/c 2/m (2) | 8A1g + 8A1u + 7A2g + 9A2u + 8B1g + 8B1u + 9B2g + 7B2u + 20Eg + 20Eu | OK; phonopy 一致 |
| 139 | I 4/m 2/m 2/m (1) | 6A1g + 6A1u + 6A2g + 9A2u + 6B1g + 6B1u + 6B2g + 7B2u + 12Eg + 16Eu | OK; phonopy 一致 |
| 140 | I 4/m 2/c 2/m (1) | 7A1g + 8A1u + 7A2g + 9A2u + 6B1g + 7B1u + 6B2g + 6B2u + 14Eg + 18Eu | OK; phonopy 一致 |
| 141 | I 4_1/a 2/m 2/d (2) | 7A1g + 6A1u + 6A2g + 8A2u + 8B1g + 6B1u + 6B2g + 7B2u + 15Eg + 15Eu | OK; phonopy 一致 |
| 142 | I 4_1/a 2/c 2/d (2) | 6A1g + 6A1u + 7A2g + 7A2u + 7B1g + 7B1u + 6B2g + 6B2u + 14Eg + 14Eu | OK; phonopy 一致 |
| 143 | P 3 (1) | 9A + 9E | OK; phonopy 一致 |
| 144 | P 3_1 (1) | 6A + 6E | OK; phonopy 一致 |
| 145 | P 3_2 (1) | 6A + 6E | OK; phonopy 一致 |
| 146 | R 3 (2) | 9A + 9E | OK; phonopy 一致 |
| 147 | P -3 (1) | 6Ag + 11Au + 6Eg + 11Eu | OK; phonopy 一致 |
| 148 | R -3 (2) | 6Ag + 11Au + 6Eg + 11Eu | OK; phonopy 一致 |
| 149 | P 3 1 2 (1) | 7A1 + 10A2 + 17E | OK; phonopy 照合不可 |
| 150 | P 3 2 1 (1) | 7A1 + 10A2 + 17E | OK; phonopy 一致 |
| 151 | P 3_1 1 2 (1) | 9A1 + 12A2 + 21E | OK; phonopy 照合不可 |
| 152 | P 3_1 2 1 (1) | 9A1 + 12A2 + 21E | OK; phonopy 一致 |
| 153 | P 3_2 1 2 (1) | 9A1 + 12A2 + 21E | OK; phonopy 照合不可 |
| 154 | P 3_2 2 1 (1) | 9A1 + 12A2 + 21E | OK; phonopy 一致 |
| 155 | R 3 2 (2) | 7A1 + 9A2 + 16E | OK; phonopy 一致 |
| 156 | P 3 m 1 (1) | 10A1 + 7A2 + 17E | OK; phonopy 一致 |
| 157 | P 3 1 m (1) | 10A1 + 7A2 + 17E | OK; phonopy 一致 |
| 158 | P 3 c 1 (1) | 9A1 + 9A2 + 18E | OK; phonopy 一致 |
| 159 | P 3 1 c (1) | 9A1 + 9A2 + 18E | OK; phonopy 一致 |
| 160 | R 3 m (2) | 9A1 + 6A2 + 15E | OK; phonopy 一致 |
| 161 | R 3 c (2) | 5A1 + 5A2 + 10E | OK; phonopy 一致 |
| 162 | P -3 1 2/m (1) | 6A1g + 7A1u + 6A2g + 10A2u + 12Eg + 17Eu | OK; phonopy 一致 |
| 163 | P -3 1 2/c (1) | 7A1g + 11A1u + 7A2g + 11A2u + 14Eg + 22Eu | OK; phonopy 一致 |
| 164 | P -3 2/m 1 (1) | 6A1g + 7A1u + 6A2g + 10A2u + 12Eg + 17Eu | OK; phonopy 一致 |
| 165 | P -3 2/c 1 (1) | 7A1g + 11A1u + 7A2g + 11A2u + 14Eg + 22Eu | OK; phonopy 一致 |
| 166 | R -3 2/m (2) | 6A1g + 7A1u + 6A2g + 10A2u + 12Eg + 17Eu | OK; phonopy 一致 |
| 167 | R -3 2/c (2) | 4A1g + 5A1u + 4A2g + 5A2u + 8Eg + 10Eu | OK; phonopy 一致 |
| 168 | P 6 (1) | 9A + 8B + 10E1 + 7E2 | OK; phonopy 一致 |
| 169 | P 6_1 (1) | 6A + 6B + 6E1 + 6E2 | OK; phonopy 一致 |
| 170 | P 6_5 (1) | 6A + 6B + 6E1 + 6E2 | OK; phonopy 一致 |
| 171 | P 6_2 (1) | 9A + 12B + 12E1 + 9E2 | OK; phonopy 一致 |
| 172 | P 6_4 (1) | 9A + 12B + 12E1 + 9E2 | OK; phonopy 一致 |
| 173 | P 6_3 (1) | 9A + 9B + 9E1 + 9E2 | OK; phonopy 一致 |
| 174 | P -6 (1) | 8A' + 9A'' + 10E' + 7E'' | OK; phonopy 一致 |
| 175 | P 6/m (1) | 6Ag + 9Au + 6Bg + 8Bu + 6E1g + 10E1u + 6E2g + 7E2u | OK; phonopy 一致 |
| 176 | P 6_3/m (1) | 7Ag + 11Au + 7Bg + 11Bu + 7E1g + 11E1u + 7E2g + 11E2u | OK; phonopy 一致 |
| 177 | P 6 2 2 (1) | 6A1 + 9A2 + 7B1 + 7B2 + 16E1 + 13E2 | OK; phonopy 一致 |
| 178 | P 6_1 2 2 (1) | 9A1 + 12A2 + 12B1 + 9B2 + 21E1 + 21E2 | OK; phonopy 照合不可 |
| 179 | P 6_5 2 2 (1) | 9A1 + 12A2 + 12B1 + 9B2 + 21E1 + 21E2 | OK; phonopy 照合不可 |
| 180 | P 6_2 2 2 (1) | 7A1 + 9A2 + 10B1 + 10B2 + 20E1 + 16E2 | OK; phonopy 照合不可 |
| 181 | P 6_4 2 2 (1) | 7A1 + 9A2 + 10B1 + 10B2 + 20E1 + 16E2 | OK; phonopy 照合不可 |
| 182 | P 6_3 2 2 (1) | 8A1 + 10A2 + 10B1 + 8B2 + 18E1 + 18E2 | OK; phonopy 照合不可 |
| 183 | P 6 m m (1) | 9A1 + 6A2 + 7B1 + 7B2 + 16E1 + 13E2 | OK; phonopy 一致 |
| 184 | P 6 c c (1) | 9A1 + 9A2 + 8B1 + 8B2 + 20E1 + 14E2 | OK; phonopy 一致 |
| 185 | P 6_3 c m (1) | 10A1 + 7A2 + 7B1 + 10B2 + 17E1 + 17E2 | OK; phonopy 一致 |
| 186 | P 6_3 m c (1) | 10A1 + 7A2 + 10B1 + 7B2 + 17E1 + 17E2 | OK; phonopy 一致 |
| 187 | P -6 m 2 (1) | 7A1' + 6A1'' + 7A2' + 9A2'' + 16E' + 13E'' | OK; phonopy 一致 |
| 188 | P -6 c 2 (1) | 8A1' + 8A1'' + 10A2' + 10A2'' + 18E' + 18E'' | OK; phonopy 一致 |
| 189 | P -6 2 m (1) | 7A1' + 6A1'' + 7A2' + 9A2'' + 16E' + 13E'' | OK; phonopy 一致 |
| 190 | P -6 2 c (1) | 8A1' + 8A1'' + 10A2' + 10A2'' + 18E' + 18E'' | OK; phonopy 照合不可 |
| 191 | P 6/m 2/m 2/m (1) | 6A1g + 6A1u + 6A2g + 9A2u + 6B1g + 7B1u + 6B2g + 7B2u + 12E1g + 16E1u + 12E2g + 13E2u | OK; phonopy 照合不可 |
| 192 | P 6/m 2/c 2/c (1) | 7A1g + 9A1u + 7A2g + 9A2u + 6B1g + 8B1u + 6B2g + 8B2u + 14E1g + 20E1u + 12E2g + 14E2u | OK; phonopy 一致 |
| 193 | P 6_3/m 2/c 2/m (1) | 7A1g + 7A1u + 6A2g + 10A2u + 7B1g + 7B1u + 6B2g + 10B2u + 13E1g + 17E1u + 13E2g + 17E2u | OK; phonopy 不一致 (7A1g + 7A1u + 6A2g + 10A2u + 6B1g + 10B1u + 7B2g + 7B2u + 13E1g + 17E1u + 13E2g + 17E2u) |
| 194 | P 6_3/m 2/m 2/c (1) | 7A1g + 7A1u + 6A2g + 10A2u + 6B1g + 10B1u + 7B2g + 7B2u + 13E1g + 17E1u + 13E2g + 17E2u | OK; phonopy 照合不可 |
| 195 | P 2 3 (1) | 3A + 3E + 14T | OK; phonopy 一致 |
| 196 | F 2 3 (1) | 4A + 4E + 16T | OK; phonopy 一致 |
| 197 | I 2 3 (1) | 4A + 4E + 18T | OK; phonopy 一致 |
| 198 | P 2_1 3 (1) | 6A + 6E + 18T | OK; phonopy 一致 |
| 199 | I 2_1 3 (1) | 6A + 6E + 20T | OK; phonopy 一致 |
| 200 | P 2/m -3 (1) | 3Ag + 3Au + 3Eg + 3Eu + 9Tg + 14Tu | OK; phonopy 一致 |
| 201 | P 2/n -3 (2) | 4Ag + 4Au + 4Eg + 4Eu + 18Tg + 18Tu | OK; phonopy 一致 |
| 202 | F 2/m -3 (1) | 4Ag + 3Au + 4Eg + 3Eu + 11Tg + 14Tu | OK; phonopy 一致 |
| 203 | F 2/d -3 (2) | 4Ag + 4Au + 4Eg + 4Eu + 15Tg + 15Tu | OK; phonopy 一致 |
| 204 | I 2/m -3 (1) | 4Ag + 3Au + 4Eg + 3Eu + 11Tg + 16Tu | OK; phonopy 一致 |
| 205 | P 2_1/a -3 (1) | 4Ag + 6Au + 4Eg + 6Eu + 12Tg + 18Tu | OK; phonopy 一致 |
| 206 | I 2_1/a -3 (1) | 5Ag + 6Au + 5Eg + 6Eu + 17Tg + 20Tu | OK; phonopy 一致 |
| 207 | P 4 3 2 (1) | 3A1 + 3A2 + 6E + 13T1 + 10T2 | OK; phonopy 一致 |
| 208 | P 4_2 3 2 (1) | 4A1 + 4A2 + 8E + 18T1 + 18T2 | OK; phonopy 一致 |
| 209 | F 4 3 2 (1) | 4A1 + 3A2 + 7E + 14T1 + 11T2 | OK; phonopy 一致 |
| 210 | F 4_1 3 2 (1) | 4A1 + 4A2 + 8E + 15T1 + 15T2 | OK; phonopy 一致 |
| 211 | I 4 3 2 (1) | 4A1 + 3A2 + 7E + 15T1 + 12T2 | OK; phonopy 一致 |
| 212 | P 4_3 3 2 (1) | 5A1 + 5A2 + 10E + 15T1 + 15T2 | OK; phonopy 一致 |
| 213 | P 4_1 3 2 (1) | 5A1 + 5A2 + 10E + 15T1 + 15T2 | OK; phonopy 一致 |
| 214 | I 4_1 3 2 (1) | 5A1 + 5A2 + 10E + 17T1 + 17T2 | OK; phonopy 一致 |
| 215 | P -4 3 m (1) | 3A1 + 3A2 + 6E + 10T1 + 13T2 | OK; phonopy 一致 |
| 216 | F -4 3 m (1) | 4A1 + 3A2 + 7E + 11T1 + 14T2 | OK; phonopy 一致 |
| 217 | I -4 3 m (1) | 4A1 + 3A2 + 7E + 12T1 + 15T2 | OK; phonopy 一致 |
| 218 | P -4 3 n (1) | 4A1 + 4A2 + 8E + 18T1 + 18T2 | OK; phonopy 一致 |
| 219 | F -4 3 c (1) | 4A1 + 4A2 + 8E + 15T1 + 15T2 | OK; phonopy 一致 |
| 220 | I -4 3 d (1) | 5A1 + 5A2 + 10E + 17T1 + 17T2 | OK; phonopy 一致 |
| 221 | P 4/m -3 2/m (1) | 3A1g + 3A1u + 3A2g + 3A2u + 6Eg + 6Eu + 9T1g + 13T1u + 9T2g + 10T2u | OK; phonopy 照合不可 |
| 222 | P 4/n -3 2/n (2) | 4A1g + 4A1u + 3A2g + 3A2u + 7Eg + 7Eu + 15T1g + 15T1u + 12T2g + 12T2u | OK; phonopy 一致 |
| 223 | P 4_2/m -3 2/n (1) | 4A1g + 3A1u + 4A2g + 3A2u + 8Eg + 6Eu + 11T1g + 16T1u + 11T2g + 16T2u | OK; phonopy 一致 |
| 224 | P 4_2/n -3 2/m (2) | 4A1g + 3A1u + 3A2g + 4A2u + 7Eg + 7Eu + 12T1g + 15T1u + 15T2g + 12T2u | OK; phonopy 一致 |
| 225 | F 4/m -3 2/m (1) | 4A1g + 3A1u + 3A2g + 3A2u + 7Eg + 6Eu + 10T1g + 13T1u + 10T2g + 10T2u | OK; phonopy 一致 |
| 226 | F 4/m -3 2/c (1) | 4A1g + 3A1u + 4A2g + 3A2u + 8Eg + 6Eu + 11T1g + 13T1u + 11T2g + 13T2u | OK; phonopy 一致 |
| 227 | F 4_1/d -3 2/m (2) | 4A1g + 3A1u + 3A2g + 4A2u + 7Eg + 7Eu + 11T1g + 13T1u + 13T2g + 11T2u | OK; phonopy 一致 |
| 228 | F 4_1/d -3 2/c (2) | 3A1g + 3A1u + 3A2g + 3A2u + 6Eg + 6Eu + 9T1g + 9T1u + 9T2g + 9T2u | OK; phonopy 一致 |
| 229 | I 4/m -3 2/m (1) | 4A1g + 3A1u + 3A2g + 3A2u + 7Eg + 6Eu + 10T1g + 14T1u + 10T2g + 11T2u | OK; phonopy 一致 |
| 230 | I 4_1/a -3 2/d (1) | 5A1g + 5A1u + 5A2g + 5A2u + 10Eg + 10Eu + 17T1g + 17T1u + 17T2g + 17T2u | OK; phonopy 一致 |
