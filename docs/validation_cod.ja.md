[English](validation_cod.md) | 日本語

# COD の実在 CIF による検査 (各空間群 1 件)

`tests/check_cod.py` の結果。出典: Crystallography Open Database (https://www.crystallography.net/cod/, パブリックドメイン)。CIF は配布物に含めず、`python tests/fetch_cod.py` で表の ID と同じものを取得して検査する (ID の一覧は `tests/data/cod/index.json`)。

| No. | COD ID | 物質 | 検出した空間群 | 原子数 | 既約分解 | 読み込み時の注記 | 結果 |
|---|---|---|---|---|---|---|---|
| 1 | [1010563](https://www.crystallography.net/cod/1010563.cif) | Al Cl3 Aluminium chloride | P-3m1 (No. 164) | 4 | A1g + 3A2u + Eg + 3Eu |  | OK |
| 2 | [2242624](https://www.crystallography.net/cod/2242624.cif) | Fe N4 Iron tetranitride | P-1 (No. 2) | 5 | 6Ag + 9Au |  | OK |
| 3 | [9006306](https://www.crystallography.net/cod/9006306.cif) | O2 Si  | P2 (No. 3) | 12 | 17A + 19B |  | OK |
| 4 | [9011672](https://www.crystallography.net/cod/9011672.cif) | Ag Cl  | P2_1/m (No. 11) | 4 | 4Ag + 2Au + 2Bg + 4Bu |  | OK |
| 5 | [1540868](https://www.crystallography.net/cod/1540868.cif) | Nb S2  | C2 (No. 5) | 6 | 4A + 5B |  | OK |
| 6 | [7707106](https://www.crystallography.net/cod/7707106.cif) | N5 Na2 sodium pentazolate | Pm (No. 6) | 7 | 14A' + 7A'' |  | OK |
| 7 | [9006305](https://www.crystallography.net/cod/9006305.cif) | O2 Si  | Pc (No. 7) | 12 | 18A' + 18A'' |  | OK |
| 8 | [1531762](https://www.crystallography.net/cod/1531762.cif) | Co O2  | Cm (No. 8) | 6 | 6A' + 3A'' |  | OK |
| 9 | [1537184](https://www.crystallography.net/cod/1537184.cif) | C K K (C N) | C2/m (No. 12) | 8 | 2Au + 4Bu |  | OK |
| 10 | [1521124](https://www.crystallography.net/cod/1521124.cif) | O2 V V O2 | P2/m (No. 10) | 6 | 4Ag + 4Au + 2Bg + 8Bu |  | OK |
| 11 | [1522310](https://www.crystallography.net/cod/1522310.cif) | Ni Ti  | P2_1/m (No. 11) | 4 | 4Ag + 2Au + 2Bg + 4Bu |  | OK |
| 12 | [1522025](https://www.crystallography.net/cod/1522025.cif) | Ni O2  | C2/m (No. 12) | 6 | 2Ag + 2Au + Bg + 4Bu |  | OK |
| 13 | [1537550](https://www.crystallography.net/cod/1537550.cif) | Cl4 Zr Zr Cl4 | P2/c (No. 13) | 10 | 7Ag + 7Au + 8Bg + 8Bu |  | OK |
| 14 | [1509319](https://www.crystallography.net/cod/1509319.cif) | Ag F2  | P2_1/c (No. 14) | 6 | 3Ag + 6Au + 3Bg + 6Bu |  | OK |
| 15 | [1011148](https://www.crystallography.net/cod/1011148.cif) | Cu O Tenorite | C2/c (No. 15) | 8 | Ag + 4Au + 2Bg + 5Bu | 1 coordinate value(s) written as rounded fractions (e.g. 0.333, 0.6667) were read as the exact fractions (1/3, 2/3, ...). | OK |
| 16 | [4002439](https://www.crystallography.net/cod/4002439.cif) | O2 Si Si O2 | P222 (No. 16) | 96 | 68A + 72B1 + 72B2 + 76B3 | 3 coordinate value(s) written as rounded fractions (e.g. 0.333, 0.6667) were read as the exact fractions (1/3, 2/3, ...). | OK |
| 17 | [1509714](https://www.crystallography.net/cod/1509714.cif) | Ag2 Se  | P222_1 (No. 17) | 12 | 8A + 10B1 + 10B2 + 8B3 |  | OK |
| 18 | [1533197](https://www.crystallography.net/cod/1533197.cif) | Bi Rb Rb Bi (S C N)4 | Pmmn (No. 59) | 4 | 2Ag + 2B1u + 2B2g + 2B2u + 2B3g + 2B3u |  | OK |
| 19 | [1520934](https://www.crystallography.net/cod/1520934.cif) | O2 Te Te O2 | P2_12_12_1 (No. 19) | 12 | 9A + 9B1 + 9B2 + 9B3 |  | OK |
| 20 | [2310733](https://www.crystallography.net/cod/2310733.cif) | K2 Zr K2 Zr F6 | Cmcm (No. 63) | 12 | 2Ag + Au + 2B1g + 4B1u + 4B2u + 2B3g + 3B3u | 2 coordinate value(s) written as rounded fractions (e.g. 0.333, 0.6667) were read as the exact fractions (1/3, 2/3, ...). | OK |
| 21 | [4344305](https://www.crystallography.net/cod/4344305.cif) | Ho Sb2  | C222 (No. 21) | 6 | A + 2B1 + 3B2 + 3B3 |  | OK |
| 22 | [1521320](https://www.crystallography.net/cod/1521320.cif) | Cu O2  | Fmmm (No. 69) | 12 | Ag + 2B1u + B2g + 2B2u + B3g + 2B3u |  | OK |
| 23 | [7209475](https://www.crystallography.net/cod/7209475.cif) | B P S4 B (P S4) | I222 (No. 23) | 12 | 3A + 5B1 + 5B2 + 5B3 |  | OK |
| 24 | [9006301](https://www.crystallography.net/cod/9006301.cif) | O2 Si  | I2_12_12_1 (No. 24) | 24 | 7A + 9B1 + 10B2 + 10B3 |  | OK |
| 25 | [1539745](https://www.crystallography.net/cod/1539745.cif) | In Sb  | Pmm2 (No. 25) | 2 | 2A1 + 2B1 + 2B2 | Sites sharing the same position were merged into one pseudo-atom: In2Sb2, Sb1In1 (mass averaged by occupancy). | OK |
| 26 | [1533907](https://www.crystallography.net/cod/1533907.cif) | O2 Se Se O2 | Pmc2_1 (No. 26) | 12 | 11A1 + 7A2 + 7B1 + 11B2 |  | OK |
| 27 | [1566654](https://www.crystallography.net/cod/1566654.cif) | D2 O Ice XIX | Pcc2 (No. 27) | 100 | 75A1 + 75A2 + 75B1 + 75B2 | Partially occupied sites were treated as fully occupied: D1a (0.43), D1b (0.45), D1c (0.52), D1d (0.6), D2a (0.57), D2b (0.55), D2c (0.48), D2d (0.4), D3a (0.5), D3b (0.5), D3c (0.5), D3d (0.5), D4a (0.34), D4b (0.78), D4c (0.71), D4d (0.32), D4e (0.22), D4f (0.66), D4g (0.68), D4h (0.29) | OK |
| 28 | [2310717](https://www.crystallography.net/cod/2310717.cif) | F Tl Tl F | Pma2 (No. 28) | 8 | 8A1 + 4A2 + 4B1 + 8B2 | 1 coordinate value(s) written as rounded fractions (e.g. 0.333, 0.6667) were read as the exact fractions (1/3, 2/3, ...). | OK |
| 29 | [1534893](https://www.crystallography.net/cod/1534893.cif) | Fe S2  | Pca2_1 (No. 29) | 12 | 9A1 + 9A2 + 9B1 + 9B2 |  | OK |
| 30 | [4001944](https://www.crystallography.net/cod/4001944.cif) | Ag I3 O8  | Pnc2 (No. 30) | 48 | 34A1 + 34A2 + 38B1 + 38B2 |  | OK |
| 31 | [2310355](https://www.crystallography.net/cod/2310355.cif) | Te2 W W Te2 | Pmn2_1 (No. 31) | 12 | 12A1 + 6A2 + 6B1 + 12B2 |  | OK |
| 32 | [7209412](https://www.crystallography.net/cod/7209412.cif) | Au2 O7 Se2 Au2 (Se O3)2 O | Pba2 (No. 32) | 22 | 16A1 + 16A2 + 17B1 + 17B2 |  | OK |
| 33 | [1528058](https://www.crystallography.net/cod/1528058.cif) | Fe P  | Pna2_1 (No. 33) | 8 | 6A1 + 6A2 + 6B1 + 6B2 |  | OK |
| 34 | [9012540](https://www.crystallography.net/cod/9012540.cif) | Fe Sb2 Seinajokite | Pnn2 (No. 34) | 6 | 4A1 + 4A2 + 5B1 + 5B2 |  | OK |
| 35 | [9008052](https://www.crystallography.net/cod/9008052.cif) | Cd Cl2 N2  | Cmmm (No. 65) | 10 | 2Ag + 2B1g + 3B1u + B2g + 3B2u + B3g + 3B3u |  | OK |
| 36 | [1538548](https://www.crystallography.net/cod/1538548.cif) | Br D D Br | Cmc2_1 (No. 36) | 8 | 4A1 + 2A2 + 2B1 + 4B2 |  | OK |
| 37 | [2003027](https://www.crystallography.net/cod/2003027.cif) | Li2 O5 Si2 Lithium Disilicate | Ccc2 (No. 37) | 36 | 13A1 + 13A2 + 14B1 + 14B2 |  | OK |
| 38 | [1520919](https://www.crystallography.net/cod/1520919.cif) | Nb S2  | Amm2 (No. 38) | 6 | 3A1 + A2 + 3B1 + 2B2 |  | OK |
| 39 | [4344349](https://www.crystallography.net/cod/4344349.cif) | S2 Ta3 Ta3 S2 | Aem2 (No. 39) | 40 | 16A1 + 14A2 + 16B1 + 14B2 |  | OK |
| 40 | [2310212](https://www.crystallography.net/cod/2310212.cif) | Pt2 U  | Cmcm (No. 63) | 12 | 2Ag + Au + 2B1g + 4B1u + 4B2u + 2B3g + 3B3u |  | OK |
| 41 | [9010972](https://www.crystallography.net/cod/9010972.cif) | O2 S  | Aea2 (No. 41) | 12 | 4A1 + 4A2 + 5B1 + 5B2 |  | OK |
| 42 | [1536584](https://www.crystallography.net/cod/1536584.cif) | S2 Ta (Sm S)1.19 Ta S2 | Fmm2 (No. 42) | 12 | 3A1 + A2 + 2B1 + 3B2 |  | OK |
| 43 | [1509692](https://www.crystallography.net/cod/1509692.cif) | Ag2 O3  | Fdd2 (No. 43) | 40 | 7A1 + 7A2 + 8B1 + 8B2 | 1 coordinate value(s) written as rounded fractions (e.g. 0.333, 0.6667) were read as the exact fractions (1/3, 2/3, ...). | OK |
| 44 | [2310656](https://www.crystallography.net/cod/2310656.cif) | Hg O  | Imm2 (No. 44) | 4 | 2A1 + 2B1 + 2B2 |  | OK |
| 45 | [1504391](https://www.crystallography.net/cod/1504391.cif) | C36 H30  | Iba2 (No. 45) | 264 | 99A1 + 99A2 + 99B1 + 99B2 |  | OK |
| 46 | [9006289](https://www.crystallography.net/cod/9006289.cif) | O2 Si  | Ima2 (No. 46) | 24 | 9A1 + 8A2 + 9B1 + 10B2 |  | OK |
| 47 | [1531802](https://www.crystallography.net/cod/1531802.cif) | Nb O  | Pmmm (No. 47) | 2 | 2B1u + 2B2u + 2B3u |  | OK |
| 48 | [7127387](https://www.crystallography.net/cod/7127387.cif) | C54 H32 O8 4,4',4'',4'''-(dibenzo[g,p]chrysene-2,7,10,15-tetrayl)tetrabenzoic acid | Pnnn (No. 48) | 188 | 70Ag + 70Au + 71B1g + 71B1u + 70B2g + 70B2u + 71B3g + 71B3u |  | OK |
| 49 | [1531068](https://www.crystallography.net/cod/1531068.cif) | O5 Ta2 Ta2 O5 | Pccm (No. 49) | 14 | 2Ag + 4Au + 4B1g + 6B1u + 3B2g + 10B2u + 3B3g + 10B3u |  | OK |
| 50 | [1000414](https://www.crystallography.net/cod/1000414.cif) | La2 Ni O4 Dilanthanum niccolate | Pban (No. 50) | 28 | 7Ag + 7Au + 13B1g + 13B1u + 9B2g + 9B2u + 13B3g + 13B3u |  | OK |
| 51 | [1510072](https://www.crystallography.net/cod/1510072.cif) | Au Cd  | Pmma (No. 51) | 4 | 2Ag + 2B1u + 2B2g + 2B2u + 2B3g + 2B3u |  | OK |
| 52 | [4002360](https://www.crystallography.net/cod/4002360.cif) | Bi3 Sr2  | Pnna (No. 52) | 20 | 6Ag + 6Au + 8B1g + 8B1u + 9B2g + 9B2u + 7B3g + 7B3u |  | OK |
| 53 | [1010271](https://www.crystallography.net/cod/1010271.cif) | H4 N4 Ammonium azide | Pmna (No. 53) | 16 | 4Ag + 5Au + 5B1g + 10B1u + 4B2g + 9B2u + 5B3g + 6B3u |  | OK |
| 54 | [7237644](https://www.crystallography.net/cod/7237644.cif) | H4 N6  | Pcca (No. 54) | 40 | 14Ag + 14Au + 15B1g + 15B1u + 15B2g + 15B2u + 16B3g + 16B3u |  | OK |
| 55 | [9003043](https://www.crystallography.net/cod/9003043.cif) | O4 Pb3 Minium | Pbam (No. 55) | 14 | 6Ag + 4Au + 6B1g + 4B1u + 3B2g + 8B2u + 3B3g + 8B3u |  | OK |
| 56 | [1534337](https://www.crystallography.net/cod/1534337.cif) | O3 Sb2 Valentinite | Pccn (No. 56) | 20 | 7Ag + 7Au + 7B1g + 7B1u + 8B2g + 8B2u + 8B3g + 8B3u |  | OK |
| 57 | [1520823](https://www.crystallography.net/cod/1520823.cif) | F Tl Tl F | Pbcm (No. 57) | 8 | 4Ag + 2Au + 4B1g + 2B1u + 2B2g + 4B2u + 2B3g + 4B3u |  | OK |
| 58 | [1008935](https://www.crystallography.net/cod/1008935.cif) | O2 Pt Platinum(IV) oxide - β | Pnnm (No. 58) | 6 | 2Ag + 2Au + 2B1g + 2B1u + B2g + 4B2u + B3g + 4B3u |  | OK |
| 59 | [1525500](https://www.crystallography.net/cod/1525500.cif) | In Sb  | Pmmn (No. 59) | 4 | 2Ag + 2B1u + 2B2g + 2B2u + 2B3g + 2B3u |  | OK |
| 60 | [1562532](https://www.crystallography.net/cod/1562532.cif) | O2 Ti 'Titanium Dioxide - II' | Pbcn (No. 60) | 12 | 4Ag + 4Au + 5B1g + 5B1u + 4B2g + 4B2u + 5B3g + 5B3u | Partially occupied sites were treated as fully occupied: Ti (0.5) | OK |
| 61 | [1509320](https://www.crystallography.net/cod/1509320.cif) | Ag F2  | Pbca (No. 61) | 12 | 3Ag + 6Au + 3B1g + 6B1u + 3B2g + 6B2u + 3B3g + 6B3u |  | OK |
| 62 | [1010473](https://www.crystallography.net/cod/1010473.cif) | Mn P Manganese phosphide (1/1) | Pnma (No. 62) | 8 | 4Ag + 2Au + 2B1g + 4B1u + 4B2g + 2B2u + 2B3g + 4B3u |  | OK |
| 63 | [1524656](https://www.crystallography.net/cod/1524656.cif) | Cu3 Ti  | Cmcm (No. 63) | 4 | Ag + B1g + B1u + B2u + B3g + B3u | Sites sharing the same position were merged into one pseudo-atom: Ti1Cu1 (mass averaged by occupancy). | OK |
| 64 | [4030189](https://www.crystallography.net/cod/4030189.cif) | D I  | Cmce (No. 64) | 16 | 4Ag + 2Au + 2B1g + 4B1u + 2B2g + 4B2u + 4B3g + 2B3u | Partially occupied sites were treated as fully occupied: D (0.5), I (0.5) | OK |
| 65 | [2310206](https://www.crystallography.net/cod/2310206.cif) | Ir V  | Cmmm (No. 65) | 8 | 2Ag + 2B1g + 2B1u + B2g + 2B2u + B3g + 2B3u |  | OK |
| 66 | [1008186](https://www.crystallography.net/cod/1008186.cif) | I3 Th THORIUM(III) IODIDE - $-beta | Cccm (No. 66) | 64 | 11Ag + 10Au + 12B1g + 11B1u + 11B2g + 15B2u + 11B3g + 15B3u |  | OK |
| 67 | [9012694](https://www.crystallography.net/cod/9012694.cif) | O Pb Litharge | Cmme (No. 67) | 8 | Ag + B1g + 2B1u + 2B2g + 2B2u + 2B3g + 2B3u |  | OK |
| 68 | [1536196](https://www.crystallography.net/cod/1536196.cif) | Pd4 Sn16  | Ccce (No. 68) | 20 | 3Ag + 3Au + 4B1g + 4B1u + 4B2g + 4B2u + 4B3g + 4B3u |  | OK |
| 69 | [1010386](https://www.crystallography.net/cod/1010386.cif) | Br H Hydrogen bromide | Fmmm (No. 69) | 8 | 2B1u + 2B2u + 2B3u |  | OK |
| 70 | [1010453](https://www.crystallography.net/cod/1010453.cif) | Si2 Ti Zhiqinite | Fddd (No. 70) | 24 | Ag + Au + 3B1g + 3B1u + 2B2g + 2B2u + 3B3g + 3B3u | 1 coordinate value(s) written as rounded fractions (e.g. 0.333, 0.6667) were read as the exact fractions (1/3, 2/3, ...). | OK |
| 71 | [1522847](https://www.crystallography.net/cod/1522847.cif) | Mo Pt2  | Immm (No. 71) | 6 | Ag + 2B1u + B2g + 2B2u + B3g + 2B3u |  | OK |
| 72 | [1010272](https://www.crystallography.net/cod/1010272.cif) | S2 Si Silicon catena-sulfide | Ibam (No. 72) | 12 | 2Ag + Au + 3B1g + 2B1u + 2B2g + 3B2u + 2B3g + 3B3u |  | OK |
| 73 | [1535987](https://www.crystallography.net/cod/1535987.cif) | Li4 N3 Ta Li4 (Ta N3) | Ibca (No. 73) | 64 | 10Ag + 10Au + 12B1g + 12B1u + 13B2g + 13B2u + 13B3g + 13B3u |  | OK |
| 74 | [1509681](https://www.crystallography.net/cod/1509681.cif) | Ag2 Ba  | Imma (No. 74) | 12 | 3Ag + Au + B1g + 3B1u + 2B2g + 3B2u + 3B3g + 2B3u | 1 coordinate value(s) written as rounded fractions (e.g. 0.333, 0.6667) were read as the exact fractions (1/3, 2/3, ...). | OK |
| 75 | [1521986](https://www.crystallography.net/cod/1521986.cif) | Ce Te2  | P4 (No. 75) | 72 | 54A + 46B + 58E | Partially occupied sites were treated as fully occupied: Te1 (0.866), Te10 (0.177), Te2 (0.134), Te3 (0.866), Te4 (0.134), Te5 (0.446), Te6 (0.377), Te7 (0.177), Te8 (0.446), Te9 (0.377) | OK |
| 76 | [1521678](https://www.crystallography.net/cod/1521678.cif) | B5 H9  | P4_1 (No. 76) | 112 | 84A + 84B + 84E |  | OK |
| 77 | [1528561](https://www.crystallography.net/cod/1528561.cif) | H2 S  | P4_2/m (No. 84) | 4 | 2Ag + Au + 2Bg + Bu + Eg + 2Eu |  | OK |
| 78 | [4002638](https://www.crystallography.net/cod/4002638.cif) | C2 Ce Rh Ce Rh C2 | P4_3 (No. 78) | 16 | 12A + 12B + 12E |  | OK |
| 79 | [2310939](https://www.crystallography.net/cod/2310939.cif) | H2 Th Th H2 | I4/mmm (No. 139) | 6 | 2A2u + B1g + Eg + 2Eu |  | OK |
| 80 | [1537320](https://www.crystallography.net/cod/1537320.cif) | As Nb Nb As | I4_1md (No. 109) | 8 | 2A1 + 2B1 + 4E | 1 coordinate value(s) written as rounded fractions (e.g. 0.333, 0.6667) were read as the exact fractions (1/3, 2/3, ...). | OK |
| 81 | [1566655](https://www.crystallography.net/cod/1566655.cif) | D2 O ice XIX | P-4 (No. 81) | 100 | 73A + 75B + 76E | Partially occupied sites were treated as fully occupied: D11 (0.5), D21 (0.5), D31 (0.253), D32 (0.747), D41 (0.5), D42 (0.723), D43 (0.5), D44 (0.277), D51 (0.343), D52 (0.299), D53 (0.747), D54 (0.614), D61 (0.5), D62 (0.785), D63 (0.5), D64 (0.215), D71 (0.657), D72 (0.386), D73 (0.253), D74 (0.701) | OK |
| 82 | [1521080](https://www.crystallography.net/cod/1521080.cif) | Ge Se2  | I-4 (No. 82) | 12 | 3A + 5B + 5E |  | OK |
| 83 | [1011329](https://www.crystallography.net/cod/1011329.cif) | O Pd Palladium oxide | P4/mmm (No. 123) | 2 | 2A2u + 2Eu |  | OK |
| 84 | [1011155](https://www.crystallography.net/cod/1011155.cif) | Pd S Vysotskite | P4_2/m (No. 84) | 16 | 5Ag + 6Au + 6Bg + 5Bu + 5Eg + 8Eu |  | OK |
| 85 | [1509387](https://www.crystallography.net/cod/1509387.cif) | Ag I  | P4/nmm (No. 129) | 4 | A1g + 2A2u + B1g + 2Eg + 2Eu |  | OK |
| 86 | [1540981](https://www.crystallography.net/cod/1540981.cif) | D Zr Zr D | P4_2/mmc (No. 131) | 8 | 2A2u + B1g + B2u + Eg + 3Eu |  | OK |
| 87 | [2310961](https://www.crystallography.net/cod/2310961.cif) | H2 Hf Hf H2 | I4/mmm (No. 139) | 6 | 2A2u + B1g + Eg + 2Eu |  | OK |
| 88 | [8103872](https://www.crystallography.net/cod/8103872.cif) | Cs S Cs S O3 F | I4_1/amd (No. 141) | 8 | 2A2u + 2B1g + 2Eg + 2Eu |  | OK |
| 89 | — | (入手不可) | | | | | |
| 90 | [1010449](https://www.crystallography.net/cod/1010449.cif) | Br6 K2 Sn Potassium hexabromostannate | P42_12 (No. 90) | 36 | 6A1 + 8A2 + 5B1 + 3B2 + 16E |  | OK |
| 91 | [1527272](https://www.crystallography.net/cod/1527272.cif) | Li2 O4 Te Li2 Te O4 | P4_122 (No. 91) | 28 | 9A1 + 12A2 + 10B1 + 11B2 + 21E |  | OK |
| 92 | [1010938](https://www.crystallography.net/cod/1010938.cif) | O2 Si Cristobalite low | P4_12_12 (No. 92) | 12 | 4A1 + 5A2 + 5B1 + 4B2 + 9E |  | OK |
| 93 | [1556349](https://www.crystallography.net/cod/1556349.cif) | C26 Cu O5  | P4_222 (No. 93) | 264 | 97A1 + 101A2 + 97B1 + 101B2 + 198E | Partially occupied sites were treated as fully occupied: C10A (0.5), C11A (0.5), C12A (0.5), C13A (0.5), C14A (0.5), C15A (0.5), C16A (0.5), C17A (0.5), C18A (0.5), C19N (0.5), C1A (0.5), C20N (0.5), C21N (0.5), C22N (0.5), C23B (0.5), C24B (0.5), C25E (0.5), C26E (0.5), C2A (0.5), C3A (0.5), C4A (0.5), C5A (0.5), C6A (0.5), C7A (0.5), C8A (0.5), C9A (0.5), Cu1 (0.5), Cu1B (0.5), O1B (0.5), O2B (0.5), O3B (0.5), O4B (0.5), O5B (0.5), O5H (0.5), O6E (0.5) | OK |
| 94 | [1511227](https://www.crystallography.net/cod/1511227.cif) | B Li3 N2 Li3 (B N2) | P4_2/mnm (No. 136) | 12 | A1g + A1u + A2g + 4A2u + 2B1g + 3B1u + 2B2g + 3Eg + 8Eu |  | OK |
| 95 | [1008093](https://www.crystallography.net/cod/1008093.cif) | Mn2 O4 Ti Dimanganese titanium oxide | P4_322 (No. 95) | 28 | 9A1 + 12A2 + 10B1 + 11B2 + 21E |  | OK |
| 96 | [2310609](https://www.crystallography.net/cod/2310609.cif) | Cu2 S  | P4_32_12 (No. 96) | 12 | 4A1 + 5A2 + 5B1 + 4B2 + 9E |  | OK |
| 97 | [1539547](https://www.crystallography.net/cod/1539547.cif) | Ca2 O20 Si8 Th Th Ca2 (Si8 O20) | I422 (No. 97) | 62 | 10A1 + 13A2 + 11B1 + 11B2 + 24E |  | OK |
| 98 | [9008973](https://www.crystallography.net/cod/9008973.cif) | Nb P  | I4_1/amd (No. 141) | 8 | 2A2u + 2B1g + 2Eg + 2Eu |  | OK |
| 99 | [1525025](https://www.crystallography.net/cod/1525025.cif) | Cu3 Pd  | P4mm (No. 99) | 28 | 21A1 + 7B1 + 28E |  | OK |
| 100 | [1526958](https://www.crystallography.net/cod/1526958.cif) | H6 O10 V3 (H3 O)2 V3 O8 | P4bm (No. 100) | 26 | 12A1 + 8A2 + 6B1 + 10B2 + 21E |  | OK |
| 101 | [1528438](https://www.crystallography.net/cod/1528438.cif) | C7 H12 Cu2 N4 S3 ((C H3)4 N) (Cu2 (N C S)3) | P4_2cm (No. 101) | 112 | 47A1 + 36A2 + 36B1 + 47B2 + 85E |  | OK |
| 102 | [4031496](https://www.crystallography.net/cod/4031496.cif) | N2 S4 S4 N2 | P4_2nm (No. 102) | 24 | 10A1 + 8A2 + 8B1 + 10B2 + 18E |  | OK |
| 103 | [1534271](https://www.crystallography.net/cod/1534271.cif) | Ta Te4  | P4cc (No. 103) | 10 | 4A1 + 4A2 + 3B1 + 3B2 + 8E |  | OK |
| 104 | [1534196](https://www.crystallography.net/cod/1534196.cif) | Hg I6 Tl4 Tl4 (Hg I6) | P4nc (No. 104) | 22 | 9A1 + 9A2 + 6B1 + 6B2 + 18E |  | OK |
| 105 | [7222662](https://www.crystallography.net/cod/7222662.cif) | Ba Ge2 P2  | P4_2mc (No. 105) | 20 | 9A1 + 5A2 + 9B1 + 5B2 + 16E |  | OK |
| 106 | — | (入手不可) | | | | | |
| 107 | [1528443](https://www.crystallography.net/cod/1528443.cif) | Ge P  | I4mm (No. 107) | 4 | 2A1 + 2E |  | OK |
| 108 | [7221364](https://www.crystallography.net/cod/7221364.cif) | Si3 Sr5 Sr5 Si3 | I4cm (No. 108) | 32 | 8A1 + 5A2 + 3B1 + 6B2 + 13E |  | OK |
| 109 | [1521083](https://www.crystallography.net/cod/1521083.cif) | Al Ge Pr Pr Al Ge | I4_1md (No. 109) | 12 | 3A1 + 3B1 + 6E |  | OK |
| 110 | [4031402](https://www.crystallography.net/cod/4031402.cif) | B20 H26 (B10 H13)2 | I4_1cd (No. 110) | 368 | 69A1 + 69A2 + 69B1 + 69B2 + 138E |  | OK |
| 111 | [1528224](https://www.crystallography.net/cod/1528224.cif) | F2 Mn Mn F2 | P-42m (No. 111) | 12 | 4A1 + 3A2 + 2B1 + 7B2 + 10E |  | OK |
| 112 | [2107218](https://www.crystallography.net/cod/2107218.cif) | C2 B2 La La B2 C2 | P4_2/mmc (No. 131) | 10 | 2A1g + 2A2g + 3A2u + 3B1g + 2B2g + 2B2u + 3Eg + 5Eu |  | OK |
| 113 | [7063440](https://www.crystallography.net/cod/7063440.cif) | Ba Se3  | P-42_1m (No. 113) | 8 | 3A1 + A2 + 2B1 + 4B2 + 7E |  | OK |
| 114 | [1545547](https://www.crystallography.net/cod/1545547.cif) | Bi2 O3  | P-42_1c (No. 114) | 20 | 7A1 + 7A2 + 7B1 + 7B2 + 16E |  | OK |
| 115 | [2106759](https://www.crystallography.net/cod/2106759.cif) | Hg I2  | P-4m2 (No. 115) | 12 | 5A1 + 3A2 + 2B1 + 6B2 + 10E |  | OK |
| 116 | [1522877](https://www.crystallography.net/cod/1522877.cif) | Ru2 Sn3  | P-4c2 (No. 116) | 20 | 5A1 + 7A2 + 6B1 + 8B2 + 17E | 1 coordinate value(s) written as rounded fractions (e.g. 0.333, 0.6667) were read as the exact fractions (1/3, 2/3, ...). | OK |
| 117 | [9012327](https://www.crystallography.net/cod/9012327.cif) | Bi2 O3 Sphaerobismoite | P-4b2 (No. 117) | 20 | 5A1 + 8A2 + 6B1 + 9B2 + 16E |  | OK |
| 118 | [1524662](https://www.crystallography.net/cod/1524662.cif) | Co In3  | P4_2/mnm (No. 136) | 16 | 3A1g + 2A1u + 2A2g + 4A2u + 2B1g + 4B1u + 3B2g + 2B2u + 4Eg + 9Eu |  | OK |
| 119 | [1011075](https://www.crystallography.net/cod/1011075.cif) | Hg2 I2 Moschelite | I4/mmm (No. 139) | 8 | 2A1g + 2A2u + 2Eg + 2Eu |  | OK |
| 120 | [1539545](https://www.crystallography.net/cod/1539545.cif) | Cr F4 Sr Sr Cr F4 | I4/mcm (No. 140) | 24 | 2A1g + 2A1u + A2g + 4A2u + B1g + 2B1u + 3B2g + B2u + 4Eg + 6Eu |  | OK |
| 121 | [1010030](https://www.crystallography.net/cod/1010030.cif) | B Fe2 Iron boride (2/1) | I-42m (No. 121) | 12 | 2A1 + A2 + B1 + 4B2 + 5E | 2 coordinate value(s) written as rounded fractions (e.g. 0.333, 0.6667) were read as the exact fractions (1/3, 2/3, ...). | OK |
| 122 | [1517802](https://www.crystallography.net/cod/1517802.cif) | C O2  | I-42d (No. 122) | 12 | A1 + 2A2 + 2B1 + 3B2 + 5E |  | OK |
| 123 | [1509184](https://www.crystallography.net/cod/1509184.cif) | Ag Ce  | P4/mmm (No. 123) | 2 | 2A2u + 2Eu |  | OK |
| 124 | [2300522](https://www.crystallography.net/cod/2300522.cif) | Nb Te4  | P4/mcc (No. 124) | 10 | 2A1g + A1u + 3A2g + 2A2u + 2B1g + B1u + 2B2g + B2u + 3Eg + 5Eu |  | OK |
| 125 | [1539051](https://www.crystallography.net/cod/1539051.cif) | Pb4 Pt  | P4/nbm (No. 125) | 10 | 2A1g + A1u + 2A2g + 3A2u + B1g + 2B1u + 2B2g + B2u + 4Eg + 4Eu |  | OK |
| 126 | [7104934](https://www.crystallography.net/cod/7104934.cif) | C108 F48  | P4/nnc (No. 126) | 320 | 60A1g + 60A1u + 60A2g + 60A2u + 60B1g + 60B1u + 60B2g + 60B2u + 120Eg + 120Eu | Partially occupied sites were treated as fully occupied: C12 (0.5) | OK |
| 127 | [2106199](https://www.crystallography.net/cod/2106199.cif) | Cs I  | P4/mmm (No. 123) | 4 | 2A2u + 2Eu |  | OK |
| 128 | [2106476](https://www.crystallography.net/cod/2106476.cif) | O20 S4 Na9 (S O4)4 Cl (H2 O2)2 | P4/mnc (No. 128) | 56 | 12A1g + 9A1u + 12A2g + 9A2u + 12B1g + 9B1u + 12B2g + 9B2u + 18Eg + 24Eu | Partially occupied sites were treated as fully occupied: O4 (0.5) | OK |
| 129 | [1010298](https://www.crystallography.net/cod/1010298.cif) | Fe Se Iron selenide | P4/nmm (No. 129) | 4 | A1g + 2A2u + B1g + 2Eg + 2Eu |  | OK |
| 130 | [1521532](https://www.crystallography.net/cod/1521532.cif) | O3 W W O3 | P4/ncc (No. 130) | 16 | 3A1g + 3A1u + 4A2g + 4A2u + 2B1g + 2B1u + B2g + B2u + 7Eg + 7Eu |  | OK |
| 131 | [1009031](https://www.crystallography.net/cod/1009031.cif) | O Pd Palladinite | P4_2/mmc (No. 131) | 4 | 2A2u + B1g + B2u + Eg + 3Eu |  | OK |
| 132 | [1509333](https://www.crystallography.net/cod/1509333.cif) | Ag F6 Ta Ag (Ta F6) | P4_2/mcm (No. 132) | 16 | 3A1g + A1u + 2A2g + 5A2u + 2B1g + 4B1u + 4B2g + B2u + 5Eg + 8Eu |  | OK |
| 133 | [1538473](https://www.crystallography.net/cod/1538473.cif) | P Ta3 Ta3 P | P4_2/nbc (No. 133) | 32 | 4A1g + 4A1u + 8A2g + 8A2u + 6B1g + 6B1u + 6B2g + 6B2u + 12Eg + 12Eu |  | OK |
| 134 | [1525095](https://www.crystallography.net/cod/1525095.cif) | B18 N2 B48 B2 N2 | I4/mmm (No. 139) | 20 | 2A1g + A2g + 4A2u + 2B1g + B2g + 2B2u + 3Eg + 6Eu |  | OK |
| 135 | [1010994](https://www.crystallography.net/cod/1010994.cif) | O2 Se Downeyite | P4_2/mbc (No. 135) | 24 | 5A1g + 3A1u + 6A2g + 4A2u + 6B1g + 4B1u + 5B2g + 3B2u + 7Eg + 11Eu |  | OK |
| 136 | [1000058](https://www.crystallography.net/cod/1000058.cif) | O2 Ru Ruthenium(IV) oxide | P4_2/mnm (No. 136) | 6 | A1g + A2g + 2A2u + B1g + 2B1u + B2g + Eg + 4Eu |  | OK |
| 137 | [1010063](https://www.crystallography.net/cod/1010063.cif) | Hg I2 MERCURY DIIODIDE | P4_2/nmc (No. 137) | 6 | A1g + 2A2u + 2B1g + B2u + 3Eg + 3Eu |  | OK |
| 138 | [1510610](https://www.crystallography.net/cod/1510610.cif) | Au Br  | P4_2/ncm (No. 138) | 8 | A1g + A1u + 3A2u + 3B1u + B2g + B2u + 2Eg + 5Eu |  | OK |
| 139 | [1523188](https://www.crystallography.net/cod/1523188.cif) | Ni Zn (Ni Zn) | I4/mmm (No. 139) | 2 | A2u + Eu | Sites sharing the same position were merged into one pseudo-atom: Zn1Ni1 (mass averaged by occupancy). | OK |
| 140 | [1009055](https://www.crystallography.net/cod/1009055.cif) | B Cr2 Chromium bromide (2/1) | I4/mcm (No. 140) | 12 | A1g + 2A2g + 2A2u + B1g + B1u + B2g + 2Eg + 3Eu | 2 coordinate value(s) written as rounded fractions (e.g. 0.333, 0.6667) were read as the exact fractions (1/3, 2/3, ...). | OK |
| 141 | [1523884](https://www.crystallography.net/cod/1523884.cif) | Hg Te (Hg Te) | I4_1/amd (No. 141) | 4 | A2u + B1g + Eg + Eu | Sites sharing the same position were merged into one pseudo-atom: Te1Hg1 (mass averaged by occupancy). | OK |
| 142 | [1521140](https://www.crystallography.net/cod/1521140.cif) | C2 Na2 Na2 C2 | I4_1/acd (No. 142) | 32 | 2A1g + 2A1u + 4A2g + 4A2u + 3B1g + 3B1u + 3B2g + 3B2u + 6Eg + 6Eu |  | OK |
| 143 | [1530161](https://www.crystallography.net/cod/1530161.cif) | H Ni2 Ni2 H | P3m1 (No. 156) | 3 | 3A1 + 3E | 4 coordinate value(s) written as rounded fractions (e.g. 0.333, 0.6667) were read as the exact fractions (1/3, 2/3, ...). | OK |
| 144 | [2105404](https://www.crystallography.net/cod/2105404.cif) | B2 O3  | P3_1 (No. 144) | 15 | 15A + 15E |  | OK |
| 145 | [4124123](https://www.crystallography.net/cod/4124123.cif) | Be Cl2  | P3_2 (No. 145) | 27 | 27A + 27E |  | OK |
| 146 | [1510571](https://www.crystallography.net/cod/1510571.cif) | Au5 Sn  | R32 (No. 155) | 18 | 2A1 + 4A2 + 6E | 2 coordinate value(s) written as rounded fractions (e.g. 0.333, 0.6667) were read as the exact fractions (1/3, 2/3, ...). | OK |
| 147 | [1536250](https://www.crystallography.net/cod/1536250.cif) | H Si Si H | P-3m1 (No. 164) | 4 | 2A1g + 2A2u + 2Eg + 2Eu | 4 coordinate value(s) written as rounded fractions (e.g. 0.333, 0.6667) were read as the exact fractions (1/3, 2/3, ...). | OK |
| 148 | [9008098](https://www.crystallography.net/cod/9008098.cif) | As I3  | R-3 (No. 148) | 8 | 4Ag + 4Au + 4Eg + 4Eu |  | OK |
| 149 | [1530197](https://www.crystallography.net/cod/1530197.cif) | U V2 U V2 O6 | P6/mmm (No. 191) | 3 | 2A2u + B2g + 2E1u + E2g | 4 coordinate value(s) written as rounded fractions (e.g. 0.333, 0.6667) were read as the exact fractions (1/3, 2/3, ...). | OK |
| 150 | [1010278](https://www.crystallography.net/cod/1010278.cif) | La2 O3 Lanthanum oxide - $-alpha | P321 (No. 150) | 5 | 2A1 + 3A2 + 5E | 2 coordinate value(s) written as rounded fractions (e.g. 0.333, 0.6667) were read as the exact fractions (1/3, 2/3, ...). | OK |
| 151 | [1538356](https://www.crystallography.net/cod/1538356.cif) | Cl3 Ti Ti Cl3 | P3_112 (No. 151) | 24 | 11A1 + 13A2 + 24E | 2 coordinate value(s) written as rounded fractions (e.g. 0.333, 0.6667) were read as the exact fractions (1/3, 2/3, ...). | OK |
| 152 | [1011369](https://www.crystallography.net/cod/1011369.cif) | Hg S Cinnabar | P3_121 (No. 152) | 6 | 2A1 + 4A2 + 6E | 3 coordinate value(s) written as rounded fractions (e.g. 0.333, 0.6667) were read as the exact fractions (1/3, 2/3, ...). | OK |
| 153 | [1010575](https://www.crystallography.net/cod/1010575.cif) | Cl3 Cr Chromium chloride | P3_212 (No. 153) | 24 | 11A1 + 13A2 + 24E |  | OK |
| 154 | [9012082](https://www.crystallography.net/cod/9012082.cif) | Hg S Cinnabar | P3_221 (No. 154) | 6 | 2A1 + 4A2 + 6E | 2 coordinate value(s) written as rounded fractions (e.g. 0.333, 0.6667) were read as the exact fractions (1/3, 2/3, ...). | OK |
| 155 | [1011250](https://www.crystallography.net/cod/1011250.cif) | Ni3 S2 Heazlewoodite | R32 (No. 155) | 5 | 2A1 + 3A2 + 5E |  | OK |
| 156 | [1530159](https://www.crystallography.net/cod/1530159.cif) | H Ni2 Ni2 H | P3m1 (No. 156) | 3 | 3A1 + 3E | 4 coordinate value(s) written as rounded fractions (e.g. 0.333, 0.6667) were read as the exact fractions (1/3, 2/3, ...). | OK |
| 157 | [1510073](https://www.crystallography.net/cod/1510073.cif) | Au Cd  | P31m (No. 157) | 18 | 12A1 + 6A2 + 18E | 4 coordinate value(s) written as rounded fractions (e.g. 0.333, 0.6667) were read as the exact fractions (1/3, 2/3, ...). | OK |
| 158 | [1520881](https://www.crystallography.net/cod/1520881.cif) | Co4 O12 Sr5 Sr5 Co4 O12 | P3c1 (No. 158) | 126 | 63A1 + 63A2 + 126E | 36 coordinate value(s) written as rounded fractions (e.g. 0.333, 0.6667) were read as the exact fractions (1/3, 2/3, ...). | OK |
| 159 | [1529015](https://www.crystallography.net/cod/1529015.cif) | O Ti6 Ti6 O | P-31c (No. 163) | 14 | 3A1g + 3A1u + 4A2g + 4A2u + 7Eg + 7Eu | 4 coordinate value(s) written as rounded fractions (e.g. 0.333, 0.6667) were read as the exact fractions (1/3, 2/3, ...). | OK |
| 160 | [1537749](https://www.crystallography.net/cod/1537749.cif) | N3 Na Na N3 | R3m (No. 160) | 4 | 4A1 + 4E | 3 coordinate value(s) written as rounded fractions (e.g. 0.333, 0.6667) were read as the exact fractions (1/3, 2/3, ...). | OK |
| 161 | [9017518](https://www.crystallography.net/cod/9017518.cif) | Fe O3 Hematite | R3c (No. 161) | 8 | 4A1 + 4A2 + 8E |  | OK |
| 162 | [7222689](https://www.crystallography.net/cod/7222689.cif) | Si Te2  | P6/mmm (No. 191) | 5 | A1g + 2A2u + B1u + B2g + E1g + 2E1u + E2g + E2u | 2 coordinate value(s) written as rounded fractions (e.g. 0.333, 0.6667) were read as the exact fractions (1/3, 2/3, ...).; Partially occupied sites were treated as fully occupied: Te1 (0.5) | OK |
| 163 | [1529955](https://www.crystallography.net/cod/1529955.cif) | O Ti6 Ti6 O | P-31m (No. 162) | 14 | 2A1g + A1u + A2g + 3A2u + 3Eg + 4Eu | 1 coordinate value(s) written as rounded fractions (e.g. 0.333, 0.6667) were read as the exact fractions (1/3, 2/3, ...). | OK |
| 164 | [1008063](https://www.crystallography.net/cod/1008063.cif) | Te2 Ti Titanium telluride | P-3m1 (No. 164) | 3 | A1g + 2A2u + Eg + 2Eu | 2 coordinate value(s) written as rounded fractions (e.g. 0.333, 0.6667) were read as the exact fractions (1/3, 2/3, ...). | OK |
| 165 | [1531032](https://www.crystallography.net/cod/1531032.cif) | K3 Rh K3 Rh (S C N)6 | P-3m1 (No. 164) | 8 | A1g + 3A2u + Eg + 3Eu | 2 coordinate value(s) written as rounded fractions (e.g. 0.333, 0.6667) were read as the exact fractions (1/3, 2/3, ...). | OK |
| 166 | [9008025](https://www.crystallography.net/cod/9008025.cif) | Rb S  | R-3m (No. 166) | 2 | 2A2u + 2Eu |  | OK |
| 167 | [1000476](https://www.crystallography.net/cod/1000476.cif) | F3 Fe Iron trifluoride | R-3c (No. 167) | 8 | A1g + 2A1u + 2A2g + 3A2u + 3Eg + 5Eu |  | OK |
| 168 | — | (入手不可) | | | | | |
| 169 | [1517832](https://www.crystallography.net/cod/1517832.cif) | Al2 S3 alpha-Al2S3 | P6_1 (No. 169) | 30 | 15A + 15B + 15E1 + 15E2 |  | OK |
| 170 | [4517545](https://www.crystallography.net/cod/4517545.cif) | In2 Se3  | P6_5 (No. 170) | 30 | 15A + 15B + 15E1 + 15E2 |  | OK |
| 171 | [4002432](https://www.crystallography.net/cod/4002432.cif) | O2 Si Si O2 | P6_2 (No. 171) | 72 | 36A + 36B + 36E1 + 36E2 |  | OK |
| 172 | [7240292](https://www.crystallography.net/cod/7240292.cif) | C19 H20 O5  | P6_4 (No. 172) | 270 | 135A + 135B + 135E1 + 135E2 | Partially occupied sites were treated as fully occupied: H4 (0.5), H5 (0.5) | OK |
| 173 | [1011053](https://www.crystallography.net/cod/1011053.cif) | C Si Moissanite 6H | P6_3mc (No. 186) | 12 | 6A1 + 6B1 + 6E1 + 6E2 | 10 coordinate value(s) written as rounded fractions (e.g. 0.333, 0.6667) were read as the exact fractions (1/3, 2/3, ...). | OK |
| 174 | [1010319](https://www.crystallography.net/cod/1010319.cif) | Ge I2 Germanium(II) iodide | P-6m2 (No. 187) | 3 | A1' + 2A2'' + 2E' + E'' | 2 coordinate value(s) written as rounded fractions (e.g. 0.333, 0.6667) were read as the exact fractions (1/3, 2/3, ...). | OK |
| 175 | [1510575](https://www.crystallography.net/cod/1510575.cif) | Au51 U14  | P6/m (No. 175) | 68 | 19Ag + 16Au + 15Bg + 18Bu + 15E1g + 20E1u + 19E2g + 14E2u | 4 coordinate value(s) written as rounded fractions (e.g. 0.333, 0.6667) were read as the exact fractions (1/3, 2/3, ...).; Partially occupied sites were treated as fully occupied: Au3 (0.5) | OK |
| 176 | [1008035](https://www.crystallography.net/cod/1008035.cif) | Br3 U Uranium(III) bromide | P6_3/m (No. 176) | 8 | 2Ag + 2Au + 2Bg + 2Bu + E1g + 3E1u + 3E2g + E2u | 2 coordinate value(s) written as rounded fractions (e.g. 0.333, 0.6667) were read as the exact fractions (1/3, 2/3, ...). | OK |
| 177 | [9017773](https://www.crystallography.net/cod/9017773.cif) | As4 Br K O6  | P6/mmm (No. 191) | 12 | 2A1g + 4A2u + B1g + 2B1u + 2B2g + B2u + 3E1g + 5E1u + 2E2g + 2E2u | 2 coordinate value(s) written as rounded fractions (e.g. 0.333, 0.6667) were read as the exact fractions (1/3, 2/3, ...). | OK |
| 178 | [1566652](https://www.crystallography.net/cod/1566652.cif) | D2 O ice XVII | P6_122 (No. 178) | 30 | 7A1 + 8A2 + 7B1 + 8B2 + 15E1 + 15E2 | Partially occupied sites were treated as fully occupied: D1 (0.5), D2 (0.5) | OK |
| 179 | [4001617](https://www.crystallography.net/cod/4001617.cif) | Ga9 Pd3 Yb2  | P6_522 (No. 179) | 84 | 19A1 + 23A2 + 20B1 + 22B2 + 42E1 + 42E2 | 3 coordinate value(s) written as rounded fractions (e.g. 0.333, 0.6667) were read as the exact fractions (1/3, 2/3, ...). | OK |
| 180 | [1010033](https://www.crystallography.net/cod/1010033.cif) | Cr Si2 Chromium silicide (1/2) | P6_222 (No. 180) | 9 | A1 + 3A2 + 2B1 + 3B2 + 5E1 + 4E2 | 4 coordinate value(s) written as rounded fractions (e.g. 0.333, 0.6667) were read as the exact fractions (1/3, 2/3, ...). | OK |
| 181 | [1528542](https://www.crystallography.net/cod/1528542.cif) | Al2 W W Al2 | P6_422 (No. 181) | 9 | A1 + 3A2 + 2B1 + 3B2 + 5E1 + 4E2 | 2 coordinate value(s) written as rounded fractions (e.g. 0.333, 0.6667) were read as the exact fractions (1/3, 2/3, ...). | OK |
| 182 | [1523551](https://www.crystallography.net/cod/1523551.cif) | Ni2 Si  | P6_3/mmc (No. 194) | 6 | 3A2u + B1u + 2B2g + 3E1u + 2E2g + E2u | 4 coordinate value(s) written as rounded fractions (e.g. 0.333, 0.6667) were read as the exact fractions (1/3, 2/3, ...). | OK |
| 183 | [1525811](https://www.crystallography.net/cod/1525811.cif) | Ta21 Te13  | P6mm (No. 183) | 136 | 46A1 + 23A2 + 37B1 + 30B2 + 70E1 + 66E2 | 6 coordinate value(s) written as rounded fractions (e.g. 0.333, 0.6667) were read as the exact fractions (1/3, 2/3, ...). | OK |
| 184 | [7110660](https://www.crystallography.net/cod/7110660.cif) | C18 H6 S3  | P6cc (No. 184) | 312 | 78A1 + 78A2 + 78B1 + 78B2 + 156E1 + 156E2 | Partially occupied sites were treated as fully occupied: C(21) (0.5), C(22) (0.5), S(4) (0.5) | OK |
| 185 | [1008751](https://www.crystallography.net/cod/1008751.cif) | F3 La Fluocerite (La) | P6_3cm (No. 185) | 24 | 8A1 + 4A2 + 4B1 + 8B2 + 12E1 + 12E2 | 2 coordinate value(s) written as rounded fractions (e.g. 0.333, 0.6667) were read as the exact fractions (1/3, 2/3, ...). | OK |
| 186 | [1010168](https://www.crystallography.net/cod/1010168.cif) | Ga N Gallium nitride | P6_3mc (No. 186) | 4 | 2A1 + 2B1 + 2E1 + 2E2 | 4 coordinate value(s) written as rounded fractions (e.g. 0.333, 0.6667) were read as the exact fractions (1/3, 2/3, ...). | OK |
| 187 | [1501516](https://www.crystallography.net/cod/1501516.cif) | C W  | P-6m2 (No. 187) | 2 | 2A2'' + 2E' | 2 coordinate value(s) written as rounded fractions (e.g. 0.333, 0.6667) were read as the exact fractions (1/3, 2/3, ...). | OK |
| 188 | [1511281](https://www.crystallography.net/cod/1511281.cif) | B O3 Tm Tm (B O3) | P6_3/mcm (No. 193) | 30 | 4A1g + A1u + 4A2g + 6A2u + 4B1g + 3B1u + 2B2g + 6B2u + 6E1g + 9E1u + 8E2g + 7E2u | 7 coordinate value(s) written as rounded fractions (e.g. 0.333, 0.6667) were read as the exact fractions (1/3, 2/3, ...). | OK |
| 189 | [2106407](https://www.crystallography.net/cod/2106407.cif) | N Ta Ta N | P-62m (No. 189) | 6 | A1' + A1'' + A2' + 3A2'' + 5E' + E'' | 2 coordinate value(s) written as rounded fractions (e.g. 0.333, 0.6667) were read as the exact fractions (1/3, 2/3, ...). | OK |
| 190 | [1530015](https://www.crystallography.net/cod/1530015.cif) | Nb S2 K0.5 (H2 O)0.4 Nb S2 | P6_3/mmc (No. 194) | 6 | A1g + 2A2u + B1u + 2B2g + E1g + 2E1u + 2E2g + E2u | 2 coordinate value(s) written as rounded fractions (e.g. 0.333, 0.6667) were read as the exact fractions (1/3, 2/3, ...). | OK |
| 191 | [1510614](https://www.crystallography.net/cod/1510614.cif) | C Au Au (C N) | P6/mmm (No. 191) | 3 | A1g + 2A2u + E1g + 2E1u | Partially occupied sites were treated as fully occupied: C1 (0.5) | OK |
| 192 | [1010541](https://www.crystallography.net/cod/1010541.cif) | Al2 Be3 O18 Si6 Triberyllium dialuminium catena-silicate | P6/mcc (No. 192) | 58 | 7A1g + 5A1u + 9A2g + 7A2u + 6B1g + 8B1u + 7B2g + 9B2u + 13E1g + 17E1u + 16E2g + 12E2u | 2 coordinate value(s) written as rounded fractions (e.g. 0.333, 0.6667) were read as the exact fractions (1/3, 2/3, ...). | OK |
| 193 | [1528149](https://www.crystallography.net/cod/1528149.cif) | Cs3 O  | P6_3/mcm (No. 193) | 8 | A1g + A2g + 2A2u + B1g + B1u + 2B2u + E1g + 3E1u + 2E2g + 2E2u |  | OK |
| 194 | [1509174](https://www.crystallography.net/cod/1509174.cif) | Ag Cd  | P6_3/mmc (No. 194) | 2 | A2u + B2g + E1u + E2g | 4 coordinate value(s) written as rounded fractions (e.g. 0.333, 0.6667) were read as the exact fractions (1/3, 2/3, ...).; Sites sharing the same position were merged into one pseudo-atom: Ag1Cd1 (mass averaged by occupancy). | OK |
| 195 | [1010043](https://www.crystallography.net/cod/1010043.cif) | I4 Sn Tin(IV) tetraiodide | P-43m (No. 215) | 5 | A1 + E + T1 + 3T2 |  | OK |
| 196 | [1523131](https://www.crystallography.net/cod/1523131.cif) | Li22 Pb5  | F-43m (No. 216) | 432 | 20A1 + 4A2 + 24E + 32T1 + 52T2 | 28 coordinate value(s) written as rounded fractions (e.g. 0.333, 0.6667) were read as the exact fractions (1/3, 2/3, ...). | OK |
| 197 | [1511475](https://www.crystallography.net/cod/1511475.cif) | B4 Li5 Li5 B4 | I-43m (No. 217) | 16 | 2A1 + 2E + 2T1 + 4T2 | Sites sharing the same position were merged into one pseudo-atom: Li1B1 (mass averaged by occupancy).; Partially occupied sites were treated as fully occupied: B2 (0.125) | OK |
| 198 | [1010034](https://www.crystallography.net/cod/1010034.cif) | Cr Si Chromium silicide (1/1) | P2_13 (No. 198) | 8 | 2A + 2E + 6T |  | OK |
| 199 | [2310270](https://www.crystallography.net/cod/2310270.cif) | Co U  | I2_13 (No. 199) | 16 | 2A + 2E + 6T |  | OK |
| 200 | [1522821](https://www.crystallography.net/cod/1522821.cif) | Mg2 Zn11  | Pm-3 (No. 200) | 39 | 6Ag + 2Au + 6Eg + 2Eu + 13Tg + 18Tu |  | OK |
| 201 | [1010382](https://www.crystallography.net/cod/1010382.cif) | O Pb2 Lead(I) oxide | Pn-3m (No. 224) | 6 | A2u + Eu + 3T1u + T2g + T2u |  | OK |
| 202 | [1528292](https://www.crystallography.net/cod/1528292.cif) | Bi2 O3  | Pm-3 (No. 200) | 248 | 4Ag + 3Au + 4Eg + 3Eu + 11Tg + 13Tu | Partially occupied sites were treated as fully occupied: Bi3 (0.5), O1 (0.25), O2 (0.25) | OK |
| 203 | [1530929](https://www.crystallography.net/cod/1530929.cif) | C Fe6 W6 Fe6 W6 C | Fd-3m (No. 227) | 104 | 2A1g + 3A2u + 2Eg + 3Eu + 3T1g + 8T1u + 6T2g + 4T2u |  | OK |
| 204 | [1538030](https://www.crystallography.net/cod/1538030.cif) | O4 Pt3 Pt3 O4 | Im-3m (No. 229) | 14 | A2u + Eu + 4T1u + 2T2u |  | OK |
| 205 | [1527512](https://www.crystallography.net/cod/1527512.cif) | Ir3 Te8  | Pa-3 (No. 205) | 12 | Ag + 2Au + Eg + 2Eu + 3Tg + 6Tu | Partially occupied sites were treated as fully occupied: Ir1 (0.75) | OK |
| 206 | [1541702](https://www.crystallography.net/cod/1541702.cif) | Ca U Ca U O3 | Ia-3 (No. 206) | 32 | Ag + 2Au + Eg + 2Eu + 5Tg + 8Tu | Sites sharing the same position were merged into one pseudo-atom: Ca1U1, U2Ca2 (mass averaged by occupancy). | OK |
| 207 | [9011053](https://www.crystallography.net/cod/9011053.cif) | Pd17 Se15 Palladseite | P432 (No. 207) | 64 | 7A1 + 7A2 + 14E + 28T1 + 22T2 |  | OK |
| 208 | [1010125](https://www.crystallography.net/cod/1010125.cif) | As2 Zn3 Zinc arsenide | Pn-3m (No. 224) | 10 | A2u + Eu + T1g + 4T1u + 2T2g + 2T2u |  | OK |
| 209 | [4302975](https://www.crystallography.net/cod/4302975.cif) | C42 H24 Cu3 O15 PCN-20 | F432 (No. 209) | 1344 | 41A1 + 43A2 + 84E + 127T1 + 125T2 |  | OK |
| 210 | [1010954](https://www.crystallography.net/cod/1010954.cif) | O2 Si Cristobalite high | Fd-3m (No. 227) | 24 | A2u + Eu + 3T1u + T2g + T2u |  | OK |
| 211 | [1537751](https://www.crystallography.net/cod/1537751.cif) | Hg4 Pt  | Im-3m (No. 229) | 10 | A2u + Eu + 3T1u + T2u |  | OK |
| 212 | [1536083](https://www.crystallography.net/cod/1536083.cif) | Si2 Sr Sr Si2 | P4_332 (No. 212) | 12 | A1 + 2A2 + 3E + 5T1 + 4T2 |  | OK |
| 213 | [1542132](https://www.crystallography.net/cod/1542132.cif) | Si2 Sr Sr Si2 | P4_132 (No. 213) | 12 | A1 + 2A2 + 3E + 5T1 + 4T2 |  | OK |
| 214 | [1509884](https://www.crystallography.net/cod/1509884.cif) | Ag3 Au Se2  | I4_132 (No. 214) | 48 | 2A1 + 3A2 + 5E + 10T1 + 9T2 |  | OK |
| 215 | [9004225](https://www.crystallography.net/cod/9004225.cif) | Fe4 N Roaldite | P-43m (No. 215) | 5 | A1 + E + T1 + 3T2 |  | OK |
| 216 | [1010076](https://www.crystallography.net/cod/1010076.cif) | Br Cu Copper(I) bromide | F-43m (No. 216) | 8 | 2T2 |  | OK |
| 217 | [1010134](https://www.crystallography.net/cod/1010134.cif) | F4 Si Silicon fluoride | I-43m (No. 217) | 10 | A1 + E + T1 + 3T2 |  | OK |
| 218 | [4124120](https://www.crystallography.net/cod/4124120.cif) | Be F2  | I-43m (No. 217) | 36 | 2A1 + 2A2 + 4E + 6T1 + 8T2 |  | OK |
| 219 | [2106953](https://www.crystallography.net/cod/2106953.cif) | Cl16 Si6 (Si5 Cl12) (Si Cl4) | F-43c (No. 219) | 176 | 5A1 + 5A2 + 10E + 17T1 + 17T2 |  | OK |
| 220 | [1010170](https://www.crystallography.net/cod/1010170.cif) | P4 Th3 Thorium phosphide | I-43d (No. 220) | 28 | A1 + 2A2 + 3E + 5T1 + 6T2 | 3 coordinate value(s) written as rounded fractions (e.g. 0.333, 0.6667) were read as the exact fractions (1/3, 2/3, ...). | OK |
| 221 | [1000050](https://www.crystallography.net/cod/1000050.cif) | Cl K Potassium chloride | Pm-3m (No. 221) | 2 | 2T1u |  | OK |
| 222 | [1001073](https://www.crystallography.net/cod/1001073.cif) | Mo3 Nd5 O16 Pentaneodymium 16-oxotrimolybdate | Pn-3n (No. 222) | 96 | 5A1g + 7A1u + 5A2g + 5A2u + 10Eg + 12Eu + 17T1g + 21T1u + 17T2g + 19T2u | 1 coordinate value(s) written as rounded fractions (e.g. 0.333, 0.6667) were read as the exact fractions (1/3, 2/3, ...). | OK |
| 223 | [1010032](https://www.crystallography.net/cod/1010032.cif) | Cr3 Si Chromium silicide (3/1) | Pm-3n (No. 223) | 8 | A2g + Eg + T1g + 3T1u + T2g + 2T2u |  | OK |
| 224 | [1000063](https://www.crystallography.net/cod/1000063.cif) | Cu2 O Cuprite | Pn-3m (No. 224) | 6 | A2u + Eu + 3T1u + T2g + T2u |  | OK |
| 225 | [1509205](https://www.crystallography.net/cod/1509205.cif) | Ag Au  | Fm-3m (No. 225) | 4 | T1u | Sites sharing the same position were merged into one pseudo-atom: Ag1Au1 (mass averaged by occupancy). | OK |
| 226 | [1524778](https://www.crystallography.net/cod/1524778.cif) | Eu Zn13  | Fm-3c (No. 226) | 112 | 2A1g + A1u + 2A2g + A2u + 4Eg + 2Eu + 5T1g + 7T1u + 4T2g + 6T2u |  | OK |
| 227 | [1537802](https://www.crystallography.net/cod/1537802.cif) | Ge Si  | Fd-3m (No. 227) | 8 | T1u + T2g | Sites sharing the same position were merged into one pseudo-atom: Si1Ge1 (mass averaged by occupancy). | OK |
| 228 | [1534373](https://www.crystallography.net/cod/1534373.cif) | H6 O6 Te Te (O H)6 | Fd-3c (No. 228) | 608 | 9A1g + 10A1u + 9A2g + 10A2u + 18Eg + 20Eu + 27T1g + 30T1u + 27T2g + 30T2u | 1 coordinate value(s) written as rounded fractions (e.g. 0.333, 0.6667) were read as the exact fractions (1/3, 2/3, ...).; Partially occupied sites were treated as fully occupied: H1 (0.5), H2 (0.5) | OK |
| 229 | [1509173](https://www.crystallography.net/cod/1509173.cif) | Ag Cd  | Im-3m (No. 229) | 2 | T1u | Sites sharing the same position were merged into one pseudo-atom: Cd1Ag1 (mass averaged by occupancy). | OK |
| 230 | [1523546](https://www.crystallography.net/cod/1523546.cif) | Ga4 Ni3  | Ia-3d (No. 230) | 112 | 2A1g + 3A1u + 3A2g + 4A2u + 5Eg + 7Eu + 10T1g + 13T1u + 9T2g + 12T2u |  | OK |
