# Track A validation (200 reps, 199 permutations, alpha = 0.05; exact-test reference at 199 perms = 0.050)

Binomial 95% band for a 0.05 rejection rate with this many reps: 
[0.020, 0.080]

| id | scenario | metric | contrast | scheme | truth shift | est. shift (sd) | truth var | est. var | rej shift | rej var | rej total | old pair-LM shift (d2) | old rej (d2) | old rej (d) |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| 1 | 2grp null | l2 | triple | block | 0.00 | -0.28 (7.02) | 0.0 | -0.9 | 0.045 | 0.050 | 0.045 | -0.28 | 0.040 | 0.030 |
| 2 | 2grp shift .3 | l2 | triple | block | 7.20 | 7.36 (8.11) | 0.0 | 1.8 | 0.175 | 0.045 | 0.070 | 7.36 | 0.175 | 0.200 |
| 3 | 2grp disp 2x | l2 | triple | block | 0.00 | -1.14 (10.77) | 160.0 | 160.0 | 0.025 | 0.925 | 0.865 | -1.14 | 0.020 | 0.075 |
| 4 | 2grp shift .3 + disp 2x | l2 | triple | block | 7.20 | 7.67 (12.27) | 160.0 | 160.9 | 0.110 | 0.930 | 0.925 | 7.67 | 0.120 | 0.270 |
| 5 | 2grp null, cor | cor | triple | block |  | -0.00 (0.04) |  | 0.0 | 0.040 | 0.035 | 0.035 |  |  |  |
| 6 | 2grp shift .3, cor | cor | triple | block |  | 0.05 (0.05) |  | 0.0 | 0.195 | 0.080 | 0.210 |  |  |  |
| 7 | 2grp disp 2x, cor | cor | triple | block |  | 0.01 (0.05) |  | 0.1 | 0.080 | 0.750 | 0.155 |  |  |  |
| 8 | 2grp null, l1 | l1 | triple | block |  | -0.03 (0.98) |  | -0.1 | 0.025 | 0.035 | 0.030 |  |  |  |
| 9 | 2grp shift .3, l1 | l1 | triple | block |  | 2.95 (1.42) |  | 0.1 | 0.750 | 0.045 | 0.255 |  |  |  |
| 10 | 2grp disp 2x, l1 | l1 | triple | block |  | 1.29 (1.29) |  | 30.5 | 0.215 | 1.000 | 1.000 |  |  |  |
| 11 | 2grp null, isotropic p200 | l2 | triple | block | 0.00 | -0.31 (5.48) | 0.0 | 0.6 | 0.050 | 0.040 | 0.035 | -0.31 | 0.055 | 0.060 |
| 12 | 2grp shift .2, iso p200 | l2 | triple | block | 8.00 | 7.98 (6.21) | 0.0 | 1.6 | 0.435 | 0.030 | 0.095 | 7.98 | 0.430 | 0.420 |
| 13 | uneq null | l2 | triple | block | 0.00 | -0.03 (10.82) | 0.0 | -2.9 | 0.040 | 0.060 | 0.050 | -0.03 | 0.050 | 0.055 |
| 14 | uneq disp 2x (B large) | l2 | triple | block | 0.00 | -2.52 (10.90) | 160.0 | 157.6 | 0.000 | 0.650 | 0.295 | -2.52 | 0.005 | 0.030 |
| 15 | uneq shift .3 | l2 | triple | block | 7.20 | 7.28 (10.84) | 0.0 | 1.0 | 0.130 | 0.035 | 0.055 | 7.28 | 0.115 | 0.125 |
| 16 | 3grp null (C shifted) | l2 | triple | block | 0.00 | -0.97 (7.31) | 0.0 | -1.2 | 0.030 | 0.060 | 0.040 | -0.97 | 0.025 | 0.025 |
| 17 | 3grp shift .3 | l2 | triple | block | 7.20 | 7.58 (9.91) | 0.0 | 3.8 | 0.150 | 0.045 | 0.080 | 7.58 | 0.150 | 0.155 |
| 18 | 3grp null, C disp | l2 | triple | block | 0.00 | 0.29 (8.78) | 0.0 | -3.7 | 0.035 | 0.040 | 0.070 | 0.29 | 0.040 | 0.050 |
| 19 | batch null, eff .4 | l2 | triple | block | 0.00 | 1.21 (8.46) | 0.0 | -3.0 | 0.055 | 0.060 | 0.035 | 0.95 | 0.025 | 0.015 |
| 20 | batch shift .3, orth | l2 | triple | block | 7.20 | 6.98 (8.84) | 0.0 | 3.6 | 0.180 | 0.035 | 0.070 | 6.85 | 0.085 | 0.035 |
| 21 | batch.conf null, eff .4 | l2 | triple | block | 0.00 | 0.19 (12.35) | 0.0 | 0.2 | 0.061 | 0.045 | 0.061 | 0.97 | 0.045 | 0.015 |
| 22 | batch.conf null, eff .8 | l2 | triple | block | 0.00 | 0.25 (12.78) | 0.0 | -0.2 | 0.035 | 0.056 | 0.040 | 0.21 | 0.025 | 0.015 |
| 23 | batch.conf shift .3, cos+.6 | l2 | triple | block | 7.20 | 7.77 (12.44) | 0.0 | 0.6 | 0.136 | 0.030 | 0.076 | 11.87 | 0.212 | 0.126 |
| 24 | batch.conf shift .3, cos-.6 | l2 | triple | block | 7.20 | 7.80 (16.03) | 0.0 | 0.9 | 0.102 | 0.056 | 0.046 | 3.41 | 0.056 | 0.015 |
| 25 | batch.conf shift .3, cor | cor | triple | block |  | 0.06 (0.06) |  | -0.0 | 0.162 | 0.010 | 0.117 |  |  |  |
| 26 | age null, eff .3 | l2 | triple | freedman-lane | 0.00 | -1.29 (9.26) | 0.0 | -1.3 | 0.020 | 0.035 | 0.035 | -2.40 | 0.025 | 0.020 |
| 27 | age shift .3, eff .3 | l2 | triple | freedman-lane | 7.20 | 7.48 (10.50) | 0.0 | -1.2 | 0.150 | 0.055 | 0.065 | 7.34 | 0.165 | 0.145 |
| 28 | age null, eff .3, iso p200 | l2 | triple | freedman-lane | 0.00 | 0.65 (7.09) | 0.0 | 1.9 | 0.020 | 0.015 | 0.010 | -3.33 | 0.000 | 0.005 |
| 29 | age null, eff .3, n40 | l2 | triple | freedman-lane | 0.00 | -0.32 (3.05) | 0.0 | -0.4 | 0.025 | 0.020 | 0.055 | -1.30 | 0.000 | 0.000 |
| 30 | age shift .2, n40 | l2 | triple | freedman-lane | 3.20 | 2.99 (3.41) | 0.0 | -0.3 | 0.180 | 0.035 | 0.055 | 2.40 | 0.095 | 0.105 |
| 31 | batch+age null | l2 | triple | freedman-lane | 0.00 | -0.31 (8.99) | 0.0 | 6.1 | 0.020 | 0.070 | 0.055 | -1.44 | 0.025 | 0.035 |
| 32 | batch+age shift .3 | l2 | triple | freedman-lane | 7.20 | 7.98 (11.94) | 0.0 | -2.3 | 0.110 | 0.060 | 0.105 | 5.11 | 0.100 | 0.085 |
| 33 | marginal null, inter .4 | l2 | marginal | block | 3.20 | 3.06 (8.57) | 0.0 | 1.4 | 0.076 | 0.061 | 0.056 |  |  |  |
| 34 | marginal shift .3 | l2 | marginal | block | 10.40 | 9.47 (10.21) | 0.0 | -0.6 | 0.215 | 0.050 | 0.095 |  |  |  |
| 35 | interaction cell null | l2 | interaction | block | 0.00 | -0.93 (17.08) | 0.0 | -2.9 | 0.035 | 0.035 | 0.070 |  |  |  |
| 36 | interaction cell shift .4 | l2 | interaction | block | 12.80 | 11.92 (15.58) | 0.0 | -2.4 | 0.201 | 0.025 | 0.070 |  |  |  |
| 37 | age slope null | l2 | numeric | block |  | -0.00 (0.03) | 0.0 |  | 0.025 |  |  |  |  |  |
| 38 | age slope eff .3 | l2 | numeric | block |  | 0.07 (0.03) | 0.0 |  | 0.655 |  |  |  |  |  |
| 39 | batch null, FL forced | l2 | triple | freedman-lane | 0.00 | 0.54 (8.77) | 0.0 | 0.1 | 0.050 | 0.075 | 0.060 | 0.52 | 0.045 | 0.020 |
| 40 | batch shift .3, FL forced | l2 | triple | freedman-lane | 7.20 | 7.87 (9.09) | 0.0 | 0.1 | 0.155 | 0.045 | 0.085 | 7.81 | 0.115 | 0.060 |
| 41 | batch null, HJ forced | l2 | triple | huh-jhun | 0.00 | 0.03 (8.85) | 0.0 | -1.9 | 0.065 |  |  | 0.15 | 0.045 | 0.020 |
| 42 | 2grp disp 2x, HJ forced | l2 | triple | huh-jhun | 0.00 | -0.51 (11.44) | 160.0 | 156.2 | 0.050 |  |  | -0.51 | 0.050 | 0.095 |

## Combination across cell types under the null (8 cell types sharing a sample-level factor)

| shared variance | global max-T rej | any cell type FWER-sig | any BH-sig | any Bonferroni-sig | any raw p<.05 |
|---|---|---|---|---|---|
| 0.0 | 0.055 | 0.055 | 0.045 | 0.045 | 0.295 |
| 0.5 | 0.035 | 0.035 | 0.020 | 0.020 | 0.220 |
| 0.9 | 0.040 | 0.040 | 0.010 | 0.010 | 0.145 |
