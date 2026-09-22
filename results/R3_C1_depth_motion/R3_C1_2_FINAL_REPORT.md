# R3-C1.2 FINAL-INTEGRITY 报告

UTC：2026-09-22T03:19:50.624847+00:00

修正：真实 endfire 波束；统一多频 D；固定 Δr 采样；FULL_WINDOW 标注。

- 波束 steered≠unsteered：**True**

## T=600 四组

| 组 | FWHM | 全窗 | D200/210 | D200/220 |
| --- | --- | --- | --- | --- |
| A | 47.8 | 0.00 | 0.00309 | 0.00672 |
| B | 47.8 | 0.00 | 0.00306 | 0.00667 |
| C | 80.0 | 1.00 | 0.00120 | 0.00931 |
| D | 80.0 | 1.00 | 0.00119 | 0.00935 |

分辨边界：**在已测10/20 m近邻内未达到可用分辨**

### `C1_QM_CZ_DEPTH_WEAK`

FINAL-INTEGRITY: real endfire B=|w^H p|^2 (differs=True), unified D=1-mean|corr_f|, dr=10.0m; 200±10/20 m still weak. T=600 A D210=0.0031 D220=0.0067 FWHM=47.8; D D210=0.0012 D220=0.0093 FWHM=80.0; hla_help=False multi_help=False. 在已测10/20 m近邻内未达到可用分辨. Known-track ceiling.

下一步：Zhu PERMANENTLY_CLOSED; Yang 2015 next (not this round)

Zhu **PERMANENTLY_CLOSED**；无 C1.x / P5 / Yang 本轮。

## 修正终判（含 T=1200 多频）

> **SUPERSEDED**：若上文仍出现 `C1_QM_CZ_DEPTH_WEAK`，以 `R3_C1_2_FINAL_DECISION.json` 与 `R3_C1_2_FINAL_CORRECTION.md` 为准。

**`C1_QM_MULTIFREQ_ONLY_CONDITIONAL`**

FINAL-INTEGRITY passed (real endfire B=|w^H p|^2 differs from unsteered; unified D=1-mean_f|corr_f|; dr=10m converged). 10 m pairs (200/210) inseparable at all T (D~0.001–0.003 << 0.05). 20 m pair (200/220) only limited multifreq gain at T=1200s (D_D=0.260 vs A~0.006) while FWHM stays 63–80 m. 10 m not reached; 20 m weak/conditional only. Not reliable 10–20 m depth discrimination under known-track ceiling.

下一步：Zhu route PERMANENTLY_CLOSED after this gate; Yang 2015 next (not this round)

Zhu **PERMANENTLY_CLOSED**。
