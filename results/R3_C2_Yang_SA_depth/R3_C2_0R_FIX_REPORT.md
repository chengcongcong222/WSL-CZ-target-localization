# R3-C2.0R-FIX 报告

UTC：2026-09-22T08:07:38.691452+00:00

## 1. Erratum RECOVERED
- δ 相位 exp(i k_m δ) 进入深度模糊；Yang 非“无需源距离”
- 不再 `C2_0_PAPER_RECOVERY_BLOCKED`（Erratum）

## 2. 本征算子
- `[D2+ω²/c²]φ=k²φ`，diag=ω²/c²−2/h²，offdiag=+1/h²
- constant-c：**PASS**（k_r 机器精度；子空间覆盖 sin(γ_n z)；近简并向量混叠不改谱）
- k_r≤ω/c_min 断言：**PASS**
- 2001/4001 收敛：**PASS**（不再用 8001 全谱）

## 3. 修正后物理结果

 f_hz  C1_corr_200_220  n_propagating  n_eff_0p2  phys_D_200_220  phys_D_200_210  max_k_r  kr_max_allowed
201.0            0.579           1206        105        0.942034        0.916184 0.841881        0.841947
235.0            0.994           1436        101        0.939057        0.934755 0.984300        0.984366
283.0            1.000           1793        104        0.914969        0.985506 1.185362        1.185428
338.0            0.067           2292         99        0.938029        0.948864 1.415745        1.415811

C1 338 Hz：**NOT_EXPLAINED**

旧 SUPERSEDED_PENDING 数值以 `*_fixed.csv` 为准。

## 4. 判定

### `C2_0R_MODE_BASELINE_VALIDATED`

Operator sign corrected; constant-c test PASS (k_r exact, subspace covers sin(γ_n z)); assert k_r≤ω/c_min PASS; 2001/4001 converged. Erratum recovered (δ). C1_338=NOT_EXPLAINED. Prior R3-C2.0R D/L_eta1/338 numbers SUPERSEDED_PENDING→_fixed.

**下一步**：no full Yang / δ-z search / Doppler / f0 / MC / P5

不做完整 Yang / δ-z / Doppler / f0 / MC / P5。
