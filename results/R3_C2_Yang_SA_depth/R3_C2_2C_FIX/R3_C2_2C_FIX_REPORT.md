# R3-C2.2C-FIX 报告

UTC: 2026-09-24T03:20:50.381015+00:00

## 状态

`C2_2C_PAPER_REPRO_CONFIRMED` → **SUPERSEDED_PENDING_INTEGRITY_FIX**

本阶段判定：**`C2_2C_FIX_IMPLEMENTATION_FAIL`**

FIELD .shd 生成并读取成功，但与 Eq.(3) 多模 modesum 在单一全局复标度下 residual~1、complex corr 低，不能互证；诊断显示 FIELD 近似单模(1/sqrt(r))而 Eq3 多模干涉不一致。按完整性规则记 IMPLEMENTATION_FAIL，不确认论文复现。

## 修复点

1. **FIELD 真交叉验证**：恢复官方 `.flp`（profile@0 km；NProf 单独行；tail `1/0.0/0.0`）；从 `.shd` 读压力；与 Eq.(3) 仅拟合 **全局复标度 c\*** 后比较 corr/residual（不再 self-compare）。
2. **ORACLE_MODE_ID**：`find_peaks` 先检测真实谱峰（prominence 主基线 3%；1/3/5/10% 敏感性）→ true k_m **一对一标注** → 仅 SELECTED 复峰进 Eq.(6)。全模态仅 `ORACLE_FULL_MODE_UPPER_BOUND`。
3. **depth gate**：`|z_hat-z_s|<=2 m` 且 `D(true)/D(max)>=0.8`；删除“4→5 最近节点”表述。
4. **Delta**：`Delta_abs = ratio * max|phi|` 全链一致。
5. **shading**：reflect padding，三窗全报。
6. **mode-count 23 vs 16**：取消 gate，仅 DIAGNOSTIC。

## 自证

- 单模 Eq1→Eq5 max_rel = 4.946e-04
- alpha 映射：alpha_m = -Im(k_complex); Im(k)<0 => alpha>0
- FIELD：FAIL/BLOCKED

## 范围

不进 E-STD / 50–60 km / 2.4 km CZ / δ 搜索 / MC / AR / P5。YANG_ROUTE_UNDECIDED。

## 停止

等待 GPT 审计。

## Peak-based Eq6 (Eq3 data; NOT paper confirm)

| zs | zr | n_peaks | n_sel | z_hat | abs_err | ratio |
| --- | --- | --- | --- | --- | --- | --- |
| 4 | 18 | 32 | 10 | 5 | 1.0 | 0.88 |
| 4 | 70 | 44 | 10 | 5 | 1.0 | 0.87 |
| 50 | 18 | 36 | 13 | 43 | 7.0 | 0.86 |
| 50 | 70 | 38 | 14 | 50 | 0.0 | 1.00 |

Depth gate would pass 3/4 but FIELD multimode cross-check failed => no PAPER_REPRO_*.
Keep: single-mode Eq1-Eq5, peak->oracle map, alpha=-Im(k), env lock.
Drop: old 4/4 confirm and field self-compare corr=1.
