# R3-C2.4A-PRIMARY-FIX 报告

UTC: 2026-09-24T09:16:56.736527+00:00

## 判定

### C2_4A_PRIMARY_METHOD_LOCKED_WITH_NOTATION_CONFLICTS_V2

锁定稿转录错误已修正；方法链不变。**≠ 方法失败。**

## 必改项（已完成）

1. **Eq.(10)**：`sin b(X_m)` 为**乘法因子**，不在分母（V1 转录错误）
2. **sin b** → `SINB_NOTATION_UNDEFINED_USE_EQ4_BEAM_FACTOR`（BF_m 自 Eq.4）
3. **FROM_EQ9_DERIVED_BM** 已写出（α→0 稳定极限）
4. **EQ19_DENOMINATOR_POWER_NOT_RESOLVED**（印刷 |A|^1 vs 常规 |A|^2）
5. **PAPER_RANGE_SAMPLE_INTERVAL_NOT_EXPLICIT**（1 s 只能标 REF7_CONSISTENT_SAMPLING_ASSUMPTION）
6. **Eq.(29)(30)** 补全：P=C/C0；90% CI = P ± 1.645√[P(1-P)/C0]
7. **HLA**：d=λ/2，`NOT_FIXED_BY_TEXT` + `PAPER_LAMBDA_REFERENCE_NOT_EXPLICIT`
8. **FIG3/4 图注面板笔误** 已记录

## 冲突登记（10 项）

见 `MMAR_IMPLEMENTATION_CONFLICT_REGISTER.csv`。

## 下一步（未执行）

R3-C2.4B：可控复指数 AR → 1990/4990 m 无噪/20dB → k̂ → Eq24 **MODE_ORDER_RECOVERY_RATE** → g(k̂) → D(z)。

本轮停止。
