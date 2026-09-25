# FORWARD_GATE_AUDIT

UTC: 2026-09-25T08:18:43.861670+00:00

## 旧 Gate 无效

单点 `c=p_fld/p_ms` ⇒ residual 恒为 0，`field_ok=True` 无信息量。

## 修正

- 复用 `r3_c2_3a_estd.py` 的 `parse_shd` / `write_env_estd` / `write_flp` / `field_pressure`
- **整条 range vector** 一个 `c* = p_MS^H p_FIELD / p_MS^H p_MS`
- 禁止逐点 scale
- 报告 complex corr / normalized residual / amp-shape / phase + **TL_SHAPE**（真正使用的 observable）
- Gate 失败 **必须停止**（不再 `if fail: pass`）

## Gate 阈值（沿用 FIELD_EQ3_MULTIMODE_VALIDATED 精神）

- complex_corr ≥ 0.95
- normalized residual ≤ 0.35
- TL_SHAPE_RMS_DB ≤ 1.5
- TL_SHAPE_CORRELATION ≥ 0.95

## 结果

ESTD_RELATIVE_TL_FORWARD_MODEL_VALIDATED

multi-range FIELD vs mode-sum: 12/12 cases pass (cc>=0.95, resid<=0.35, TL_rms<=1.5 dB, TL_corr>=0.95); one global c* per (f,zs); no per-point scale; REANCHOR-1 90 cases not recomputed

REANCHOR-1 的 90 case **不重算**；仅解除 forward Gate。
