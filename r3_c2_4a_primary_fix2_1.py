#!/usr/bin/env python3
"""R3-C2.4A-PRIMARY-FIX2.1: minimal source-fidelity patch. No AR, no 4B."""
from __future__ import annotations

import hashlib
import json
from datetime import datetime, timezone
from pathlib import Path

ROOT = Path(__file__).resolve().parent
PDF = ROOT / "literature" / "liang2018_mmar_7824671.pdf"
OUT = ROOT / "results" / "R3_C2_AR_MMAR" / "R3_C2_4A_PRIMARY_FIX2"
V2 = ROOT / "results" / "R3_C2_AR_MMAR" / "R3_C2_4A_PRIMARY_FIX"
NOW = datetime.now(timezone.utc).isoformat()
SHA = hashlib.sha256(PDF.read_bytes()).hexdigest()


def main() -> int:
    # 3) V2 superseded marker
    (V2 / "SUPERSEDED_BY_R3_C2_4A_PRIMARY_FIX2_1.md").write_text(
        f"""# SUPERSEDED_BY_R3_C2_4A_PRIMARY_FIX2_1

UTC: {NOW}

`results/R3_C2_AR_MMAR/R3_C2_4A_PRIMARY_FIX/`（V2）不得作为实现输入。
唯一入口：`R3_C2_4A_PRIMARY_FIX2/MMAR_EQUATION_LOCK_FINAL.md`（FIX2.1 修订版）。
""",
        encoding="utf-8",
    )

    # 1+4) rewrite EQUATION_LOCK_FINAL with strict PRINTED_PRIMARY for Eq.9/19/24
    lock = OUT / "MMAR_EQUATION_LOCK_FINAL.md"
    t = lock.read_text(encoding="utf-8")
    # replace Eq.(9) block to restore printed sin b
    old9_start = "## Eq.(9)"
    old10_start = "## Eq.(10)"
    i9 = t.find(old9_start)
    i10 = t.find(old10_start)
    if i9 < 0 or i10 < 0 or i10 <= i9:
        raise RuntimeError("Eq block markers missing")
    new9 = r"""## Eq.(9) `PRINTED_PRIMARY`（严格原文印刷）

```latex
g(k_r,z_r)\sim\sum_{m=1}^{M}
\frac{\phi_m(z_s)\phi_m(z_r)}{\sqrt{k_r k_m}}\,\sin b(X_m)
\int_{r_0}^{r_0+R} e^{j(k_r-k_m)r-\alpha_m r}\,dr
=\sum_{m=1}^{M} a_m\frac{\phi_m(z_s)\phi_m(z_r)}{k_r-k_m+j\alpha_m}
```

积分核：**exp[j(k_r−k_m)r − α_m r]**。此处 **sin b(X_m) 为论文印刷原样**，未替换为 BF_m / mathcal B_m。

---

## Eq.(9) `DERIVED_IMPLEMENTATION`（单独解释层）

```latex
\sin b(X_m)\ \Rightarrow\ BF_m=\frac{1}{2L+1}\cdot\frac{\sin[(L+1/2)d X_m]}{\sin[(d/2) X_m]}
```

标记：**`DERIVED_EQ4_BEAM_FACTOR_INTERPRETATION`**（由 Eq.(4)→Eq.(7) 代数，非论文明文定义 sin b）。

实现公式（`FROM_EQ9_DERIVED_BM`，与 PRINTED_PRIMARY 分离）：

```latex
b_m^{\mathrm{impl}}=\frac{2e^{-\alpha_m r'}}{\alpha_m k_m}\sinh\left(\frac{\alpha_m R}{2}\right)\phi_m(z_s)\,BF_m,
\quad r'=r_0+R/2
```

"""
    t = t[:i9] + new9 + t[i10:]

    # Eq.(19) status: not PRIMARY_CONFLICT
    t = t.replace(
        "## Eq.(19) `PRINTED_PRIMARY` + `PRIMARY_CONFLICT`（分母幂次）",
        "## Eq.(19) `PRINTED_PRIMARY` + `STANDARD_AR_PSD_MODULUS_POWER_2_CONTROL`（分母幂次对照，非原文冲突）",
    )
    t = t.replace(
        "冲突：`EQ19_DENOMINATOR_POWER_NOT_RESOLVED`（印刷幂次 1；DSP 常规 PSD 幂次 2 仅作 CONTROL）。",
        "对照：`EQ19_DENOMINATOR_POWER_NOT_RESOLVED` → 记为需控制实验的实现歧义；**非 PRIMARY_CONFLICT**。\n"
        "- `PRINTED_EQ19_MODULUS_POWER_1`\n"
        "- `STANDARD_AR_PSD_MODULUS_POWER_2_CONTROL`（DSP 常规，非论文）",
    )
    # Eq.(24) status
    t = t.replace(
        "## Eq.(24) `PRINTED_PRIMARY` + `PRIMARY_CONFLICT`（求解器未指定）",
        "## Eq.(24) `PRINTED_PRIMARY` + `PAPER_DOES_NOT_SPECIFY_NUMERICAL_SOLVER`（实现细节未指定，非公式冲突）",
    )
    lock.write_text(t, encoding="utf-8")

    # 2) AR order preregistration
    (OUT / "AR_ORDER_PREREGISTRATION.md").write_text(
        f"""# AR_ORDER_PREREGISTRATION

UTC: {NOW}
基线 commit：4c3dbc9ae539bbb2e732c107097baaa6489f2a62

在任何 AR 数值结果出现前冻结。**禁止**根据 D(z) 或 MODE_ORDER_RECOVERY_RATE 调整 p。

## 打印（PRIMARY）

Liang Eq.(18)：p often set to (2/3)(2L+1)；仿真 L=5 → 2L+1=11。
记号冲突：`AR_SAMPLE_COUNT_NOTATION_CONFLICT` / `PRIMARY_NOTATION_CONFLICT`（2L+1 同时为阵元数与 y[i] 长度）。

## 预注册两个解释分支（4B 都跑）

### A. `PRINTED_LITERAL_CONTROL`

```latex
p=\\left\\lfloor\\frac{{2}}{{3}}(2L+1)\\right\\rfloor=\\left\\lfloor\\frac{{2\\times 11}}{{3}}\\right\\rfloor=7
```

仅“原文字面控制”，**不得**声称正确物理解释。

### B. `OUR_MOVING_SAMPLE_INTERPRETATION`

将 Eq.(17)/(18) 碰撞的 2L+1 解释为移动距离样本数 **N_r**：

```latex
p=\\left\\lfloor\\frac{{2}}{{3}}N_r\\right\\rfloor
```

明确：**我们的实现解释**，不得反写成论文原文。

## 规则

- 两分支均可测试
- **禁止**按最终深度 / mode-order 结果再选 p
- 结果报告必须分列 A/B
""",
        encoding="utf-8",
    )

    # 5) paper config final
    (OUT / "MMAR_PAPER_CONFIG_FINAL.md").write_text(
        f"""# MMAR_PAPER_CONFIG_FINAL

UTC: {NOW}
PRIMARY SHA256: `{SHA}`
唯一配置入口（4B 不得从 superseded V1/V2 继承参数）

| 项 | 值 | 状态 |
| --- | --- | --- |
| frequency | 350 Hz | PRINTED_PRIMARY |
| HLA N_elem | 2L+1=**11**（L=5） | PRINTED_PRIMARY |
| d | **λ/2**；参考声速 **主文未明确** | `PAPER_LAMBDA_REFERENCE_NOT_EXPLICIT` |
| HLA depth | 70 m | PRINTED_PRIMARY |
| z_s | 4 m / 50 m | PRINTED_PRIMARY |
| r0 | 5010 m（算法当未知） | PRINTED_PRIMARY |
| v | 2.5 m/s | PRINTED_PRIMARY |
| motion | 沿 beam 方向远离 | PRINTED_PRIMARY |
| look angle | θ̂=θ=θ | PRINTED_PRIMARY（ORACLE_LOOK_DIRECTION） |
| depth during aperture | fixed | PRINTED_PRIMARY |
| range span | **1990 m**（insufficient） / **4990 m**（sufficient） | PRINTED_PRIMARY |
| SNR | 20 / 5 / −5 dB | PRINTED_PRIMARY |
| field | KRAKEN | PRINTED_PRIMARY |
| bottom | referred to Yang2015 Ref.[7] | REFERENCED |
| MC | C0=500；correct \\|ẑ−z\\|≤5 m；90% CI Eq.30 | PRINTED_PRIMARY |
| δ | **0** = `ORACLE_OFFSET_ALIGNMENT` | 项目约定 |
| Δ regularizer | ~½ max\\|φ\\| | PRINTED_PRIMARY（Eq.23） |

## 非 Liang 主文明示

| 项 | 说明 |
| --- | --- |
| **Δr / Δt** | 主文本节**未给**时间/距离采样间隔 |
| 若用 Δt=1 s | 只能标 **`REF7_CONSISTENT_SAMPLING_ASSUMPTION`** |
| 空间无折叠 | Δr < π/k_max；见 `SPATIAL_SAMPLING_ALIAS_AUDIT.md` |
| AR p | 见 `AR_ORDER_PREREGISTRATION.md`（两分支） |
""",
        encoding="utf-8",
    )

    decision = "C2_4A_PRIMARY_LOCK_FINAL_PENDING_GPT_SPOTCHECK"
    (OUT / "FIX2_1_DECISION.json").write_text(
        json.dumps(
            {
                "stage": "R3-C2.4A-PRIMARY-FIX2.1",
                "decision": decision,
                "gpt_audit": "MINIMAL_SOURCE_FIDELITY_PATCH_REQUIRED → applied",
                "patches": [
                    "Eq9 PRINTED_PRIMARY restores sin b(X_m); BF_m only in DERIVED_IMPLEMENTATION",
                    "AR_ORDER_PREREGISTRATION.md p=7 CONTROL + floor(2 N_r/3) OUR_MOVING_SAMPLE_INTERPRETATION",
                    "V2 SUPERSEDED_BY_R3_C2_4A_PRIMARY_FIX2_1.md",
                    "Eq19/24 not PRIMARY_CONFLICT",
                    "MMAR_PAPER_CONFIG_FINAL.md",
                ],
                "not_done": ["AR", "MMAR", "FIELD", "4B", "E-STD", "MC", "P5"],
                "created_utc": NOW,
            },
            ensure_ascii=False,
            indent=2,
        ),
        encoding="utf-8",
    )
    (OUT / "FIX2_1_REPORT.md").write_text(
        f"""# R3-C2.4A-PRIMARY-FIX2.1 报告

UTC: {NOW}

## 判定

### `{decision}`

## 五项补丁

1. Eq.(9) `PRINTED_PRIMARY` 恢复原文 **sin b(X_m)**；BF_m 仅在 `DERIVED_EQ4_BEAM_FACTOR_INTERPRETATION`
2. **`AR_ORDER_PREREGISTRATION.md`**：`PRINTED_LITERAL_CONTROL` p=7；`OUR_MOVING_SAMPLE_INTERPRETATION` p=⌊2N_r/3⌋；禁止按 D(z) 调 p
3. V2 目录 **`SUPERSEDED_BY_R3_C2_4A_PRIMARY_FIX2_1.md`**
4. Eq.(19)→`STANDARD_AR_PSD_MODULUS_POWER_2_CONTROL`；Eq.(24)→`PAPER_DOES_NOT_SPECIFY_NUMERICAL_SOLVER`（均非 PRIMARY_CONFLICT）
5. **`MMAR_PAPER_CONFIG_FINAL.md`**（Δr/Δt 非主文明示）

## 停止

不写 AR、不进 R3-C2.4B。
""",
        encoding="utf-8",
    )
    (OUT / "GPT_SYNC.md").write_text(
        f"# R3-C2.4A-PRIMARY-FIX2.1\\n\\n**{decision}**\\n\\nEq9 原文锁+AR p 预注册+V2 superseded+状态分类+PAPER_CONFIG_FINAL。\\n",
        encoding="utf-8",
    )
    print("DECISION", decision)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
