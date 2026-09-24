#!/usr/bin/env python3
"""R3-C2.4A-PRIMARY-FIX: correct lock artifacts only. No AR code, no repro."""
from __future__ import annotations

import hashlib
import json
from datetime import datetime, timezone
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parent
PDF = ROOT / "literature" / "liang2018_mmar_7824671.pdf"
OUT = ROOT / "results" / "R3_C2_AR_MMAR" / "R3_C2_4A_PRIMARY_FIX"
OUT.mkdir(parents=True, exist_ok=True)
NOW = datetime.now(timezone.utc).isoformat()
SHA = hashlib.sha256(PDF.read_bytes()).hexdigest()

CONFLICTS = [
    "AR_SAMPLE_COUNT_NOTATION_CONFLICT",
    "PAPER_RANGE_SAMPLE_INTERVAL_NOT_EXPLICIT",
    "EQ10_PRINTED_FORM_INTERNAL_INCONSISTENCY",
    "EQ12_ALPHA_SUBSCRIPT_CONFLICT",
    "SINB_NOTATION_UNDEFINED_USE_EQ4_BEAM_FACTOR",
    "EQ19_DENOMINATOR_POWER_NOT_RESOLVED",
    "AR_NORMALIZED_SPECTRAL_AXIS_MAPPING_NOT_EXPLICIT",
    "PAPER_DOES_NOT_SPECIFY_NUMERICAL_SOLVER",
    "PAPER_LAMBDA_REFERENCE_NOT_EXPLICIT",
    "MMAR_INHERITS_OFFSET_RANGE_CAVEAT",
]


def T(s: str) -> str:
    return s.replace("__NOW__", NOW).replace("__SHA__", SHA)


def main() -> int:
    # ---- Eq lock V2 ----
    (OUT / "MMAR_EQUATION_LOCK_V2.md").write_text(
        T(
            r"""# MMAR_EQUATION_LOCK_V2

UTC: __NOW__
PRIMARY SHA256: `__SHA__`

## 修正记录（相对 V1）

| 项 | V1 错误 | V2 正确 |
| --- | --- | --- |
| Eq.(10) | sin b(X_m) 误放在**分母** | **乘在分式之后** |
| sin b(X_m) | 似为未定义函数 | `SINB_NOTATION_UNDEFINED_USE_EQ4_BEAM_FACTOR` |
| Eq.(29)(30) | 文字占位 | 按 PDF 补全 |
| HLA d | d_m=0.5 m | d=λ/2，数值未由正文钉死 |

## Eq.(1)–(4) HLA

```latex
(1) p_l(r_i,z_r)=\sqrt{2\pi}e^{-j\pi/4}\sum_m \phi_m(z_s)\phi_m(z_r)\frac{\exp\{-j(k_m-j\alpha_m)[r_i+ld\sin\theta_i]\}}{\sqrt{k_m[r_i+ld\sin\theta_i]}}
(2) p_l=\sum_m A_m\exp\{-j(k_m-j\alpha_m)[r_i+ld\sin\theta_i]\},\ A_m=\sqrt{2\pi}e^{-j\pi/4}\phi_m(z_s)\phi_m(z_r)/\sqrt{k_m r_i}
(3) B(\hat\theta_i)=\frac{1}{2L+1}\sum_{l=-L}^{L}e^{j k l d\sin\hat\theta_i}p_l(r_i,z_r)
(4) B=\frac{1}{2L+1}\sum_m A_m e^{-jk_m r_i-\alpha_m r_i}\underbrace{\frac{\sin[(L+1/2)dX_m]}{\sin((d/2)X_m)}}_{\mathrm{array\ factor\ (Eq.4)}},\
X_m=-(k_m-j\alpha_m)\sin\theta_i+k\sin\hat\theta_i
```

**BF_m（实现用，源自 Eq.4）**：

```latex
BF_m=\frac{1}{2L+1}\cdot\frac{\sin[(L+1/2)dX_m]}{\sin((d/2)X_m]}
```

论文后文的 `sin b(X_m)` **无独立定义** → 使用 BF_m，禁止自造 `sinb()`。

## Eq.(5)–(8)

```latex
(5) Hankel pair p \leftrightarrow g(k_r)
(6) g(k_r)\sim\frac{e^{i\pi/4}}{\sqrt{2\pi k_r}}\int p e^{ik_r r}\sqrt{r}\,dr
(7) g(k_r,z_r)=\frac{e^{i\pi/4}}{\sqrt{2\pi k_r}}\int_{r_0}^{r_0+R}B(r)e^{ik_r r}S(r)\,dr
(8) S(r)=\langle|B(r)|^2\rangle^{-1/2}
```

## Eq.(9)（积分形式 — 实现主依据）

```latex
g\sim\sum_m \frac{\phi_m(z_s)\phi_m(z_r)}{\sqrt{k_r k_m}}\,BF_m
\int_{r_0}^{r_0+R} e^{j(k_r-k_m)r-\alpha_m r}\,dr
=\sum_m a_m\frac{\phi_m(z_s)\phi_m(z_r)}{k_r-k_m+j\alpha_m}
```

积分核：**exp[j(k_r−k_m)r − α_m r]**（α 在实部衰减项，不在 j[…] 内）。

## Eq.(10) — **印刷原样（乘法因子，非分母）**

```latex
a_m=\frac{ e^{j[(k_r-k_m)-\alpha_m](r_0+R)}-e^{j[(k_r-k_m)-\alpha_m]r_0} }{j\sqrt{k_r k_m}}\;\sin b(X_m)
```

即：

```latex
a_m=\Big[\cdots\Big]/\big(j\sqrt{k_rk_m}\big)\ \times\ \sin b(X_m)
```

**不是** `/ [j√(k_r k_m) sin b(X_m)]`。

仍保留 `EQ10_PRINTED_FORM_INTERNAL_INCONSISTENCY`：α_m 被印进 j[(k_r−k_m)−α_m]。

## Eq.(11)–(16)

```latex
(11) g(k_m,z_r)\sim b_m\phi_m(z_r)
(12) PRINTED: b_m=\frac{2e^{-\alpha_0 r'}}{\alpha_m k_m}\sinh(\alpha_m R/2)\phi_m(z_s)\sin b(X_m),\ r'=r_0+R/2
(13)–(16) g=\Phi b,\ \Phi=\mathrm{diag}(\phi_m(z_r)),\ b=[b_m]
```

Eq.(12) `e^{-α_0 r'}` → `EQ12_ALPHA_SUBSCRIPT_CONFLICT`。

## FROM_EQ9_DERIVED_BM（实现闭式）

令 Δ_k=k_r−k_m，

```latex
I_m=\int_{r_0}^{r_0+R}e^{(j\Delta_k-\alpha_m)r}dr
=\frac{e^{(j\Delta_k-\alpha_m)(r_0+R)}-e^{(j\Delta_k-\alpha_m)r_0}}{j\Delta_k-\alpha_m}
```

k_r=k_m 时：

```latex
I_m=\frac{2e^{-\alpha_m r'}}{\alpha_m}\sinh\left(\frac{\alpha_m R}{2}\right),\quad r'=r_0+R/2
```

故：

```latex
b_m^{\mathrm{impl}}=\frac{2e^{-\alpha_m r'}}{\alpha_m k_m}\sinh\left(\frac{\alpha_m R}{2}\right)\phi_m(z_s)\,BF_m
```

标 **`FROM_EQ9_DERIVED_BM`**；**不用 printed α_0**。实现 α→0 稳定极限：sinh(x)/x 形式。

## Eq.(17)–(19) AR

```latex
(17) y[i]=B(r_i)S(r_i),\ i=1,\ldots,2L+1
(18) y[i]=-\sum_{k=1}^{p}a[k]y[i-k]+u[i],\ p\ \mathrm{often}\ (2/3)(2L+1)
(19) PRINTED: P_{AR}(l)=\frac{\sigma^2}{\left|1+\sum_{k=1}^{p}a[k]e^{-ilk}\right|}
```

`EQ19_DENOMINATOR_POWER_NOT_RESOLVED`：印刷 **|A(l)|^1**；DSP 常规 **|A(l)|^{-2}** 仅作 STANDARD_AR_PSD_CONTROL，非原文。

## Eq.(20)–(23) 深度函数

```latex
(20) D(z)=\varphi(z)\,b b^H\,\varphi^H(z)
(21) \varphi(z)=[\phi_1(z),\ldots,\phi_M(z)]
(22) b=(\Phi+U)^{-1}g
(23) U=\mathrm{diag}([\Delta/\phi_1(z_r),\ldots,\Delta/\phi_M(z_r)]),\
\Delta\ \text{on the order of one-half of max mode function}
```

**LIANG_PAPER_REGULARIZER Δ≈½ max|φ|**（≠ Yang 0.1）。

## Eq.(24)–(27) ordered matching

```latex
(24) \min (k-k_0)^H(k-k_0)\ \mathrm{s.t.}\ k_0(1)<\cdots<k_0(M_0)
(25) k=[k_1,\ldots,k_{M_0}]^T
(26) k'=[k'_1,\ldots,k'_M]^T
(27) k_0\subseteq k'
```

`PAPER_DOES_NOT_SPECIFY_NUMERICAL_SOLVER`；若自写 DP：**`OUR_EXACT_SOLVER_FOR_EQ24`**。

## Eq.(28)–(30)

```latex
(28) SNR=10\lg(P_s/P_n)|_{r=r_0}
(29) P=C/C_0,\quad C_0=500,\quad \mathrm{correct}:|\hat z-z_\mathrm{true}|\le 5\ \mathrm{m}
(30) 90\%\ \mathrm{CI}: \ P \pm 1.645\sqrt{P(1-P)/C_0}
```
"""
        ),
        encoding="utf-8",
    )

    # equation table V2
    rows = []
    for i in range(1, 31):
        rows.append(
            {
                "equation_id": f"Eq.({i})",
                "status": "LOCKED_V2",
                "sinb_note": "USE_EQ4_BEAM_FACTOR" if i in (9, 10, 12) else "",
                "implementation": (
                    "FROM_EQ9_DERIVED_BM"
                    if i in (11, 12)
                    else (
                        "PRINTED_EQ19 + STANDARD_AR_PSD_CONTROL"
                        if i == 19
                        else (
                            "ORDERED_SUBSET + OUR_EXACT_SOLVER_FOR_EQ24"
                            if i in (24, 25, 26, 27)
                            else ("P=C/C0" if i == 29 else ("90%CI" if i == 30 else "AS_LOCKED"))
                        )
                    )
                ),
            }
        )
    pd.DataFrame(rows).to_csv(OUT / "MMAR_EQUATION_TABLE_V2.csv", index=False, encoding="utf-8-sig")

    # algebra audit V2
    (OUT / "EQ7_TO_EQ12_ALGEBRA_AUDIT_V2.md").write_text(
        T(
            r"""# EQ7_TO_EQ12_ALGEBRA_AUDIT_V2

UTC: __NOW__

## Eq.(9) 积分核

exp[j(k_r−k_m)r − α_m r] → 分母 j(k_r−k_m)−α_m = j(k_r−k_m+jα_m)

与印刷第二形式 (k_r−k_m+jα_m) **一致**。

## Eq.(10) 转录修正

**正确印刷**：sin b(X_m) **乘在** /[j√(k_r k_m)] **之后**，**不在分母**。

```latex
a_m=\frac{e^{j[(k_r-k_m)-\alpha_m](r_0+R)}-e^{j[(k_r-k_m)-\alpha_m]r_0}}{j\sqrt{k_r k_m}}\;\sin b(X_m)
```

## 仍保留内部不一致

`EQ10_PRINTED_FORM_INTERNAL_INCONSISTENCY`：指数里 α 进 j[(k_r−k_m)−α_m]=j(k_r−k_m)−jα_m，
与积分核系数 j(k_r−k_m)−α_m 不符。

→ 实现：**FROM_EQ9_DERIVED_BM** / Eq.(9) 积分。

## Eq.(12)

e^{−α_0 r'} → `EQ12_ALPHA_SUBSCRIPT_CONFLICT`；实现用 e^{−α_m r'}。

## sin b(X_m)

`SINB_NOTATION_UNDEFINED_USE_EQ4_BEAM_FACTOR`：

```latex
BF_m=\frac{1}{2L+1}\frac{\sin[(L+1/2)dX_m]}{\sin((d/2)X_m]}
```
"""
        ),
        encoding="utf-8",
    )

    # conflict register
    pd.DataFrame(
        [
            ("AR_SAMPLE_COUNT_NOTATION_CONFLICT", "Eq17 i=1..2L+1 vs motion samples; L=5→11 elements", "OPEN"),
            ("PAPER_RANGE_SAMPLE_INTERVAL_NOT_EXPLICIT", "v=2.5, R=1990/4990 given; Δt not in this section", "OPEN; if use 1s label REF7_CONSISTENT_SAMPLING_ASSUMPTION"),
            ("EQ10_PRINTED_FORM_INTERNAL_INCONSISTENCY", "α inside j[(kr-km)-α]", "OPEN; impl FROM_EQ9_DERIVED_BM"),
            ("EQ12_ALPHA_SUBSCRIPT_CONFLICT", "e^{-α0 r'} printed", "OPEN; impl α_m"),
            ("SINB_NOTATION_UNDEFINED_USE_EQ4_BEAM_FACTOR", "sin b(X_m) undefined; use Eq.4 array factor BF_m", "OPEN"),
            ("EQ19_DENOMINATOR_POWER_NOT_RESOLVED", "printed |A|^{-1}; standard AR PSD |A|^{-2}", "OPEN; PRINTED + STANDARD_AR_PSD_CONTROL; peak-position unit test first"),
            ("AR_NORMALIZED_SPECTRAL_AXIS_MAPPING_NOT_EXPLICIT", "l→k_r not stated", "OPEN; test l=±k_r Δr via controlled exp + paper scale"),
            ("PAPER_DOES_NOT_SPECIFY_NUMERICAL_SOLVER", "Eq24 solver not named", "OPEN; OUR_EXACT_SOLVER_FOR_EQ24 if DP"),
            ("PAPER_LAMBDA_REFERENCE_NOT_EXPLICIT", "d=λ/2 but c_ref for λ not stated", "OPEN; d=NOT_FIXED_BY_TEXT; later d=c_ref/(2f) + sensitivity"),
            ("MMAR_INHERITS_OFFSET_RANGE_CAVEAT", "no absolute range claim vs Yang2018 Erratum", "OPEN; δ=0 ORACLE_OFFSET_ALIGNMENT"),
        ],
        columns=["id", "evidence", "status"],
    ).to_csv(OUT / "MMAR_IMPLEMENTATION_CONFLICT_REGISTER.csv", index=False, encoding="utf-8-sig")

    (OUT / "MMAR_AR_IMPLEMENTATION_CONTRACT.md").write_text(
        T(
            r"""# MMAR_AR_IMPLEMENTATION_CONTRACT

UTC: __NOW__

## AR 输入

y[i] = B(r_i) S(r_i)  （HLA beam × Eq.8 shading）

不是 raw hydrophone / AR peak / Hankel spectrum。

## 样本与阶数

- `AR_SAMPLE_COUNT_NOTATION_CONFLICT`：印刷 2L+1 与运动序列冲突
- p often (2/3)(2L+1)
- `PAPER_RANGE_SAMPLE_INTERVAL_NOT_EXPLICIT`：Δt 未在 Liang 本节写明
  - 若用 1 s → **`REF7_CONSISTENT_SAMPLING_ASSUMPTION`**（不得写 LIANG_EXPLICIT_DT）

## 谱

- PRINTED_EQ19_MODULUS_POWER_1
- STANDARD_AR_PSD_CONTROL_MODULUS_POWER_2（对照，非原文）
- 先做 controlled complex-exp：两式 **峰位置** 是否一致（不得按深度择优）

## 峰 → k̂ → 幅度

- AR 峰位 → k̂_m
- 幅度 **只** 从 Hankel g(k̂_m)
- 禁止 AR peak height → amplitude

## l → k_r

`AR_NORMALIZED_SPECTRAL_AXIS_MAPPING_NOT_EXPLICIT`

允许待验证假设：l = ± k_r Δr (mod 2π)，须 unit test + paper scale 决定符号/尺度。

## Eq.24 求解

`PAPER_DOES_NOT_SPECIFY_NUMERICAL_SOLVER`；DP → `OUR_EXACT_SOLVER_FOR_EQ24`。

## b_m

`FROM_EQ9_DERIVED_BM` + BF_m（Eq.4）。

## δ

`ORACLE_OFFSET_ALIGNMENT`（δ=0）。
"""
        ),
        encoding="utf-8",
    )

    pd.DataFrame(
        [
            ("A_PAPER_EQUIV", 11, "lambda/2", "NOT_FIXED_BY_TEXT", "PAPER_LAMBDA_REFERENCE_NOT_EXPLICIT; later d=c_ref/(2f)"),
            ("B_PROJECT_14M", 8, 2.0, 14.0, "8 el × 2 m = 14 m aperture"),
            ("C_4EL_SCAN_L5", 4, "L/3", "scan", "CONDITION_ENHANCEMENT_SCAN L=5m"),
            ("C_4EL_SCAN_L10", 4, "L/3", "scan", "L=10m"),
            ("C_4EL_SCAN_L20", 4, "L/3", "scan", "L=20m"),
            ("C_4EL_SCAN_L40", 4, "L/3", "scan", "L=40m"),
            ("C_4EL_SCAN_L80", 4, "L/3", "scan", "L=80m"),
        ],
        columns=["config", "N_elem", "d_expression", "d_m", "note"],
    ).to_csv(OUT / "MMAR_HLA_MIGRATION_TABLE_V2.csv", index=False, encoding="utf-8-sig")

    (OUT / "MMAR_PAPER_FIDELITY_TARGETS_V2.md").write_text(
        T(
            """# MMAR_PAPER_FIDELITY_TARGETS_V2

UTC: __NOW__

## FIG3_FIG4_CAPTION_PANEL_LABEL_TYPO

图注后半句重复 (a,c,e) 与正文矛盾。

**以正文+图布局为准**：

- SAB：(a), (c), (e) = SNR 20 / 5 / −5
- MMAR：(b), (d), (f) = SNR 20 / 5 / −5

## Targets（下阶段复现用，本轮不跑）

| ID | 内容 |
| --- | --- |
| Fig3 | sufficient 4990 m spectra SAB vs MMAR |
| Fig4 | insufficient 1990 m spectra |
| Fig5 | sufficient span depth |
| Fig6 | insufficient span depth |
| Fig7 | P=C/C0 vs span; correct if |ẑ−z|≤5 m; C0=500; 90% CI Eq.30 |

分层：R0 无噪/高SNR → R1 20dB → R2 5dB → R3 −5dB；MC 后置。
"""
        ),
        encoding="utf-8",
    )

    decision = "C2_4A_PRIMARY_METHOD_LOCKED_WITH_NOTATION_CONFLICTS_V2"
    (OUT / "R3_C2_4A_PRIMARY_FIX_DECISION.json").write_text(
        json.dumps(
            {
                "stage": "R3-C2.4A-PRIMARY-FIX",
                "rc3_c2_4a_primary_fix_decision": decision,
                "why": "锁定稿转录错误已修正（Eq10 乘法因子；Eq29/30 补全；sinb→BF_m；HLA d=λ/2）；冲突/未决项登记完整；≠方法失败",
                "fixes": [
                    "Eq10 sinb multiplier not denominator",
                    "SINB_NOTATION_UNDEFINED_USE_EQ4_BEAM_FACTOR",
                    "FROM_EQ9_DERIVED_BM written",
                    "EQ19_DENOMINATOR_POWER_NOT_RESOLVED",
                    "PAPER_RANGE_SAMPLE_INTERVAL_NOT_EXPLICIT",
                    "Eq29 P=C/C0; Eq30 90% CI",
                    "FIG3_FIG4_CAPTION_PANEL_LABEL_TYPO",
                    "HLA d=lambda/2 NOT_FIXED_BY_TEXT",
                ],
                "conflicts": CONFLICTS,
                "yang_route": "YANG_ROUTE_UNDECIDED",
                "mmar_route": "RC3-C2-AR/MMAR INDEPENDENT CANDIDATE",
                "not_done": ["AR code", "MMAR", "FIELD", "paper repro", "E-STD", "MC", "P5"],
                "created_utc": NOW,
            },
            ensure_ascii=False,
            indent=2,
        ),
        encoding="utf-8",
    )

    (OUT / "R3_C2_4A_PRIMARY_FIX_REPORT.md").write_text(
        T(
            """# R3-C2.4A-PRIMARY-FIX 报告

UTC: __NOW__

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
"""
        ),
        encoding="utf-8",
    )
    (OUT / "R3_C2_4A_PRIMARY_FIX_GPT_SYNC.md").write_text(
        T(f"# R3-C2.4A-PRIMARY-FIX\\n\\n**{decision}**\\n\\nEq10 乘法因子已改；10 项冲突已登记。\\n"),
        encoding="utf-8",
    )
    print("DECISION", decision)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
