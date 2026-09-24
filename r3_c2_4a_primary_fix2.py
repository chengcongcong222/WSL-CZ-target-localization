#!/usr/bin/env python3
"""R3-C2.4A-PRIMARY-FIX2: final lock integrity. No AR algorithm, no MMAR repro."""
from __future__ import annotations

import hashlib
import json
import re
from datetime import datetime, timezone
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parent
PDF = ROOT / "literature" / "liang2018_mmar_7824671.pdf"
OUT = ROOT / "results" / "R3_C2_AR_MMAR" / "R3_C2_4A_PRIMARY_FIX2"
OUT.mkdir(parents=True, exist_ok=True)
NOW = datetime.now(timezone.utc).isoformat()
SHA = hashlib.sha256(PDF.read_bytes()).hexdigest()

# Yang/Liang sim ranges from paper texts
K_FIG = (1.32, 1.50)  # rad/m from Liang Fig.3/4 axis
DR_YANG_REF7 = 2.5  # m — Yang/Ref.7 sampling (NOT Liang explicit)


def main() -> int:
    # ---- 1) SUPERSESSION ----
    v1 = ROOT / "results" / "R3_C2_AR_MMAR" / "R3_C2_4A_PRIMARY"
    (v1 / "SUPERSEDED_BY_R3_C2_4A_PRIMARY_FIX2.md").write_text(
        f"""# SUPERSEDED

UTC: {NOW}

`results/R3_C2_AR_MMAR/R3_C2_4A_PRIMARY/` 与 V1 锁定稿
→ **SUPERSEDED_BY_R3_C2_4A_PRIMARY_FIX2**

历史保留供审计，**不得作为后续实现输入**。
唯一公式锁：`R3_C2_4A_PRIMARY_FIX2/MMAR_EQUATION_LOCK_FINAL.md`
""",
        encoding="utf-8",
    )
    lock_py = ROOT / "r3_c2_4a_primary_lock.py"
    if lock_py.exists():
        t = lock_py.read_text(encoding="utf-8")
        if "SUPERSEDED_BY_R3_C2_4A_PRIMARY_FIX2" not in t:
            lock_py.write_text(
                '"""SUPERSEDED_BY_R3_C2_4A_PRIMARY_FIX2 — do not run; see r3_c2_4a_primary_fix2.py"""\n'
                "raise SystemExit('SUPERSEDED_BY_R3_C2_4A_PRIMARY_FIX2 — use MMAR_EQUATION_LOCK_FINAL.md')\n",
                encoding="utf-8",
            )
    fix_py = ROOT / "r3_c2_4a_primary_fix.py"
    if fix_py.exists():
        t = fix_py.read_text(encoding="utf-8")
        if "SUPERSEDED_BY_R3_C2_4A_PRIMARY_FIX2" not in t:
            fix_py.write_text(
                '"""SUPERSEDED_BY_R3_C2_4A_PRIMARY_FIX2 — do not run; see r3_c2_4a_primary_fix2.py"""\n'
                "raise SystemExit('SUPERSEDED_BY_R3_C2_4A_PRIMARY_FIX2 — use MMAR_EQUATION_LOCK_FINAL.md')\n",
                encoding="utf-8",
            )

    # ---- 2) FINAL EQUATION LOCK ----
    (OUT / "MMAR_EQUATION_LOCK_FINAL.md").write_text(
        T := f"""# MMAR_EQUATION_LOCK_FINAL

UTC: {NOW}
PRIMARY SHA256: `{SHA}`
SUPERSEDES: R3_C2_4A_PRIMARY (V1), R3_C2_4A_PRIMARY_FIX (V2)

状态标记：`PRINTED_PRIMARY` / `PRIMARY_CONFLICT` / `DERIVED_IMPLEMENTATION`

---

## Eq.(1) `PRINTED_PRIMARY`

```latex
p_l(r_i,z_r)=\\sqrt{{2\\pi}}\\,e^{{-j(\\pi/4)}}\\sum_{{m=1}}^{{M}}
\\phi_m(z_s)\\phi_m(z_r)
\\frac{{\\exp\\{{-j(k_m-j\\alpha_m)[r_i+ld\\sin\\theta_i]\\}}}}{{\\sqrt{{k_m[r_i+ld\\sin\\theta_i]}}}}
```

## Eq.(2) `PRINTED_PRIMARY`

```latex
p_l(r_i,z_r)=\\sum_{{m=1}}^{{M}} A_m\\exp\\{{-j(k_m-j\\alpha_m)[r_i+ld\\sin\\theta_i]\\}},
\\quad A_m=\\sqrt{{2\\pi}}e^{{-j\\pi/4}}\\frac{{\\phi_m(z_s)\\phi_m(z_r)}}{{\\sqrt{{k_m r_i}}}}
```

（far-field：r_i+ld\\sin\\theta_i\\approx r_i）

## Eq.(3) `PRINTED_PRIMARY`

```latex
B(\\hat\\theta_i)=\\frac{{1}}{{2L+1}}\\sum_{{l=-L}}^{{L}} e^{{j k l d \\sin\\hat\\theta_i}}\\,p_l(r_i,z_r)
```

## Eq.(4) `PRINTED_PRIMARY`

```latex
B(\\hat\\theta_i)=\\frac{{1}}{{2L+1}}\\sum_{{m=1}}^{{M}} A_m e^{{-jk_m r_i-\\alpha_m r_i}}
\\frac{{\\sin[(L+1/2)d X_m]}}{{\\sin[(d/2) X_m]}},
\\quad X_m=-(k_m-j\\alpha_m)\\sin\\theta_i+k\\sin\\hat\\theta_i
```

**BF_m（阵因子，Eq.4 原式）**：

```latex
BF_m=\\frac{{1}}{{2L+1}}\\cdot\\frac{{\\sin[(L+1/2)d X_m]}}{{\\sin[(d/2) X_m]}}
```

（分母书写必须为 `\\sin[(d/2)X_m]`，不得 `sin((d/2)X_m]`。）

## Eq.(5) `PRINTED_PRIMARY`

```latex
p(r;z_s,z_r)=\\int_0^\\infty g(k_r;z_s,z_r)J_0(k_r r)k_r\\,dk_r,
\\quad g(k_r;z_s,z_r)=\\int_0^\\infty p(r;z_s,z_r)J_0(k_r r)r\\,dr
```

## Eq.(6) `PRINTED_PRIMARY`

```latex
g(k_r;z_s,z_r)\\sim\\frac{{e^{{i\\pi/4}}}}{{\\sqrt{{2\\pi k_r}}}}\\int_{{-\\infty}}^{{+\\infty}} p(r;z_s,z_r)e^{{i k_r r}}\\sqrt{{r}}\\,dr,
\\quad k_r r\\gg 1
```

## Eq.(7) `PRINTED_PRIMARY`

```latex
g(k_r,z_r)=\\frac{{e^{{i\\pi/4}}}}{{\\sqrt{{2\\pi k_r}}}}\\int_{{r_0}}^{{r_0+R}} B(r)\\,e^{{i k_r r}}\\,S(r)\\,dr,
\\quad k_r r_0\\gg 1
```

## Eq.(8) `PRINTED_PRIMARY`

```latex
S(r)=\\langle |B(r)|^2\\rangle^{{-1/2}}
```

## Eq.(9) `PRINTED_PRIMARY`（积分核）+ `DERIVED_IMPLEMENTATION`（第二形式用 Eq.9 重推）

```latex
g(k_r,z_r)\\sim\\sum_{{m=1}}^{{M}}
\\frac{{\\phi_m(z_s)\\phi_m(z_r)}}{{\\sqrt{{k_r k_m}}}}\\,\\mathcal{{B}}_m
\\int_{{r_0}}^{{r_0+R}} e^{{j(k_r-k_m)r-\\alpha_m r}}\\,dr
=\\sum_{{m=1}}^{{M}} a_m\\frac{{\\phi_m(z_s)\\phi_m(z_r)}}{{k_r-k_m+j\\alpha_m}}
```

积分核 **exp[j(k_r−k_m)r − α_m r]**。

其中 \\(\\mathcal{{B}}_m\\) 在论文印刷中写作 **sin b(X_m)**（见下）。

## Eq.(10) `PRIMARY_CONFLICT`（印刷原样）

```latex
a_m=\\frac{{ e^{{j[(k_r-k_m)-\\alpha_m](r_0+R)}}-e^{{j[(k_r-k_m)-\\alpha_m]r_0}} }}{{j\\sqrt{{k_r k_m}}}}\\;\\sin b(X_m)
```

- **sin b(X_m) 为乘法因子，不在分母**（相对 V1 转录已修正）
- 仍冲突：α_m 印入 j[(k_r−k_m)−α_m]，与 Eq.(9) 积分核不符

## Eq.(11) `PRINTED_PRIMARY`

```latex
g(k_m,z_r)\\sim b_m\\phi_m(z_r)
```

## Eq.(12) `PRIMARY_CONFLICT`（印刷原样）

```latex
b_m=\\frac{{2 e^{{-\\alpha_0 r'}}}}{{\\alpha_m k_m}}\\sinh\\left(\\frac{{\\alpha_m R}}{{2}}\\right)\\phi_m(z_s)\\sin b(X_m),
\\quad r'=r_0+R/2
```

冲突：`e^{{-α_0 r'}}` 下标 α_0（应为 mode 的 α_m）。

## Eq.(13) `PRINTED_PRIMARY`

```latex
\\mathbf{{g}}=\\Phi\\cdot\\mathbf{{b}}
```

## Eq.(14) `PRINTED_PRIMARY`

```latex
\\mathbf{{g}}=[g(k_1,z_r),g(k_2,z_r),\\ldots,g(k_M,z_r)]^T
```

## Eq.(15) `PRINTED_PRIMARY`

```latex
\\Phi=\\mathrm{{diag}}([\\phi_1(z_r),\\phi_2(z_r),\\ldots,\\phi_M(z_r)])
```

## Eq.(16) `PRINTED_PRIMARY`

```latex
\\mathbf{{b}}=[b_1,b_2,\\ldots,b_M]^T
```

## Eq.(17) `PRINTED_PRIMARY` + `PRIMARY_CONFLICT`（样本记号）

```latex
y[i]=B(r_i)S(r_i),\\quad i=1,2,\\ldots,2L+1
```

冲突：`2L+1` 同时用作 HLA 阵元数与 AR 数据长度（见 AR_INTERPRETATION_GATE）。

## Eq.(18) `PRINTED_PRIMARY` + `PRIMARY_CONFLICT`（阶数）

```latex
y[i]=-\\sum_{{k=1}}^{{p}} a[k]y[i-k]+u[i],
\\quad p\\ \\mathrm{{often\\ set\\ to}}\\ (2/3)(2L+1)
```

## Eq.(19) `PRINTED_PRIMARY` + `PRIMARY_CONFLICT`（分母幂次）

```latex
P_{{AR}}(l)=\\frac{{\\sigma^2}}{{\\left|1+\\sum_{{k=1}}^{{p}}a[k]\\exp[-ilk]\\right|}}
```

冲突：`EQ19_DENOMINATOR_POWER_NOT_RESOLVED`（印刷幂次 1；DSP 常规 PSD 幂次 2 仅作 CONTROL）。

## Eq.(20) `PRINTED_PRIMARY`

```latex
D(z)=\\varphi(z)\\,\\mathbf{{b}}\\mathbf{{b}}^H\\,\\varphi^H(z)
```

## Eq.(21) `PRINTED_PRIMARY`

```latex
\\varphi(z)=[\\phi_1(z),\\phi_2(z),\\ldots,\\phi_M(z)]
```

## Eq.(22) `PRINTED_PRIMARY`

```latex
\\mathbf{{b}}=(\\Phi+U)^{{-1}}\\mathbf{{g}}
```

## Eq.(23) `PRINTED_PRIMARY`

```latex
U=\\mathrm{{diag}}\\left(\\left[\\frac{{\\Delta}}{{\\phi_1(z_r)}},\\frac{{\\Delta}}{{\\phi_2(z_r)}},\\ldots,\\frac{{\\Delta}}{{\\phi_M(z_r)}}\\right]\\right)
```

Δ：**on the order of one-half of the maximum value of the mode function**（`LIANG_PAPER_REGULARIZER`）。

## Eq.(24) `PRINTED_PRIMARY` + `PRIMARY_CONFLICT`（求解器未指定）

```latex
\\min_{{k_0}}(\\mathbf{{k}}-\\mathbf{{k}}_0)^H(\\mathbf{{k}}-\\mathbf{{k}}_0)
\\quad \\mathrm{{s.t.}}\\quad k_0(1)<k_0(2)<\\cdots<k_0(M_0)
```

## Eq.(25) `PRINTED_PRIMARY`

```latex
\\mathbf{{k}}=[k_1,k_2,\\ldots,k_{{M_0}}]^T
```

## Eq.(26) `PRINTED_PRIMARY`

```latex
\\mathbf{{k}}'=[k'_1,k'_2,\\ldots,k'_M]^T
```

## Eq.(27) `PRINTED_PRIMARY`

```latex
\\mathbf{{k}}_0\\subseteq\\mathbf{{k}}'
```

## Eq.(28) `PRINTED_PRIMARY`

```latex
\\mathrm{{SNR}}=10\\lg\\frac{{P_s}}{{P_n}}\\Big|_{{r=r_0}}
```

## Eq.(29) `PRINTED_PRIMARY`

```latex
P=\\frac{{C}}{{C_0}},\\quad C_0=500,
\\quad \\mathrm{{correct}}:\\ |\\hat{{z}}-z_{{\\mathrm{{true}}}}|\\le 5\\ \\mathrm{{m}}
```

## Eq.(30) `PRINTED_PRIMARY`

```latex
P\\pm 1.645\\sqrt{{\\frac{{P(1-P)}}{{C_0}}}}\\quad (90\\%\\ \\mathrm{{CI}})
```

---

## DERIVED_IMPLEMENTATION 汇总

| 代号 | 内容 | 来源 |
| --- | --- | --- |
| `DERIVED_EQ4_BEAM_FACTOR_INTERPRETATION` | 后文 `sin b(X_m)` 解释为 Eq.(4) 阵因子 BF_m；**论文未独立定义 sin b** | Eq.(4)→Eq.(7) 代数 |
| `FROM_EQ9_DERIVED_BM` | `b_m=\\frac{{2e^{{-\\alpha_m r'}}}}{{\\alpha_m k_m}}\\sinh(\\alpha_m R/2)\\phi_m(z_s)\\,BF_m` | Eq.(9) 积分 |
| `OUR_EXACT_SOLVER_FOR_EQ24` | （契约，本轮不跑）DP 有序子序列 | 我方实现 |
""",
        encoding="utf-8",
    )

    # static check for typo
    typo_hits = []
    for p in ROOT.rglob("*.md"):
        if "R3_C2_4A" not in str(p) and "FIX2" not in str(p):
            # still scan MMAR-related
            if "MMAR" not in p.name and "liang" not in p.name.lower():
                continue
        try:
            t = p.read_text(encoding="utf-8", errors="replace")
        except Exception:
            continue
        if "sin((d/2)X_m]" in t or "sin((d/2)X_m]" in t:
            typo_hits.append(str(p))
        if "sin((d/2)X_m]" in t:
            pass
    # also fix V2 lock if typo present in our generated files
    for p in [
        ROOT / "results/R3_C2_AR_MMAR/R3_C2_4A_PRIMARY_FIX/MMAR_EQUATION_LOCK_V2.md",
        ROOT / "results/R3_C2_AR_MMAR/R3_C2_4A_PRIMARY_FIX/EQ7_TO_EQ12_ALGEBRA_AUDIT_V2.md",
        ROOT / "r3_c2_4a_primary_fix.py",
    ]:
        if p.exists():
            t = p.read_text(encoding="utf-8", errors="replace")
            t2 = t.replace("sin((d/2)X_m]", "sin[(d/2)X_m]")
            if t2 != t:
                p.write_text(t2, encoding="utf-8")
                typo_hits.append(f"FIXED:{p}")

    (OUT / "SUPERSESSION_AUDIT.md").write_text(
        f"""# SUPERSESSION_AUDIT

UTC: {NOW}

| 工件 | 状态 |
| --- | --- |
| r3_c2_4a_primary_lock.py | SUPERSEDED_BY_R3_C2_4A_PRIMARY_FIX2（重跑即退出） |
| r3_c2_4a_primary_fix.py | SUPERSEDED_BY_R3_C2_4A_PRIMARY_FIX2（重跑即退出） |
| results/.../R3_C2_4A_PRIMARY/ | SUPERSEDED；仅审计 |
| results/.../R3_C2_4A_PRIMARY_FIX/ | SUPERSEDED_BY_FIX2；仅审计 |
| **MMAR_EQUATION_LOCK_FINAL.md** | **唯一实现输入** |

BF_m 括号静态检查命中/修复：{typo_hits if typo_hits else "none（生成物中已用 \\sin[(d/2)X_m]）"}
""",
        encoding="utf-8",
    )

    # ---- 5) spatial sampling / alias audit ----
    rows = []
    for kmin, kmax, tag in [(K_FIG[0], K_FIG[1], "LIANG_FIG_K"), (1.32, 1.50, "dup")]:
        if tag == "dup":
            continue
        for dr, dr_tag in [
            (DR_YANG_REF7, "REF7_CONSISTENT_SAMPLING_ASSUMPTION dr=2.5m"),
            (0.5, "control dr=0.5m"),
            (0.3, "control dr=0.3m"),
        ]:
            for k in (kmin, (kmin + kmax) / 2, kmax):
                w = k * dr
                folds = bool(abs(w) > np.pi)
                w_wrapped = (w + np.pi) % (2 * np.pi) - np.pi
                rows.append(
                    {
                        "k_rad_per_m": k,
                        "dr_m": dr,
                        "k_dr_rad": w,
                        "k_dr_over_pi": w / np.pi,
                        "folds_if_spatial_Nyquist": folds,
                        "omega_s_wrapped": w_wrapped,
                        "dr_tag": dr_tag,
                    }
                )
        # no-fold control: k_max * dr < pi
        dr_nf = 0.9 * np.pi / kmax
        rows.append(
            {
                "k_rad_per_m": kmax,
                "dr_m": dr_nf,
                "k_dr_rad": kmax * dr_nf,
                "k_dr_over_pi": kmax * dr_nf / np.pi,
                "folds_if_spatial_Nyquist": False,
                "omega_s_wrapped": np.nan,
                "dr_tag": f"NO_FOLD_CONTROL dr=0.9*pi/kmax={dr_nf:.4f}m",
            }
        )
    import pandas as pd

    pd.DataFrame(rows).to_csv(OUT / "spatial_sampling_alias_table.csv", index=False)

    w_ref = 1.50 * DR_YANG_REF7
    folds_ref = abs(w_ref) > np.pi
    (OUT / "SPATIAL_SAMPLING_ALIAS_AUDIT.md").write_text(
        f"""# SPATIAL_SAMPLING_ALIAS_AUDIT

UTC: {NOW}

## 推导（只做解析，不评 MMAR 成败）

运动距离样本：r_i = r_0 + i Δr

模态相位因子（Eq.2 空间部分）：y_i ∝ exp[−j k_m r_i]

⇒ 相邻样本：**y_i ∝ exp[−j i k_m Δr]**

AR 归一化角频率与物理水平波数：

$$
\\omega_s \\equiv \\pm k_m\\,\\Delta r \\pmod{{2\\pi}}
$$

存在 **符号** 与 **分支** 选择；须由可控复指数单元测试 + 论文谱横轴（Fig.3/4 约 **k=1.32–1.50 rad/m**）共同确定。

## 历史 Yang/Ref.7 采样 Δr≈2.5 m 解析检查

| k (rad/m) | kΔr | (kΔr)/π | 空间折叠？ |
| --- | --- | --- | --- |
| 1.32 | {1.32*2.5:.3f} | {1.32*2.5/np.pi:.3f} | {'YES' if 1.32*2.5>np.pi else 'NO'} |
| 1.50 | {w_ref:.3f} | {w_ref/np.pi:.3f} | {'YES' if folds_ref else 'NO'} |

**k_max Δr = {w_ref:.3f} > π** ⇒ 按空间 Nyquist 会发生 **波数折叠/混叠**（不能用 Δr=2.5 m 无折叠覆盖 1.5 rad/m）。

## 无折叠 control

要求 k_max Δr < π ⇒ Δr < π/k_max ≈ **{np.pi/1.50:.4f} m**

预注册 control：**Δr = 0.9 π / k_max ≈ {0.9*np.pi/1.50:.4f} m**（表中 NO_FOLD_CONTROL）。

## 边界

本审计 **不** 评价 MMAR 成败；只冻结实现解释与采样约束。
`PAPER_RANGE_SAMPLE_INTERVAL_NOT_EXPLICIT` 仍成立；Δt=1 s 仅 `REF7_CONSISTENT_SAMPLING_ASSUMPTION`。
""",
        encoding="utf-8",
    )

    # ---- 6) AR interpretation gate ----
    (OUT / "AR_INTERPRETATION_GATE.md").write_text(
        f"""# AR_INTERPRETATION_GATE

UTC: {NOW}

## 记号

| 量 | 符号 | 规则 |
| --- | --- | --- |
| HLA 阵元数 | **2L+1** | 仅阵元 |
| 移动距离样本数 | **N_r** | 独立符号；禁止再用 2L+1 |
| AR 阶数 | p | 打印 (2/3)(2L+1) → `PRIMARY_NOTATION_CONFLICT`；若用 N_r 规则 → **我们的实现解释**，不得反写成论文 |

## 采样

- Δt=1 s **不是** Liang 明文
- 借 Yang/Ref.7 → **`REF7_CONSISTENT_SAMPLING_ASSUMPTION`**
- 空间无折叠 control：Δr < π/k_max（见 SPATIAL_SAMPLING_ALIAS_AUDIT）

## Eq.(19)

- `PRINTED_EQ19_MODULUS_POWER_1`
- `STANDARD_AR_PSD_MODULUS_POWER_2_CONTROL`
- **先** 可控复指数验证峰**位置**；**禁止**按 D(z) 深度结果择优

## l → k_r

`AR_NORMALIZED_SPECTRAL_AXIS_MAPPING_NOT_EXPLICIT`

待验证假设：ω_s ≡ ± k_m Δr (mod 2π)；符号/分支由 unit test + paper scale 决定，**本轮不冻结为实现公式**。

## Eq.(24)

- `PAPER_DOES_NOT_SPECIFY_NUMERICAL_SOLVER`
- DP → **`OUR_EXACT_SOLVER_FOR_EQ24`**（本轮只写契约）
- 核心指标：**MODE_ORDER_RECOVERY_RATE**（非 Rayleigh unique）
""",
        encoding="utf-8",
    )

    # ---- contract final ----
    (OUT / "MMAR_IMPLEMENTATION_CONTRACT_FINAL.md").write_text(
        f"""# MMAR_IMPLEMENTATION_CONTRACT_FINAL

UTC: {NOW}

唯一公式锁：`MMAR_EQUATION_LOCK_FINAL.md`

## 输入链

HLA B(θ̂) [Eq3–4, BF_m] → S(r) [Eq8] → y[i]=B S [Eq17, **N_r** 样本]
→ AR P_AR [Eq19 + 幂次对照] → k̂ → Hankel g(k̂) [Eq7]
→ ORDERED_SUBSET Eq24（`OUR_EXACT_SOLVER_FOR_EQ24`）
→ b=(Φ+U)^{-1}g [Eq22–23, Δ~½max|φ|] → D(z) [Eq20]

## 实现闭式

`FROM_EQ9_DERIVED_BM` + `DERIVED_EQ4_BEAM_FACTOR_INTERPRETATION`（非“论文定义 sinb”）

## 采样

`REF7_CONSISTENT_SAMPLING_ASSUMPTION` 或无折叠 control；须报 `SPATIAL_SAMPLING_ALIAS_AUDIT` 结果。

## δ

`ORACLE_OFFSET_ALIGNMENT`（δ=0）；`MMAR_INHERITS_OFFSET_RANGE_CAVEAT`

## 本轮不做

AR 正式算法 / MMAR 深度复现 / FIELD / E-STD / MC / P5
""",
        encoding="utf-8",
    )

    decision = "C2_4A_PRIMARY_LOCK_V2_INTEGRITY_FIXED_PENDING_GPT_AUDIT"
    (OUT / "FIX2_DECISION.json").write_text(
        json.dumps(
            {
                "stage": "R3-C2.4A-PRIMARY-FIX2",
                "fix2_decision": decision,
                "why": "V1/V2 已 SUPERSEDED；唯一 Eq.(1)–(30) FINAL 锁；BF_m 括号/ sinb 解释、AR 采样与波数折叠审计、实现契约已清；待 GPT 审计",
                "superseded": [
                    "r3_c2_4a_primary_lock.py",
                    "r3_c2_4a_primary_fix.py",
                    "R3_C2_4A_PRIMARY/",
                    "R3_C2_4A_PRIMARY_FIX/",
                ],
                "unique_lock": "MMAR_EQUATION_LOCK_FINAL.md",
                "alias_finding": "k_max*dr(2.5m)>pi at k=1.5 rad/m; NO_FOLD_CONTROL dr~1.88m",
                "not_done": ["AR algorithm", "MMAR depth repro", "FIELD", "E-STD", "MC", "P5", "R3-C2.4B"],
                "created_utc": NOW,
            },
            ensure_ascii=False,
            indent=2,
        ),
        encoding="utf-8",
    )
    (OUT / "FIX2_REPORT.md").write_text(
        f"""# R3-C2.4A-PRIMARY-FIX2 报告

UTC: {NOW}

## 判定

### `{decision}`

## 完成

1. V1/V2 **SUPERSEDED_BY_R3_C2_4A_PRIMARY_FIX2**（旧脚本重跑即退出）
2. **MMAR_EQUATION_LOCK_FINAL.md**：Eq.(1)–(30) 全式；状态 PRINTED_PRIMARY / PRIMARY_CONFLICT / DERIVED_IMPLEMENTATION
3. BF_m 分母 `\\sin[(d/2)X_m]`；sin b → `DERIVED_EQ4_BEAM_FACTOR_INTERPRETATION`
4. **FROM_EQ9_DERIVED_BM** 保留；Eq10/Eq12 冲突保留
5. **SPATIAL_SAMPLING_ALIAS_AUDIT**：y_i∝e^{{-ji k_m Δr}}；ω_s≡±k_m Δr (mod 2π)；Δr=2.5 m @ k=1.5 **折叠**；NO_FOLD_CONTROL Δr≈{0.9*np.pi/1.5:.3f} m
6. **AR_INTERPRETATION_GATE**：N_r 独立；Δt 仅 REF7 假设；Eq19 峰位测试先行；Eq24 求解器边界
7. 实现契约 FINAL

## 停止

不写 AR、不跑 MMAR/FIELD/E-STD/MC/P5；**不进入 R3-C2.4B**。
""",
        encoding="utf-8",
    )
    (OUT / "GPT_SYNC.md").write_text(
        f"# R3-C2.4A-PRIMARY-FIX2\\n\\n**{decision}**\\n\\n唯一锁=MMAR_EQUATION_LOCK_FINAL.md；kΔr@2.5m 折叠已审计。\\n待 GPT 审计后再开 4B。\\n",
        encoding="utf-8",
    )
    print("DECISION", decision)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
