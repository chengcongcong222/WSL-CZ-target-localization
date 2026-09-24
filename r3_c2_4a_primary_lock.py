#!/usr/bin/env python3
"""R3-C2.4A-PRIMARY: lock Liang2018 MMAR formulas + internal consistency audit. No AR code."""
from __future__ import annotations

import hashlib
import json
from datetime import datetime, timezone
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parent
PDF = ROOT / "literature" / "liang2018_mmar_7824671.pdf"
OUT = ROOT / "results" / "R3_C2_AR_MMAR" / "R3_C2_4A_PRIMARY"
OUT.mkdir(parents=True, exist_ok=True)
NOW = datetime.now(timezone.utc).isoformat()
SHA = hashlib.sha256(PDF.read_bytes()).hexdigest()
NPAGES = 15


def T(s: str) -> str:
    return s.replace("__NOW__", NOW).replace("__SHA__", SHA)


EQS = [
    ("Eq.(1)", r"p_l(r_i,z_r)=\sqrt{2\pi}e^{-j(\pi/4)}\sum_{m=1}^M \phi_m(z_s)\phi_m(z_r)\frac{\exp\{-j(k_m-j\alpha_m)[r_i+ld\sin\theta_i]\}}{\sqrt{k_m[r_i+ld\sin\theta_i]}}", "HLA element normal-mode field with path r_i+ld sinθ_i", "P2"),
    ("Eq.(2)", r"p_l(r_i,z_r)=\sum_m A_m\exp\{-j(k_m-j\alpha_m)[r_i+ld\sin\theta_i]\},\quad A_m=\sqrt{2\pi}e^{-j\pi/4}\frac{\phi_m(z_s)\phi_m(z_r)}{\sqrt{k_m r_i}}", "far-field r_i+ld sinθ≈r_i; modal amp A_m", "P3"),
    ("Eq.(3)", r"B(\hat\theta_i)=\frac{1}{2L+1}\sum_{l=-L}^{L} e^{j k l d \sin\hat\theta_i}\, p_l(r_i,z_r)", "conventional HLA beamforming weights", "P3"),
    ("Eq.(4)", r"B(\hat\theta_i)=\frac{1}{2L+1}\sum_m A_m e^{-jk_m r_i-\alpha_m r_i}\frac{\sin[(L+1/2)d X_m]}{\sin((d/2)X_m)},\ X_m=-(k_m-j\alpha_m)\sin\theta_i+k\sin\hat\theta_i", "array factor; BEAM_FACTOR_FROM_EQ4 = sin[(L+1/2)dX]/sin((d/2)X)", "P3"),
    ("Eq.(5)", r"p(r;z_s,z_r)=\int_0^\infty g(k_r)J_0(k_r r)k_r dk_r,\quad g(k_r)=\int_0^\infty p(r)J_0(k_r r) r dr", "Hankel pair", "P3"),
    ("Eq.(6)", r"g(k_r;z_s,z_r)\sim\frac{e^{i\pi/4}}{\sqrt{2\pi k_r}}\int_{-\infty}^{+\infty}p(r;z_s,z_r)e^{i k_r r}\sqrt{r}\,dr,\ k_r r\gg1", "FT approx of Hankel", "P3"),
    ("Eq.(7)", r"g(k_r,z_r)=\frac{e^{i\pi/4}}{\sqrt{2\pi k_r}}\int_{r_0}^{r_0+R} B(r)\,e^{i k_r r}\,S(r)\,dr,\ k_r r_0\gg1", "generalized Hankel on beam output", "P3"),
    ("Eq.(8)", r"S(r)=\langle |B(r)|^2\rangle^{-1/2}", "data spreading compensation", "P3"),
    ("Eq.(9)", r"g\sim\sum_m \frac{\phi_m(z_s)\phi_m(z_r)}{\sqrt{k_r k_m}}\sin b(X_m)\int_{r_0}^{r_0+R}e^{j(k_r-k_m)r-\alpha_m r}dr=\sum_m a_m\frac{\phi_m(z_s)\phi_m(z_r)}{k_r-k_m+j\alpha_m}", "closed form after Eq4+Eq7; integral kernel exp[j(kr-km)r - αm r]", "P3"),
    ("Eq.(10)", r"a_m=\frac{e^{j[(k_r-k_m)-\alpha_m](r_0+R)}-e^{j[(k_r-k_m)-\alpha_m]r_0}}{j\sqrt{k_r k_m}\,\sin b(X_m)}", "PRINTED a_m — algebra audit flagged", "P3"),
    ("Eq.(11)", r"g(k_m,z_r)\sim b_m\phi_m(z_r)", "peak at k_r=k_m", "P3"),
    ("Eq.(12)", r"b_m=\frac{2e^{-\alpha_0 r'}}{\alpha_m k_m}\sinh\left(\frac{\alpha_m R}{2}\right)\phi_m(z_s)\sin b(X_m),\ r'=r_0+R/2", "PRINTED b_m with e^{-α0 r'} — α0 subscript flagged", "P3"),
    ("Eq.(13)", r"\mathbf{g}=\Phi\cdot\mathbf{b}", "matrix form", "P3"),
    ("Eq.(14)", r"\mathbf{g}=[g(k_1,z_r),\ldots,g(k_M,z_r)]^T", "spectrum vector", "P3"),
    ("Eq.(15)", r"\Phi=\mathrm{diag}([\phi_1(z_r),\ldots,\phi_M(z_r)])", "diag mode functions at zr", "P3"),
    ("Eq.(16)", r"\mathbf{b}=[b_1,\ldots,b_M]^T", "modal coefficients", "P3"),
    ("Eq.(17)", r"y[i]=B(r_i)S(r_i),\quad i=1,2,\ldots,2L+1", "AR input = beamformed × S; SAMPLE COUNT CONFLICT", "P4"),
    ("Eq.(18)", r"y[i]=-\sum_{k=1}^{p}a[k]y[i-k]+u[i],\quad p\ \mathrm{often}=(2/3)(2L+1)", "AR model; order notation CONFLICT", "P4"),
    ("Eq.(19)", r"P_{AR}(l)=\frac{\sigma^2}{\left|1+\sum_{k=1}^{p}a[k]\exp[-ilk]\right|}", "AR spectrum; l→k_r mapping NOT_EXPLICIT", "P4"),
    ("Eq.(20)", r"D(z)=\varphi(z)\,\mathbf{b}\mathbf{b}^H\varphi^H(z)", "matched-mode depth function", "P4"),
    ("Eq.(21)", r"\varphi(z)=[\phi_1(z),\phi_2(z),\ldots,\phi_M(z)]", "replica row", "P4"),
    ("Eq.(22)", r"\mathbf{b}=(\Phi+U)^{-1}\mathbf{g}", "regularized solve", "P4"),
    ("Eq.(23)", r"U=\mathrm{diag}\left(\left[\frac{\Delta}{\phi_1(z_r)},\frac{\Delta}{\phi_2(z_r)},\ldots,\frac{\Delta}{\phi_M(z_r)}\right]\right),\ \Delta\sim\tfrac12\max|\phi|", "Liang regularization — NOT Yang 0.1", "P4"),
    ("Eq.(24)", r"\min_{k_0}(\mathbf{k}-\mathbf{k}_0)^H(\mathbf{k}-\mathbf{k}_0)\ \mathrm{s.t.}\ k_0(1)<k_0(2)<\cdots<k_0(M_0)", "ORDERED_SUBSET_GLOBAL_MATCH", "P4"),
    ("Eq.(25)", r"\mathbf{k}=[k_1,\ldots,k_{M_0}]^T", "estimated wavenumbers", "P4"),
    ("Eq.(26)", r"\mathbf{k}'=[k'_1,\ldots,k'_M]^T", "KRAKEN/model wavenumbers", "P4"),
    ("Eq.(27)", r"\mathbf{k}_0\subseteq\mathbf{k}'", "ordered subset of model k'", "P4"),
    ("Eq.(28)", r"\mathrm{SNR}=10\lg\frac{P_s}{P_n}\Big|_{r=r_0}", "SNR at r0", "P5"),
    ("Eq.(29)", r"(\text{correct-depth probability — see paper Sec.4 MC})\ \ | \hat z-z_\mathrm{true}|\le 5\ \mathrm{m}", "MC success criterion |err|<=5m; C0=500", "P5+"),
    ("Eq.(30)", r"(\text{CI for probability — binomial CI as printed in paper})", "confidence interval of correct-depth rate", "P5+"),
]


def main() -> int:
    # primary lock
    (OUT / "MMAR_PRIMARY_SOURCE_LOCK.md").write_text(
        T(
            """# MMAR_PRIMARY_SOURCE_LOCK

UTC: __NOW__

## PRIMARY_SOURCE_AVAILABLE = TRUE

| 项 | 值 |
| --- | --- |
| filename | literature/liang2018_mmar_7824671.pdf |
| user path | C:\\Users\\ccc\\Desktop\\liang.pdf |
| SHA256 | `__SHA__` |
| pages | %d |
| title | Match-Mode Autoregressive Method for Moving Source Depth Estimation in Shallow Water Waveguides |
| authors | Liang Guo-Long; Zhang Yi-Feng; Zou Nan; Wang Jin-Jin |
| journal | Mathematical Problems in Engineering 2018, Article ID 7824671 |
| DOI | 10.1155/2018/7824671 |
| published | 31 December 2018 |
| received / accepted | 23 Nov 2018 / 16 Dec 2018 |

作者以主文为准（勿用二手错误拼写）。

状态：`C2_4A_PRIMARY_PDF_REQUIRED` → **SUPERSEDED_BY_USER_PRIMARY_PDF**
新增：**PRIMARY_SOURCE_AVAILABLE**

保持：`C2_3A_FOURIER_MODE_IDENTITY_LIMITED`、`ORACLE_EQ5_DEPTH_SIGNATURE_SURVIVES`、`YANG_ROUTE_UNDECIDED`
MMAR = **RC3-C2-AR/MMAR INDEPENDENT CANDIDATE**（不是 YANG_AR_VERSION）
"""
            % NPAGES
        ),
        encoding="utf-8",
    )

    # equation lock + tables
    rows = []
    parts = [T("# MMAR_EQUATION_LOCK\n\nUTC: __NOW__\n\n视觉+文本双通道；冲突以 PDF 视觉为准。\n")]
    for eid, tex, role, page in EQS:
        rows.append(
            {
                "equation_id": eid,
                "exact_formula_latex": tex,
                "paper_page": page,
                "symbols": "see SYMBOL_TABLE / lock body",
                "input": "",
                "output": "",
                "purpose": role,
                "implementation_required": "YES" if eid.split(".")[1].strip("()") in "1234789101112131718192021222324" else "DOC",
                "confidence": "VISUAL_TEXT",
            }
        )
        parts.append(f"\n## {eid} ({page})\n\n```latex\n{tex}\n```\n\n{role}\n")
    pd.DataFrame(rows).to_csv(OUT / "MMAR_EQUATION_TABLE.csv", index=False, encoding="utf-8-sig")
    (OUT / "MMAR_EQUATION_LOCK.md").write_text("".join(parts), encoding="utf-8")

    syms = [
        ("p_l(r_i,z_r)", "pressure at HLA element l, snapshot i", "Eq1"),
        ("r_i", "source-to-HLA-center range at i-th moment (r_i=r_0+i V Δt)", "Eq1"),
        ("l", "HLA element index -L..L; 2L+1 elements", "Eq1-4"),
        ("d", "element spacing (paper: λ/2)", "Eq1"),
        ("θ_i, θ̂_i", "true/look angle; paper assumes θ̂=θ=θ", "Eq3"),
        ("k_m, α_m, φ_m", "mode wavenumber, attenuation, depth function", "Eq1"),
        ("X_m", "-(k_m-j α_m) sinθ_i + k sin θ̂_i", "Eq4"),
        ("sin b(X_m)", "BEAM_FACTOR_FROM_EQ4 = sin[(L+1/2)dX]/sin((d/2)X)  [notation in later eqs]", "Eq4/9/12"),
        ("B(r)", "HLA beamformed complex output", "Eq3/7"),
        ("S(r)", "⟨|B|²⟩^{-1/2} spreading compensation", "Eq8"),
        ("R", "range span during observation", "Eq7"),
        ("r0", "unknown initial range (paper sim 5010 m)", "Eq7"),
        ("y[i]", "AR input = B(r_i) S(r_i)", "Eq17"),
        ("p", "AR order, often (2/3)(2L+1) — NOTATION CONFLICT", "Eq18"),
        ("a[k], σ², u[i]", "AR coefficients, noise var, white noise", "Eq18-19"),
        ("P_AR(l)", "AR spectrum; l→k_r NOT_EXPLICIT", "Eq19"),
        ("D(z)", "matched-mode depth function φ b b^H φ^H", "Eq20"),
        ("b", "(Φ+U)^{-1} g", "Eq22"),
        ("U, Δ", "diag(Δ/φ_m(z_r)); Δ ~ (1/2) max|φ|  (NOT Yang 0.1)", "Eq23"),
        ("k, k', k0", "estimated / model / ordered subset wavenumbers", "Eq24-27"),
    ]
    pd.DataFrame(syms, columns=["symbol", "definition", "where"]).to_csv(
        OUT / "MMAR_SYMBOL_TABLE.csv", index=False, encoding="utf-8-sig"
    )

    (OUT / "MMAR_METHOD_IDENTITY.md").write_text(
        T(
            """# MMAR_METHOD_IDENTITY

UTC: __NOW__

方法链（主文确认）：

**HLA beamforming (Eq1–4) → AR 估 k̂ (Eq17–19) → generalized Hankel 在 k̂ 处取 g(k̂) (Eq7) → ordered-subset mode matching (Eq24–27) → matched-mode D(z) (Eq20–23)**

| 模块 | 来源 | 角色 |
| --- | --- | --- |
| HLA beamforming | Liang2018 | 相对 Yang 单听器 SAB **新增** |
| moving-range sequence | 同 Yang SAB 类 | 运动合成距离 |
| AR wavenumber | Liang2018 核心 | **只定 k̂** |
| Hankel amplitudes | Yang/Ref 传统 | 在 k̂ 处取幅度 |
| ordered mode matching | Liang2018 Eq24–27 | **不是 nearest / Hungarian** |
| matched-mode depth | Liang2018 Eq20–23 | 矩阵形式，非 Yang Eq6 原文 |

**不得称为 Yang2015 的 AR 版。**
AR peak height **不得**作 modal amplitude。
"""
        ),
        encoding="utf-8",
    )

    (OUT / "MMAR_OBSERVABLE_CHAIN.md").write_text(
        T(
            """# MMAR_OBSERVABLE_CHAIN

UTC: __NOW__

```
raw / FIELD pressure on HLA
  → Eq.(3) B(θ̂_i)          # HLA beam
  → Eq.(8) S(r)             # data spreading
  → Eq.(17) y[i]=B S        # AR input
  → Eq.(18–19) P_AR(l)      # ONLY k̂_m (peak locations)
  → Eq.(7)  g(k̂_m, z_r)    # complex amplitudes from Hankel
  → Eq.(22–23) b = (Φ+U)^{-1} g
  → Eq.(20) D(z)=φ b b^H φ^H
  → peak of D(z)
```

禁止：AR peak height → amplitude。
禁止：把 y[i] 当成 raw single hydrophone。
"""
        ),
        encoding="utf-8",
    )

    (OUT / "MMAR_AR_ORDER_AUDIT.md").write_text(
        T(
            """# MMAR_AR_ORDER_AUDIT

UTC: __NOW__

## AR_SAMPLE_COUNT_NOTATION_CONFLICT

| 证据 | 内容 |
| --- | --- |
| 正文 | HLA 有 **2L+1** receivers |
| Eq.(17) | y[i], i=1…**2L+1** |
| Eq.(18) | p often **(2/3)(2L+1)** |
| 仿真 | **L=5 → 11 elements** |

字面解释：N_data=11, p≈7。与 4990/1990 m 运动孔径谱估计物理不一致。

| 解释 | 状态 |
| --- | --- |
| PRINTED_LITERAL: N_data=11, p≈7 | 记录 |
| PHYSICALLY_REQUIRED_MOVING_SEQUENCE: 由 range samples 组成，**论文无独立样本数符号** | 记录 |
| 本轮选定 | **禁止决定** |

## AR estimator

- **modified covariance**（避免 spectral line splitting）
- 峰位置 → k̂_m；**峰高不得作幅度**（large variance）

## AR_NORMALIZED_SPECTRAL_AXIS_MAPPING_NOT_EXPLICIT

P_AR(l) 的 l→物理 k_r 映射**原文未明确**。不得假定 k=ω/c 或 k=l/Δr。
"""
        ),
        encoding="utf-8",
    )

    (OUT / "EQ7_TO_EQ12_ALGEBRA_AUDIT.md").write_text(
        T(
            """# EQ7_TO_EQ12_ALGEBRA_AUDIT

UTC: __NOW__

## Eq.(9) 积分核（印刷）

exp[ j(k_r−k_m) r − α_m r ]

解析积分：

∫ exp[j(k_r−k_m)r − α_m r] dr = Δexp / [ j(k_r−k_m) − α_m ]
= Δexp / [ j( k_r−k_m + j α_m ) ]

故闭式分母应与 **(k_r − k_m + j α_m)** 一致 — Eq.(9) 第二形式分母印刷为 **(k_r−k_m+jα_m)** ✓

## Eq.(10) 印刷形式

a_m = [ e^{ j[(k_r−k_m)−α_m](r0+R) } − e^{ j[(k_r−k_m)−α_m] r0 } ] / ( j√(k_r k_m) sin b(X_m) )

印刷把 α_m 放进 **j[(k_r−k_m)−α_m] = j(k_r−k_m) − j α_m**，
与积分核指数 **j(k_r−k_m)r − α_m r** 的系数 **j(k_r−k_m) − α_m** 不一致。

### 状态：`EQ10_PRINTED_FORM_INTERNAL_INCONSISTENCY`

实现只允许 **Eq.(9) 积分形式** 或由 Eq.(9) 重推的闭式。

## Eq.(12) α 下标

印刷：**e^{−α_0 r'}**（α_0 而非 α_m）。
正文/Eq9/Eq11 要求 mode-specific α_m。

### 状态：`EQ12_ALPHA_SUBSCRIPT_CONFLICT`

实现：**FROM_EQ9_DERIVED_BM**，禁止 PRINTED_EQ12_BLIND_COPY。
（对比 Yang2015：b_m 用 e^{−α_m r0}，分母 α_m k_m，sinh(α_m R/2)，另乘 sin b(X_m)。）
"""
        ),
        encoding="utf-8",
    )

    (OUT / "MMAR_MODE_ORDER_MATCHING.md").write_text(
        T(
            """# MMAR_MODE_ORDER_MATCHING

UTC: __NOW__

Eq.(24)–(27)：

- 估计向量 k=[k_1…k_{M0}]^T
- 模型 KRAKEN k'=[k'_1…k'_M]^T
- 从 k' 中选 **有序子集** k0，满足 k0(1)<…<k0(M0)
- 最小化 (k−k0)^H (k−k0)

类型：**ORDERED_SUBSET_GLOBAL_MATCH**

不是：independent nearest neighbor / Hungarian / C2.3A Rayleigh unique-only。

## 对 C2.3A 的含义

MMAR **不要求** 每个 AR 峰在 Rayleigh cell 内唯一 mode。
它用高分辨 k̂ + 全局有序子集匹配定 mode order。

后续指标应是 **MODE_ORDER_RECOVERY_RATE**，不是 N_UNIQUE_RAYLEIGH_MODES。
"""
        ),
        encoding="utf-8",
    )

    (OUT / "MMAR_OFFSET_RANGE_AUDIT.md").write_text(
        T(
            """# MMAR_OFFSET_RANGE_AUDIT

UTC: __NOW__

论文主张 without knowing absolute source range；仿真 r0=5010 m 当作未知。

## MMAR_INHERITS_OFFSET_RANGE_CAVEAT

Yang 2018 Erratum：合成距离积分与数据距离存在 offset-range δ，
g = b_m φ_m(z_r) exp(i k_m δ)；真实数据需搜索 δ。

因此：

- 不得写 “MMAR 不需要距离”
- paper reproduction：**δ=0 = ORACLE_OFFSET_ALIGNMENT**
- 以后：RC2-conditioned δ search（本轮不做）
"""
        ),
        encoding="utf-8",
    )

    (OUT / "MMAR_PAPER_CONFIG.md").write_text(
        T(
            """# MMAR_PAPER_CONFIG

UTC: __NOW__

| 项 | 值 | 来源 |
| --- | --- | --- |
| frequency | 350 Hz | §4 |
| HLA | 2L+1, L=5 → **11** elements | §2/§4 |
| d | **λ/2** | §4 |
| HLA depth | **70 m** | §4 |
| z_s | 4 m / 50 m | §4 |
| r0 | 5010 m（算法当未知） | §4 |
| v | 2.5 m/s | §4 |
| motion | away, 沿 beam 方向（对比用） | §4 |
| SNR | 20 / 5 / −5 dB | Eq.28, §4 |
| sufficient span | **4990 m** | §4 |
| insufficient span | **1990 m** | §4 |
| field | KRAKEN | §4 |
| bottom | referred to [7]=Yang2015 Ref.7 | §4 |
| MC | C0=**500**；correct: \|ẑ−z_true\|≤**5 m** | §4 |
| Fig7 span scan | 1990…4990 m | §4 |
| look angle | θ̂=θ=θ, source along beam | §2 |
| depth during aperture | **fixed**（快变则失效） | §1/§5 |

### ORACLE_LOOK_DIRECTION + PAPER_RADIAL_MOTION

E-STD 不得免费继承。
"""
        ),
        encoding="utf-8",
    )

    (OUT / "MMAR_PAPER_FIDELITY_TARGETS.md").write_text(
        T(
            """# MMAR_PAPER_FIDELITY_TARGETS

UTC: __NOW__

| Target | 内容 | 复现阶段 |
| --- | --- | --- |
| Fig.3 | sufficient 4990 m, SAB vs MMAR spectra | PAPER-R1 |
| Fig.4 | insufficient 1990 m, SAB vs MMAR | PAPER-R1 |
| Fig.5 | sufficient span depth | PAPER-R1 |
| Fig.6 | insufficient span depth | PAPER-R1 |
| Fig.7 | P(correct) vs span, SNR 20/5/−5 | PAPER-R1–R3 + MC |

分层：PAPER-R0 无噪/高SNR机制 → R1 20dB → R2 5dB → R3 −5dB；MC 暂缓。

Fig.7 描述（target only）：
- 20 dB：4490–4990 m 两法都高；1990–3990 m MMAR>SAB
- 5 dB：MMAR ≳2490 m 较好
- −5 dB：MMAR ≳2990 m 较好

本轮不运行数值复现。
"""
        ),
        encoding="utf-8",
    )

    (OUT / "MMAR_ESTD_TRANSFER_AUDIT.md").write_text(
        T(
            """# MMAR_ESTD_TRANSFER_AUDIT

UTC: __NOW__

| paper assumption | E-STD | 判定 |
| --- | --- | --- |
| range-independent | 近 CZ 近似 | APPROX |
| source depth fixed in aperture | UUV 短窗理论可 | YES_THEORETICAL；**FIXED_DEPTH_DURING_APERTURE_REQUIRED** |
| known radial increments | 理论上限 | YES_UPPER_BOUND |
| known f0 | S2 上限 | YES_S2 |
| HLA | 项目有拖曳 HLA | YES_HARDWARE_CLASS |
| θ̂=θ | 需 RC2 bearing | **NOT_FREE**；本轮 ORACLE_LOOK_DIRECTION |
| source along beam | 一般不成立 | NO_GENERAL |
| far-field r_i+ld sinθ≈r_i | 14 m 阵 @50 km | YES_NEAR |
| shallow modal density | CZ 更密 | **DIFFERENT** |
| Δ≈½ max\|φ\| | Liang 正则 | LIANG_PAPER_REGULARIZER |
| δ | Erratum | ORACLE_OFFSET_ALIGNMENT |
"""
        ),
        encoding="utf-8",
    )

    (OUT / "FOURIER_MODE_SET_SCOPE_CORRECTION.md").write_text(
        T(
            """# FOURIER_MODE_SET_SCOPE_CORRECTION

UTC: __NOW__

SET-D_REFRACTED_TRAPPED 由 zgrid 仅 150–250 m 判 turning：

→ **INVALID_WITH_CURRENT_ZGRID**（无法判 surface / deep turning / bottom）

Fourier 稳健冻结仅基于：

- SET-A ALL
- SET-B receiver-observable
- SET-C search-band supported

三者仍 n_rayleigh_unique=0。

标签作用域：

**FOURIER_MODE_IDENTITY_LIMIT_ROBUST_TO_MODE_SET**
→ **ROBUST_ACROSS_NONORACLE_OBSERVABILITY_SETS_A_B_C**
"""
        ),
        encoding="utf-8",
    )

    pd.DataFrame(
        [
            ("A_PAPER_EQUIV", 11, 0.5, "d=lambda/2 @350Hz≈0.5*1500/350≈2.14m — use actual lambda/2 in implementation; paper L=5", "WAIT exact lambda"),
            ("B_PROJECT_14M", 8, 2.0, "8 elements × d=2 m → physical aperture 14 m (corrected)", "PROJECT_REPRESENTATIVE"),
            ("C_4EL_SCAN_L5", 4, 5.0 / 3, "4-element uniform, d=L/3 for historical L=5m", "CONDITION_ENHANCEMENT_SCAN"),
            ("C_4EL_SCAN_L10", 4, 10.0 / 3, "d=L/3, L=10 m", "CONDITION_ENHANCEMENT_SCAN"),
            ("C_4EL_SCAN_L20", 4, 20.0 / 3, "d=L/3, L=20 m", "CONDITION_ENHANCEMENT_SCAN"),
            ("C_4EL_SCAN_L40", 4, 40.0 / 3, "d=L/3, L=40 m", "CONDITION_ENHANCEMENT_SCAN"),
            ("C_4EL_SCAN_L80", 4, 80.0 / 3, "d=L/3, L=80 m", "CONDITION_ENHANCEMENT_SCAN"),
        ],
        columns=["config", "N_elem", "d_m", "note", "status"],
    ).to_csv(OUT / "MMAR_HLA_MIGRATION_TABLE_CORRECTED.csv", index=False, encoding="utf-8-sig")

    # AR resolution wording correction (recompute from previous csv if present)
    res_note = T(
        """# AR resolution wording

E-STD 相邻 mode Δk ~ 1e-4–3e-4 rad/m

→ 冻结为 **ESTD_ADJACENT_MODE_SPACING_SCALE**

Fourier 等效孔径 ~20–65 km（2π/Δk）。

**不得**称 AR_HARD_REQUIRED_RESOLUTION：Eq.(24) 全局 ordered matching 不要求所有相邻 mode 逐条分开。

真指标：k̂ 误差 vs ordered-subset **MODE_ORDER_RECOVERY_RATE**。
"""
    )
    (OUT / "AR_RESOLUTION_WORDING.md").write_text(res_note, encoding="utf-8")

    decision = "C2_4A_PRIMARY_METHOD_LOCKED_WITH_NOTATION_CONFLICTS"
    conflicts = [
        "AR_SAMPLE_COUNT_NOTATION_CONFLICT",
        "EQ10_PRINTED_FORM_INTERNAL_INCONSISTENCY",
        "EQ12_ALPHA_SUBSCRIPT_CONFLICT",
        "AR_NORMALIZED_SPECTRAL_AXIS_MAPPING_NOT_EXPLICIT",
        "MMAR_INHERITS_OFFSET_RANGE_CAVEAT",
    ]
    why = (
        "主文 PDF 已锁定；Eq.(1)–(28)+方法链/mode-order/正则化/仿真配置已恢复。"
        "存在上述原文内部记号/印刷冲突，故 WITH_NOTATION_CONFLICTS；实现须用 Eq9 积分重推闭式、"
        "AR 只定 k̂、ordered-subset matching、Liang Δ~½max|φ|、δ=0。冲突≠方法失败。"
    )
    dec = {
        "stage": "R3-C2.4A-PRIMARY",
        "rc3_c2_4a_primary_decision": decision,
        "why": why,
        "conflicts": conflicts,
        "primary_source_available": True,
        "sha256": SHA,
        "method_chain": "HLA beam → AR k̂ → Hankel g(k̂) → ordered-subset match → D(z)",
        "regularizer": "LIANG_PAPER_REGULARIZER Delta~0.5 max|phi|",
        "offset_range": "delta=0 ORACLE_OFFSET_ALIGNMENT",
        "mode_id_rule": "ORDERED_SUBSET_GLOBAL_MATCH Eq24-27",
        "yang_route": "YANG_ROUTE_UNDECIDED",
        "mmar_route": "RC3-C2-AR/MMAR INDEPENDENT CANDIDATE",
        "not_done": ["AR code", "MMAR run", "FIELD", "numerical repro", "E-STD AR", "MC", "P5"],
        "created_utc": NOW,
    }
    (OUT / "R3_C2_4A_PRIMARY_DECISION.json").write_text(json.dumps(dec, ensure_ascii=False, indent=2), encoding="utf-8")

    (OUT / "R3_C2_4A_PRIMARY_REPORT.md").write_text(
        T(
            """# R3-C2.4A-PRIMARY 报告

UTC: __NOW__

## 状态

`PRIMARY_SOURCE_AVAILABLE`；判定：

### C2_4A_PRIMARY_METHOD_LOCKED_WITH_NOTATION_CONFLICTS

## 方法链（主文）

HLA Eq1–4 → AR k̂ Eq17–19 → Hankel g(k̂) Eq7 → **ordered-subset** Eq24–27 → D(z) Eq20–23

AR **只定波数**；幅度来自 Hankel。正则化 Δ~½ max|φ|（非 Yang 0.1）。

## 冲突清单（≠方法失败）

1. AR_SAMPLE_COUNT_NOTATION_CONFLICT（2L+1 阵元 vs 运动样本）
2. EQ10_PRINTED_FORM_INTERNAL_INCONSISTENCY（α 进 j[·]）
3. EQ12_ALPHA_SUBSCRIPT_CONFLICT（e^{-α0 r'}）
4. AR_NORMALIZED_SPECTRAL_AXIS_MAPPING_NOT_EXPLICIT
5. MMAR_INHERITS_OFFSET_RANGE_CAVEAT（δ=0=ORACLE_OFFSET_ALIGNMENT）

实现约束：Eq9 积分重推；FROM_EQ9_DERIVED_BM。

## 其它修正

- SET-D → INVALID_WITH_CURRENT_ZGRID；Fourier 稳健范围改为 SET-A/B/C
- HLA：8×2 m=14 m（非 14×1 m）；4 元条件增强为 scan
- AR 分辨率措辞：ESTD_ADJACENT_MODE_SPACING_SCALE（非 HARD requirement）

## 停止

不编码 AR、不跑 MMAR/FIELD、不进 E-STD AR、不 MC、不进 P5。
"""
        ),
        encoding="utf-8",
    )
    (OUT / "R3_C2_4A_PRIMARY_GPT_SYNC.md").write_text(
        T(f"# R3-C2.4A-PRIMARY\\n\\n**{decision}**\\n\\nconflicts={conflicts}\\nδ=0 ORACLE；Δ~½max|φ|；ORDERED_SUBSET_MATCH\\n"),
        encoding="utf-8",
    )
    print("DECISION", decision)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
