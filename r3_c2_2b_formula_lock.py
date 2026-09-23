#!/usr/bin/env python3
"""R3-C2.2B: lock Yang 2015 formulas + paper reproduction protocol. No algorithm run."""
from __future__ import annotations

import hashlib
import json
from datetime import datetime, timezone
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parent
OUT = ROOT / "results" / "R3_C2_Yang_SA_depth" / "R3_C2_2B"
PDF = ROOT / "literature" / "yang2015_Yang2015_JASA138_1678.pdf"
OUT.mkdir(parents=True, exist_ok=True)
NOW = datetime.now(timezone.utc).isoformat()
SHA = hashlib.sha256(PDF.read_bytes()).hexdigest()
NPAGES = 10


def T(s: str) -> str:
    return s.replace("__NOW__", NOW).replace("__SHA__", SHA)


EQ = [
    {
        "equation_id": "Eq.(1)",
        "exact_formula_latex": r"g(k_r,z_r)=\frac{e^{i\pi/4}}{\sqrt{2\pi k_r}}\int_{r_1}^{r_2} p(r,z_r)\,e^{i k_r r}\,S(r)\,dr,\quad k_r r_1\gg 1",
        "paper_page": "1679 (PDF p.3)",
        "all_symbols": r"g(k_r,z_r): SA modal beam output / wavenumber spectrum; k_r: steering wavenumber; z_r: receiver depth; p(r,z_r): block pressure at range r; S(r): spreading-loss shading; r1,r2: aperture ends; e^{i\pi/4}: cylindrical-wave prefactor",
        "input": "p(r_n,z_r) on radial samples, S(r), k_r grid",
        "output": "complex spectrum g(k_r,z_r)",
        "physical_role": "synthetic-aperture modal beamforming (generalized Hankel if S=sqrt(r))",
        "used_in_single_receiver": "YES (core)",
        "used_in_VLA": "YES (per hydrophone, then summed)",
        "implementation_required": "YES",
    },
    {
        "equation_id": "Eq.(2)",
        "exact_formula_latex": r"S(r)=\langle |p(r,z_r)|^2\rangle^{-1/2}",
        "paper_page": "1679 (PDF p.3)",
        "all_symbols": r"\langle\cdot\rangle: range averaging/smoothing; |p|^2: intensity",
        "input": "pressure intensity vs range",
        "output": "shading S(r) compensating spreading loss (~sqrt(r))",
        "physical_role": "weighting in Eq.(1) integrand (NUMERATOR of integrand, multiplies p*exp(i k_r r))",
        "used_in_single_receiver": "YES",
        "used_in_VLA": "replaced by Eq.(8)",
        "implementation_required": "YES",
    },
    {
        "equation_id": "Eq.(3)",
        "exact_formula_latex": r"p(r,z_r)=\sum_{m=1}^{M}\sqrt{\frac{2\pi}{k_m r}}\,\phi_m(z_s)\phi_m(z_r)\,e^{-i k_m r-\alpha_m r-i\pi/4}",
        "paper_page": "1679 (PDF p.3)",
        "all_symbols": r"phi_m(z): mode depth function; k_m: mode wavenumber; alpha_m: mode attenuation (>0); z_s: source depth; M: number of modes; spatial phase e^{-i k_m r} (as printed)",
        "input": "modal parameters, r, z_s, z_r",
        "output": "normal-mode pressure (range-independent)",
        "physical_role": "propagator used to derive Eq.(4)/(5)",
        "used_in_single_receiver": "derivation / oracle peaks",
        "used_in_VLA": "same",
        "implementation_required": "for algebra check and ORACLE mode peaks",
    },
    {
        "equation_id": "Eq.(4)",
        "exact_formula_latex": (
            r"g(k_r,z_r)\approx\sum_{m=1}^{M}\frac{\phi_m(z_s)\phi_m(z_r)}{\sqrt{k_r k_m}}"
            r"\int_{r_1}^{r_2} e^{i(k_r-k_m)r-\alpha_m r}\,dr"
            r"=\sum_{m=1}^{M} a_m\frac{\phi_m(z_s)\phi_m(z_r)}{k_r-k_m-i\alpha_m}"
        ),
        "paper_page": "1679 (PDF p.3)",
        "all_symbols": r"a_m=[e^{i(k_r-k_m)r_2-\alpha_m r_2}-e^{i(k_r-k_m)r_1-\alpha_m r_1}]/(i k_m) (AS PRINTED)",
        "input": "Eq.(3) into Eq.(1) with S~sqrt(r)",
        "output": "spectrum as sum of modal Lorentzians",
        "physical_role": "closed form after substitution; peaks at k_r=k_m",
        "used_in_single_receiver": "YES (equivalent form)",
        "used_in_VLA": "YES",
        "implementation_required": "YES - implement FIRST FORM (integral) or Eq.(5) peak form; see algebra check",
    },
    {
        "equation_id": "Eq.(5)",
        "exact_formula_latex": r"g(k_m,z_r)=b_m\,\phi_m(z_r)",
        "paper_page": "1679 (PDF p.3)",
        "all_symbols": r"b_m=(2 e^{-\alpha_m r_0})/(\alpha_m k_m)\sinh[(\alpha_m \Delta R)/2]\,\phi_m(z_s); r_0=(r_2+r_1)/2; \Delta R=r_2-r_1=vT",
        "input": "spectrum value at modal peak",
        "output": "modal peak amplitude b_m phi_m(z_r)",
        "physical_role": "peak relation; b_m carries z_s through phi_m(z_s) and range factors",
        "used_in_single_receiver": "YES (feeds Eq.6)",
        "used_in_VLA": "YES (Eq.9 with delta_m)",
        "implementation_required": "YES",
    },
    {
        "equation_id": "Eq.(6)",
        "exact_formula_latex": (
            r"D(z)=\left|\sum_{m=1}^{M}\phi_m(z)\,\frac{g(k_m,z_r)}{\bar{\phi}_m(z_r)}\right|^2,"
            r"\quad \bar{\phi}_m^{-1}(z_r)=\frac{\phi_m(z_r)}{\phi_m^2(z_r)+\Delta^2}"
        ),
        "paper_page": "1679 (PDF p.3)",
        "all_symbols": r"D(z): source-depth ambiguity; \bar{\phi}_m^{-1}: regularized inverse mode depth function at z_r; \Delta: empirical regularizer (order 0.1*max|phi_m|)",
        "input": "g(k_m,z_r) at identified modal peaks; nominal phi_m(.)",
        "output": "depth ambiguity D(z) on search grid",
        "physical_role": "single-hydrophone source-depth estimator (Bartlett class)",
        "used_in_single_receiver": "YES (core)",
        "used_in_VLA": "NO (uses Eq.9 instead)",
        "implementation_required": "YES",
    },
    {
        "equation_id": "Eq.(7)",
        "exact_formula_latex": r"\tilde{\theta}(n+1)=\tilde{\theta}(n)+K_1 u(n)+K_2\sum_{m=0}^{n} u(m)",
        "paper_page": "1681 (PDF p.5)",
        "all_symbols": r"theta~(n): PLL output phase; u(m): phase error data-oscillator; K1=0.001; K2=0.0001 (2nd-order PLL)",
        "input": "narrowband complex time series",
        "output": "tracked instantaneous phase theta~(t)",
        "physical_role": "remove fast random phase (jitter + medium)",
        "used_in_single_receiver": "YES for real data",
        "used_in_VLA": "YES for real data",
        "implementation_required": "YES for real-data path; optional for pure simulation",
    },
    {
        "equation_id": "Eq.(8)",
        "exact_formula_latex": r"S(r)=\left\langle \sum_{j=1}^{N}|p(r,z_j)|^2\,\Delta z\right\rangle^{-1/2}",
        "paper_page": "1683 (PDF p.7)",
        "all_symbols": r"z_j: VLA hydrophone depths; \Delta z: phone spacing; N: #phones",
        "input": "VLA pressure vs range",
        "output": "depth-integrated intensity shading",
        "physical_role": "VLA replacement of Eq.(2); S~sqrt(r) ignoring alpha_m differences",
        "used_in_single_receiver": "NO",
        "used_in_VLA": "YES",
        "implementation_required": "for VLA path only",
    },
    {
        "equation_id": "Eq.(9)",
        "exact_formula_latex": (
            r"D(z_j)=\left|\sum_{m=1}^{M} g(k_m,z_j)\right|^2"
            r"=\left|\sum_{m=1}^{M}\phi_m(z_j)\phi_m(z_s)\delta_m\right|^2,\quad j=1,\ldots,N"
        ),
        "paper_page": "1683 (PDF p.7)",
        "all_symbols": r"delta_m=(2 e^{-\alpha_m r_0})/(\alpha_m k_m)\sinh[(\alpha_m \Delta R)/2] (mode shading coefficient); note delta_m = b_m/phi_m(z_s)",
        "input": "per-phone SA spectra at modal peaks g(k_m,z_j)",
        "output": "beam power vs hydrophone depth = depth ambiguity",
        "physical_role": "VLA model-independent depth estimation (no phi_m needed experimentally)",
        "used_in_single_receiver": "NO",
        "used_in_VLA": "YES",
        "implementation_required": "for VLA path only",
    },
    {
        "equation_id": "Eq.(10)",
        "exact_formula_latex": r"D_{\mathrm{MMP}}(r=r_s,z)=\left|\sum_{m=1}^{M}\phi_m(z)\phi_m(z_s)\delta_m'\right|^2,\quad \delta_m'\equiv e^{-2\alpha_m r_s}/k_m",
        "paper_page": "1684 (PDF p.8)",
        "all_symbols": r"r_s: true source range; delta_m': MMP shading; related D_MMP(r,z)=(r/2pi)|sum[A_m^{rplc}]^* A_m^{data}|^2",
        "input": "matched-mode amplitudes",
        "output": "MMP range-depth ambiguity at r=r_s",
        "physical_role": "relation to matched-mode processing (comparison only)",
        "used_in_single_receiver": "NO",
        "used_in_VLA": "comparison",
        "implementation_required": "NO for Yang SA depth (documentation only)",
    },
]

APP = [
    {
        "equation_id": "Eq.(A1)",
        "exact_formula_latex": (
            r"p(r,z_0,t)=b(t)\sum_{m=1}^{M}\sqrt{\frac{2\pi}{k_m r}}\phi_m(z_s)\phi_m(z_0)"
            r"\times e^{-i(2\pi f_c t - k_m r)-\alpha_m r-i\pi/4}+n(t)"
        ),
        "physical_role": "moving-source narrowband time series (as printed; spatial phase +i k_m r after distributing -i)",
        "note": "SPATIAL PHASE SIGN OPPOSITE to Eq.(3) e^{-i k_m r}",
    },
    {
        "equation_id": "Eq.(A2)",
        "exact_formula_latex": (
            r"p(r,z_0,t_{n,l})\simeq e^{i\theta(t_{n,l})}\sum_{m=1}^{M}\sqrt{\frac{2\pi}{k_m r_n}}\phi_m(z_s)\phi_m(z_0)"
            r"\times e^{-i(2\pi f_c T_n - k_m r_n)-\alpha_m r_n-i\pi/4}\,e^{-i 2\pi f_D (l\Delta t)}+n(t_{n,l})"
        ),
        "physical_role": "Doppler-shifted block model; f_D=f_c(1+v/c); r_n=v T_n; dr=v dt",
        "note": "assumes constant Doppler for dominant modes; |b|~1; e^{-alpha_m (l dr)}~1",
    },
    {
        "equation_id": "Eq.(A3)",
        "exact_formula_latex": (
            r"\tilde{p}(r_n,z_0,f)=\sum_{l=1}^{L} e^{i 2\pi f t_{n,l}}\left(e^{-i\hat{\theta}(t_{n,l})}p(r_{n,l},t_{n,l})\right)"
            r"= e^{i 2\pi (f-f_c)T_n}\sum_{l=1}^{L} e^{i 2\pi f (l\Delta t)} e^{-i 2\pi f_D (l\Delta t)} e^{i\Delta\theta(t_{n,l})}"
            r"\times \sum_{m=1}^{M}\sqrt{\frac{2\pi}{k_m r}}\phi_m(z_s)\phi_m(z_0)e^{i k_m r_n-\alpha_m r_n-i\pi/4}"
        ),
        "physical_role": "PLL/FFT block pressure estimator (printed modal factor has e^{+i k_m r_n})",
        "note": "SPATIAL PHASE +i k_m r_n as printed; conflicts with Eq.(3)",
    },
    {
        "equation_id": "Eq.(A4)",
        "exact_formula_latex": (
            r"\tilde{p}(r_n,z_0,f_D)= L\, e^{i 2\pi (\Delta f) T_n}\, e^{i\overline{\Delta\theta}(T_n)}\, p(r_n,z_0)"
        ),
        "physical_role": "block-by-block pressure recovery; df=f_D-f_c; mean dtheta; p(r_n,z_0) stated to be Eq.(3)",
        "note": "if dtheta=0, correct Doppler phase by e^{-i 2 pi (df) T_n} then use p in Eq.(1)",
    },
]


def write_pdf_note() -> None:
    body = T(
        """# YANG2015_PRIMARY_PDF_NOTE

UTC: __NOW__

## PRIMARY_SOURCE_AVAILABLE = TRUE

| 项 | 值 |
| --- | --- |
| filename | `literature/yang2015_Yang2015_JASA138_1678.pdf` |
| user path | `C:\\Users\\ccc\\Desktop\\yang2015.pdf` |
| SHA256 | `__SHA__` |
| pages | %d |
| DOI | 10.1121/1.4929748 |
| journal pages | 1678-1686 |
| title | Source depth estimation based on synthetic aperture beamfoming for a moving source |
| author | T. C. Yang (Zhejiang University) |

## 状态更新

- `C2_2A_PRIMARY_FORMULA_NOT_RECOVERED` -> **SUPERSEDED_BY_USER_PRIMARY_PDF**
- `YANG_ROUTE_UNDECIDED` **保持不变**

## 审计原则执行记录

- 未依赖旧 Yang-like estimator
- 未从旧代码反推公式
- 公式以 **PDF 页面视觉** 为准，正文文本仅作辅助
- 关键符号（e^{i pi/4}、steering/propagation 相位、alpha_m、S(r) 位置、正则化逆）均对照原页

Erratum PDF：**本轮未提供** -> `CARRY_FORWARD_FROM_PRIOR_ERRATUM_RECOVERY`；论文复现基线 **delta=0**。
"""
        % NPAGES
    )
    (OUT / "YANG2015_PRIMARY_PDF_NOTE.md").write_text(body, encoding="utf-8")


def write_equation_lock() -> None:
    rows = []
    for e in EQ:
        rows.append(
            {
                "equation_id": e["equation_id"],
                "exact_formula_latex": e["exact_formula_latex"],
                "paper_page": e["paper_page"],
                "all_symbols": e["all_symbols"],
                "input": e["input"],
                "output": e["output"],
                "physical_role": e["physical_role"],
                "used_in_single_receiver": e["used_in_single_receiver"],
                "used_in_VLA": e["used_in_VLA"],
                "implementation_required": e["implementation_required"],
            }
        )
    pd.DataFrame(rows).to_csv(
        OUT / "YANG2015_EQUATION_TABLE_LOCKED.csv", index=False, encoding="utf-8-sig"
    )

    parts = [
        T(
            """# YANG2015_EQUATION_LOCK

UTC: __NOW__
PRIMARY SHA256: `__SHA__`

视觉核对页：PDF p.3 (JASA 1679) Eq.(1)-(6)；p.5 Eq.(7)；p.7 Eq.(8)-(9)；p.8 Eq.(10)。
"""
        )
    ]
    for e in EQ:
        parts.append(
            "\n## "
            + e["equation_id"]
            + "  ·  "
            + e["paper_page"]
            + "\n\n```latex\n"
            + e["exact_formula_latex"]
            + "\n```\n\n"
            + "- 符号："
            + e["all_symbols"]
            + "\n- 输入："
            + e["input"]
            + "\n- 输出："
            + e["output"]
            + "\n- 角色："
            + e["physical_role"]
            + "\n- 单水听器："
            + e["used_in_single_receiver"]
            + " · VLA："
            + e["used_in_VLA"]
            + "\n- 实现："
            + e["implementation_required"]
            + "\n"
        )
    (OUT / "YANG2015_EQUATION_LOCK.md").write_text("".join(parts), encoding="utf-8")


def write_appendix() -> None:
    parts = [
        T(
            """# YANG2015_APPENDIX_EQUATIONS

UTC: __NOW__

视觉核对：PDF p.9 / JASA 1685。

## A1-A4
"""
        )
    ]
    for a in APP:
        parts.append(
            "\n### "
            + a["equation_id"]
            + "\n\n```latex\n"
            + a["exact_formula_latex"]
            + "\n```\n\n"
            + "- 角色："
            + a["physical_role"]
            + "\n- 注记："
            + a["note"]
            + "\n"
        )
    parts.append(
        T(
            """
## 完整观测链

```
原始时域 p(r,z0,t) [A1]
  -> Doppler f_D / dr 分块 [A2]
  -> PLL 相位 theta_hat + 低通 [Sec.III / Eq.7]
  -> 块 FFT 压力 p~(r_n,z0,f) [A3]
  -> Doppler 相位校正 -> p(r_n,z0) [A4]
  -> Eq.(1) 波数谱 g(k_r,z_r)
  -> 模态峰识别 (ORACLE / PRACTICAL)
  -> Eq.(6) D(z) -> 归一化 -> argmax
```

## 相位约定冲突（必读）

- Eq.(3) 印刷：空间相位 **e^{-i k_m r}**
- A1/A2/A3 印刷：空间相位 **e^{+i k_m r}**（由 e^{-i(2 pi f_c t - k_m r)} 与 e^{+i k_m r_n} 可见）
- A4 声称 p(r_n,z0) 即 Eq.(3)

**主方法实现必须采用 Eq.(3) 约定**（与 Eq.(1)(4)(5) 代数闭环）。
若从 A3 直接得到 +i k_m r 的压力，须与 Eq.(1) steering 约定一致（或取共轭）后再进谱估计。
此为 **PAPER_APPENDIX_VS_EQ3_PHASE_CONVENTION**，不阻塞主链锁定。
"""
        )
    )
    (OUT / "YANG2015_APPENDIX_EQUATIONS.md").write_text("".join(parts), encoding="utf-8")


def write_algebra_check() -> None:
    body = T(
        r"""# YANG2015_EQ1_EQ5_ALGEBRA_CHECK

UTC: __NOW__

## 目标

验证 **Eq.(3) -> Eq.(1) -> Eq.(4) ->（k_r=k_m）-> Eq.(5)** 代数闭环，
捕捉 steering 相位、传播相位、alpha_m、sqrt(k_r k_m) 因子错误。

## 代入（S(r)~sqrt(r)）

Eq.(1)xEq.(3)：

p*S*e^{i k_r r} ~ sum_m sqrt(2 pi / k_m) phi_m(z_s) phi_m(z_r) e^{i(k_r-k_m)r - alpha_m r - i pi/4}

再乘 e^{i pi/4}/sqrt(2 pi k_r) 并积分：

g ~ sum_m [phi_s phi_r / sqrt(k_r k_m)] int_{r1}^{r2} e^{i(k_r-k_m)r - alpha_m r} dr

**= Eq.(4) 第一形式。** 闭环 PASS

相位核对：
- steering **e^{+i k_r r}**（Eq.1）
- 传播 **e^{-i k_m r}**（Eq.3）
- 合成 **e^{i(k_r-k_m)r}** -> 峰在 k_r=k_m
- 前因子 e^{i pi/4}*e^{-i pi/4}=1
- sqrt(2 pi / k_m)/sqrt(2 pi k_r)=1/sqrt(k_r k_m)
- S~sqrt(r) 消掉 1/sqrt(r)

## k_r=k_m -> Eq.(5)

int_{r1}^{r2} e^{-alpha_m r} dr = (e^{-alpha_m r1}-e^{-alpha_m r2})/alpha_m

= e^{-alpha_m r0} * 2 sinh(alpha_m dR/2)/alpha_m ,  r0=(r1+r2)/2, dR=r2-r1

=> g(k_m,z_r)= [2 e^{-alpha_m r0} sinh(alpha_m dR/2)/(alpha_m k_m)] phi_s phi_r

= b_m phi_m(z_r),

**b_m = (2 e^{-alpha_m r0})/(alpha_m k_m) * sinh[(alpha_m dR)/2] * phi_m(z_s)**

**= Eq.(5) 印刷形式。** 闭环 PASS

## Eq.(4) 第二形式（印刷）内部问题

印刷：

a_m = [e^{i(k_r-k_m)r2 - alpha_m r2} - e^{i(k_r-k_m)r1 - alpha_m r1}]/(i k_m)

g = sum a_m phi_s phi_r / (k_r - k_m - i alpha_m)

与第一形式逐点比较：

1. 积分闭式分母应为 i(k_r-k_m)-alpha_m = i(k_r-k_m+i alpha_m)，对应 **(k_r-k_m+i alpha_m)**，
   印刷为 **(k_r-k_m-i alpha_m)**（alpha_m 虚部符号相反）。
2. a_m 分子印刷为 r2 项减 r1 项；在 k_r=k_m 处给出 **-b_m**（与 Eq.5 差符号）。
   若分子改为 r1-r2，峰值符号才与 Eq.5 一致。

**结论：**
- Eq3->Eq4(第一形式)->Eq5：**PASS / CLOSED**
- Eq.(4) 第二形式与第一形式：**PAPER_INTERNAL_INCONSISTENT**（疑似排版笔误）
- **实现必须使用第一形式或 Eq.(5)**，不得按第二形式未校正符号直接编码
- 本项 **不是** 我方转录冲突（视觉转录与印刷一致）

## 检查清单映射

| 项 | 结果 |
| --- | --- |
| Eq1 steering 符号 | PASS（e^{+i k_r r}） |
| Eq1 前因子 e^{i pi/4}/sqrt(2 pi k_r) | PASS |
| Eq3 传播相位/衰减 | PASS（e^{-i k_m r - alpha_m r - i pi/4}） |
| S(r) 在积分分子 | PASS |
| sqrt(k_r k_m) 因子 | PASS |
| Eq5 b_m 与峰值积分 | PASS |
| Eq4 第二形式自洽 | FAIL（论文内部，见上） |
"""
    )
    (OUT / "YANG2015_EQ1_EQ5_ALGEBRA_CHECK.md").write_text(body, encoding="utf-8")


def write_bm() -> None:
    body = T(
        r"""# YANG2015_BM_DEFINITION

UTC: __NOW__

## 主文 Eq.(5) 印刷定义（视觉锁定）

$$
b_m = \frac{2 e^{-\alpha_m r_0}}{\alpha_m k_m} \sinh\left[\frac{\alpha_m \Delta R}{2}\right] \phi_m(z_s)
$$

- $r_0 \equiv (r_2+r_1)/2$
- $\Delta R = r_2-r_1 = vT$（合成孔径径向跨度）
- $b_m$ 为 **实数**（论文原文：is a real number）

## 因子分解核对

| 因子 | 是否在 b_m 中 |
| --- | --- |
| $\phi_m(z_s)$ | **YES** |
| $\alpha_m$ | **YES**（$e^{-\alpha_m r_0}$ 与 $\sinh(\alpha_m \Delta R/2)$） |
| $k_m$ | **YES**（分母 $\alpha_m k_m$） |
| $r_0$ | **YES** |
| $\Delta R$ | **YES** |
| source-level / 源幅度 | **NO** |
| $1/\sqrt{k_m}$ 独立因子 | **NO**（并入 $1/k_m$，来自 $\sqrt{k_r k_m}|_{k_r=k_m}$） |
| $\phi_m(z_r)$ | **NO**（在 $g=b_m\phi_m(z_r)$ 外侧） |
| range spreading $\sqrt{r}$ | 已在 Eq.(1) 的 S(r)~sqrt(r) 中处理，不在 b_m |

## 与 delta_m（Eq.9）关系

$$
\delta_m = \frac{2 e^{-\alpha_m r_0}}{\alpha_m k_m}\sinh[(\alpha_m \Delta R)/2] = b_m / \phi_m(z_s)
$$

即 **mode shading coefficient**（Ref.13）。

## 禁止

不得再用 $A_m=\phi_m(z_s)\phi_m(z_r)$ 代替 Yang 的 $b_m$。
$A_m$ 仅是去掉 $k_m,\alpha_m,r_0,\Delta R$ 因子后的简化深度乘积，**不是** Eq.(5) 系数。

## Erratum（carry-forward）

$\delta$ 初始距离偏移时：$g(k_m,z_r)=b_m\phi_m(z_r)e^{i k_m \delta}$。

**PAPER_REPRODUCTION_BASELINE: delta = 0**（本轮不重核 Erratum 全文）。
"""
    )
    (OUT / "YANG2015_BM_DEFINITION.md").write_text(body, encoding="utf-8")


def write_eq6_observable() -> None:
    body = T(
        r"""# YANG2015_EQ6_OBSERVABLE_NOTE

UTC: __NOW__

## Eq.(6) 印刷形式

$$
D(z)=\left|\sum_{m=1}^{M} \phi_m(z)\,\frac{g(k_m,z_r)}{\bar{\phi}_m(z_r)}\right|^2
$$

正则化逆（**视觉锁定**）：

$$
\bar{\phi}_m^{-1}(z_r)=\frac{\phi_m(z_r)}{\phi_m^2(z_r)+\Delta^2}
$$

故

$$
\frac{g(k_m,z_r)}{\bar{\phi}_m(z_r)} = g(k_m,z_r)\cdot\bar{\phi}_m^{-1}(z_r)
= g(k_m,z_r)\,\frac{\phi_m(z_r)}{\phi_m^2(z_r)+\Delta^2}
$$

## observable 到底是什么

| 问题 | 锁定答案 |
| --- | --- |
| 输入是复值还是幅值？ | **复值谱峰** $g(k_m,z_r)$ 进入求和后再取模平方（Bartlett） |
| 是否只用 spectrum amplitudes？ | 附录原文：*Eq.(6) uses only the spectrum amplitudes, not the corresponding mode wavenumbers* |
| $k_m$ 数值是否进深度 score？ | **NO**。$k_m$ 只用于 **mode identification**（峰<->模态编号） |
| 印刷记号 $g(k_m,z_r)$ | 表示在模态波数峰处的谱值；实现取该峰处的复谱 |

**禁止** `sqrt(member_energy)` 或任何自造 score。

## 模态编号错误的影响

论文 Sec.III：

- 编号 $n$ 与 $m$ 错配 -> 深度误差量级 **$H/n - H/m$**（$H$ 为水深）
- 实测：用模拟谱 + 相邻模态波数差及其增长率辅助；作者承认 **educated guess**
- 可通过微调 residual Doppler 平移谱横轴使数据谱与模型谱对齐

## 正则化 Delta

| 项 | 值 |
| --- | --- |
| 符号 | $\Delta$ |
| 精确定义 | $\bar{\phi}_m^{-1}=\phi_m/(\phi_m^2+\Delta^2)$ |
| 论文推荐量级 | *on the order of one tenth of the maximum value of the mode depth functions* |
| 是否唯一数值 | **否** -> 状态 **`EMPIRICAL_NOT_UNIQUE`** |
| 作用 | 防止 $\phi_m(z_r)$ 零交叉处逆模态函数爆炸 |

**不得**把 0.1 写成论文固定参数。

## 深度分布归一化（Sec.II 原文）

normalized depth distribution = Eq.(6) 的 $D(z)$ **除以搜索深度网格上 $D(z)$ 的总和**。

这是论文展示用归一化分布，**不是**另造 correlation score。
"""
    )
    (OUT / "YANG2015_EQ6_OBSERVABLE_NOTE.md").write_text(body, encoding="utf-8")


def write_repro_and_ref7() -> None:
    body = T(
        r"""# YANG2015_PAPER_REPRO_CONFIG

UTC: __NOW__

## DIRECTLY_SPECIFIED_IN_YANG2015（论文模拟）

| 参数 | 值 | 出处 |
| --- | --- | --- |
| source depth z_s | 4 m（shallow）, 50 m（deep） | Sec.II |
| receiver depth z_r | 18 m, 70 m | Sec.II |
| frequency | 350 Hz | Sec.II |
| speed | 5 kn | Sec.II |
| motion | moving **away** from receiver | Sec.II |
| initial range r1 | 5010 m（模拟生成用；**接收端算法不把 r1 当已知输入**） | Sec.II |
| sampling | 1 s | Sec.II |
| synthetic range span | **约 4990 m** | Sec.II |
| VLA demo | 20 elements, dz=4 m, covering 2-80 m | Sec.IV B |
| SWellEx96 real data | 127 Hz zs=9 m; 130 Hz zs=54 m; z_r=bottom phone; ~35 min | Sec.III |

## REFERENCED_TO_REF7

Ref.7 = T. C. Yang, *Data-based matched-mode source localization for a moving source*, JASA **135**, 1218-1230 (2014).

| 项目 | 状态 |
| --- | --- |
| downward refractive SSP（具体剖面） | **NEEDS_REF7** |
| water depth | **NEEDS_REF7** |
| bottom properties | **NEEDS_REF7** |
| mode attenuation / model details | **NEEDS_REF7** |
| 其它传播配置 | **NEEDS_REF7** |
| Eq.(1) 在 Ref.7 中称 generalized Hankel transform | 已记录 |

**PAPER_REPRO_ENV = `PAPER_REPRO_NEEDS_REF7`**

## 模态峰识别（必须分层）

| 标签 | 定义 | 论文依据 |
| --- | --- | --- |
| **ORACLE_MODE_ID** | 用模型 true k_m 标注谱峰 | Fig.1 "+" true modal wavenumbers |
| **PRACTICAL_MODE_ID** | mode spacing + 相邻 spacing 增长率 + 模拟谱辅助；含 educated guess | Sec.III 实测流程 |

**禁止混为一谈。** 忠实复现论文仿真用 ORACLE_MODE_ID。

## 孔径条件（Sec.IV C）

- 分辨 n/m 阶模态：L > lambda_nm = 2 pi / |k_n-k_m|
- 示例：lambda_12 约 2.3 km（**论文环境示例**）
- 1 km 短孔径 -> 论文建议高分辨算法（Ref.12 AR 等）

**2.3 km 不得直接迁移为 E-STD 固定阈值**；E-STD 须用 KRAKEN k_m 重算。

## PAPER_REPRODUCTION_BASELINE

- **delta = 0**（Erratum 修正 carry-forward；无 Erratum PDF 本轮重核）
- 第一轮忠实复现 **不得**把 2.4 km E-STD 场景塞进论文条件
- 论文孔径 **约 4.99 km** 优先
"""
    )
    (OUT / "YANG2015_PAPER_REPRO_CONFIG.md").write_text(body, encoding="utf-8")

    body2 = T(
        r"""# YANG2015_REF7_DEPENDENCY

UTC: __NOW__

Yang2015 明确：*simulated data given in Ref.7 for a downward refractive SSP. (The readers are referred to Ref.7 for details.)*

## 完全复现论文模拟仍缺（Yang2015 未给，不得猜）

1. 具体 SSP 数值剖面（Fig.1 of Ref.7）
2. 水深 H
3. 底质声学参数（密度、声速、衰减）
4. 模态求解/衰减模型细节
5. 产生 p(r,z) 的传播配置（采样几何、源级归一等）
6. Fig.1 所用 true k_m 表（ORACLE 标注）

## 允许状态

- `PAPER_REPRO_NEEDS_REF7` <- **当前**
- `PAPER_REPRO_ENV_COMPLETE`（取得 Ref.7 并锁定后）

## 边界

- Eq.(1)-(10)、A1-A4 **不依赖** Ref.7 即可锁定（已完成）
- 论文条件复现数值结果 **依赖** Ref.7
- 即使 NEEDS_REF7，**YANG_ROUTE 仍为 UNDECIDED**，不得 pass
"""
    )
    (OUT / "YANG2015_REF7_DEPENDENCY.md").write_text(body2, encoding="utf-8")


def write_impl_spec() -> None:
    body = T(
        r"""# YANG2015_IMPLEMENTATION_SPEC

UTC: __NOW__

**本轮禁止运行。** 仅伪代码规范。

## INPUT

- moving CW pressure blocks p(r_n,z_r)（或原始时域走 A1-A4）
- 可选：nominal 环境 -> phi_m(z)
- 可选：true k_m（ORACLE_MODE_ID）

## 单水听器主路径

```
# 0) 压力块（若从时域）
raw time series
  -> [A2] Doppler / range increments r_n = r1 + n * dr
  -> [Eq.7] PLL + low-pass  (real data)
  -> [A3/A4] block pressure p(r_n, z_r)

# 1) shading
S(r_n) = <|p|^2>^{-1/2}          # Eq.(2)

# 2) SA spectrum
g(k_r, z_r) = e^{i pi/4}/sqrt(2 pi k_r)
              * sum_n p(r_n,z_r) e^{i k_r r_n} S(r_n) dr   # Eq.(1)
# 实现等价：Eq.(4) 第一形式积分，禁止未校正的 Eq.(4) 第二形式

# 3) modal peaks
detect peaks in |g(k_r,.)|
mode identification:
  ORACLE_MODE_ID     = match true k_m          # paper simulation
  PRACTICAL_MODE_ID  = spacing + sim guide     # real-data style
g_m <- complex peak value g(k_m, z_r)

# 4) regularized inverse at receiver
phi_m(z_r);  phi_bar_inv = phi_m(z_r) / (phi_m(z_r)^2 + Delta^2)
# Delta = EMPIRICAL_NOT_UNIQUE ~ 0.1 * max|phi_m|  (order-of-magnitude only)

# 5) depth ambiguity Eq.(6)
for z in search_grid:
    D(z) = | sum_m phi_m(z) * g_m * phi_bar_inv_m |^2

# 6) normalized distribution (paper display convention)
Dn(z) = D(z) / sum_z D(z)

# 7) estimate
z_hat = argmax_z D(z)    # or Dn; same argmax
```

## VLA 路径（Sec.IV B）

```
S(r) = < sum_j |p(r,z_j)|^2 dz >^{-1/2}     # Eq.(8)
per phone j: g_j = Eq.(1)
D(z_j) = | sum_m g(k_m, z_j) |^2            # Eq.(9)
P(z_j|z_s) = K D(z_j), K = [sum_j D(z_j)]^{-1}
z_hat = argmax_j D(z_j)
# 不需要 phi_m / SSP
```

## 明确禁止（本轮）

- 运行任何 Yang / E-STD 仿真
- 输出 z_hat / FWHM / PSL / 性能 pass-fail
- delta 扫描、Doppler 误差、MC、AR、P5
- 关闭任何候选算法
- 使用 A_m=phi_s phi_r 或 sqrt(member_energy) 作为 Yang score
"""
    )
    (OUT / "YANG2015_IMPLEMENTATION_SPEC.md").write_text(body, encoding="utf-8")


def write_checklist() -> None:
    rows = [
        ("Eq1_sign_checked", "PASS", "steering e^{+i k_r r}, prefactor e^{i pi/4}/sqrt(2 pi k_r) visual ok"),
        ("Eq1_weight_checked", "PASS", "S(r) in integrand numerator; S=<|p|^2>^{-1/2}; S~sqrt(r)->Hankel"),
        ("Eq3_phase_checked", "PASS", "e^{-i k_m r - alpha_m r - i pi/4} and sqrt(2 pi/(k_m r)) visual ok"),
        ("Eq4_algebra_checked", "PASS", "Eq3->Eq1->Eq4 first form->Eq5 closed; second form see note"),
        ("Eq4_second_form_consistent", "FAIL", "printed second form vs first form (alpha sign / a_m numerator order)"),
        ("Eq5_bm_checked", "PASS", "b_m=(2e^{-alpha r0})/(alpha k_m)sinh[(alpha dR)/2]phi_s matches peak integral"),
        ("Eq6_regularization_checked", "PASS", "phi_bar_inv=phi_m/(phi_m^2+Delta^2); Delta=EMPIRICAL_NOT_UNIQUE"),
        ("Eq6_observable_checked", "PASS", "complex peaks g(k_m,z_r) in Bartlett sum; k_m not in depth score"),
        ("depth_normalization_checked", "PASS", "D(z)/sum_grid D(z) is paper normalized distribution"),
        ("mode_identification_role_checked", "PASS", "ORACLE_MODE_ID vs PRACTICAL_MODE_ID separated"),
        ("delta_zero_baseline_checked", "PASS", "PAPER_REPRODUCTION_BASELINE delta=0; erratum carry-forward"),
        ("paper_span_4990m_checked", "PASS", "paper aperture ~4990 m; 2.4 km excluded from first repro"),
        ("Ref7_dependency_checked", "PASS", "gap listed; status PAPER_REPRO_NEEDS_REF7"),
        ("appendix_vs_eq3_phase", "UNRESOLVED", "A1-A3 spatial phase +i k_m r vs Eq.(3) -i k_m r; main chain uses Eq.(3)"),
        ("eq7_pll_coefficients", "PASS", "K1=0.001, K2=0.0001 confirmed"),
    ]
    pd.DataFrame(
        [{"item": a, "status": b, "note": c} for a, b, c in rows]
    ).to_csv(OUT / "YANG2015_IMPLEMENTATION_CHECKLIST.csv", index=False, encoding="utf-8-sig")

    unresolved = [r for r in rows if r[1] == "UNRESOLVED"]
    failed = [r for r in rows if r[1] == "FAIL"]
    (OUT / "YANG2015_IMPLEMENTATION_CHECKLIST_STATUS.md").write_text(
        "# checklist status\n\n"
        f"- FAIL: {len(failed)} -> {[r[0] for r in failed]}\n"
        f"- UNRESOLVED: {len(unresolved)} -> {[r[0] for r in unresolved]}（**禁止自动忽略**）\n"
        "- 其余 PASS\n\n"
        "UNRESOLVED 处理：记入实现规范相位约定决策；**不**静默丢弃。\n"
        "FAIL 处理：实现禁用 Eq.(4) 第二形式未校正编码；使用第一形式/Eq.(5)。\n",
        encoding="utf-8",
    )


def write_report_decision() -> None:
    decision = "C2_2B_METHOD_SPEC_LOCKED"
    why = (
        "Yang 2015 主文 PDF 已由用户提供并完成视觉+文本双通道锁定：Eq.(1)-(10) 与附录 A1-A4 全部写入；"
        "Eq3->Eq4(第一形式)->Eq5 代数闭环 PASS；Eq.(6) 复谱峰 observable 与正则化逆 "
        "phi_bar_inv=phi_m/(phi_m^2+Delta^2) 已锁定，Delta=EMPIRICAL_NOT_UNIQUE；归一化 D/sum D 已锁定。"
        "论文复现条件已分 DIRECTLY_SPECIFIED / REFERENCED_TO_REF7；环境完整复现状态 PAPER_REPRO_NEEDS_REF7。"
        "Eq.(4) 第二形式与附录空间相位存在论文内部不一致，已记录且实现规避。"
        "未运行任何算法。YANG_ROUTE_UNDECIDED。"
    )

    dec = {
        "stage": "R3-C2.2B",
        "rc3c2_2b_decision": decision,
        "why": why,
        "three_state": {
            "kraken_data_baseline": "C2_1_PARSER_VALIDATED",
            "yang_formula_recovery": "PRIMARY_SOURCE_LOCKED (user PDF)",
            "yang_route_in_E_STD": "YANG_ROUTE_UNDECIDED",
        },
        "supersedes": {
            "C2_2A_PRIMARY_FORMULA_NOT_RECOVERED": "SUPERSEDED_BY_USER_PRIMARY_PDF",
        },
        "paper_repro_env": "PAPER_REPRO_NEEDS_REF7",
        "paper_reproduction_baseline": {
            "delta": 0,
            "range_span_m": 4990,
            "freq_hz": 350,
            "speed_kn": 5,
        },
        "mode_id_layers": ["ORACLE_MODE_ID", "PRACTICAL_MODE_ID"],
        "regularization": {
            "symbol": "Delta",
            "status": "EMPIRICAL_NOT_UNIQUE",
            "order": "0.1 * max|phi_m|",
        },
        "eq6_observable": "complex modal peaks g(k_m,z_r); wavenumbers NOT in depth score",
        "known_paper_internal_issues": [
            "Eq.(4) second form vs first form (alpha sign / a_m numerator order)",
            "Appendix A1-A3 spatial phase +i k_m r vs Eq.(3) e^{-i k_m r}",
        ],
        "checklist_unresolved": ["appendix_vs_eq3_phase"],
        "checklist_failed": ["Eq4_second_form_consistent"],
        "not_done": [
            "run Yang",
            "run E-STD",
            "z_hat",
            "FWHM",
            "PSL",
            "delta scan",
            "Doppler error",
            "MC",
            "AR",
            "P5",
            "close any candidate",
        ],
        "primary_pdf": {
            "path": "literature/yang2015_Yang2015_JASA138_1678.pdf",
            "sha256": SHA,
            "pages": NPAGES,
            "doi": "10.1121/1.4929748",
        },
        "created_utc": NOW,
    }
    (OUT / "R3_C2_2B_DECISION.json").write_text(
        json.dumps(dec, ensure_ascii=False, indent=2), encoding="utf-8"
    )

    report = T(
        r"""# R3-C2.2B 报告：Yang 2015 主文公式锁定 + 复现协议锁定

UTC：__NOW__

## 0. 历史状态

- `C2_2A_PRIMARY_FORMULA_NOT_RECOVERED` -> **SUPERSEDED_BY_USER_PRIMARY_PDF**
- KRAKEN：`C2_1_PARSER_VALIDATED`
- 公式恢复：**PRIMARY_SOURCE_LOCKED**
- `YANG_ROUTE_UNDECIDED`（**不变**）
- 本轮 **禁止** Yang 性能 pass/fail

## 1. PDF 归档

见 `YANG2015_PRIMARY_PDF_NOTE.md`。SHA256 `__SHA__`，10 页，JASA 1678-1686。

## 2-3. Eq.(1)-(10) 与 A1-A4

见 `YANG2015_EQUATION_LOCK.md`、`YANG2015_EQUATION_TABLE_LOCKED.csv`、`YANG2015_APPENDIX_EQUATIONS.md`。

全部乘法因子按原页锁定（含 e^{i pi/4}、sqrt(2 pi/(k_m r))、alpha_m、sinh、正则化逆）。

## 4. 代数闭环

见 `YANG2015_EQ1_EQ5_ALGEBRA_CHECK.md`。

**Eq3->Eq4(第一形式)->Eq5：PASS。**
Eq.(4) 第二形式与第一形式不一致（论文内部），实现禁用未校正第二形式。

## 5-8. b_m / Eq.(6) / 正则化 / 归一化

见 `YANG2015_BM_DEFINITION.md`、`YANG2015_EQ6_OBSERVABLE_NOTE.md`。

- b_m=(2 e^{-alpha_m r0})/(alpha_m k_m) sinh[(alpha_m dR)/2] phi_m(z_s)
- Eq.(6) 用 **复谱峰**；**不用**波数数值进深度 score
- phi_bar_inv = phi_m/(phi_m^2+Delta^2)，Delta **EMPIRICAL_NOT_UNIQUE**
- 归一化 D(z)/sum_grid D(z)

## 9. Erratum

`CARRY_FORWARD_FROM_PRIOR_ERRATUM_RECOVERY`；**delta=0** 基线。本轮无 Erratum PDF，未声称重核全文。

## 10-13. 复现配置 / Ref.7 / 峰识别 / 孔径

见 `YANG2015_PAPER_REPRO_CONFIG.md`、`YANG2015_REF7_DEPENDENCY.md`。

- 论文模拟：350 Hz、5 kn、z_s=4/50、z_r=18/70、span~**4990 m**、r1=5010（仅生成）
- **ORACLE_MODE_ID** vs **PRACTICAL_MODE_ID** 分层
- L > 2 pi/|k_n-k_m|，示例 lambda_12~2.3 km（示例，非 E-STD 阈值）
- 环境：**PAPER_REPRO_NEEDS_REF7**

## 14-15. 实现规范与自检

见 `YANG2015_IMPLEMENTATION_SPEC.md`、`YANG2015_IMPLEMENTATION_CHECKLIST.csv`。

FAIL：`Eq4_second_form_consistent`（已规避）。
UNRESOLVED：`appendix_vs_eq3_phase`（已记录，主链按 Eq.(3)）——**未忽略**。

## 16. 判定

### __DECISION__

__WHY__

即使 `PAPER_REPRO_NEEDS_REF7`，**YANG_ROUTE 仍为 UNDECIDED，不得 pass**。

## 17. 停止声明

不跑 Yang / E-STD / z_hat / FWHM / PSL / delta 扫描 / Doppler 误差 / MC / AR / P5；不关闭任何候选。
"""
    )
    report = report.replace("__DECISION__", decision).replace("__WHY__", why)
    (OUT / "R3_C2_2B_REPORT.md").write_text(report, encoding="utf-8")

    sync = T(
        """# R3-C2.2B GPT SYNC

UTC：__NOW__

## 判定

**__DECISION__**

## 三态

- KRAKEN：VALIDATED
- Yang 公式：PRIMARY_SOURCE_LOCKED（用户 PDF）
- Yang 路线：UNDECIDED

## 关键锁定

1. Eq.(1) SA 波束：e^{i pi/4}/sqrt(2 pi k_r) * int p e^{i k_r r} S dr
2. Eq.(5) b_m=(2e^{-alpha r0})/(alpha k_m)sinh(alpha dR/2)phi_s
3. Eq.(6) D(z)=|sum phi_m(z) g phi_m(z_r)/(phi_m^2+Delta^2)|^2，Delta 非唯一
4. observable=复谱峰；k_m 不进深度 score
5. 论文孔径~4990 m；Ref.7 缺口 -> PAPER_REPRO_NEEDS_REF7
6. delta=0 基线

## 已知论文内部问题

- Eq.(4) 第二形式符号
- 附录 vs Eq.(3) 空间相位

## 未做

算法运行、性能判定、关闭候选。

## 下一阶段（未执行）

论文原条件复现 -> 单因素自检 -> E-STD 迁移（需另授权）。
"""
    ).replace("__DECISION__", decision)
    (OUT / "R3_C2_2B_GPT_SYNC.md").write_text(sync, encoding="utf-8")
    print("DECISION", decision)


def main() -> None:
    write_pdf_note()
    write_equation_lock()
    write_appendix()
    write_algebra_check()
    write_bm()
    write_eq6_observable()
    write_repro_and_ref7()
    write_impl_spec()
    write_checklist()
    write_report_decision()


if __name__ == "__main__":
    main()
