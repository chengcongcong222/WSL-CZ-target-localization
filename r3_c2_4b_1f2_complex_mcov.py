#!/usr/bin/env python3
"""R3-C2.4B-1F2: complex modified-covariance correction + p=1 identity gate."""
from __future__ import annotations

import json
from datetime import datetime, timezone
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parent
OUT = ROOT / "results" / "R3_C2_AR_MMAR" / "R3_C2_4B_1F2_COMPLEX_MCOV"
OUT.mkdir(parents=True, exist_ok=True)
NOW = datetime.now(timezone.utc).isoformat()

K_BAND = (1.32, 1.50)
K_PHYS = [1.32, 1.41, 1.50]
N_R = 200
P_LIT = 7
P_MOVING = int(np.floor(2 * N_R / 3))  # 133
TWO_PI = 2 * np.pi


def complex_mcov_ar(y: np.ndarray, p: int) -> dict:
    """Complex modified covariance.

    e_f[n] = y[n] + sum_k a[k] y[n-k]
    e_b*[n] = conj(y[n-p]) + sum_k a[k] conj(y[n-p+k])
    """
    y = np.asarray(y, dtype=complex)
    N = y.size
    Af, bf, Ab, bb = [], [], [], []
    for n in range(p, N):
        Af.append(np.array([y[n - k] for k in range(1, p + 1)]))
        bf.append(-y[n])
        # backward conjugated block
        Ab.append(np.array([np.conj(y[n - p + k]) for k in range(1, p + 1)]))
        bb.append(-np.conj(y[n - p]))
    Af = np.asarray(Af)
    bf = np.asarray(bf)
    Ab = np.asarray(Ab)
    bb = np.asarray(bb)
    A = np.vstack([Af, Ab])
    b = np.concatenate([bf, bb])
    rank = int(np.linalg.matrix_rank(A))
    cond = float(np.linalg.cond(A))
    a, *_ = np.linalg.lstsq(A, b, rcond=None)
    u = np.array(
        [y[n] + np.dot(a, np.array([y[n - k] for k in range(1, p + 1)])) for n in range(p, N)]
    )
    sigma2 = float(np.mean(np.abs(u) ** 2))
    return {
        "a": a,
        "sigma2": sigma2,
        "n_forward": N - p,
        "n_backward": N - p,
        "n_equations": int(A.shape[0]),
        "n_params": p,
        "rank": rank,
        "cond": cond,
        "backward_is_conj": True,
        "data_model_rank_limit": rank < p,
    }


def ar_spectrum(a, sigma2, omega, power=1):
    a = np.asarray(a, dtype=complex)
    acc = np.ones_like(omega, dtype=complex)
    for k in range(1, a.size + 1):
        acc = acc + a[k - 1] * np.exp(-1j * omega * k)
    return sigma2 / np.maximum(np.abs(acc) ** power, 1e-30)


def detect_peaks(P, omega, rel_prom=0.05):
    n = P.size
    prom = rel_prom * float(P.max())
    cand = [i for i in range(1, n - 1) if P[i] >= P[i - 1] and P[i] >= P[i + 1] and P[i] > 0]
    scored = []
    for i in cand:
        l = i
        while l > 0 and P[l - 1] <= P[i]:
            l -= 1
        r = i
        while r < n - 1 and P[r + 1] <= P[i]:
            r += 1
        pr = P[i] - max(P[l : i + 1].min(), P[i : r + 1].min())
        if pr >= prom:
            scored.append(i)
    scored.sort(key=lambda i: -P[i])
    keep = []
    for i in scored:
        if all(abs(i - j) > 3 for j, _ in keep):
            keep.append((i, float(P[i])))
    keep.sort()
    idx = [i for i, _ in keep]
    return idx, np.array([float(omega[i]) for i in idx])


def wrap(x):
    return (x + np.pi) % TWO_PI - np.pi


def map_k_from_omega(w, dr, k_band):
    """k = -(omega + 2π m)/dr  (frozen analytic); pick branches in physical band."""
    cands = [-(w + TWO_PI * m) / dr for m in range(-8, 9)]
    in_band = sorted({round(c, 9) for c in cands if k_band[0] - 1e-9 <= c <= k_band[1] + 1e-9})
    return in_band


def gen_tone(ks, dr, N, amps=None, noise=0.0, seed=0):
    n = np.arange(N)
    y = np.zeros(N, dtype=complex)
    amps = amps or [1.0] * len(ks)
    for k, A in zip(ks, amps):
        y = y + A * np.exp(-1j * k * n * dr)
    if noise:
        rng = np.random.default_rng(seed)
        y = y + (noise / np.sqrt(2)) * (rng.standard_normal(N) + 1j * rng.standard_normal(N))
    return y


def peak_scores(ws, dr, k_true_list, tol=1e-3):
    all_k = []
    for w in ws:
        all_k.extend(map_k_from_omega(float(w), dr, K_BAND))
    all_k = sorted(set(all_k))
    hits = sum(1 for kt in k_true_list if any(abs(k - kt) < tol for k in all_k))
    false = sum(1 for k in all_k if all(abs(k - kt) > tol for kt in k_true_list))
    recall = hits / max(len(k_true_list), 1)
    precision = hits / max(len(all_k), 1) if all_k else 0.0
    f1 = 2 * precision * recall / max(precision + recall, 1e-30)
    sel = None
    return {
        "n_detected_peaks_in_band": len(all_k),
        "true_peak_recalled": int(hits),
        "false_peak_count": int(false),
        "recall": float(recall),
        "precision": float(precision),
        "f1": float(f1),
        "peak_count_err": int(len(all_k) - len(k_true_list)),
        "recovered_ks": ";".join(f"{k:.6f}" for k in all_k),
    }


def main() -> int:
    # ---- p=1 identity gate ----
    p1_rows = []
    gate_ok = True
    for k in K_PHYS:
        for dr in (0.5, 1.885, 2.5):
            y = gen_tone([k], dr, 80)
            fit = complex_mcov_ar(y, 1)
            a1 = fit["a"][0]
            a1_exact = -np.exp(-1j * k * dr)
            coef_err = float(abs(a1 - a1_exact))
            omega = np.linspace(-np.pi, np.pi, 8192, endpoint=False)
            P = ar_spectrum(fit["a"], fit["sigma2"], omega, power=1)
            idx, ws = detect_peaks(P, omega, rel_prom=0.05)
            w_peak = float(ws[np.argmax(P[idx])]) if idx else np.nan
            w_expect = wrap(-k * dr)
            peak_err = float(abs(w_peak - w_expect)) if w_peak == w_peak else np.nan
            ok = bool(coef_err < 1e-6 and peak_err == peak_err and peak_err < 1e-3)
            gate_ok = gate_ok and ok
            p1_rows.append(
                {
                    "k_true": k,
                    "dr": dr,
                    "a1_re": float(a1.real),
                    "a1_im": float(a1.imag),
                    "a1_exact_re": float(a1_exact.real),
                    "a1_exact_im": float(a1_exact.imag),
                    "complex_coef_err": coef_err,
                    "omega_peak": w_peak,
                    "omega_expected_wrap_minus_kdr": w_expect,
                    "peak_err_rad": peak_err,
                    "P1_OK": ok,
                    "rank": fit["rank"],
                    "cond": fit["cond"],
                }
            )
    pd.DataFrame(p1_rows).to_csv(OUT / "P1_COMPLEX_EXP_IDENTITY.csv", index=False)
    p1_status = "P1_COMPLEX_EXP_IDENTITY_VALIDATED" if gate_ok else "P1_COMPLEX_EXP_IDENTITY_FAILED"

    (OUT / "COMPLEX_MCOV_FORMULA_AUDIT.md").write_text(
        f"""# COMPLEX_MCOV_FORMULA_AUDIT

UTC: {NOW}

`COMPLEX_BACKWARD_CONJUGATION_REQUIRED`

```text
Af[n,:] = [y[n-1], ..., y[n-p]]
bf[n]   = -y[n]
Ab[n,:] = [conj(y[n-p+1]), ..., conj(y[n])]
bb[n]   = -conj(y[n-p])
```

对应 e_b*[n] = y*[n-p] + Σ a[k] y*[n-p+k]（backward 用 a* 后整体共轭）。

p=1 恒等式：y=e^{{-jknΔr}} ⇒ a1=−e^{{-jkΔr}}；ω_peak=wrap(−kΔr)。

Gate：**{p1_status}**
""",
        encoding="utf-8",
    )

    (OUT / "SIGN_MAPPING_FINAL.md").write_text(
        f"""# SIGN_MAPPING_FINAL

UTC: {NOW}

**不再 majority vote。** 由 p=1 解析恒等式 + 数值验证冻结：

$$
\\omega_{{\\rm peak}}=\\mathrm{{wrap}}(-k\\Delta r)
\\quad\\Rightarrow\\quad
k=-\\frac{{\\omega+2\\pi m}}{{\\Delta r}}
$$

分支：物理带 [1.32,1.50] 选取合法 m。

p=1 gate：**{p1_status}**
""",
        encoding="utf-8",
    )

    # ---- performance after gate ----
    if not gate_ok:
        decision = "AR_K_RECOVERY_UNIT_FAILED"
        why = "P1_COMPLEX_EXP_IDENTITY_FAILED — 不跑 p=7/133 性能"
        sdf = pd.DataFrame()
        mdf = pd.DataFrame()
    else:
        p_branches = [(7, "PRINTED_LITERAL_CONTROL"), (P_MOVING, "OUR_MOVING_SAMPLE_INTERPRETATION")]
        single_rows = []
        for k in K_PHYS:
            for dr, dr_tag in [(1.885, "NO_FOLD_CONTROL"), (0.5, "AUX"), (2.5, "REF7_2p5")]:
                y = gen_tone([k], dr, N_R)
                for p, p_tag in p_branches:
                    fit = complex_mcov_ar(y, p)
                    omega = np.linspace(-np.pi, np.pi, 16384, endpoint=False)
                    for power, pow_tag in [
                        (1, "PRINTED_EQ19_MODULUS_POWER_1"),
                        (2, "STANDARD_AR_PSD_MODULUS_POWER_2_CONTROL"),
                    ]:
                        P = ar_spectrum(fit["a"], fit["sigma2"], omega, power=power)
                        idx, ws = detect_peaks(P, omega, rel_prom=0.05)
                        met = peak_scores(ws, dr, [k])
                        sel_k = np.nan
                        sel_err = np.nan
                        if idx:
                            bi = int(np.argmax(P[idx]))
                            in_b = map_k_from_omega(ws[bi], dr, K_BAND)
                            if in_b:
                                sel_k = in_b[0]
                                sel_err = abs(sel_k - k)
                        folds = abs(k * dr) > np.pi
                        in_true = map_k_from_omega(wrap(-k * dr), dr, K_BAND)
                        if not folds:
                            alias = "NO_FOLD_DIRECT_RECOVERY"
                        elif len(in_true) == 1:
                            alias = "FOLDED_BUT_BRANCH_RECOVERABLE"
                        else:
                            alias = "FOLDED_BRANCH_AMBIGUOUS"
                        met.update(
                            {
                                "k_true": k,
                                "dr": dr,
                                "dr_tag": dr_tag,
                                "p": p,
                                "p_tag": p_tag,
                                "power": power,
                                "power_tag": pow_tag,
                                "selected_peak_k": sel_k,
                                "selected_abs_err": sel_err,
                                "alias_class": alias,
                                "rank": fit["rank"],
                                "cond": fit["cond"],
                                "data_model_rank_limit": fit["data_model_rank_limit"],
                            }
                        )
                        single_rows.append(met)
        sdf = pd.DataFrame(single_rows)
        sdf.to_csv(OUT / "SINGLE_TONE_RECOVERY_FINAL.csv", index=False)

        multi_defs = [
            ("2_peak_wide", [1.35, 1.48]),
            ("3_peak_wide", [1.33, 1.41, 1.49]),
            ("5_peak_wide", [1.32, 1.36, 1.41, 1.46, 1.50]),
            ("2_peak_close", [1.40, 1.405]),
            ("3_peak_close", [1.40, 1.404, 1.408]),
            ("5_peak_close", [1.38, 1.39, 1.40, 1.41, 1.42]),
        ]
        multi_rows = []
        for name, ks in multi_defs:
            for dr, dr_tag in [(1.885, "NO_FOLD"), (2.5, "REF7_2p5")]:
                y = gen_tone(ks, dr, N_R)
                for p, p_tag in p_branches:
                    fit = complex_mcov_ar(y, p)
                    for power, pow_tag in [(1, "PRINTED_EQ19_POWER_1"), (2, "STANDARD_PSD_POWER_2")]:
                        omega = np.linspace(-np.pi, np.pi, 16384, endpoint=False)
                        P = ar_spectrum(fit["a"], fit["sigma2"], omega, power=power)
                        idx, ws = detect_peaks(P, omega, rel_prom=0.05)
                        met = peak_scores(ws, dr, ks)
                        met.update(
                            {
                                "case": name,
                                "n_true": len(ks),
                                "dr_tag": dr_tag,
                                "p_tag": p_tag,
                                "power": power,
                                "power_tag": pow_tag,
                                "complete_success": bool(met["recall"] == 1.0 and met["false_peak_count"] == 0),
                                "rank": fit["rank"],
                                "cond": fit["cond"],
                            }
                        )
                        multi_rows.append(met)
        mdf = pd.DataFrame(multi_rows)
        mdf.to_csv(OUT / "MULTITONE_RECOVERY_FINAL.csv", index=False)

        prim = sdf[sdf["power"] == 1]
        f1_st = float(prim["f1"].mean())
        rec_st = float(prim["recall"].mean())
        prec_st = float(prim["precision"].mean())
        multi_ok = float(mdf[mdf["power"] == 1]["complete_success"].mean())
        multi_f1 = float(mdf[mdf["power"] == 1]["f1"].mean())

        # numerical conditioning diagnostic: tiny complex gaussian noise, fixed seed
        cond_rows = []
        for p, p_tag in p_branches:
            for dr in (1.885, 2.5):
                y0 = gen_tone([1.41], dr, N_R)
                fit0 = complex_mcov_ar(y0, p)
                w0 = np.linspace(-np.pi, np.pi, 4096, endpoint=False)
                P0 = ar_spectrum(fit0["a"], fit0["sigma2"], w0, 1)
                i0, ws0 = detect_peaks(P0, w0, rel_prom=0.05)
                k0s = map_k_from_omega(ws0[int(np.argmax(P0[i0]))], dr, K_BAND) if i0 else []
                a0 = fit0["a"]
                for eps in (1e-6, 1e-4):
                    y1 = y0 + eps / np.sqrt(2) * (
                        np.random.default_rng(42).standard_normal(N_R)
                        + 1j * np.random.default_rng(43).standard_normal(N_R)
                    )
                    fit1 = complex_mcov_ar(y1, p)
                    P1 = ar_spectrum(fit1["a"], fit1["sigma2"], w0, 1)
                    i1, ws1 = detect_peaks(P1, w0, rel_prom=0.05)
                    k1s = map_k_from_omega(ws1[int(np.argmax(P1[i1]))], dr, K_BAND) if i1 else []
                    cond_rows.append(
                        {
                            "p_tag": p_tag,
                            "p": p,
                            "dr": dr,
                            "eps": eps,
                            "coef_rel_change": float(np.linalg.norm(fit1["a"] - a0) / max(np.linalg.norm(a0), 1e-30)),
                            "peak_k0": k0s[0] if k0s else np.nan,
                            "peak_k1": k1s[0] if k1s else np.nan,
                            "peak_shift": abs((k0s[0] if k0s else np.nan) - (k1s[0] if k1s else np.nan)),
                            "rank": fit1["rank"],
                            "cond": fit1["cond"],
                            "data_model_rank_limit": fit1["data_model_rank_limit"],
                        }
                    )
        pd.DataFrame(cond_rows).to_csv(OUT / "numerical_conditioning_diagnostic.csv", index=False)

        (OUT / "RANK_VS_CONDITIONING_AUDIT.md").write_text(
            f"""# RANK_VS_CONDITIONING_AUDIT

UTC: {NOW}

## DATA_MODEL_RANK_LIMIT

无噪声 q 个复指数是**低秩模型**：单音 rank≪p 不等于 p=133 在真实多模上不稳。
单音 rank deficiency 记为 `DATA_MODEL_RANK_LIMIT`。

## NUMERICAL_CONDITIONING

独立 diagnostic：固定 seed、ε=1e-6/1e-4 复高斯扰动，看系数/峰位敏感度
（`numerical_conditioning_diagnostic.csv`）。**不计入** MMAR 成功率，**不**用来调 p。

## 性能分支（未改）

p=7 `PRINTED_LITERAL_CONTROL`；p=133 `OUR_MOVING_SAMPLE_INTERPRETATION`。
""",
            encoding="utf-8",
        )

        if f1_st >= 0.9 and multi_f1 >= 0.75 and multi_ok >= 0.5:
            decision = "AR_K_RECOVERY_UNIT_VALIDATED"
        elif f1_st >= 0.5:
            decision = "AR_K_RECOVERY_PARTIAL"
        else:
            # check if conditioning (not just rank) blows up
            cdf = pd.DataFrame(cond_rows)
            if bool((cdf["coef_rel_change"] > 10).any() and (cdf["peak_shift"] > 1e-3).any()):
                decision = "AR_K_RECOVERY_NUMERICALLY_UNSTABLE"
            elif f1_st >= 0.3:
                decision = "AR_K_RECOVERY_PARTIAL"
            else:
                decision = "AR_K_RECOVERY_UNIT_FAILED"
        why = f"P1={p1_status}; F1_single={f1_st:.3f}; multi_complete={multi_ok:.2f}; multi_F1={multi_f1:.3f}"

    (OUT / "R3_C2_4B_1F2_DECISION.json").write_text(
        json.dumps(
            {
                "stage": "R3-C2.4B-1F2",
                "decision": decision,
                "why": why,
                "p1_identity_gate": p1_status,
                "complex_backward_conjugation": True,
                "sign_mapping": "k=-(omega+2pi m)/dr from p=1 analytic wrap(-k dr)",
                "eq24": "OUR_EXACT_SOLVER_FOR_EQ24_VALIDATED (not re-run)",
                "not_done": ["4B-2", "Hankel", "D(z)", "FIELD", "E-STD", "MC", "P5"],
                "created_utc": NOW,
            },
            ensure_ascii=False,
            indent=2,
        ),
        encoding="utf-8",
    )
    (OUT / "R3_C2_4B_1F2_REPORT.md").write_text(
        f"""# R3-C2.4B-1F2 complex modified covariance

UTC: {NOW}

## 判定

### `{decision}`

{why}

## 修复

1. **`COMPLEX_BACKWARD_CONJUGATION_REQUIRED`**：Ab=conj(y[...])，bb=−conj(y[n−p])
2. **p=1 硬门** `{p1_status}`：a1=−e^{{−jkΔr}}；ω=wrap(−kΔr)
3. sign mapping 由解析恒等式冻结（无 majority vote）
4. rank：`DATA_MODEL_RANK_LIMIT` ≠ 数值不稳；另做 ε 扰动 conditioning diagnostic

## 指标（power=1）

见 `SINGLE_TONE_RECOVERY_FINAL.csv` / `MULTITONE_RECOVERY_FINAL.csv`。

## 停止

不进 4B-2 / Hankel / D(z) / FIELD / E-STD / MC / P5。
""",
        encoding="utf-8",
    )
    (OUT / "GPT_SYNC.md").write_text(
        f"# R3-C2.4B-1F2\\n\\n**{decision}**\\n\\n{p1_status}\\n{why}\\n",
        encoding="utf-8",
    )
    print("DECISION", decision)
    print("P1", p1_status, "why", why)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
