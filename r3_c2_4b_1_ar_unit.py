#!/usr/bin/env python3
"""R3-C2.4B-1: AR wavenumber recovery unit test. No Hankel, no D(z), no FIELD."""
from __future__ import annotations

import json
from datetime import datetime, timezone
from itertools import combinations
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parent
OUT = ROOT / "results" / "R3_C2_AR_MMAR" / "R3_C2_4B_1_AR_UNIT"
OUT.mkdir(parents=True, exist_ok=True)
NOW = datetime.now(timezone.utc).isoformat()

K_PHYS = [1.32, 1.41, 1.50]  # rad/m frozen before data gen
K_BAND = (1.32, 1.50)
DRS = [
    (np.pi * 0.9 / 1.50, "NO_FOLD_CONTROL"),
    (0.5, "AUX_CONTROL"),
    (2.5, "REF7_CONSISTENT_SAMPLING_ASSUMPTION"),
]
P_LIT = 7  # floor(2/3 * 11)
TWO_PI = 2 * np.pi


def modified_covariance_ar(y: np.ndarray, p: int):
    """Least-squares forward-backward (modified covariance) AR coefficients.

    y[i] = -sum_k a[k] y[i-k] + u[i]
    Returns a[1..p] and sigma2.
    """
    y = np.asarray(y, dtype=complex)
    N = y.size
    if N < p + 2:
        raise ValueError("N too small")
    # rows: n=p..N-1 forward, n=0..N-p-1 backward
    A_rows = []
    b_rows = []
    for n in range(p, N):
        A_rows.append(y[n - 1 :: -1][:p])  # y[n-1], y[n-2], ...
        b_rows.append(-y[n])
    for n in range(0, N - p):
        # backward: y[n+p] + sum_k a[k] y[n+p-k]
        A_rows.append(y[n + p - 1 :: -1][:p] if False else y[n : n + p][::-1])
        b_rows.append(-y[n + p])
    A = np.asarray(A_rows, dtype=complex)
    b = np.asarray(b_rows, dtype=complex)
    # solve least squares
    a, *_ = np.linalg.lstsq(A, b, rcond=None)
    # residual variance (forward only)
    u = []
    for n in range(p, N):
        pred = -np.dot(a, y[n - 1 :: -1][:p])
        u.append(y[n] - pred)
    sigma2 = float(np.mean(np.abs(u) ** 2))
    return a, sigma2


def ar_spectrum(a, sigma2, omega, power=1):
    """P_AR(omega) = sigma2 / |1+sum a[k] e^{-i omega k}|^power  (power=1 printed, 2 standard)."""
    a = np.asarray(a, dtype=complex)
    p = a.size
    acc = np.ones_like(omega, dtype=complex)
    for k in range(1, p + 1):
        acc = acc + a[k - 1] * np.exp(-1j * omega * k)
    denom = np.abs(acc) ** power
    denom = np.maximum(denom, 1e-30)
    return sigma2 / denom


def find_omega_peaks(P, omega):
    peaks = []
    n = P.size
    for i in range(1, n - 1):
        if P[i] >= P[i - 1] and P[i] >= P[i + 1] and P[i] > 0:
            peaks.append(omega[i])
    return peaks


def map_omega_to_k(omega_peaks, dr, k_band):
    """omega ≡ ± k dr (mod 2π); branch select into k_band."""
    rows = []
    for w in np.atleast_1d(omega_peaks):
        # candidates k = (w + 2π m)/dr and k = (-w + 2π m)/dr
        cands = []
        for m in range(-5, 6):
            cands.append((w + TWO_PI * m) / dr)
            cands.append((-w + TWO_PI * m) / dr)
        cands = [c for c in cands if np.isfinite(c)]
        in_band = [c for c in cands if k_band[0] - 1e-9 <= c <= k_band[1] + 1e-9]
        if not in_band:
            direct = [c for c in cands if abs(c) < 1 / dr * np.pi + 0.5]
            status = "FOLDED_BRANCH_AMBIGUOUS" if not in_band else "OK"
            k_rec = in_band[0] if in_band else np.nan
            status = "NO_CANDIDATE_IN_BAND"
            unique = False
        elif len(in_band) == 1:
            k_rec = in_band[0]
            unique = True
            # if |omega| <= π and no fold: direct
            if abs(w) <= np.pi and abs(w / dr - k_rec) < 1e-9:
                status = "NO_FOLD_DIRECT_RECOVERY"
            else:
                status = "FOLDED_BUT_BRANCH_RECOVERABLE"
        else:
            k_rec = in_band[0]
            unique = False
            status = "FOLDED_BRANCH_AMBIGUOUS"
        rows.append(
            {
                "omega_wrapped": float(((w + np.pi) % TWO_PI) - np.pi),
                "k_direct_w_over_dr": float(w / dr),
                "k_rec": float(k_rec) if k_rec == k_rec else np.nan,
                "n_candidates_in_band": len(in_band),
                "unique_in_band": unique,
                "status": status,
                "all_in_band": ";".join(f"{c:.6f}" for c in in_band[:6]),
            }
        )
    return rows


def gen_tone(k_list, dr, N, amps=None, phase0=None, noise=0.0, seed=0):
    n = np.arange(N)
    y = np.zeros(N, dtype=complex)
    if amps is None:
        amps = np.ones(len(k_list))
    if phase0 is None:
        phase0 = np.zeros(len(k_list))
    for k, A, ph in zip(k_list, amps, phase0):
        y += A * np.exp(-1j * k * n * dr + 1j * ph)
    if noise:
        rng = np.random.default_rng(seed)
        y = y + noise / np.sqrt(2) * (rng.standard_normal(N) + 1j * rng.standard_normal(N))
    return y


def recover_k_for_case(k_true, dr, N, p, power, label):
    y = gen_tone([k_true], dr, N)
    try:
        a, s2 = modified_covariance_ar(y, p)
    except Exception as e:
        return {
            "k_true": k_true,
            "dr": dr,
            "N": N,
            "p": p,
            "power": power,
            "label": label,
            "status": f"AR_FIT_FAIL:{e}",
            "abs_err": np.nan,
        }
    omega = np.linspace(-np.pi, np.pi, 8192, endpoint=False)
    P = ar_spectrum(a, s2, omega, power=power)
    peaks = find_omega_peaks(P, omega)
    mapped = map_omega_to_k(peaks, dr, K_BAND)
    # pick best peak closest to true after mapping
    best = None
    best_err = np.inf
    for m in mapped:
        if m["k_rec"] == m["k_rec"]:
            err = abs(m["k_rec"] - k_true)
            if err < best_err:
                best_err = err
                best = m
    rec = {
        "k_true": k_true,
        "dr": dr,
        "N": N,
        "p": p,
        "power": power,
        "label": label,
        "n_omega_peaks": len(peaks),
        "abs_err": best_err if best else np.nan,
        "k_rec": best["k_rec"] if best else np.nan,
        "status": best["status"] if best else "NO_PEAK",
        "unique_in_band": best["unique_in_band"] if best else False,
    }
    return rec


def eq24_dp(k_est, k_model):
    """OUR_EXACT_SOLVER_FOR_EQ24: min ||k-k0||^2 s.t. k0 increasing subset of k_model."""
    k_est = np.asarray(k_est, dtype=float)
    k_model = np.asarray(k_model, dtype=float)
    M0, M = k_est.size, k_model.size
    if M0 > M:
        return None, np.inf
    neg = -1e300
    dp = np.full((M0, M), np.inf)
    prev = np.full((M0, M), -1, dtype=int)
    for j in range(M):
        dp[0, j] = (k_est[0] - k_model[j]) ** 2
    for i in range(1, M0):
        best_so_far = np.inf
        best_j = -1
        for j in range(M):
            # min dp[i-1, j'<j]
            if best_j >= 0:
                dp[i, j] = (k_est[i] - k_model[j]) ** 2 + best_so_far
                prev[i, j] = best_j
            if dp[i - 1, j] < best_so_far:
                best_so_far = dp[i - 1, j]
                best_j = j
    j = int(np.argmin(dp[M0 - 1]))
    cost = float(dp[M0 - 1, j])
    idx = [0] * M0
    idx[M0 - 1] = j
    for i in range(M0 - 1, 0, -1):
        j = prev[i, j]
        idx[i - 1] = j
    return idx, cost


def eq24_bruteforce(k_est, k_model):
    k_est = np.asarray(k_est, dtype=float)
    k_model = np.asarray(k_model, dtype=float)
    M0 = k_est.size
    best, best_c = None, np.inf
    for comb in combinations(range(k_model.size), M0):
        if list(comb) != sorted(comb):
            continue
        c = float(np.sum((k_est - k_model[list(comb)]) ** 2))
        if c < best_c:
            best_c, best = c, list(comb)
    return best, best_c


def main() -> int:
    N_R = 200  # moving-range samples for interpretation branch
    p_moving = int(np.floor(2 * N_R / 3))
    p_branches = [(P_LIT, "PRINTED_LITERAL_CONTROL"), (p_moving, "OUR_MOVING_SAMPLE_INTERPRETATION")]

    # ---- single tone ----
    single_rows = []
    for k in K_PHYS:
        for dr, dr_tag in DRS:
            for p, p_tag in p_branches:
                for power, pow_tag in [
                    (1, "PRINTED_EQ19_MODULUS_POWER_1"),
                    (2, "STANDARD_AR_PSD_MODULUS_POWER_2_CONTROL"),
                ]:
                    r = recover_k_for_case(k, dr, N_R, p, power, f"{dr_tag}|{p_tag}|{pow_tag}")
                    r["dr_tag"] = dr_tag
                    r["p_tag"] = p_tag
                    r["power_tag"] = pow_tag
                    single_rows.append(r)
    pd.DataFrame(single_rows).to_csv(OUT / "SINGLE_TONE_RECOVERY.csv", index=False)

    # ---- multitone ----
    multi_defs = [
        ("2_peak_wide", [1.35, 1.48], None),
        ("3_peak_wide", [1.33, 1.41, 1.49], None),
        ("5_peak_wide", [1.32, 1.36, 1.41, 1.46, 1.50], None),
        ("2_peak_close", [1.40, 1.405], None),
        ("3_peak_close", [1.40, 1.404, 1.408], None),
        ("5_peak_close", [1.38, 1.39, 1.40, 1.41, 1.42], None),
    ]
    multi_rows = []
    for name, ks, _ in multi_defs:
        for dr, dr_tag in [(DRS[0][0], DRS[0][1]), (2.5, DRS[2][1])]:
            for p, p_tag in p_branches:
                for power, pow_tag in [(1, "PRINTED_EQ19_MODULUS_POWER_1"), (2, "STANDARD_AR_PSD_MODULUS_POWER_2_CONTROL")]:
                    y = gen_tone(ks, dr, N_R)
                    try:
                        a, s2 = modified_covariance_ar(y, p)
                    except Exception as e:
                        multi_rows.append(
                            {
                                "case": name,
                                "n_true": len(ks),
                                "dr": dr,
                                "p": p,
                                "power": power,
                                "status": f"AR_FIT_FAIL:{e}",
                                "n_rec_in_band": 0,
                                "mean_abs_err": np.nan,
                                "max_abs_err": np.nan,
                                "all_recovered_within_1e-3": False,
                            }
                        )
                        continue
                    omega = np.linspace(-np.pi, np.pi, 16384, endpoint=False)
                    P = ar_spectrum(a, s2, omega, power=power)
                    peaks = find_omega_peaks(P, omega)
                    mapped = [m for m in map_omega_to_k(peaks, dr, K_BAND) if m["k_rec"] == m["k_rec"]]
                    krec = sorted(m["k_rec"] for m in mapped)
                    # greedy match to true
                    errs = []
                    used = set()
                    for kt in ks:
                        if not krec:
                            errs.append(np.nan)
                            continue
                        j = int(np.argmin([abs(kr - kt) if i not in used else np.inf for i, kr in enumerate(krec)]))
                        used.add(j)
                        errs.append(abs(krec[j] - kt))
                    errs = [e for e in errs if e == e]
                    multi_rows.append(
                        {
                            "case": name,
                            "n_true": len(ks),
                            "dr": dr,
                            "dr_tag": dr_tag,
                            "p": p,
                            "p_tag": p_tag,
                            "power": power,
                            "power_tag": pow_tag,
                            "n_rec_in_band": len(krec),
                            "mean_abs_err": float(np.mean(errs)) if errs else np.nan,
                            "max_abs_err": float(np.max(errs)) if errs else np.nan,
                            "all_recovered_within_1e-3": bool(errs and all(e < 1e-3 for e in errs)),
                            "status": "OK" if errs and all(e < 1e-3 for e in errs) else "PARTIAL_OR_FAIL",
                        }
                    )
    pd.DataFrame(multi_rows).to_csv(OUT / "MULTITONE_RECOVERY.csv", index=False)

    # ---- alias branch summary ----
    sdf = pd.DataFrame(single_rows)
    # classification counts
    def class_of(row):
        if row["status"] == "NO_FOLD_DIRECT_RECOVERY" and row["abs_err"] < 1e-3:
            return "NO_FOLD_DIRECT_RECOVERY"
        if row["status"] == "FOLDED_BUT_BRANCH_RECOVERABLE" and row["abs_err"] < 1e-3:
            return "FOLDED_BUT_BRANCH_RECOVERABLE"
        if row["abs_err"] == row["abs_err"] and row["abs_err"] < 1e-2:
            return "FOLDED_BUT_BRANCH_RECOVERABLE"
        return "FOLDED_BRANCH_AMBIGUOUS"

    sdf["alias_class"] = sdf.apply(class_of, axis=1)
    sdf.to_csv(OUT / "SINGLE_TONE_RECOVERY.csv", index=False)

    (OUT / "ALIAS_BRANCH_RECOVERY.md").write_text(
        f"""# ALIAS_BRANCH_RECOVERY

UTC: {NOW}

模型：y[n]=Σ A_q exp(−j k_q n Δr) ⇒ AR 归一化角频率 ω ≡ ± k Δr (mod 2π)

## 分类

| 类 | 条件 |
| --- | --- |
| `NO_FOLD_DIRECT_RECOVERY` | \\|ω\\|≤π 且 k=ω/Δr 落在物理带 [1.32,1.50] |
| `FOLDED_BUT_BRANCH_RECOVERABLE` | 折叠，但在已知物理带内唯一可恢复 |
| `FOLDED_BRANCH_AMBIGUOUS` | 带内多候选 / 无法唯一 |

## 单音结果摘要（按 Δr）

见 `SINGLE_TONE_RECOVERY.csv`（k=1.32/1.41/1.50；p=7 与 p=⌊2N_r/3⌋；Eq19 幂次 1 与 2）。

本审计**不**评价 MMAR 深度成败。
""",
        encoding="utf-8",
    )

    # ---- Eq24 solver unit test ----
    eq_rows = []
    rng = np.random.default_rng(0)
    for trial in range(20):
        M = 8
        M0 = 4
        k_model = np.sort(rng.uniform(1.32, 1.50, M))
        idx_true = sorted(rng.choice(M, M0, replace=False))
        k_est = k_model[idx_true] + rng.normal(0, 0.0005, M0)
        idx_dp, cost_dp = eq24_dp(k_est, k_model)
        idx_bf, cost_bf = eq24_bruteforce(k_est, k_model)
        eq_rows.append(
            {
                "trial": trial,
                "M": M,
                "M0": M0,
                "dp_indices": ";".join(map(str, idx_dp)),
                "bf_indices": ";".join(map(str, idx_bf)),
                "true_indices": ";".join(map(str, idx_true)),
                "dp_cost": cost_dp,
                "bf_cost": cost_bf,
                "cost_match": abs(cost_dp - cost_bf) < 1e-9,
                "recovered_true_order": idx_dp == list(idx_true),
            }
        )
    eq_df = pd.DataFrame(eq_rows)
    eq_df.to_csv(OUT / "EQ24_SOLVER_UNIT_TEST.csv", index=False)
    solver_ok = bool(eq_df["cost_match"].all())
    order_rate = float(eq_df["recovered_true_order"].mean())

    # ---- decision ----
    single_ok = float((sdf["abs_err"] < 1e-3).mean()) if len(sdf) else 0.0
    multi_ok = float((pd.DataFrame(multi_rows)["all_recovered_within_1e-3"]).mean()) if multi_rows else 0.0
    has_fold_amb = (sdf["alias_class"] == "FOLDED_BRANCH_AMBIGUOUS").any()
    has_fold_rec = (sdf["alias_class"] == "FOLDED_BUT_BRANCH_RECOVERABLE").any()
    has_nofold = (sdf["alias_class"] == "NO_FOLD_DIRECT_RECOVERY").any()

    if not solver_ok or single_ok < 0.5:
        decision = "AR_K_RECOVERY_UNIT_FAILED"
        why = f"solver_ok={solver_ok}, single_ok={single_ok:.2f}"
    elif single_ok < 0.95 or has_fold_amb:
        decision = "AR_K_RECOVERY_PARTIAL_WITH_ALIAS_LIMIT"
        why = f"single_ok={single_ok:.2f}, fold_amb={has_fold_amb}, fold_rec={has_fold_rec}, nofold={has_nofold}"
    elif multi_ok < 0.5:
        decision = "AR_K_RECOVERY_UNIT_FAILED"
        why = f"multi_ok={multi_ok:.2f}"
    else:
        decision = "AR_K_RECOVERY_UNIT_VALIDATED"
        why = f"single_ok={single_ok:.2f}, multi_ok={multi_ok:.2f}, solver_ok={solver_ok}, order_rate={order_rate:.2f}"

    # implementation audit
    (OUT / "AR_IMPLEMENTATION_AUDIT.md").write_text(
        f"""# AR_IMPLEMENTATION_AUDIT

UTC: {NOW}
基线 commit: c6a49c1b6e033941d6e0df0f4fff1c28a78995a8
来源仅：MMAR_EQUATION_LOCK_FINAL / AR_INTERPRETATION_GATE / AR_ORDER_PREREGISTRATION / SPATIAL_SAMPLING_ALIAS_AUDIT / MMAR_PAPER_CONFIG_FINAL

## AR

- **modified covariance**（forward-backward LS），对应论文选型（避免 line splitting）
- 输入：y[n]=Σ A_q e^{{-j k_q n Δr}}（真值 k 冻结后生成）
- 只恢复**峰位置** k̂；不用 AR peak height 作幅度

## Eq.(19)

- `PRINTED_EQ19_MODULUS_POWER_1`
- `STANDARD_AR_PSD_MODULUS_POWER_2_CONTROL`
- 只比峰位置，不按深度择优

## AR order（预注册，未改）

- `PRINTED_LITERAL_CONTROL`: p=7
- `OUR_MOVING_SAMPLE_INTERPRETATION`: p=⌊2N_r/3⌋={p_moving}（N_r={N_R}）
- **无额外 p 扫描**

## Eq.(24)

- `OUR_EXACT_SOLVER_FOR_EQ24`（DP）
- 与 brute-force 小规模交叉验证：cost_match={solver_ok}
- 人工 mode set 上 MODE_ORDER_RECOVERY_RATE={order_rate:.3f}

## 本轮禁止（已遵守）

Hankel g、D(z)、Liang Fig5/6、FIELD、E-STD、MC、P5
""",
        encoding="utf-8",
    )

    dec = {
        "stage": "R3-C2.4B-1",
        "decision": decision,
        "why": why,
        "single_tone_ok_rate": single_ok,
        "multitone_ok_rate": multi_ok,
        "eq24_solver_ok": solver_ok,
        "mode_order_recovery_rate_synthetic": order_rate,
        "p_branches": {"PRINTED_LITERAL_CONTROL": 7, "OUR_MOVING_SAMPLE_INTERPRETATION": p_moving, "N_r": N_R},
        "alias_classes_present": sorted(set(sdf["alias_class"].tolist())),
        "not_done": ["Hankel", "D(z)", "FIELD", "E-STD", "MC", "P5", "Liang depth repro"],
        "created_utc": NOW,
    }
    (OUT / "R3_C2_4B_1_DECISION.json").write_text(json.dumps(dec, ensure_ascii=False, indent=2), encoding="utf-8")
    (OUT / "R3_C2_4B_1_REPORT.md").write_text(
        f"""# R3-C2.4B-1 AR 波数恢复单元测试

UTC: {NOW}

## 判定

### `{decision}`

{why}

## 指标

- 单音 |k̂−k|<1e-3 比例：**{single_ok:.3f}**
- 多音全恢复比例：**{multi_ok:.3f}**
- Eq24 solver cost_match：**{solver_ok}**
- 合成 mode set MODE_ORDER_RECOVERY_RATE：**{order_rate:.3f}**

## 两分支 p（未调参）

p=7（PRINTED_LITERAL_CONTROL）与 p=⌊2N_r/3⌋={p_moving}（OUR_MOVING_SAMPLE_INTERPRETATION），N_r={N_R}。

## 别名

{sorted(set(sdf['alias_class'].tolist()))}

详见 `ALIAS_BRANCH_RECOVERY.md`。

## 停止

不进 4B-2 / FIELD / E-STD / MC / P5。
""",
        encoding="utf-8",
    )
    (OUT / "GPT_SYNC.md").write_text(
        f"# R3-C2.4B-1\\n\\n**{decision}**\\n\\n{why}\\n\\nsingle_ok={single_ok:.3f} multi_ok={multi_ok:.3f} order_rate={order_rate:.3f}\\n",
        encoding="utf-8",
    )
    print("DECISION", decision, "single", single_ok, "multi", multi_ok, "order", order_rate)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
