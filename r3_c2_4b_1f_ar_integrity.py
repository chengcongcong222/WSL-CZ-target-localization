#!/usr/bin/env python3
"""R3-C2.4B-1F: AR implementation integrity fix. No Hankel/D(z)/FIELD/E-STD."""
from __future__ import annotations

import json
from datetime import datetime, timezone
from itertools import combinations
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parent
OUT = ROOT / "results" / "R3_C2_AR_MMAR" / "R3_C2_4B_1F_AR_INTEGRITY"
OUT.mkdir(parents=True, exist_ok=True)
NOW = datetime.now(timezone.utc).isoformat()

K_BAND = (1.32, 1.50)
K_PHYS = [1.32, 1.41, 1.50]
N_R = 200
P_LIT = 7
P_MOVING = int(np.floor(2 * N_R / 3))  # 133
TWO_PI = 2 * np.pi


def modified_covariance_ar(y: np.ndarray, p: int):
    """True modified-covariance (forward + backward LS).

    Forward  (n=p..N-1):  y[n]   = -sum_k a[k] y[n-k]
    Backward (n=p..N-1):  y[n-p] = -sum_k a[k] y[n-p+k]
    """
    y = np.asarray(y, dtype=complex)
    N = y.size
    Af, bf, Ab, bb = [], [], [], []
    for n in range(p, N):
        Af.append(np.array([y[n - k] for k in range(1, p + 1)]))
        bf.append(-y[n])
        Ab.append(np.array([y[n - p + k] for k in range(1, p + 1)]))
        bb.append(-y[n - p])
    A = np.vstack([np.asarray(Af), np.asarray(Ab)])
    b = np.concatenate([np.asarray(bf), np.asarray(bb)])
    # rank / condition
    rank = int(np.linalg.matrix_rank(A))
    cond = float(np.linalg.cond(A))
    a, *_ = np.linalg.lstsq(A, b, rcond=None)
    # residual on forward
    u = []
    for n in range(p, N):
        u.append(y[n] + np.dot(a, np.array([y[n - k] for k in range(1, p + 1)])))
    sigma2 = float(np.mean(np.abs(u) ** 2))
    n_f, n_b = N - p, N - p
    dup = n_f == n_b and n_f > 0 and np.allclose(np.asarray(Af), np.asarray(Ab)[:, ::-1] if False else np.asarray(Ab))
    # check if Ab rows equal Af rows (the bug)
    same_as_fwd = False
    if Af and Ab:
        same_as_fwd = bool(np.allclose(np.asarray(Af), np.asarray(Ab)))
    return {
        "a": a,
        "sigma2": sigma2,
        "n_forward": n_f,
        "n_backward": n_b,
        "n_params": p,
        "n_equations": int(A.shape[0]),
        "rank": rank,
        "cond": cond,
        "backward_dup_forward": same_as_fwd,
        "sufficient_independent": rank >= p and n_f + n_b > p,
    }


def ar_spectrum(a, sigma2, omega, power=1):
    a = np.asarray(a, dtype=complex)
    p = a.size
    acc = np.ones_like(omega, dtype=complex)
    for k in range(1, p + 1):
        acc = acc + a[k - 1] * np.exp(-1j * omega * k)
    denom = np.maximum(np.abs(acc) ** power, 1e-30)
    return sigma2 / denom


def detect_peaks(P, omega, rel_prom=0.05):
    """Truth-independent peak detector (fixed threshold)."""
    n = P.size
    prom = rel_prom * float(P.max())
    cand = []
    for i in range(1, n - 1):
        if P[i] >= P[i - 1] and P[i] >= P[i + 1] and P[i] > 0:
            cand.append(i)
    peaks = []
    for i in cand:
        # simple prominence: height - max(min-left, min-right until higher)
        l = i
        while l > 0 and P[l - 1] <= P[i]:
            l -= 1
        r = i
        while r < n - 1 and P[r + 1] <= P[i]:
            r += 1
        pr = P[i] - max(P[l : i + 1].min(), P[i : r + 1].min())
        if pr >= prom:
            peaks.append(i)
    # distance: merge within 3 bins keep higher
    peaks.sort(key=lambda i: -P[i])
    keep = []
    for i in peaks:
        if all(abs(i - j) > 3 for j, _ in keep):
            keep.append((i, P[i]))
    keep.sort()
    idx = [i for i, _ in keep]
    return idx, np.array([omega[i] for i in idx])


def freeze_sign_mapping(p=1, dr=1.0, k=1.41, N=100):
    """SIGN_MAPPING_DIAGNOSTIC_ONLY: y=e^{-j k n dr} → which omega peak vs k."""
    n = np.arange(N)
    y = np.exp(-1j * k * n * dr)
    fit = modified_covariance_ar(y, p)
    omega = np.linspace(-np.pi, np.pi, 4096, endpoint=False)
    P = ar_spectrum(fit["a"], fit["sigma2"], omega, power=1)
    idx, ws = detect_peaks(P, omega, rel_prom=0.05)
    if ws.size == 0:
        return {"ok": False}
    w = float(ws[np.argmax(P[idx])])
    k_plus = w / dr
    k_minus = -w / dr
    err_plus = abs(k_plus - k)
    err_minus = abs(k_minus - k)
    mapping = "k=+omega/dr" if err_plus < err_minus else "k=-omega/dr"
    return {
        "ok": True,
        "k_true": k,
        "omega_peak": w,
        "k_if_plus": k_plus,
        "k_if_minus": k_minus,
        "err_plus": err_plus,
        "err_minus": err_minus,
        "frozen_mapping": mapping,
    }


def omega_to_k(w, dr, sign, k_band):
    """Apply frozen sign: k = sign * omega / dr, then fold branches into k_band."""
    k0 = sign * w / dr
    cands = [k0 + sign * 0 + m * TWO_PI / dr for m in range(-8, 9)]
    # also exact: k such that sign*k*dr ≡ w (mod 2π) ⇒ k = sign*(w+2πm)/dr
    cands = [sign * (w + TWO_PI * m) / dr for m in range(-8, 9)]
    in_band = sorted({round(c, 9) for c in cands if k_band[0] - 1e-9 <= c <= k_band[1] + 1e-9})
    return in_band


def classify_alias(k_true, dr, in_band, mapping_sign):
    w_true = mapping_sign * k_true * dr
    w_true_wrapped = (w_true + np.pi) % TWO_PI - np.pi
    folds = abs(k_true * dr) > np.pi
    if not folds:
        if in_band and abs(in_band[0] - k_true) < 1e-6:
            return "NO_FOLD_DIRECT_RECOVERY"
        return "NO_FOLD_DIRECT_RECOVERY" if not folds else "FOLDED_BUT_BRANCH_RECOVERABLE"
    if len(in_band) == 1:
        return "FOLDED_BUT_BRANCH_RECOVERABLE"
    if len(in_band) > 1:
        return "FOLDED_BRANCH_AMBIGUOUS"
    return "FOLDED_BRANCH_AMBIGUOUS"


def gen_tone(ks, dr, N, amps=None):
    n = np.arange(N)
    y = np.zeros(N, complex)
    amps = amps or [1.0] * len(ks)
    for k, A in zip(ks, amps):
        y = y + A * np.exp(-1j * k * n * dr)
    return y


def peak_metrics(omega_peaks, dr, sign, k_band, k_true_list, tol=1e-3):
    """Truth-independent peaks; score vs truth after freeze."""
    rec = []
    for w in omega_peaks:
        in_b = omega_to_k(float(w), dr, sign, k_band)
        rec.append((float(w), in_b))
    all_k = sorted({k for _, ib in rec for k in ib})
    # selected: highest not available here; caller passes order
    hits = 0
    false_peaks = 0
    for kt in k_true_list:
        if any(abs(k - kt) < tol for k in all_k):
            hits += 1
    # false: in_band candidates not near any true
    for k in all_k:
        if all(abs(k - kt) > tol for kt in k_true_list):
            false_peaks += 1
    recall = hits / max(len(k_true_list), 1)
    precision = hits / max(len(all_k), 1) if all_k else 0.0
    f1 = 2 * precision * recall / max(precision + recall, 1e-30)
    return {
        "n_detected_peaks_in_band": len(all_k),
        "true_peak_recalled": int(hits),
        "false_peak_count": int(false_peaks),
        "recall": float(recall),
        "precision": float(precision),
        "f1": float(f1),
        "peak_count_err": int(len(all_k) - len(k_true_list)),
        "recovered_ks": ";".join(f"{k:.6f}" for k in all_k),
        "abs_err_list": ";".join(
            f"{min(abs(k-kt) for k in all_k):.3e}" if all_k else "nan" for kt in k_true_list
        ),
        "max_abs_err": float(
            max((min(abs(k - kt) for k in all_k) for kt in k_true_list), default=np.nan)
        )
        if all_k
        else np.nan,
    }


def eq24_dp(k_est, k_model):
    k_est = np.asarray(k_est, float)
    k_model = np.asarray(k_model, float)
    M0, M = k_est.size, k_model.size
    if M0 > M:
        return None, np.inf
    dp = np.full((M0, M), np.inf)
    prev = np.full((M0, M), -1, int)
    for j in range(M):
        dp[0, j] = (k_est[0] - k_model[j]) ** 2
    for i in range(1, M0):
        best = np.inf
        bj = -1
        for j in range(M):
            if bj >= 0:
                dp[i, j] = (k_est[i] - k_model[j]) ** 2 + best
                prev[i, j] = bj
            if dp[i - 1, j] < best:
                best = dp[i - 1, j]
                bj = j
    j = int(np.argmin(dp[M0 - 1]))
    idx = [j]
    for i in range(M0 - 1, 0, -1):
        j = int(prev[i, j])
        idx.append(j)
    idx = idx[::-1]
    return idx, float(dp[M0 - 1, idx[-1]])


def eq24_bf(k_est, k_model):
    k_est = np.asarray(k_est, float)
    k_model = np.asarray(k_model, float)
    M0 = k_est.size
    best, bc = None, np.inf
    for comb in combinations(range(k_model.size), M0):
        if list(comb) != sorted(comb):
            continue
        c = float(np.sum((k_est - k_model[list(comb)]) ** 2))
        if c < bc:
            bc, best = c, list(comb)
    return best, bc


def main() -> int:
    # ---- sign mapping diagnostic (p=1 only) ----
    sign_rows = []
    for k in K_PHYS:
        for dr in [1.885, 0.5, 2.5]:
            r = freeze_sign_mapping(p=1, dr=dr, k=k, N=80)
            r["dr"] = dr
            sign_rows.append(r)
    # majority vote freeze
    maps = [r["frozen_mapping"] for r in sign_rows if r.get("ok")]
    frozen_map = max(set(maps), key=maps.count) if maps else "k=-omega/dr"
    sign = 1.0 if frozen_map == "k=+omega/dr" else -1.0
    pd.DataFrame(sign_rows).to_csv(OUT / "SIGN_MAPPING_DIAGNOSTIC.csv", index=False)

    # ---- matrix audit for real MC ----
    audits = []
    for p in (1, 7, 133):
        y = gen_tone([1.41], 1.885, N_R)
        fit = modified_covariance_ar(y, p)
        audits.append(
            {
                "p": p,
                "N_r": N_R,
                "n_forward": fit["n_forward"],
                "n_backward": fit["n_backward"],
                "n_equations": fit["n_equations"],
                "expected_2_Nr_minus_p": 2 * (N_R - p),
                "n_params": fit["n_params"],
                "rank": fit["rank"],
                "cond": fit["cond"],
                "backward_dup_forward": fit["backward_dup_forward"],
                "sufficient_independent": fit["sufficient_independent"],
            }
        )
    pd.DataFrame(audits).to_csv(OUT / "modified_covariance_matrix_audit.csv", index=False)
    audit_md = [
        "# MODIFIED_COVARIANCE_MATRIX_AUDIT",
        "",
        f"UTC: {NOW}",
        "",
        "Forward (n=p..N-1): A=[y[n-1]…y[n-p]] → −y[n]",
        "Backward (n=p..N-1): A=[y[n-p+1]…y[n]] → −y[n-p]  （与 forward **独立**）",
        "",
        f"2(N_r−p) at N_r=200,p=133: **{2*(200-133)}** 方程 / **133** 参数 → 仍偏紧（rank 见 csv）。",
        "",
        pd.DataFrame(audits).to_markdown(index=False) if False else pd.DataFrame(audits).to_string(index=False),
        "",
        "backward_dup_forward 必须为 False（修复后）。",
    ]
    (OUT / "MODIFIED_COVARIANCE_MATRIX_AUDIT.md").write_text("\n".join(audit_md), encoding="utf-8")

    # ---- single / multi with truth-independent peaks ----
    p_branches = [
        (7, "PRINTED_LITERAL_CONTROL"),
        (P_MOVING, "OUR_MOVING_SAMPLE_INTERPRETATION"),
    ]
    powers = [(1, "PRINTED_EQ19_MODULUS_POWER_1"), (2, "STANDARD_AR_PSD_MODULUS_POWER_2_CONTROL")]

    single_rows = []
    for k in K_PHYS:
        for dr, dr_tag in [(1.885, "NO_FOLD_CONTROL"), (0.5, "AUX_CONTROL"), (2.5, "REF7_2p5")]:
            y = gen_tone([k], dr, N_R)
            for p, p_tag in p_branches:
                fit = modified_covariance_ar(y, p)
                omega = np.linspace(-np.pi, np.pi, 16384, endpoint=False)
                peak_pos_1 = peak_pos_2 = None
                met1 = met2 = {}
                for power, pow_tag in powers:
                    P = ar_spectrum(fit["a"], fit["sigma2"], omega, power=power)
                    idx, ws = detect_peaks(P, omega, rel_prom=0.05)
                    met = peak_metrics(ws, dr, sign, K_BAND, [k])
                    # selected peak = strongest mapped peak (truth-independent)
                    sel_k = np.nan
                    if idx:
                        best = int(np.argmax(P[idx]))
                        in_b = omega_to_k(float(ws[best]), dr, sign, K_BAND)
                        sel_k = in_b[0] if in_b else np.nan
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
                            "selected_abs_err": abs(sel_k - k) if sel_k == sel_k else np.nan,
                            "alias_class": classify_alias(k, dr, omega_to_k(float(ws[0]), dr, sign, K_BAND) if len(ws) else [], sign)
                            if True
                            else "",
                            "backward_dup_forward": fit["backward_dup_forward"],
                            "rank": fit["rank"],
                            "cond": fit["cond"],
                        }
                    )
                    # alias class from geometry only
                    folds = abs(k * dr) > np.pi
                    in_b = omega_to_k(sign * k * dr, dr, sign, K_BAND)
                    if not folds:
                        met["alias_class"] = "NO_FOLD_DIRECT_RECOVERY"
                    elif len(in_b) == 1:
                        met["alias_class"] = "FOLDED_BUT_BRANCH_RECOVERABLE"
                    else:
                        met["alias_class"] = "FOLDED_BRANCH_AMBIGUOUS"
                    single_rows.append(met)
                    if power == 1:
                        peak_pos_1 = list(ws)
                    else:
                        peak_pos_2 = list(ws)
                # peak position consistency power1 vs 2
                if peak_pos_1 is not None and peak_pos_2 is not None:
                    same = (
                        len(peak_pos_1) == len(peak_pos_2)
                        and all(abs(a - b) < 1e-6 for a, b in zip(sorted(peak_pos_1), sorted(peak_pos_2)))
                    )
                    single_rows[-1]["eq19_peakpos_same_as_power2"] = same
                    single_rows[-2]["eq19_peakpos_same_as_power2"] = same

    sdf = pd.DataFrame(single_rows)
    sdf.to_csv(OUT / "SINGLE_TONE_RECOVERY_CORRECTED.csv", index=False)

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
        for dr, dr_tag in [(1.885, "NO_FOLD_CONTROL"), (2.5, "REF7_2p5")]:
            y = gen_tone(ks, dr, N_R)
            for p, p_tag in p_branches:
                fit = modified_covariance_ar(y, p)
                for power, pow_tag in powers:
                    omega = np.linspace(-np.pi, np.pi, 16384, endpoint=False)
                    P = ar_spectrum(fit["a"], fit["sigma2"], omega, power=power)
                    idx, ws = detect_peaks(P, omega, rel_prom=0.05)
                    met = peak_metrics(ws, dr, sign, K_BAND, ks)
                    met.update(
                        {
                            "case": name,
                            "n_true": len(ks),
                            "dr": dr,
                            "dr_tag": dr_tag,
                            "p": p,
                            "p_tag": p_tag,
                            "power": power,
                            "power_tag": pow_tag,
                            "backward_dup_forward": fit["backward_dup_forward"],
                            "rank": fit["rank"],
                            "complete_success": bool(met["recall"] == 1.0 and met["false_peak_count"] == 0),
                        }
                    )
                    multi_rows.append(met)
    mdf = pd.DataFrame(multi_rows)
    mdf.to_csv(OUT / "MULTITONE_RECOVERY_CORRECTED.csv", index=False)

    (OUT / "PEAK_SELECTOR_CONTRACT.md").write_text(
        f"""# PEAK_SELECTOR_CONTRACT

UTC: {NOW}

## 禁止

`choose peak closest to k_true` — **oracle peak selection 已删除**。

## Truth-independent 流程

1. 用固定 prominence=5% global max 检测 **全部** AR 谱峰（不知 k_true）
2. **冻结** 峰集合
3. 用冻结映射 k=sign·ω/Δr + 物理带分支 → k̂ 集合
4. `selected_peak` = **谱峰高度最大**的映射峰（与真值无关）
5. 真值**只做事后评分**：recall / precision / F1 / abs err / false_peak_count / peak_count_err

## p=1 SIGN_MAPPING_DIAGNOSTIC_ONLY

只用于冻结 ω↔k 符号；**不进**性能统计、**不是**新性能分支。

性能分支仍仅：

- `PRINTED_LITERAL_CONTROL` p=7
- `OUR_MOVING_SAMPLE_INTERPRETATION` p={P_MOVING}

## Eq.(19)

幂次 1 与 2 同时跑；报告峰位置是否一致（单调变换，不重复计成功样本）。
""",
        encoding="utf-8",
    )

    (OUT / "ALIAS_CLASSIFICATION_CORRECTED.md").write_text(
        f"""# ALIAS_CLASSIFICATION_CORRECTED

UTC: {NOW}

符号映射（p=1 diagnostic 冻结）：**{frozen_map}**（sign={sign:+.0f}）

分类仅由几何+物理带决定（**不**看恢复误差）：

| 类 | 条件 |
| --- | --- |
| `NO_FOLD_DIRECT_RECOVERY` | kΔr < π |
| `FOLDED_BUT_BRANCH_RECOVERABLE` | kΔr > π 且物理带 [1.32,1.50] 内**唯一**分支 |
| `FOLDED_BRANCH_AMBIGUOUS` | 物理带内多分支 |

修复前把 no-fold 标成 folded 的错误已废除。

sign_rows / single 表见同目录 csv。
""",
        encoding="utf-8",
    )

    # Eq24 recheck (solver unchanged)
    rng = np.random.default_rng(1)
    eq_rows = []
    for trial in range(15):
        M, M0 = 8, 4
        k_model = np.sort(rng.uniform(1.32, 1.50, M))
        k_est = np.sort(rng.uniform(1.32, 1.50, M0))
        idx_dp, c_dp = eq24_dp(k_est, k_model)
        idx_bf, c_bf = eq24_bf(k_est, k_model)
        eq_rows.append(
            {
                "trial": trial,
                "cost_match": abs(c_dp - c_bf) < 1e-9,
                "dp_cost": c_dp,
                "bf_cost": c_bf,
            }
        )
    eq_df = pd.DataFrame(eq_rows)
    eq_df.to_csv(OUT / "EQ24_SOLVER_UNIT_TEST.csv", index=False)
    solver_ok = bool(eq_df["cost_match"].all())

    # metrics (power=1 primary; do not double count power=2)
    prim = sdf[(sdf["power"] == 1)]
    recall_st = float(prim["true_peak_recalled"].mean()) if len(prim) else 0
    prec_st = float(prim["precision"].mean()) if len(prim) else 0
    f1_st = float(prim["f1"].mean()) if len(prim) else 0
    fp_st = float(prim["false_peak_count"].mean()) if len(prim) else 0
    sel_err = float(prim["selected_abs_err"].mean()) if len(prim) else np.nan
    multi_ok = float(mdf[mdf["power"] == 1]["complete_success"].mean()) if len(mdf) else 0
    multi_f1 = float(mdf[mdf["power"] == 1]["f1"].mean()) if len(mdf) else 0
    dup_fwd = bool(sdf["backward_dup_forward"].any())
    p133 = prim[prim["p"] == P_MOVING]
    p133_unstable = bool((p133["rank"] < p133["p"]).any() or (p133["cond"] > 1e12).any()) if len(p133) else True

    if dup_fwd or not solver_ok:
        decision = "AR_K_RECOVERY_UNIT_FAILED"
        why = f"dup_fwd={dup_fwd} solver_ok={solver_ok}"
    elif p133_unstable and f1_st < 0.8:
        decision = "AR_K_RECOVERY_NUMERICALLY_UNSTABLE"
        why = f"p=133 rank/cond unstable; single F1={f1_st:.3f}"
    elif f1_st >= 0.9 and multi_f1 >= 0.8 and multi_ok >= 0.5:
        decision = "AR_K_RECOVERY_UNIT_VALIDATED"
        why = f"F1_single={f1_st:.3f} F1_multi={multi_f1:.3f} multi_complete={multi_ok:.2f}"
    elif f1_st >= 0.5:
        decision = "AR_K_RECOVERY_PARTIAL"
        why = f"F1_single={f1_st:.3f} F1_multi={multi_f1:.3f} multi_complete={multi_ok:.2f} p133_unstable={p133_unstable}"
    else:
        decision = "AR_K_RECOVERY_UNIT_FAILED"
        why = f"F1_single={f1_st:.3f}"

    (OUT / "R3_C2_4B_1F_DECISION.json").write_text(
        json.dumps(
            {
                "stage": "R3-C2.4B-1F",
                "decision": decision,
                "why": why,
                "frozen_sign_mapping": frozen_map,
                "modified_cov_backward_dup_forward": dup_fwd,
                "eq24_solver_ok": solver_ok,
                "single_f1_power1": f1_st,
                "single_precision": prec_st,
                "single_recall": recall_st,
                "single_fp_mean": fp_st,
                "multi_complete_rate_power1": multi_ok,
                "multi_f1_power1": multi_f1,
                "p133_unstable": p133_unstable,
                "p_branches": {"PRINTED_LITERAL_CONTROL": 7, "OUR_MOVING_SAMPLE_INTERPRETATION": P_MOVING},
                "not_done": ["Hankel", "D(z)", "FIELD", "E-STD", "MC", "P5", "4B-2"],
                "created_utc": NOW,
            },
            ensure_ascii=False,
            indent=2,
        ),
        encoding="utf-8",
    )
    (OUT / "R3_C2_4B_1F_REPORT.md").write_text(
        f"""# R3-C2.4B-1F AR implementation integrity

UTC: {NOW}

## 判定

### `{decision}`

{why}

## 修复

1. **真 modified covariance**：backward 为 y[n−p]=−Σa[k]y[n−p+k]（与 forward 独立）；`backward_dup_forward={dup_fwd}`
2. **oracle 选峰删除**：峰集先冻结；selected=最强峰；真值只评分
3. **符号冻结**（p=1 diagnostic）：**{frozen_map}**
4. alias 类按 **kΔr 与物理带** 分类，不按误差改标签

## 指标（power=1，不与 power=2 重复计数）

| 项 | 值 |
| --- | --- |
| single F1 / recall / precision | {f1_st:.3f} / {recall_st:.3f} / {prec_st:.3f} |
| single mean FP | {fp_st:.2f} |
| multi complete / F1 | {multi_ok:.2f} / {multi_f1:.3f} |
| p=133 unstable | {p133_unstable} |
| Eq24 solver | {solver_ok} |

N_r=200,p=133：2(N_r−p)={2*(200-133)} 约束 / 133 参数。

## 停止

不进 4B-2 / Hankel / D(z) / FIELD / E-STD / MC / P5。
""",
        encoding="utf-8",
    )
    (OUT / "GPT_SYNC.md").write_text(
        f"# R3-C2.4B-1F\\n\\n**{decision}**\\n\\n{why}\\n\\nsign={frozen_map}; F1_st={f1_st:.3f}; multi_complete={multi_ok:.2f}\\n",
        encoding="utf-8",
    )
    print("DECISION", decision)
    print("sign", frozen_map, "F1", f1_st, "multi_ok", multi_ok, "p133_unstable", p133_unstable, "dup_fwd", dup_fwd)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
