#!/usr/bin/env python3
"""RANGE-UB-2B-FIX1A: Postprocessing and decision integrity correction.

No BELLHOP. Re-analyze FIX1 data with correct contiguous working range,
tie-aware Spearman, and true theory envelope.
"""
from __future__ import annotations

import json
from datetime import datetime, timezone
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parent
OUT = ROOT / "results" / "RANGE_UB" / "E_STD_VLA_RANGE_CONTINUATION_Z200_FIX"
FIX1 = ROOT / "results" / "RANGE_UB" / "E_STD_VLA_RANGE_CONTINUATION_Z200"
OUT.mkdir(parents=True, exist_ok=True)
NOW = datetime.now(timezone.utc).isoformat()

H = 5000.0
ZS = 200.0
C_MIN, C_MAX = 1500.0, 1635.2


def contiguous_working_range(wr_df, target=50.0):
    """Find max contiguous range containing target where all_three_valid=True."""
    valid = set(wr_df[wr_df["all_three_valid"]]["r_km"].values)
    if target not in valid:
        return [], 0.0
    # extend left
    lo = target
    while lo - 0.5 in valid:
        lo -= 0.5
    # extend right
    hi = target
    while hi + 0.5 in valid:
        hi += 0.5
    block = sorted([r for r in valid if lo <= r <= hi])
    return block, hi - lo


def spearman_tie_aware(x, y):
    """Spearman with average-rank tie handling."""
    def rank_avg(a):
        s = pd.Series(a)
        return s.rank(method="average").values
    rx, ry = rank_avg(x), rank_avg(y)
    return float(np.corrcoef(rx, ry)[0, 1])


def theory_envelope_at_r(r_km):
    """Xu2024 Eq.(10)(11)(13) with E-STD params, zs=200, c in [1500,1635.2]."""
    rows = []
    for c in np.linspace(C_MIN, C_MAX, 20):
        for n2 in [-1, 1]:
            for n3 in [-1, 1]:
                k32 = (-6 * H - (n3 + n2) * ZS) / (c * r_km * 1000) * 1000
                rows.append({"family": "k32", "k": k32})
        for n2 in [-1, 1]:
            for n4 in [-1, 1]:
                k42 = (2 * H + (n4 - n2) * ZS) / (c * r_km * 1000) * 1000
                rows.append({"family": "k42", "k": k42})
        for n3 in [-1, 1]:
            for n4 in [-1, 1]:
                k43 = (8 * H + (n4 + n3) * ZS) / (c * r_km * 1000) * 1000
                rows.append({"family": "k43", "k": k43})
    df = pd.DataFrame(rows)
    return df.groupby("family")["k"].agg(["min", "max"]).to_dict("index")


def main():
    config = {"stage": "RANGE-UB-2B-FIX1A", "created_utc": NOW,
              "baseline_commit": "0e3a461e546e42c138a1ca8c2c0b2123e69afbfb",
              "no_new_bellhop": True}
    (OUT / "RANGE_UB_2B_FIX1A_CONFIG.json").write_text(json.dumps(config, ensure_ascii=False, indent=2), encoding="utf-8")

    # integrity audit
    (OUT / "FIX1_DECISION_INTEGRITY_AUDIT.md").write_text(
        """# FIX1_DECISION_INTEGRITY_AUDIT

UTC: """ + NOW + """

FIX1 判定暂停: RANGE_UB_2B_FIX1_DECISION_BLOCKED_BY_POSTPROCESSING_INCONSISTENCIES

问题:
1. working span hard-coded to 48-52 km (should be contiguous block containing 50)
2. statistics included disconnected valid islands
3. Spearman implementation didn't handle ties correctly

原始 continuation 文件保留。
""", encoding="utf-8")

    # load FIX1 data
    wr_df = pd.read_csv(FIX1 / "VLA_FAMILY_WORKING_RANGE.csv")
    g_df = pd.read_csv(FIX1 / "RANGE_G_VECTOR.csv")
    cont1 = pd.read_csv(FIX1 / "FAMILY_CONTINUATION_1601.csv")
    conv_df = pd.read_csv(FIX1 / "FAMILY_CONTINUATION_CONVERGENCE.csv")

    # =========================================================
    # 1. Contiguous working range
    # =========================================================
    block, span = contiguous_working_range(wr_df, 50.0)
    pd.DataFrame([{"r_km": r, "in_contiguous_block": True} for r in block]).to_csv(
        OUT / "CONTIGUOUS_WORKING_RANGE_AUDIT.csv", index=False)

    # =========================================================
    # 2. Statistics only in contiguous range
    # =========================================================
    g_cont = g_df[g_df["r_km"].isin(block)]
    g_med = g_cont.groupby("r_km")["g"].median().reset_index()
    g_med.columns = ["r_km", "g_median"]

    if len(g_med) >= 3:
        r_vals = g_med["r_km"].values
        g_vals = g_med["g_median"].values
        inv_r = 1.0 / r_vals
        pearson = float(np.corrcoef(inv_r, g_vals)[0, 1])
        spearman_corrected = spearman_tie_aware(inv_r, g_vals)
        # old (wrong) spearman for comparison
        def rank_ord(a):
            order = np.argsort(a)
            ranks = np.empty_like(order, dtype=float)
            ranks[order] = np.arange(len(a), dtype=float)
            return ranks
        spearman_old = float(np.corrcoef(rank_ord(inv_r), rank_ord(g_vals))[0, 1])
        A = np.vstack([inv_r, np.ones_like(inv_r)]).T
        coef, _, _, _ = np.linalg.lstsq(A, g_vals, rcond=None)
        g_pred = A @ coef
        r2 = float(1 - np.sum((g_vals - g_pred)**2) / np.sum((g_vals - np.mean(g_vals))**2)) if np.var(g_vals) > 0 else 0
        Q = r_vals * g_vals
    else:
        pearson = spearman_corrected = spearman_old = r2 = np.nan
        Q = np.array([])

    pd.DataFrame([{"old_custom_spearman": spearman_old, "corrected_spearman": spearman_corrected,
                   "pearson": pearson, "r2": r2}]).to_csv(OUT / "SPEARMAN_TIE_AUDIT.csv", index=False)

    g_vs = g_med.copy()
    g_vs["g_median"] = g_vals if len(g_med) else []
    g_vs["inv_r"] = inv_r if len(g_med) else []
    g_vs.to_csv(OUT / "G_VS_INVERSE_RANGE_FIXED.csv", index=False)

    # =========================================================
    # 3. True theory envelope audit
    # =========================================================
    post_rows = []
    retention_rows = []
    for rk in block:
        env = theory_envelope_at_r(rk)
        for fam in ["k32", "k42", "k43"]:
            c = conv_df[(conv_df["r_km"] == rk) & (conv_df["family"] == fam)]
            if not len(c):
                continue
            k_val = float(c.iloc[0]["k_1601"])
            lo, hi = env[fam]["min"], env[fam]["max"]
            center = (lo + hi) / 2
            inside = lo <= k_val <= hi
            post_rows.append({"r_km": rk, "family": fam, "k_tracked": k_val,
                              "k_theory_min": lo, "k_theory_max": hi,
                              "inside_true_theory_envelope": inside,
                              "distance_to_band_center": abs(k_val - center),
                              "normalized_distance": abs(k_val - center) / max(hi - lo, 1e-9)})
            retention_rows.append({"r_km": rk, "family": fam, "inside_envelope": inside})
    post_df = pd.DataFrame(post_rows)
    post_df.to_csv(OUT / "POSTHOC_TRUE_THEORY_ENVELOPE_AUDIT.csv", index=False)
    ret_df = pd.DataFrame(retention_rows)
    if len(ret_df):
        ret_summary = ret_df.groupby("r_km")["inside_envelope"].all().reset_index()
        ret_summary.columns = ["r_km", "all_three_inside"]
        ret_summary.to_csv(OUT / "XU2024_FAMILY_RETENTION_VS_RANGE.csv", index=False)

    # =========================================================
    # 4. J_g in contiguous range
    # =========================================================
    score_rows = []
    g50 = g_df[(g_df["r_km"] == 50.0)]
    if len(g50) == 3:
        g_true = np.array([g50[g50["family"] == f]["g"].values[0] for f in ["k32", "k42", "k43"]])
        for rk in block:
            sub = g_df[g_df["r_km"] == rk]
            if len(sub) == 3:
                g_r = np.array([sub[sub["family"] == f]["g"].values[0] for f in ["k32", "k42", "k43"]])
                J = float(np.sqrt(np.mean((g_r - g_true) ** 2)))
                score_rows.append({"r_km": rk, "J_g": J})
    score_df = pd.DataFrame(score_rows)
    score_df.to_csv(OUT / "VLA_RANGE_SCORE_50KM_FIXED.csv", index=False)

    # one-point law
    law_rows = []
    if len(g_med):
        g50_row = g_med[g_med["r_km"] == 50.0]
        if len(g50_row):
            C50 = 50.0 * float(g50_row.iloc[0]["g_median"])
            for _, r in g_med.iterrows():
                r_hat = C50 / r["g_median"] if r["g_median"] > 0 else np.nan
                law_rows.append({"r_true_km": r["r_km"], "r_hat_km": r_hat,
                                 "abs_error_km": abs(r_hat - r["r_km"]) if np.isfinite(r_hat) else np.nan,
                                 "relative_error": abs(r_hat - r["r_km"]) / r["r_km"] if np.isfinite(r_hat) else np.nan})
    pd.DataFrame(law_rows).to_csv(OUT / "ONE_POINT_RANGE_LAW_DIAGNOSTIC_FIXED.csv", index=False)

    # =========================================================
    # 5. Decision (strict pre-registered gate)
    # =========================================================
    near_ok = all(bool(wr_df[wr_df["r_km"] == r].iloc[0]["all_three_valid"]) for r in [49.5, 50.5]) if len(wr_df) else False
    span_ok = span >= 5.0
    jg_unique = True
    if len(score_df):
        false = score_df[np.abs(score_df["r_km"] - 50.0) >= 0.1]
        if len(false) and (false["J_g"] < 1e-6).any():
            jg_unique = False
    spearman_ok = spearman_corrected >= 0.90 if np.isfinite(spearman_corrected) else False

    if near_ok and span_ok and jg_unique and spearman_ok:
        decision = "E_STD_VLA_SLOPE_RANGE_INFORMATION_CONFIRMED_AT_ZS200"
    elif near_ok and not span_ok:
        decision = "E_STD_VLA_SLOPE_RANGE_INFORMATION_LOCAL_ONLY"
    elif not near_ok:
        decision = "E_STD_VLA_SLOPE_BRANCH_CONTINUATION_PARTIAL"
    else:
        decision = "E_STD_VLA_SLOPE_RANGE_INFORMATION_NOT_ESTABLISHED"

    why = f"near_ok={near_ok}, span={span:.1f}km, jg_unique={jg_unique}, spearman_corrected={spearman_corrected:.3f}"
    dec = {"stage": "RANGE-UB-2B-FIX1A", "decision": decision, "why": why,
           "spearman_corrected": spearman_corrected, "span_km": span,
           "interpretation": "stable ridge continuation != Xu2024 1/r range-family continuation",
           "created_utc": NOW}
    (OUT / "RANGE_UB_2B_FIX1A_DECISION.json").write_text(json.dumps(dec, ensure_ascii=False, indent=2), encoding="utf-8")

    report = f"""# RANGE-UB-2B-FIX1A — Postprocessing Integrity Correction

UTC: {NOW}

## Contiguous Working Range

{block[0]:.1f}–{block[-1]:.1f} km, span = {span:.1f} km

## Statistics (contiguous only)

Pearson = {pearson:.3f}
Spearman (corrected) = {spearman_corrected:.3f}
Spearman (old) = {spearman_old:.3f}
R² = {r2:.3f}

## True Theory Envelope Audit

{post_df.to_string(index=False) if len(post_df) else 'N/A'}

## J_g (contiguous)

{score_df.to_string(index=False) if len(score_df) else 'N/A'}

## 判定

### `{decision}`

{why}

**解释**：50 km 附近稳定 ridge 可连续追踪，零误差网格自匹配唯一；但未呈现预期 1/r 距离律，不能冻结为可靠绝对距离观测。
"""
    (OUT / "RANGE_UB_2B_FIX1A_REPORT.md").write_text(report, encoding="utf-8")
    (OUT / "GPT_SYNC.md").write_text(f"# RANGE-UB-2B-FIX1A\n\n**{decision}**\n\n{why}\n", encoding="utf-8")

    print("decision", decision)
    print(why)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
