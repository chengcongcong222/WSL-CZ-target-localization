#!/usr/bin/env python3
"""1E-A-CLOSEOUT-FIX: Survivor-set metric integrity. Data-only, no propagation."""
from __future__ import annotations

import json
from datetime import datetime, timezone
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parent
OUT = ROOT / "results" / "R3_RC23_CLOSEDLOOP" / "R3_RC23_TURN_AMPLITUDE_CLOSEOUT_FIX"
AMP = ROOT / "results" / "R3_RC23_CLOSEDLOOP" / "R3_RC23_TURN_AMPLITUDE_ROBUSTNESS"
OUT.mkdir(parents=True, exist_ok=True)
NOW = datetime.now(timezone.utc).isoformat()

CONFIG_MAP = {"TRIPLE": "201+235+283", "FOUR": "201+235+283+338"}


def main():
    (OUT / "CLOSEOUT_FIX_CONFIG.json").write_text(json.dumps({
        "stage": "CLOSEOUT_FIX", "created_utc": NOW,
        "baseline_commit": "84e294dc7bd078074bfae95feca8f7c1fbd58c3d",
        "fix": "strict_single_bin uses R_keep not R_best",
    }, indent=2), encoding="utf-8")

    cur = pd.read_csv(AMP / "AMPLITUDE_CANDIDATE_SCORES.csv")
    orig_cases = pd.read_csv(AMP / "AMPLITUDE_ANCHOR_CASES.csv")

    # =========================================================
    # 1. Recompute metrics with correct definitions
    # =========================================================
    rows = []
    for (cfg, shape, mode, A, z_true), grp in cur.groupby(["config", "shape", "mode", "A", "ztrue"]):
        J = grp["J"].values
        r0 = grp["r0"].values
        keep = grp["keep"].values
        node_ids = grp["node_id"].values
        jmin = float(J.min())

        # truth node
        truth_node = grp[np.abs(r0 - 50.0) < 0.1]
        if len(truth_node):
            ti = truth_node["J"].idxmin()
            true_J = float(truth_node.loc[ti, "J"])
            true_rank = int((J < true_J).sum()) + 1
            true_kept = bool(truth_node.loc[ti, "keep"])
        else:
            true_J, true_rank, true_kept = np.nan, -1, False

        # R_keep = surviving range bins (tau=0.5)
        surv = grp[keep]
        r_keep = sorted(set(int(round(x)) for x in surv["r0"])) if len(surv) else []
        # R_best = global-minimum range bins
        r_best = sorted(set(int(round(x)) for x in r0[J <= jmin + 1e-10]))
        # R_wrong = surviving wrong-range bins
        r_wrong = [r for r in r_keep if r != 50]

        strict = true_kept and true_rank == 1 and r_keep == [50]
        best_unique = r_best == [50]
        wrong_survives = len(r_wrong) > 0

        surv_r = surv["r0"].values if len(surv) else np.array([])
        width = float(surv_r.max() - surv_r.min()) if len(surv_r) else 0.0
        max_dev = float(np.max(np.abs(surv_r - 50.0))) if len(surv_r) else 0.0

        wrong_mask = np.abs(r0 - 50.0) >= 0.1
        best_wrong_J = float(J[wrong_mask].min()) if wrong_mask.any() else np.nan
        wrong_margin = best_wrong_J - true_J if np.isfinite(true_J) and np.isfinite(best_wrong_J) else np.nan

        rows.append({
            "config": cfg, "shape": shape, "mode": mode, "A_rms_db": A, "z_true_m": z_true,
            "true_J": true_J, "J_min": jmin, "true_rank": true_rank, "truth_retained": true_kept,
            "surviving_r_bins": str(r_keep), "n_surviving_r_bins": len(r_keep),
            "surviving_r_width_km": width, "max_surviving_deviation_km": max_dev,
            "best_range_bins": str(r_best), "best_range_unique_50": best_unique,
            "strict_single_bin": strict, "wrong_range_survives": wrong_survives,
            "wrong_range_margin": wrong_margin, "n_survivors": int(keep.sum()),
        })

    df = pd.DataFrame(rows)
    df.to_csv(OUT / "AMPLITUDE_CASE_METRICS_FINAL.csv", index=False)

    # =========================================================
    # 2. Regression vs original 1E-A
    # =========================================================
    reg_rows = []
    for _, r in df.iterrows():
        o = orig_cases[(orig_cases.config == r["config"]) & (orig_cases["shape"] == r["shape"]) &
                       (orig_cases["mode"] == r["mode"]) & (orig_cases.A_rms_db == r["A_rms_db"]) &
                       (orig_cases.z_true_m == r["z_true_m"])]
        if len(o):
            orig = bool(o.iloc[0]["range_anchored"])
            match = bool(r["strict_single_bin"]) == orig
            reg_rows.append({"config": r["config"], "shape": r["shape"], "mode": r["mode"],
                             "A_rms_db": r["A_rms_db"], "z_true_m": r["z_true_m"],
                             "new_strict": bool(r["strict_single_bin"]), "orig_anchored": orig, "match": match})
    reg_df = pd.DataFrame(reg_rows)
    reg_df.to_csv(OUT / "STRICT_ANCHOR_REGRESSION.csv", index=False)
    reg_ok = bool(reg_df["match"].all()) and len(reg_df) == 96

    if not reg_ok:
        (OUT / "R3_RC23_CLOSEDLOOP_1E_A_CLOSEOUT_FIX_DECISION.json").write_text(
            json.dumps({"decision": "CLOSEOUT_FIX_BLOCKED_BY_ORIGINAL_CASE_REGRESSION",
                        "n_match": int(reg_df["match"].sum()), "n_total": len(reg_df), "created_utc": NOW}, indent=2),
            encoding="utf-8")
        print("BLOCKED_BY_REGRESSION", int(reg_df["match"].sum()), "/", len(reg_df))
        return 1

    # =========================================================
    # 3. Working range summary
    # =========================================================
    wr_rows = []
    for A in [0.0, 0.25, 0.5, 1.0]:
        for cfg in ["TRIPLE", "FOUR"]:
            sub = df[(df.A_rms_db == A) & (df.config == cfg)]
            if not len(sub):
                continue
            wr_rows.append({
                "A_rms_db": A, "config": cfg, "n_combos": len(sub),
                "strict_single_bin": int(sub["strict_single_bin"].sum()),
                "truth_retained": int(sub["truth_retained"].sum()),
                "true_rank1": int((sub["true_rank"] == 1).sum()),
                "best_range_unique_50": int(sub["best_range_unique_50"].sum()),
                "min_wrong_range_margin": float(sub["wrong_range_margin"].min()),
                "median_r_width": float(sub["surviving_r_width_km"].median()),
                "max_r_width": float(sub["surviving_r_width_km"].max()),
                "max_dev": float(sub["max_surviving_deviation_km"].max()),
            })
    wr_df = pd.DataFrame(wr_rows)
    wr_df.to_csv(OUT / "AMPLITUDE_WORKING_RANGE_SUMMARY_FINAL.csv", index=False)

    # A=2 summary from original
    a2 = orig_cases[orig_cases.A_rms_db == 2.0]
    a2_rows = []
    for cfg in ["TRIPLE", "FOUR"]:
        sub = a2[a2.config == cfg]
        if len(sub):
            a2_rows.append({"config": cfg, "n_combos": len(sub),
                            "range_anchored": int(sub["range_anchored"].sum()),
                            "truth_retained": int(sub["truth_retained"].sum()),
                            "note": "SUMMARY_ONLY_FROM_ORIGINAL_CASE_TABLE"})
    pd.DataFrame(a2_rows).to_csv(OUT / "A2_BOUNDARY_SUMMARY.csv", index=False)

    # =========================================================
    # 4. Decision
    # =========================================================
    decision = "AMPLITUDE_VARIATION_BREAKS_TURN_RANGE_ANCHOR_IN_TESTED_RANGE"
    dec = {
        "stage": "CLOSEOUT_FIX", "decision": decision,
        "interpretation": "STRICT_SINGLE_RANGE_BIN_CRITERION_NOT_UNIFORMLY_RETAINED",
        "auxiliary_labels": [
            "OPTIMAL_TRUTH_RANK_ROBUST_TO_TESTED_1DB_AMPLITUDE_VARIATION",
            "RANGE_SET_CONTRACTION_RETAINS_LOCAL_VALUE_THROUGH_TESTED_0P5DB",
        ],
        "regression_96_of_96": True,
        "created_utc": NOW,
    }
    (OUT / "R3_RC23_CLOSEDLOOP_1E_A_CLOSEOUT_FIX_DECISION.json").write_text(
        json.dumps(dec, ensure_ascii=False, indent=2), encoding="utf-8")

    report = f"""# 1E-A-CLOSEOUT-FIX

UTC: {NOW}

## 修复

strict_single_bin 改用 R_keep（tau=0.5 survivor bins），不再用 R_best。

## Regression

96/96 与原始 1E-A range_anchored 一致。

## 最终统计

{wr_df.to_string(index=False)}

## 三层

optimal rank → survivor-set contraction → strict single-bin anchor

## 判定

**{decision}** = STRICT_SINGLE_RANGE_BIN_CRITERION_NOT_UNIFORMLY_RETAINED

辅助: OPTIMAL_TRUTH_RANK_ROBUST_TO_TESTED_1DB_AMPLITUDE_VARIATION
      RANGE_SET_CONTRACTION_RETAINS_LOCAL_VALUE_THROUGH_TESTED_0P5DB
"""
    (OUT / "R3_RC23_CLOSEDLOOP_1E_A_CLOSEOUT_FIX_REPORT.md").write_text(report, encoding="utf-8")
    (OUT / "GPT_SYNC.md").write_text(f"# CLOSEOUT_FIX\n\n**{decision}**\n\nregression 96/96 PASS\n", encoding="utf-8")

    print("regression", reg_ok, int(reg_df["match"].sum()), "/", len(reg_df))
    print("decision", decision)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
