#!/usr/bin/env python3
"""RANGE-UB-2B-FIX0B: Theory-envelope and family-peak integrity fix.

Fixes:
1. k43 enumeration (independent n3, n4 loops)
2. Local max check on full P(k) (not just within band)
3. Per-zs theory envelope
4. Boundary diagnostics + prominence
"""
from __future__ import annotations

import json
from datetime import datetime, timezone
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parent
OUT = ROOT / "results" / "RANGE_UB" / "E_STD_VLA_SLOPE_FAMILY_ID_FIX"
FIX0 = ROOT / "results" / "RANGE_UB" / "E_STD_VLA_DELAY_DEPTH_RIDGE"
OUT.mkdir(parents=True, exist_ok=True)
NOW = datetime.now(timezone.utc).isoformat()

H = 5000.0
R = 50000.0
ZS_LIST = [180.0, 200.0, 220.0]
C_MIN, C_MAX = 1500.0, 1635.2
NBEAMS = [1601, 3201]
DK = 0.002  # ms/m


def theory_envelope_by_zs():
    """Enumerate k32/k42/k43 with correct independent loops."""
    rows = []
    for zs in ZS_LIST:
        for c in np.linspace(C_MIN, C_MAX, 20):
            for n2 in [-1, 1]:
                for n3 in [-1, 1]:
                    # k32 = [-6H - (n3+n2)*zs] / (cR)
                    k32 = (-6 * H - (n3 + n2) * zs) / (c * R) * 1000
                    rows.append({"zs": zs, "c": c, "family": "k32", "n_vals": f"n2={n2},n3={n3}", "k_ms_per_m": k32})
            for n2 in [-1, 1]:
                for n4 in [-1, 1]:
                    # k42 = [2H + (n4-n2)*zs] / (cR)
                    k42 = (2 * H + (n4 - n2) * zs) / (c * R) * 1000
                    rows.append({"zs": zs, "c": c, "family": "k42", "n_vals": f"n2={n2},n4={n4}", "k_ms_per_m": k42})
            # k43: INDEPENDENT n3, n4 loops
            for n3 in [-1, 1]:
                for n4 in [-1, 1]:
                    k43 = (8 * H + (n4 + n3) * zs) / (c * R) * 1000
                    rows.append({"zs": zs, "c": c, "family": "k43", "n_vals": f"n3={n3},n4={n4}", "k_ms_per_m": k43})
    return pd.DataFrame(rows)


def load_radon_map(nb, zs):
    fp = FIX0 / f"RADON_SLOPE_MAP_nb{nb}_zs{int(zs)}.csv"
    return pd.read_csv(fp) if fp.exists() else None


def slope_profile(df):
    return df.groupby("k_ms_per_m")["score"].max().reset_index().sort_values("k_ms_per_m").reset_index(drop=True)


def find_family_peak_full(prof, k_lo, k_hi, full_df=None):
    """Find max in band, but check local max on FULL profile."""
    mask = (prof["k_ms_per_m"] >= k_lo) & (prof["k_ms_per_m"] <= k_hi)
    sub = prof[mask]
    if len(sub) == 0:
        return None
    ks_all = prof["k_ms_per_m"].values
    ps_all = prof["score"].values
    best_i = int(np.argmax(sub["score"].values))
    k_star = float(sub["k_ms_per_m"].values[best_i])
    p_star = float(sub["score"].values[best_i])

    # find k_star position in FULL profile
    full_idx = int(np.argmin(np.abs(ks_all - k_star)))
    # check neighbors in FULL profile
    p_left = float(ps_all[full_idx - 1]) if full_idx > 0 else -np.inf
    p_right = float(ps_all[full_idx + 1]) if full_idx < len(ps_all) - 1 else -np.inf
    is_local = p_star > p_left and p_star >= p_right

    # boundary distance in bins
    bins_lo = int(np.sum(ks_all < k_star))
    bins_hi = int(np.sum(ks_all > k_star))
    at_edge = (k_star - k_lo < DK * 1.5) or (k_hi - k_star < DK * 1.5)

    b_star = np.nan
    if full_df is not None:
        rows = full_df[full_df["k_ms_per_m"] == k_star]
        if len(rows) > 0:
            b_star = float(rows.loc[rows["score"].idxmax(), "b_ms"])

    prominence = p_star / (0.5 * (p_left + p_right)) if (p_left > -np.inf and p_right > -np.inf and (p_left + p_right) > 0) else np.nan

    return {
        "k_star": k_star, "b_star": b_star, "score": p_star,
        "full_profile_local_max": is_local,
        "p_left": p_left, "p_right": p_right,
        "at_band_edge": at_edge,
        "distance_to_lower_edge_bins": bins_lo,
        "distance_to_upper_edge_bins": bins_hi,
        "prominence_ratio": prominence,
    }


def main():
    config = {"stage": "RANGE-UB-2B-FIX0B", "created_utc": NOW,
              "baseline_commit": "1772f34628b71f801ec512836ca7eb2fa8a1bd55",
              "fixes": ["k43 independent n3/n4 loops", "full P(k) local-max check", "per-zs envelope"]}
    (OUT / "RANGE_UB_2B_FIX0B_CONFIG.json").write_text(json.dumps(config, ensure_ascii=False, indent=2), encoding="utf-8")

    # =========================================================
    # 1. Theory enumeration
    # =========================================================
    env_df = theory_envelope_by_zs()
    env_df.to_csv(OUT / "THEORY_ENUMERATION_AUDIT.csv", index=False)
    env_by_zs = env_df.groupby(["zs", "family"])["k_ms_per_m"].agg(["min", "max"]).reset_index()
    env_by_zs.columns = ["zs", "family", "k_min_ms_per_m", "k_max_ms_per_m"]
    env_by_zs.to_csv(OUT / "E_STD_XU2024_SLOPE_ENVELOPE_BY_ZS.csv", index=False)

    # =========================================================
    # 2. Family peak detection with full-profile local-max
    # =========================================================
    peak_rows = []
    boundary_rows = []
    prominence_rows = []
    conv_rows = []
    g_rows = []

    for zs in ZS_LIST:
        maps = {}
        full_maps = {}
        for nb in NBEAMS:
            df = load_radon_map(nb, zs)
            if df is not None:
                maps[nb] = slope_profile(df)
                full_maps[nb] = df

        for fam in ["k32", "k42", "k43"]:
            env_row = env_by_zs[(env_by_zs["zs"] == zs) & (env_by_zs["family"] == fam)]
            if not len(env_row):
                continue
            k_lo = float(env_row.iloc[0]["k_min_ms_per_m"])
            k_hi = float(env_row.iloc[0]["k_max_ms_per_m"])

            results = {}
            for nb in NBEAMS:
                prof = maps.get(nb)
                if prof is None:
                    continue
                peak = find_family_peak_full(prof, k_lo, k_hi, full_maps.get(nb))
                results[nb] = peak
                if peak:
                    peak_rows.append({"zs": zs, "family": fam, "nbeams": nb, "k_lo": k_lo, "k_hi": k_hi,
                                      "k_star": peak["k_star"], "score": peak["score"],
                                      "full_profile_local_max": peak["full_profile_local_max"],
                                      "at_band_edge": peak["at_band_edge"], "detected": peak["full_profile_local_max"]})
                    boundary_rows.append({"zs": zs, "family": fam, "nbeams": nb,
                                          "k_lo": k_lo, "k_hi": k_hi, "k_star": peak["k_star"],
                                          "distance_to_lower_edge_bins": peak["distance_to_lower_edge_bins"],
                                          "distance_to_upper_edge_bins": peak["distance_to_upper_edge_bins"],
                                          "p_left": peak["p_left"], "p_star": peak["score"], "p_right": peak["p_right"],
                                          "full_profile_local_max": peak["full_profile_local_max"]})
                    prominence_rows.append({"zs": zs, "family": fam, "nbeams": nb,
                                            "k_star": peak["k_star"], "prominence_ratio": peak["prominence_ratio"]})

            # convergence (only if both have full_profile_local_max)
            r1, r2 = results.get(1601), results.get(3201)
            if r1 and r2 and r1["full_profile_local_max"] and r2["full_profile_local_max"]:
                sign_ok = np.sign(r1["k_star"]) == np.sign(r2["k_star"])
                rel_diff = abs(r1["k_star"] - r2["k_star"]) / max(abs(r2["k_star"]), 0.02)
                b_diff = abs(r1["b_star"] - r2["b_star"]) if np.isfinite(r1["b_star"]) and np.isfinite(r2["b_star"]) else 999
                ok = sign_ok and rel_diff <= 0.05 and b_diff <= 50.0
                conv_rows.append({"zs": zs, "family": fam, "k_1601": r1["k_star"], "k_3201": r2["k_star"],
                                  "rel_diff": rel_diff, "b_diff_ms": b_diff, "sign_ok": sign_ok, "pass": ok})
            else:
                conv_rows.append({"zs": zs, "family": fam, "pass": False, "note": "NOT_LOCAL_MAX_OR_MISSING"})

            # g-vector (only if local max detected at 1601)
            if r1 and r1["full_profile_local_max"]:
                if fam == "k32":
                    g_rows.append({"zs": zs, "family": fam, "k": r1["k_star"], "g": -r1["k_star"] / 6})
                elif fam == "k42":
                    g_rows.append({"zs": zs, "family": fam, "k": r1["k_star"], "g": r1["k_star"] / 2})
                else:
                    g_rows.append({"zs": zs, "family": fam, "k": r1["k_star"], "g": r1["k_star"] / 8})
            else:
                g_rows.append({"zs": zs, "family": fam, "k": np.nan, "g": np.nan})

    peak_df = pd.DataFrame(peak_rows)
    peak_df.to_csv(OUT / "FAMILY_PEAK_DETECTION.csv", index=False)
    pd.DataFrame(boundary_rows).to_csv(OUT / "FAMILY_PEAK_BOUNDARY_AUDIT.csv", index=False)
    pd.DataFrame(prominence_rows).to_csv(OUT / "FAMILY_RIDGE_PROMINENCE.csv", index=False)
    conv_df = pd.DataFrame(conv_rows)
    conv_df.to_csv(OUT / "XU2024_FAMILY_RIDGE_CONVERGENCE_FIX.csv", index=False)
    g_df = pd.DataFrame(g_rows)
    g_df.to_csv(OUT / "THEORY_FAMILY_G_VECTOR_FIX.csv", index=False)

    # =========================================================
    # 3. Decision
    # =========================================================
    n_local = int(peak_df["detected"].sum()) if len(peak_df) else 0
    all_9 = bool(conv_df["pass"].all()) and len(conv_df) == 9
    # check zs=200 all three families
    z200_ok = all(
        (conv_df[(conv_df["zs"] == 200) & (conv_df["family"] == fam)]["pass"].any())
        for fam in ["k32", "k42", "k43"]
    ) if len(conv_df) else False

    if all_9 and n_local == 18:
        decision = "XU2024_SECONDARY_SLOPE_FAMILIES_IDENTIFIED_AT_50KM"
    elif n_local >= 12 and z200_ok:
        decision = "XU2024_SECONDARY_SLOPE_FAMILIES_PARTIAL_AT_50KM"
    elif n_local >= 6:
        decision = "THEORY_BAND_ENERGY_PRESENT_BUT_FAMILY_RIDGE_IDENTITY_PARTIAL"
    else:
        decision = "GENERIC_RIDGES_EXIST_XU2024_FAMILY_IDENTITY_NOT_ESTABLISHED"

    why = f"n_local_max={n_local}/18, all_9_conv={all_9}, zs200_ok={z200_ok}"
    dec = {"stage": "RANGE-UB-2B-FIX0B", "decision": decision, "why": why, "created_utc": NOW}
    (OUT / "RANGE_UB_2B_FIX0B_DECISION.json").write_text(json.dumps(dec, ensure_ascii=False, indent=2), encoding="utf-8")

    report = f"""# RANGE-UB-2B-FIX0B — Theory-Envelope & Family-Peak Integrity

UTC: {NOW}

## 修复

1. k43 枚举：独立 n3/n4 循环（原 n3 固定 +1）
2. 局部极大：完整 P(k) 判断（含 band 外邻点）
3. Per-zs 理论包络

## 理论包络（修正后）

{env_by_zs.to_string(index=False)}

## Family Peak Detection（完整 P(k) 局部极大）

{peak_df.to_string(index=False) if len(peak_df) else 'N/A'}

## Boundary Audit

{pd.DataFrame(boundary_rows).to_string(index=False) if boundary_rows else 'N/A'}

## Convergence

{conv_df.to_string(index=False)}

## g-vector

{g_df.to_string(index=False) if len(g_df) else 'N/A'}

## 判定

### `{decision}`

{why}
"""
    (OUT / "RANGE_UB_2B_FIX0B_REPORT.md").write_text(report, encoding="utf-8")
    (OUT / "GPT_SYNC.md").write_text(f"# RANGE-UB-2B-FIX0B\n\n**{decision}**\n\n{why}\n", encoding="utf-8")

    print("decision", decision)
    print(why)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
