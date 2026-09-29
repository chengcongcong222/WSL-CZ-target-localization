#!/usr/bin/env python3
"""RANGE-UB-2B: E-STD VLA secondary-delay-slope reference upper bound.

BELLHOP arrivals at VLA depths -> path class mapping -> cluster delays
-> slope fitting for tau32, tau42, tau43 -> range profile.
"""
from __future__ import annotations

import json
import subprocess
from datetime import datetime, timezone
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parent
OUT = ROOT / "results" / "RANGE_UB" / "E_STD_VLA_SECONDARY_SLOPE"
ZGRID = ROOT / "results" / "R3_C2_Yang_SA_depth" / "R3_C2_1" / "_kraken_zgrid"
AT_BIN = ROOT / "tools" / "acoustics_toolbox" / "atWin10" / "at" / "bin"
WORK = OUT / "_bellhop"
OUT.mkdir(parents=True, exist_ok=True)
WORK.mkdir(parents=True, exist_ok=True)
NOW = datetime.now(timezone.utc).isoformat()

FREQ = 235.0
H_DEPTH = 5000.0
ZBOX_M = 5001.0
RBOX_KM = 70.0
Z_S_ALL = np.arange(150.0, 250.0 + 1e-9, 5.0)  # 21
R_ALL_KM = np.arange(45.0, 60.0 + 1e-9, 0.5)  # 31
ZR_VLA = np.arange(20.0, 1620.0 + 1e-9, 10.0)  # 161
Z_TRUE_LIST = [180.0, 200.0, 220.0]
NBEAMS_CONV = [1601, 3201]
CONV_RS = [45.0, 50.0, 56.0, 60.0]
CONV_ZS = [180.0, 200.0, 220.0]
NBEAMS_MAIN = 1601
SLOPE_REL_TOL = 0.05
MIN_DEPTH_FRAC = 0.70
MIN_DEPTH_SPAN = 1200.0
MIN_R2 = 0.95


def write_env(path, zs, nbeams, tag):
    src = ZGRID / "zgrid_f235.env"
    lines = src.read_text(encoding="utf-8").splitlines()
    ssp = []
    for ln in lines[5:]:
        p = ln.split()
        if len(p) >= 2:
            try:
                ssp.append(f"  {float(p[0]):.1f} {float(p[1]):.4f} /")
            except ValueError:
                break
        else:
            break
    n_rz = len(ZR_VLA)
    out = [f"'{tag}'", f"{FREQ:.1f}", "1", "'CVN'", f"  51 0.0 {H_DEPTH:.1f}", *ssp,
           "'R' 0.0", "1", f"{zs:.1f} /",
           str(n_rz),
           f"{ZR_VLA[0]:.1f} {ZR_VLA[-1]:.1f} /",
           f"{len(R_ALL_KM)}", f"{R_ALL_KM[0]:.1f} {R_ALL_KM[-1]:.1f} /",
           "'A'", str(nbeams), "-85.0 85.0 /",
           f"0.0 {ZBOX_M:.1f} {RBOX_KM:.1f}"]
    path.write_text("\n".join(out) + "\n", encoding="utf-8")


def parse_arr(fp):
    if not fp.exists():
        return []
    lines = [l.strip() for l in fp.read_text(encoding="utf-8", errors="replace").splitlines() if l.strip()]
    if not lines:
        return []
    h = lines[0].split()
    Nsd, Nrd, Nrr = int(float(h[1])), int(float(h[2])), int(float(h[3]))
    idx = 1
    sd = [float(x) for x in lines[idx].split()[:Nsd]]; idx += 1
    rd = [float(x) for x in lines[idx].split()[:Nrd]]; idx += 1
    rr = [float(x) for x in lines[idx].split()[:Nrr]]; idx += 1
    out = []
    for isd in range(Nsd):
        if idx >= len(lines): break
        idx += 1  # Narrmx
        for ird in range(Nrd):
            for ir in range(Nrr):
                if idx >= len(lines): break
                narr = int(float(lines[idx])); idx += 1
                for _ in range(narr):
                    if idx >= len(lines): break
                    p = lines[idx].split(); idx += 1
                    if len(p) >= 8:
                        out.append({
                            "zs": sd[isd] if isd < len(sd) else None,
                            "zr": rd[ird] if ird < len(rd) else None,
                            "r_km": rr[ir] / 1000 if ir < len(rr) else None,
                            "amp": float(p[0]), "delay_s": float(p[2]),
                            "rcv_ang": float(p[5]),
                            "n_top": int(float(p[6])), "n_bot": int(float(p[7])),
                        })
    return out


def classify_path(a):
    """Oracle path class: reflection_order = n_top+n_bot + sign(rcv_ang)."""
    ro = a["n_top"] + a["n_bot"]
    pos = a["rcv_ang"] > 0
    if ro == 1 and pos:
        return "C2"
    elif ro == 2 and not pos:
        return "C3"
    elif ro == 2 and pos:
        return "C4"
    elif ro == 1 and not pos:
        return "C1_neg"  # first-order negative (not in secondary set)
    elif ro == 0:
        return "C1_rev"  # reversed/refracted
    else:
        return f"OTHER_ro{ro}"


def cluster_delay(arrs, zr_val, cluster):
    """Energy-weighted delay for a cluster at given zr."""
    sel = [a for a in arrs if abs(a["zr"] - zr_val) < 0.5 and classify_path(a) == cluster]
    if not sel:
        return np.nan, 0
    w = np.array([a["amp"] ** 2 for a in sel])
    t = np.array([a["delay_s"] for a in sel])
    if w.sum() <= 0:
        return np.nan, len(sel)
    return float(np.sum(w * t) / np.sum(w)), len(sel)


def fit_slope(zr_vals, tau_vals):
    """OLS fit tau = k*zr + b. Returns k (ms/m), intercept, R2, n_valid, span."""
    mask = np.isfinite(zr_vals) & np.isfinite(tau_vals)
    z, t = zr_vals[mask], tau_vals[mask]
    n = len(z)
    if n < 3:
        return np.nan, np.nan, 0.0, n, 0.0
    span = float(z.max() - z.min())
    # OLS
    zmean, tmean = z.mean(), t.mean()
    dz = z - zmean
    dt = t - tmean
    k = float(np.sum(dz * dt) / np.sum(dz ** 2))  # ms/m
    b = float(tmean - k * zmean)
    t_pred = k * z + b
    ss_res = np.sum((t - t_pred) ** 2)
    ss_tot = np.sum((t - tmean) ** 2)
    r2 = float(1 - ss_res / ss_tot) if ss_tot > 0 else 0.0
    return k, b, r2, n, span


def main():
    config = {
        "stage": "RANGE-UB-2B", "created_utc": NOW,
        "baseline_commit": "e6e932e24a3a4345b2fee4a09bfb5256068ef29d",
        "observable": "XU2024_SECONDARY_DELAY_SLOPE_VECTOR",
        "clusters": "C2(ro=1,+), C3(ro=2,-), C4(ro=2,+)",
        "slope_method": "ORACLE_LINEAR_SLOPE_FIT",
        "vla": "PAPER_SCALE_VLA_REFERENCE_APERTURE zr=20:10:1620",
        "nbeams_conv": NBEAMS_CONV, "nbeams_main": NBEAMS_MAIN,
    }
    (OUT / "RANGE_UB_2B_CONFIG.json").write_text(json.dumps(config, ensure_ascii=False, indent=2), encoding="utf-8")

    # 2A correction
    (OUT / "2A_TABLE2_CONFLICT_CORRECTION.md").write_text(
        """# 2A_TABLE2_CONFLICT_CORRECTION

UTC: """ + NOW + """

Table 2 k42 行存在原文内部算术不一致：
- k42_est = 0.1088 ms/m
- R_printed = 53.84 km
- error_printed = 4.0%
- R_formula = 2H/(c*k42) = 52.007 km
- error_formula = 4.014%

至少一个印刷值有误，不能断言哪个。
正确标签：TABLE2_K42_INTERNAL_ARITHMETIC_INCONSISTENCY
主状态：XU2024_VLA_RANGE_FORMULA_NUMERICALLY_CLOSED_WITH_TABLE2_K42_CONFLICT
""", encoding="utf-8")

    # =========================================================
    # 1. Ray density convergence (12 points)
    # =========================================================
    conv_slopes = {}  # (nb, zs, rk) -> {k32, k42, k43}
    path_audit_rows = []
    for nb in NBEAMS_CONV:
        for zs in CONV_ZS:
            tag = f"vla_nb{nb}_zs{int(zs)}"
            write_env(WORK / f"{tag}.env", zs, nb, tag)
            subprocess.run([str(AT_BIN / "bellhop.exe"), tag], cwd=str(WORK), capture_output=True, timeout=600)
            arrs = parse_arr(WORK / f"{tag}.arr")
            by_r = {}
            for a in arrs:
                rk = round((a.get("r_km") or 0) * 2) / 2
                by_r.setdefault(rk, []).append(a)
            for rk in R_ALL_KM:
                ra = by_r.get(round(rk * 2) / 2, [])
                # cluster delays at each zr
                zr_list = sorted(set(a["zr"] for a in ra if a["zr"] is not None))
                t2_map, t3_map, t4_map = {}, {}, {}
                for zr in zr_list:
                    t2, n2 = cluster_delay(ra, zr, "C2")
                    t3, n3 = cluster_delay(ra, zr, "C3")
                    t4, n4 = cluster_delay(ra, zr, "C4")
                    t2_map[zr] = t2; t3_map[zr] = t3; t4_map[zr] = t4
                # slopes
                zr_arr = np.array(sorted(zr_list))
                tau32 = np.array([(t3_map.get(z, np.nan) - t2_map.get(z, np.nan)) * 1000 if np.isfinite(t3_map.get(z, np.nan)) and np.isfinite(t2_map.get(z, np.nan)) else np.nan for z in zr_arr])
                tau42 = np.array([(t4_map.get(z, np.nan) - t2_map.get(z, np.nan)) * 1000 if np.isfinite(t4_map.get(z, np.nan)) and np.isfinite(t2_map.get(z, np.nan)) else np.nan for z in zr_arr])
                tau43 = np.array([(t4_map.get(z, np.nan) - t3_map.get(z, np.nan)) * 1000 if np.isfinite(t4_map.get(z, np.nan)) and np.isfinite(t3_map.get(z, np.nan)) else np.nan for z in zr_arr])
                k32, _, r2_32, n32, sp32 = fit_slope(zr_arr, tau32)
                k42, _, r2_42, n42, sp42 = fit_slope(zr_arr, tau42)
                k43, _, r2_43, n43, sp43 = fit_slope(zr_arr, tau43)
                conv_slopes[(nb, zs, rk)] = {"k32": k32, "k42": k42, "k43": k43,
                                              "r2_32": r2_32, "r2_42": r2_42, "r2_43": r2_43,
                                              "n32": n32, "n42": n42, "n43": n43}
                # path audit at representative zr
                if rk in CONV_RS and zs in CONV_ZS and nb == 1601:
                    for zr in [200.0, 600.0, 1000.0, 1400.0]:
                        for a in ra:
                            if abs(a["zr"] - zr) < 0.5:
                                path_audit_rows.append({"zs": zs, "r_km": rk, "zr": zr,
                                                        "n_top": a["n_top"], "n_bot": a["n_bot"],
                                                        "rcv_ang": a["rcv_ang"], "delay_s": a["delay_s"],
                                                        "amp": a["amp"], "class": classify_path(a)})

    # convergence check
    conv_rows = []
    conv_pass = True
    for zs in CONV_ZS:
        for rk in CONV_RS:
            s1 = conv_slopes.get((1601, zs, rk), {})
            s2 = conv_slopes.get((3201, zs, rk), {})
            results = {}
            ok_all = True
            for key in ["k32", "k42", "k43"]:
                k1, k2 = s1.get(key, np.nan), s2.get(key, np.nan)
                if np.isfinite(k1) and np.isfinite(k2):
                    sign_ok = np.sign(k1) == np.sign(k2)
                    denom = max(abs(k2), 0.02)
                    rel_diff = abs(k1 - k2) / denom
                    ok = sign_ok and rel_diff <= SLOPE_REL_TOL
                    results[key] = {"k1601": k1, "k3201": k2, "rel_diff": rel_diff, "sign_ok": sign_ok, "ok": ok}
                    if not ok: ok_all = False
                else:
                    results[key] = {"ok": False}
                    ok_all = False
            conv_rows.append({"zs": zs, "r_km": rk,
                              "k32_1601": s1.get("k32", np.nan), "k32_3201": s2.get("k32", np.nan),
                              "k42_1601": s1.get("k42", np.nan), "k42_3201": s2.get("k42", np.nan),
                              "k43_1601": s1.get("k43", np.nan), "k43_3201": s2.get("k43", np.nan),
                              "pass": ok_all})
            if not ok_all: conv_pass = False

    conv_df = pd.DataFrame(conv_rows)
    conv_df.to_csv(OUT / "VLA_SLOPE_CONVERGENCE.csv", index=False)
    pd.DataFrame(path_audit_rows).to_csv(OUT / "VLA_CLUSTER_DELAY_SPOTCHECK.csv", index=False)

    # path class mapping audit
    (OUT / "PATH_CLASS_MAPPING_AUDIT.md").write_text(
        f"""# PATH_CLASS_MAPPING_AUDIT

UTC: {NOW}

## 分类规则 (ORACLE_PATH_CLASS_VLA_UPPER_BOUND)

| 类 | reflection_order | rcv_ang | 含义 |
|---|---|---|---|
| C2 | n_top+n_bot=1 | >0 | 一阶反射正到达角 |
| C3 | n_top+n_bot=2 | <0 | 二阶反射负到达角 |
| C4 | n_top+n_bot=2 | >0 | 二阶反射正到达角 |

## Spotcheck

共 {len(path_audit_rows)} 条记录于代表点。

## Convergence

{'PASS' if conv_pass else 'FAIL'} (1601 vs 3201, 12 points)
""", encoding="utf-8")

    if not conv_pass:
        (OUT / "RANGE_UB_2B_DECISION.json").write_text(
            json.dumps({"stage": "RANGE-UB-2B", "decision": "RANGE_UB_2B_BLOCKED_BY_PATH_CLASS_MAPPING",
                        "created_utc": NOW}, indent=2), encoding="utf-8")
        print("NOT_NUMERICALLY_STABLE")
        return 1

    # =========================================================
    # 2. Full slope library at 1601
    # =========================================================
    slope_rows = []
    for zs in Z_S_ALL:
        tag = f"vla_zs{int(zs)}"
        write_env(WORK / f"{tag}.env", zs, NBEAMS_MAIN, tag)
        subprocess.run([str(AT_BIN / "bellhop.exe"), tag], cwd=str(WORK), capture_output=True, timeout=600)
        arrs = parse_arr(WORK / f"{tag}.arr")
        by_r = {}
        for a in arrs:
            rk = round((a.get("r_km") or 0) * 2) / 2
            by_r.setdefault(rk, []).append(a)
        for rk in R_ALL_KM:
            ra = by_r.get(round(rk * 2) / 2, [])
            zr_list = sorted(set(a["zr"] for a in ra if a["zr"] is not None))
            t2_map, t3_map, t4_map = {}, {}, {}
            for zr in zr_list:
                t2, _ = cluster_delay(ra, zr, "C2")
                t3, _ = cluster_delay(ra, zr, "C3")
                t4, _ = cluster_delay(ra, zr, "C4")
                t2_map[zr] = t2; t3_map[zr] = t3; t4_map[zr] = t4
            zr_arr = np.array(sorted(zr_list))
            tau32 = np.array([(t3_map.get(z, np.nan) - t2_map.get(z, np.nan)) * 1000 if np.isfinite(t3_map.get(z, np.nan)) and np.isfinite(t2_map.get(z, np.nan)) else np.nan for z in zr_arr])
            tau42 = np.array([(t4_map.get(z, np.nan) - t2_map.get(z, np.nan)) * 1000 if np.isfinite(t4_map.get(z, np.nan)) and np.isfinite(t2_map.get(z, np.nan)) else np.nan for z in zr_arr])
            tau43 = np.array([(t4_map.get(z, np.nan) - t3_map.get(z, np.nan)) * 1000 if np.isfinite(t4_map.get(z, np.nan)) and np.isfinite(t3_map.get(z, np.nan)) else np.nan for z in zr_arr])
            k32, _, r2_32, n32, sp32 = fit_slope(zr_arr, tau32)
            k42, _, r2_42, n42, sp42 = fit_slope(zr_arr, tau42)
            k43, _, r2_43, n43, sp43 = fit_slope(zr_arr, tau43)

            # admission
            n_zr = len(zr_list)
            frac32 = n32 / n_zr if n_zr > 0 else 0
            frac42 = n42 / n_zr if n_zr > 0 else 0
            frac43 = n43 / n_zr if n_zr > 0 else 0
            valid32 = frac32 >= MIN_DEPTH_FRAC and sp32 >= MIN_DEPTH_SPAN and r2_32 >= MIN_R2 and np.isfinite(k32) and k32 < 0
            valid42 = frac42 >= MIN_DEPTH_FRAC and sp42 >= MIN_DEPTH_SPAN and r2_42 >= MIN_R2 and np.isfinite(k42) and k42 > 0
            valid43 = frac43 >= MIN_DEPTH_FRAC and sp43 >= MIN_DEPTH_SPAN and r2_43 >= MIN_R2 and np.isfinite(k43) and k43 > 0

            # g-vector
            g32 = -k32 / 6 if np.isfinite(k32) else np.nan
            g42 = k42 / 2 if np.isfinite(k42) else np.nan
            g43 = k43 / 8 if np.isfinite(k43) else np.nan

            slope_rows.append({
                "zs_m": zs, "r_km": rk,
                "k32_ms_per_m": k32, "k42_ms_per_m": k42, "k43_ms_per_m": k43,
                "r2_32": r2_32, "r2_42": r2_42, "r2_43": r2_43,
                "n32": n32, "n42": n42, "n43": n43,
                "span32": sp32, "span42": sp42, "span43": sp43,
                "valid32": valid32, "valid42": valid42, "valid43": valid43,
                "g32": g32, "g42": g42, "g43": g43,
                "g_valid": valid32 and valid42 and valid43,
            })

    slope_df = pd.DataFrame(slope_rows)
    slope_df.to_csv(OUT / "VLA_SLOPE_LIBRARY.csv", index=False)

    # =========================================================
    # 3. Range profile
    # =========================================================
    profile_rows = []
    local_rows = []
    gap_rows = []
    rank_rows = []

    for z_true in Z_TRUE_LIST:
        truth_row = slope_df[(slope_df["zs_m"] == z_true) & (np.abs(slope_df["r_km"] - 50.0) < 0.1)]
        if not len(truth_row) or not truth_row.iloc[0]["g_valid"]:
            rank_rows.append({"z_true_m": z_true, "status": "TRUTH_G_INVALID"})
            continue
        g_true = np.array([truth_row.iloc[0]["g32"], truth_row.iloc[0]["g42"], truth_row.iloc[0]["g43"]])

        scores = []
        for _, r in slope_df.iterrows():
            if not r["g_valid"]:
                scores.append({"zs": r["zs_m"], "r_km": r["r_km"], "J": np.nan})
                continue
            g = np.array([r["g32"], r["g42"], r["g43"]])
            J = float(np.sqrt(np.mean((g - g_true) ** 2)))
            scores.append({"zs": r["zs_m"], "r_km": r["r_km"], "J": J})

        scores_df = pd.DataFrame(scores)
        valid = scores_df[np.isfinite(scores_df["J"])]
        if len(valid) < 10:
            rank_rows.append({"z_true_m": z_true, "status": "INSUFFICIENT_VALID"})
            continue

        Js = valid["J"].values
        rs = valid["r_km"].values
        zs_arr = valid["zs"].values
        order = np.argsort(Js)
        tm = (np.abs(rs - 50.0) < 0.1) & (np.abs(zs_arr - z_true) < 0.1)
        if tm.any():
            ti = int(np.where(tm)[0][int(np.argmin(Js[tm]))])
            true_rank = int(np.where(order == ti)[0][0]) + 1
            true_J = float(Js[ti])
        else:
            true_rank, true_J = -1, np.nan
        rank_rows.append({"z_true_m": z_true, "true_rank": true_rank, "true_J": true_J, "n_valid": len(valid)})

        for rk in R_ALL_KM:
            m = np.abs(rs - rk) < 0.1
            if m.any():
                J_star = float(Js[m].min())
                z_star = float(zs_arr[m][int(np.argmin(Js[m]))])
            else:
                J_star = z_star = np.nan
            profile_rows.append({"z_true_m": z_true, "r_km": rk, "J_VLA_star": J_star, "z_star_m": z_star})
            if 48.0 <= rk <= 52.0:
                local_rows.append({"z_true_m": z_true, "r_km": rk, "J_VLA_star": J_star, "z_star_m": z_star})

        # global gap
        false_gaps = []
        for rk in R_ALL_KM:
            if abs(rk - 50.0) < 0.1: continue
            m = np.abs(rs - rk) < 0.1
            if m.any():
                J_star = float(Js[m].min())
                z_star = float(zs_arr[m][int(np.argmin(Js[m]))])
                false_gaps.append({"r_km": rk, "J": J_star, "z_star": z_star})
        false_gaps.sort(key=lambda x: x["J"])
        if false_gaps:
            gap_rows.append({
                "z_true_m": z_true,
                "nearest_false_r_km": false_gaps[0]["r_km"],
                "global_min_false_gap": false_gaps[0]["J"],
                "nearest_false_z_star": false_gaps[0]["z_star"],
            })

    prof_df = pd.DataFrame(profile_rows)
    prof_df.to_csv(OUT / "VLA_RANGE_PROFILE.csv", index=False)
    local_df = pd.DataFrame(local_rows)
    local_df.to_csv(OUT / "VLA_LOCAL_RANGE_PROFILE.csv", index=False)
    gap_df = pd.DataFrame(gap_rows)
    gap_df.to_csv(OUT / "VLA_GLOBAL_FALSE_RANGE_GAP.csv", index=False)
    rank_df = pd.DataFrame(rank_rows)

    # working range map
    wr_rows = []
    for rk in R_ALL_KM:
        sub = slope_df[np.abs(slope_df["r_km"] - rk) < 0.1]
        n_valid = int(sub["g_valid"].sum())
        wr_rows.append({"r_km": rk, "n_zs_valid": n_valid, "n_zs_total": len(sub),
                        "slope_availability": n_valid / max(len(sub), 1)})
    pd.DataFrame(wr_rows).to_csv(OUT / "VLA_SLOPE_WORKING_RANGE_MAP.csv", index=False)

    # source depth invariance
    inv_rows = []
    for rk in R_ALL_KM:
        sub = slope_df[(np.abs(slope_df["r_km"] - rk) < 0.1) & slope_df["g_valid"]]
        if len(sub) >= 3:
            inv_rows.append({
                "r_km": rk,
                "g32_mean": sub["g32"].mean(), "g32_std": sub["g32"].std(),
                "g42_mean": sub["g42"].mean(), "g42_std": sub["g42"].std(),
                "g43_mean": sub["g43"].mean(), "g43_std": sub["g43"].std(),
                "n_zs": len(sub),
            })
    pd.DataFrame(inv_rows).to_csv(OUT / "SOURCE_DEPTH_INVARIANCE.csv", index=False)

    # =========================================================
    # 4. Decision
    # =========================================================
    def range_unique(z_true):
        prof = [p for p in profile_rows if p["z_true_m"] == z_true and np.isfinite(p.get("J_VLA_star", np.nan))]
        if not prof: return False
        k = int(np.argmin([p["J_VLA_star"] for p in prof]))
        if abs(prof[k]["r_km"] - 50.0) > 0.1: return False
        for p in prof:
            if abs(p["r_km"] - 50.0) > 0.1 and p["J_VLA_star"] < 1e-6: return False
        return True

    all_unique = all(range_unique(z) for z in Z_TRUE_LIST)
    any_unique = any(range_unique(z) for z in Z_TRUE_LIST)
    all_rank1 = bool((rank_df["true_rank"] == 1).all()) if len(rank_df) else False

    if all_unique and all_rank1:
        decision = "E_STD_VLA_SECONDARY_SLOPE_RANGE_INFORMATION_CONFIRMED"
    elif any_unique:
        decision = "E_STD_VLA_SECONDARY_SLOPE_RANGE_INFORMATION_PARTIAL"
    else:
        decision = "E_STD_VLA_SLOPE_TRANSFER_NOT_ESTABLISHED"

    why = f"all_unique={all_unique}, any_unique={any_unique}, all_rank1={all_rank1}, conv={conv_pass}"
    dec = {
        "stage": "RANGE-UB-2B", "decision": decision, "why": why,
        "oracle_label": "ORACLE_PATH_CLASS_VLA_SLOPE_UPPER_BOUND",
        "not_HLA_capability": True,
        "created_utc": NOW,
    }
    (OUT / "RANGE_UB_2B_DECISION.json").write_text(json.dumps(dec, ensure_ascii=False, indent=2), encoding="utf-8")

    report = f"""# RANGE-UB-2B — E-STD VLA Secondary Delay Slope Upper Bound

UTC: {NOW}

基线：`e6e932e24a3a4345b2fee4a09bfb5256068ef29d`

## Observable

`XU2024_SECONDARY_DELAY_SLOPE_VECTOR`: k32, k42, k43 from BELLHOP cluster delays at VLA depths.

Oracle: `ORACLE_PATH_CLASS_VLA_SLOPE_UPPER_BOUND`（非 HLA 可实现）。

## Convergence (1601 vs 3201)

{'PASS' if conv_pass else 'FAIL'}

## 判定

### `{decision}`

{why}

## 排序诊断

{rank_df.to_string(index=False)}

## 全局假距离 gap

{gap_df.to_string(index=False) if len(gap_df) else 'N/A'}

## 局部距离 profile (48–52 km)

{local_df.to_string(index=False) if len(local_df) else 'N/A'}

## Working range map

{pd.DataFrame(wr_rows).to_string(index=False) if wr_rows else 'N/A'}

## 未做

噪声、SNR、自相关提取器、Radon、HLA 融合、RC2/RC3、P5。
"""
    (OUT / "RANGE_UB_2B_REPORT.md").write_text(report, encoding="utf-8")
    (OUT / "GPT_SYNC.md").write_text(f"# RANGE-UB-2B\n\n**{decision}**\n\n{why}\n", encoding="utf-8")

    print("conv_pass", conv_pass)
    print("decision", decision)
    print(why)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
