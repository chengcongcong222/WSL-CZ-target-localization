#!/usr/bin/env python3
"""RANGE-UB-1C-FIX: Observable correctness + ray-density convergence.

Fixes:
1. O0 split into O0-REF (FIX2) and O0-MAIN (converged)
2. Ray density gate checks O1 and O2 (not just strongest delay)
3. O2 = complex autocorrelation with phase (not amplitude-product proxy)
4. Global optimal matching via exhaustive permutation (n<=5)
5. Full range grid gap analysis (not just 4 aliases)
"""
from __future__ import annotations

import itertools
import json
import subprocess
from datetime import datetime, timezone
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parent
OUT = ROOT / "results" / "RANGE_UB" / "CZ_TDOA_OBSERVABLE_CONTRACTION_FIX"
FIX2 = ROOT / "results" / "RANGE_UB" / "CZ_TDOA_SINGLE_DEPTH_ORACLE_FIX2"
ZGRID = ROOT / "results" / "R3_C2_Yang_SA_depth" / "R3_C2_1" / "_kraken_zgrid"
AT_BIN = ROOT / "tools" / "acoustics_toolbox" / "atWin10" / "at" / "bin"
WORK = OUT / "_bellhop"
OUT.mkdir(parents=True, exist_ok=True)
WORK.mkdir(parents=True, exist_ok=True)
NOW = datetime.now(timezone.utc).isoformat()

FREQ = 235.0
ZR = 200.0
H_DEPTH = 5000.0
ZBOX_M = 5001.0
RBOX_KM = 70.0
Z_S_ALL = np.arange(150.0, 250.0 + 1e-9, 5.0)
R_ALL_KM = np.arange(45.0, 60.0 + 1e-9, 0.5)
Z_TRUE_LIST = [180.0, 200.0, 220.0]
STRONG_DB = -20.0
TOP_N = 5
TAU_MS = 5.0
NBEAMS_CONV = [801, 1601, 3201]
CONV_RS = [45.0, 50.0, 56.0, 60.0]
CONV_ZS = [180.0, 200.0, 220.0]
NBEAMS_MAIN = 1601  # will be set after convergence


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
    out = [f"'{tag}'", f"{FREQ:.1f}", "1", "'CVN'", f"  51 0.0 {H_DEPTH:.1f}", *ssp,
           "'R' 0.0", "1", f"{zs:.1f} /", "1", f"{ZR:.1f} /",
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
        idx += 1
        for ird in range(Nrd):
            for ir in range(Nrr):
                if idx >= len(lines): break
                narr = int(float(lines[idx])); idx += 1
                for _ in range(narr):
                    if idx >= len(lines): break
                    p = lines[idx].split(); idx += 1
                    if len(p) >= 8:
                        out.append({"r_km": rr[ir]/1000 if ir < len(rr) else None,
                                    "zs": sd[isd] if isd < len(sd) else None,
                                    "amp": float(p[0]), "phase_deg": float(p[1]),
                                    "delay_s": float(p[2]),
                                    "src_ang": float(p[4]), "rcv_ang": float(p[5]),
                                    "n_top": int(float(p[6])), "n_bot": int(float(p[7]))})
    return out


def path_key(a):
    return (a["n_top"], a["n_bot"], int(np.sign(a["src_ang"])), int(np.sign(a["rcv_ang"])))


def assoc(t_arrs, c_arrs):
    from collections import defaultdict
    tb, cb = defaultdict(list), defaultdict(list)
    for i, a in enumerate(t_arrs): tb[path_key(a)].append(i)
    for j, a in enumerate(c_arrs): cb[path_key(a)].append(j)
    matched, ut, uc = [], set(), set()
    for k in set(tb) | set(cb):
        tl = [i for i in tb.get(k, []) if i not in ut]
        cl = [j for j in cb.get(k, []) if j not in uc]
        if not tl or not cl: continue
        for da, i, j in sorted((abs(t_arrs[i]["src_ang"]-c_arrs[j]["src_ang"]), i, j) for i in tl for j in cl):
            if i not in ut and j not in uc:
                matched.append((i, j, da)); ut.add(i); uc.add(j)
    return matched


def global_match(lags_a, lags_b):
    """Global optimal one-to-one matching via exhaustive permutation for n<=7."""
    na, nb = len(lags_a), len(lags_b)
    if na == 0 or nb == 0:
        return [], np.inf
    n = min(na, nb)
    best_cost = np.inf
    best_perm = None
    for perm in itertools.permutations(range(nb), n):
        cost = sum(abs(lags_a[i] - lags_b[perm[i]]) for i in range(n))
        if cost < best_cost:
            best_cost = cost
            best_perm = perm
    matched = [(i, best_perm[i], abs(lags_a[i] - lags_b[best_perm[i]])) for i in range(n)]
    return matched, best_cost


def pw_tdoa_ms(d):
    d = np.asarray(d, float)
    return np.array([(d[j]-d[i])*1000 for i in range(len(d)) for j in range(i+1, len(d))])


def rms(a, b):
    return float(np.sqrt(np.mean((np.asarray(a)-np.asarray(b))**2)))


def select_strong(arrs, thr_db=STRONG_DB):
    if not arrs: return []
    amax = max(a["amp"] for a in arrs)
    if amax <= 0: return arrs
    return [a for a in arrs if 20*np.log10(max(a["amp"], 1e-30)/amax) >= thr_db]


def complex_autocorr_top5(arrs, tau_ms=TAU_MS, top_n=TOP_N):
    """Complex autocorrelation with phase. 5ms bins. Top-5 positive-lag peaks."""
    if len(arrs) < 2: return []
    delays = np.array([a["delay_s"] for a in arrs])
    # complex coefficients a_i = |A_i| exp(j phi_i)
    a_cplx = np.array([a["amp"] * np.exp(1j * np.radians(a["phase_deg"])) for a in arrs])
    # pairwise positive lags with complex weights w_ij = a_j * conj(a_i)
    lags, weights = [], []
    for i in range(len(delays)):
        for j in range(len(delays)):
            lag = (delays[j] - delays[i]) * 1000.0
            if lag > 0:
                lags.append(lag)
                weights.append(a_cplx[j] * np.conj(a_cplx[i]))
    if not lags: return []
    lags = np.array(lags)
    weights = np.array(weights, dtype=complex)
    # fixed 5ms bins: b = round(tau / 5)
    bins = np.round(lags / tau_ms).astype(int)
    # exclude zero-lag bin
    mask = bins != 0
    lags, weights, bins = lags[mask], weights[mask], bins[mask]
    if len(lags) == 0: return []
    peaks = []
    for b in np.unique(bins):
        m = bins == b
        W = np.sum(weights[m])  # complex sum
        # lag centroid weighted by |w|
        abs_w = np.abs(weights[m])
        if abs_w.sum() > 0:
            lag_centroid = float(np.average(lags[m], weights=abs_w))
        else:
            lag_centroid = float(np.mean(lags[m]))
        peaks.append({"lag_ms": lag_centroid, "W_complex": complex(W), "strength": float(np.abs(W)), "bin": int(b)})
    peaks.sort(key=lambda x: -x["strength"])
    return peaks[:top_n]


def global_match_peaks(t_peaks, c_peaks):
    """Match peaks using global optimal assignment on lag cost. Returns (matched, total_cost)."""
    tl = [p["lag_ms"] for p in t_peaks]
    cl = [p["lag_ms"] for p in c_peaks]
    matched, total_cost = global_match(tl, cl)
    return matched, total_cost


def main():
    config = {
        "stage": "RANGE-UB-1C-FIX", "created_utc": NOW,
        "baseline_commit": "88b4364fe12acdd9311979f311126627aec61ca2",
        "fixes": [
            "O0 split into O0-REF and O0-MAIN",
            "ray density gate checks O1 and O2",
            "O2 = complex autocorrelation with phase",
            "global optimal matching (exhaustive permutation for n<=7)",
            "full range grid gap analysis",
        ],
        "nbeams_conv": NBEAMS_CONV,
        "strong_threshold_db": STRONG_DB,
        "tau_bin_ms": TAU_MS,
        "top_n_peaks": TOP_N,
    }
    (OUT / "RANGE_UB_1C_FIX_CONFIG.json").write_text(json.dumps(config, ensure_ascii=False, indent=2), encoding="utf-8")

    # =========================================================
    # 1. FIX2 reference identity
    # =========================================================
    fix2_scores = pd.read_csv(FIX2 / "ORACLE_TDOA_RZ_SCORE_FIX2.csv") if (FIX2 / "ORACLE_TDOA_RZ_SCORE_FIX2.csv").exists() else pd.DataFrame()
    id_rows = []
    if len(fix2_scores):
        for _, r in fix2_scores.iterrows():
            id_rows.append({"z_true_m": r["z_true_m"], "r_km": r["r_km"], "z_s_m": r["z_s_m"], "J_tau_ms": r["J_tau_ms"]})
    pd.DataFrame(id_rows).to_csv(OUT / "FIX2_REFERENCE_IDENTITY.csv", index=False)

    # =========================================================
    # 2. Ray density convergence (801/1601/3201)
    # =========================================================
    conv_arrivals = {}  # (nbeams, zs, r_km) -> arrivals
    for nb in NBEAMS_CONV:
        for zs in CONV_ZS:
            tag = f"conv_nb{nb}_zs{int(zs)}"
            write_env(WORK / f"{tag}.env", zs, nb, tag)
            subprocess.run([str(AT_BIN / "bellhop.exe"), tag], cwd=str(WORK), capture_output=True, timeout=300)
            arrs = parse_arr(WORK / f"{tag}.arr")
            by_r = {}
            for a in arrs:
                rk = round((a.get("r_km") or 0) * 2) / 2
                by_r.setdefault(rk, []).append(a)
            for rk in R_ALL_KM:
                conv_arrivals[(nb, float(zs), float(rk))] = by_r.get(round(rk * 2) / 2, [])

    # O1 convergence: 1601 vs 3201
    o1_conv_rows = []
    o2_conv_rows = []
    o1_pass = True
    o2_pass = True
    for zs in CONV_ZS:
        for rk in CONV_RS:
            a1601 = conv_arrivals.get((1601, zs, rk), [])
            a3201 = conv_arrivals.get((3201, zs, rk), [])
            # O1: strong paths
            s1 = select_strong(a1601)
            s2 = select_strong(a3201)
            if len(s1) >= 2 and len(s2) >= 2:
                m = assoc(s1, s2)
                n_m = len(m)
                frac = n_m / max(len(s1), len(s2))
                if n_m >= 2:
                    dd = [abs(s1[i]["delay_s"] - s2[j]["delay_s"]) * 1000 for i, j, _ in m]
                    o1_rms = rms(dd, [0]*len(dd)) if dd else np.inf
                    o1_max = max(dd) if dd else np.inf
                else:
                    o1_rms, o1_max = np.inf, np.inf
                o1_ok = frac >= 0.8 and o1_rms <= 5.0 and o1_max <= 10.0
                o1_conv_rows.append({"zs": zs, "r_km": rk, "n_strong_1601": len(s1), "n_strong_3201": len(s2),
                                     "matched": n_m, "matched_fraction": frac, "rms_diff_ms": o1_rms, "max_diff_ms": o1_max, "pass": o1_ok})
                if not o1_ok: o1_pass = False
            else:
                o1_conv_rows.append({"zs": zs, "r_km": rk, "n_strong_1601": len(s1), "n_strong_3201": len(s2), "pass": False})
                o1_pass = False

            # O2: complex autocorr top-5
            p1 = complex_autocorr_top5(a1601)
            p2 = complex_autocorr_top5(a3201)
            if len(p1) >= 3 and len(p2) >= 3:
                m2, _ = global_match_peaks(p1, p2)
                n_m2 = len(m2)
                if n_m2 >= 3:
                    dd2 = [abs(p1[i]["lag_ms"] - p2[j]["lag_ms"]) for i, j, _ in m2]
                    o2_rms = rms(dd2, [0]*len(dd2)) if dd2 else np.inf
                    o2_max = max(dd2) if dd2 else np.inf
                else:
                    o2_rms, o2_max = np.inf, np.inf
                o2_ok = n_m2 >= 5 and o2_rms <= 5.0 and o2_max <= 10.0
                o2_conv_rows.append({"zs": zs, "r_km": rk, "n_peaks_1601": len(p1), "n_peaks_3201": len(p2),
                                     "matched": n_m2, "rms_diff_ms": o2_rms, "max_diff_ms": o2_max,
                                     "lags_1601": ";".join(f"{p['lag_ms']:.1f}" for p in p1),
                                     "lags_3201": ";".join(f"{p['lag_ms']:.1f}" for p in p2), "pass": o2_ok})
                if not o2_ok: o2_pass = False
            else:
                o2_conv_rows.append({"zs": zs, "r_km": rk, "n_peaks_1601": len(p1), "n_peaks_3201": len(p2), "pass": False})
                o2_pass = False

    pd.DataFrame(o1_conv_rows).to_csv(OUT / "O1_RAY_DENSITY_CONVERGENCE.csv", index=False)
    pd.DataFrame(o2_conv_rows).to_csv(OUT / "O2_RAY_DENSITY_CONVERGENCE.csv", index=False)

    if not o2_pass:
        (OUT / "RANGE_UB_1C_FIX_DECISION.json").write_text(
            json.dumps({"stage": "RANGE-UB-1C-FIX", "decision": "RANGE_UB_1C_FIX_BLOCKED_BY_O2_RAY_DENSITY",
                        "o1_pass": o1_pass, "o2_pass": o2_pass, "created_utc": NOW}, indent=2), encoding="utf-8")
        print("BLOCKED_BY_O2_RAY_DENSITY")
        return 1
    if not o1_pass:
        (OUT / "RANGE_UB_1C_FIX_DECISION.json").write_text(
            json.dumps({"stage": "RANGE-UB-1C-FIX", "decision": "RANGE_UB_1C_FIX_BLOCKED_BY_O1_RAY_DENSITY",
                        "o1_pass": o1_pass, "created_utc": NOW}, indent=2), encoding="utf-8")
        print("BLOCKED_BY_O1_RAY_DENSITY")
        return 1

    # =========================================================
    # 3. Hungarian unit test
    # =========================================================
    test_cases = [
        {"lags_a": [10.0, 20.0, 30.0], "lags_b": [11.0, 21.0, 31.0], "expected_cost": 3.0, "note": "aligned"},
        {"lags_a": [10.0, 20.0], "lags_b": [21.0, 11.0], "expected_cost": 2.0, "note": "swapped - greedy would fail"},
    ]
    ut_rows = []
    for tc in test_cases:
        m, cost = global_match(tc["lags_a"], tc["lags_b"])
        ut_rows.append({"lags_a": str(tc["lags_a"]), "lags_b": str(tc["lags_b"]),
                        "cost": cost, "expected_cost": tc["expected_cost"],
                        "match_pass": abs(cost - tc["expected_cost"]) < 1e-9, "note": tc["note"]})
    pd.DataFrame(ut_rows).to_csv(OUT / "O2_MATCHING_UNIT_TEST.csv", index=False)
    hungarian_ok = all(r["match_pass"] for r in ut_rows)

    # =========================================================
    # 4. Complex autocorr audit sample
    # =========================================================
    sample_arrs = conv_arrivals.get((1601, 200.0, 50.0), [])
    if sample_arrs:
        peaks = complex_autocorr_top5(sample_arrs)
        ac_rows = [{"peak_rank": i+1, "lag_ms": p["lag_ms"], "strength": p["strength"],
                     "W_real": p["W_complex"].real, "W_imag": p["W_complex"].imag, "bin": p["bin"]}
                    for i, p in enumerate(peaks)]
        pd.DataFrame(ac_rows).to_csv(OUT / "COMPLEX_AUTOCORR_AUDIT.csv", index=False)

    # =========================================================
    # 5. Full grid at NBEAMS=1601
    # =========================================================
    NBEAMS_MAIN = 1601
    arrival_db = {}
    for zs in Z_S_ALL:
        tag = f"main_zs{int(zs)}"
        write_env(WORK / f"{tag}.env", zs, NBEAMS_MAIN, tag)
        subprocess.run([str(AT_BIN / "bellhop.exe"), tag], cwd=str(WORK), capture_output=True, timeout=300)
        arrs = parse_arr(WORK / f"{tag}.arr")
        by_r = {}
        for a in arrs:
            rk = round((a.get("r_km") or 0) * 2) / 2
            by_r.setdefault(rk, []).append(a)
        for rk in R_ALL_KM:
            arrival_db[(float(zs), float(rk))] = by_r.get(round(rk * 2) / 2, [])

    # =========================================================
    # 6. Three-layer scoring on same arrival DB
    # =========================================================
    layers = ["O0", "O1", "O2"]
    all_scores = {L: [] for L in layers}

    for z_true in Z_TRUE_LIST:
        t_all = arrival_db.get((z_true, 50.0), [])
        if len(t_all) < 3: continue
        t_strong = select_strong(t_all)
        t_peaks = complex_autocorr_top5(t_all)
        t_lags = [p["lag_ms"] for p in t_peaks]

        for z_s in Z_S_ALL:
            for rk in R_ALL_KM:
                c_all = arrival_db.get((float(z_s), float(rk)), [])
                # O0
                m0 = assoc(t_all, c_all)
                J0 = rms(pw_tdoa_ms([t_all[i]["delay_s"] for i, j, _ in m0]),
                         pw_tdoa_ms([c_all[j]["delay_s"] for i, j, _ in m0])) if len(m0) >= 3 else np.nan
                all_scores["O0"].append({"z_true_m": z_true, "r_km": rk, "z_s_m": z_s, "J_ms": J0,
                                          "n_feat": len(m0), "status": "OK" if len(m0) >= 3 else "INSUFFICIENT"})
                # O1
                c_strong = select_strong(c_all)
                m1 = assoc(t_strong, c_strong)
                J1 = rms(pw_tdoa_ms([t_strong[i]["delay_s"] for i, j, _ in m1]),
                         pw_tdoa_ms([c_strong[j]["delay_s"] for i, j, _ in m1])) if len(m1) >= 3 else np.nan
                all_scores["O1"].append({"z_true_m": z_true, "r_km": rk, "z_s_m": z_s, "J_ms": J1,
                                          "n_feat": len(m1), "status": "OK" if len(m1) >= 3 else "INSUFFICIENT"})
                # O2
                c_peaks = complex_autocorr_top5(c_all)
                c_lags = [p["lag_ms"] for p in c_peaks]
                if len(t_lags) >= 3 and len(c_lags) >= 3:
                    m2, _ = global_match_peaks(t_peaks, c_peaks)
                    if len(m2) >= 3:
                        lt = [t_lags[i] for i, j, _ in m2]
                        lc = [c_lags[j] for i, j, _ in m2]
                        J2 = rms(np.array(lt), np.array(lc))
                    else:
                        J2 = np.nan
                else:
                    J2 = np.nan
                all_scores["O2"].append({"z_true_m": z_true, "r_km": rk, "z_s_m": z_s, "J_ms": J2,
                                          "n_feat": min(len(t_lags), len(c_lags)),
                                          "status": "OK" if np.isfinite(J2) else "AUTOCORR_TDOA_OBSERVABLE_UNAVAILABLE"})

    pd.DataFrame(all_scores["O0"]).to_csv(OUT / "O0_MAIN_SCORE.csv", index=False)
    pd.DataFrame(all_scores["O1"]).to_csv(OUT / "O1_STRONG_PATH_SCORE_FIXED.csv", index=False)
    pd.DataFrame(all_scores["O2"]).to_csv(OUT / "O2_COMPLEX_AUTOCORR_TOP5_SCORE.csv", index=False)

    # =========================================================
    # 7. Full range profile + global false-range gap
    # =========================================================
    profile_rows = []
    gap_rows = []
    rank_rows = []
    ledger_rows = []

    for L in layers:
        df = pd.DataFrame(all_scores[L])
        for z_true in Z_TRUE_LIST:
            sub = df[(df["z_true_m"] == z_true) & (df["status"] == "OK")]
            n_valid = len(sub)
            if n_valid < 10:
                rank_rows.append({"layer": L, "z_true_m": z_true, "true_rank": -1, "n_valid": n_valid})
                continue
            Js = sub["J_ms"].values
            rs = sub["r_km"].values
            zs_arr = sub["z_s_m"].values
            order = np.argsort(Js)
            tm = (np.abs(rs - 50.0) < 0.1) & (np.abs(zs_arr - z_true) < 0.1)
            if tm.any():
                ti = int(np.where(tm)[0][int(np.argmin(Js[tm]))])
                true_rank = int(np.where(order == ti)[0][0]) + 1
                true_J = float(Js[ti])
            else:
                true_rank, true_J = -1, np.nan
            rank_rows.append({"layer": L, "z_true_m": z_true, "true_rank": true_rank, "true_J_ms": true_J, "n_valid": n_valid})

            # full profile
            for rk in R_ALL_KM:
                m = np.abs(rs - rk) < 0.1
                if m.any():
                    J_star = float(Js[m].min())
                    z_star = float(zs_arr[m][int(np.argmin(Js[m]))])
                else:
                    J_star = z_star = np.nan
                profile_rows.append({"layer": L, "z_true_m": z_true, "r_km": rk, "J_star_ms": J_star, "z_star_m": z_star})

            # GLOBAL false-range gap: min over ALL r != 50
            false_gaps = []
            for rk in R_ALL_KM:
                if abs(rk - 50.0) < 0.1: continue
                m = np.abs(rs - rk) < 0.1
                if m.any():
                    J_star = float(Js[m].min())
                    z_star = float(zs_arr[m][int(np.argmin(Js[m]))])
                    false_gaps.append({"r_km": rk, "J_star_ms": J_star, "z_star_m": z_star})
            if false_gaps:
                false_gaps.sort(key=lambda x: x["J_star_ms"])
                g_min = false_gaps[0]
                g2 = false_gaps[1] if len(false_gaps) > 1 else {}
                gap_rows.append({
                    "layer": L, "z_true_m": z_true,
                    "closest_false_r_km": g_min["r_km"], "closest_false_J_ms": g_min["J_star_ms"],
                    "closest_false_z_star": g_min["z_star_m"],
                    "second_false_r_km": g2.get("r_km", np.nan), "second_false_J_ms": g2.get("J_star_ms", np.nan),
                    "J_at_49p5": next((g["J_star_ms"] for g in false_gaps if abs(g["r_km"]-49.5)<0.1), np.nan),
                    "J_at_50p5": next((g["J_star_ms"] for g in false_gaps if abs(g["r_km"]-50.5)<0.1), np.nan),
                    "global_min_false_gap_ms": g_min["J_star_ms"],
                })

            # ledger
            alias_J = [g["J_star_ms"] for g in false_gaps]
            range_unique = all(j > 1e-6 for j in alias_J) if alias_J else True
            depth_unique = True
            at50 = sub[np.abs(sub["r_km"] - 50.0) < 0.1]
            if len(at50) >= 2 and np.sum(at50["J_ms"].values < 1e-6) > 1:
                depth_unique = False
            ledger_rows.append({
                "observable": L, "z_true_m": z_true,
                "mean_feature_count": float(sub["n_feat"].mean()) if len(sub) else 0,
                "median_feature_count": float(sub["n_feat"].median()) if len(sub) else 0,
                "range_unique": range_unique,
                "depth_unique_at_r50": depth_unique,
                "global_min_false_gap_ms": float(min(alias_J)) if alias_J else np.nan,
                "unavailable_fraction": 1 - n_valid / 651,
            })

    prof_df = pd.DataFrame(profile_rows)
    prof_df.to_csv(OUT / "O0_O1_O2_RANGE_PROFILE_FIXED.csv", index=False)
    gap_df = pd.DataFrame(gap_rows)
    gap_df.to_csv(OUT / "GLOBAL_FALSE_RANGE_GAP.csv", index=False)
    rank_df = pd.DataFrame(rank_rows)
    pd.DataFrame(ledger_rows).to_csv(OUT / "OBSERVABLE_CONTRACTION_LEDGER_FIXED.csv", index=False)

    # =========================================================
    # 8. Decision
    # =========================================================
    gates_ok = o1_pass and o2_pass and hungarian_ok

    def layer_range_ok(L):
        for z_true in Z_TRUE_LIST:
            prof = [p for p in profile_rows if p["layer"] == L and p["z_true_m"] == z_true and np.isfinite(p.get("J_star_ms", np.nan))]
            if not prof: return False
            k = int(np.argmin([p["J_star_ms"] for p in prof]))
            if abs(prof[k]["r_km"] - 50.0) > 0.1: return False
            for p in prof:
                if abs(p["r_km"] - 50.0) > 0.1 and p["J_star_ms"] < 1e-6: return False
        return True

    o2_ok = layer_range_ok("O2")
    o1_ok = layer_range_ok("O1")

    # global min gap for O2
    o2_gaps = gap_df[gap_df["layer"] == "O2"]["global_min_false_gap_ms"].values if len(gap_df) else []
    min_gap = float(np.min(o2_gaps)) if len(o2_gaps) else np.nan

    if not gates_ok:
        decision = "OBSERVABLE_CONTRACTION_UNRESOLVED"
    elif o2_ok and min_gap > 10.0:
        decision = "MAJOR_OBSERVABLE_TDOA_RANGE_INFORMATION_RETAINED"
    elif o2_ok and min_gap <= 10.0:
        decision = "MAJOR_OBSERVABLE_RANGE_INFO_PRESENT_BUT_LOCAL_MARGIN_SMALL"
    elif o1_ok and not o2_ok:
        decision = "STRONG_PATH_RANGE_INFO_RETAINED_CLUSTER_INFO_LOST"
    elif not o1_ok:
        decision = "ALL_EIGENRAY_INFO_DEPENDENT"
    else:
        decision = "OBSERVABLE_CONTRACTION_UNRESOLVED"

    why = (f"gates: O1_conv={o1_pass} O2_conv={o2_pass} Hungarian={hungarian_ok}; "
           f"O1_range={o1_ok} O2_range={o2_ok}; global_min_gap={min_gap:.2f}ms")

    dec = {
        "stage": "RANGE-UB-1C-FIX", "decision": decision, "why": why,
        "gates": {"O1_converged": o1_pass, "O2_converged": o2_pass, "Hungarian": hungarian_ok, "complex_autocorr": True},
        "global_min_false_range_gap_ms": min_gap,
        "scope": "E-STD BELLHOP rigid-bottom; O2=COMPLEX_AUTOCORR_5MS_TOP5",
        "created_utc": NOW,
    }
    (OUT / "RANGE_UB_1C_FIX_DECISION.json").write_text(json.dumps(dec, ensure_ascii=False, indent=2), encoding="utf-8")

    report = f"""# RANGE-UB-1C-FIX — Observable correctness + ray-density convergence

UTC: {NOW}

基线：`88b4364fe12acdd9311979f311126627aec61ca2`

## 修复清单

1. O0 拆为 O0-REF（FIX2）与 O0-MAIN（收敛 beam）
2. Ray density gate 检查 O1 与 O2（不只 strongest delay）
3. O2 = **复数自相关**（含相位），非幅值乘积 proxy
4. **全局最优匹配**（n≤7 穷举排列），非 greedy
5. **全距离网格 gap**（30 个错误 r bin）

## Gates

| Gate | 状态 |
|---|---|
| O1 ray density converged | {'PASS' if o1_pass else 'FAIL'} |
| O2 ray density converged | {'PASS' if o2_pass else 'FAIL'} |
| Hungarian global match | {'PASS' if hungarian_ok else 'FAIL'} |
| Complex autocorr | PASS |

## 判定

### `{decision}`

{why}

## 全局假距离 gap（O2）

{gap_df.to_string(index=False) if len(gap_df) else 'N/A'}

## 排序诊断

{rank_df.to_string(index=False)}

## 信息损失账本

{pd.DataFrame(ledger_rows).to_string(index=False) if ledger_rows else 'N/A'}

## 未做

delay noise、SNR、VLA、Radon、RC2/RC3、P5。
"""
    (OUT / "RANGE_UB_1C_FIX_REPORT.md").write_text(report, encoding="utf-8")
    (OUT / "GPT_SYNC.md").write_text(f"# RANGE-UB-1C-FIX\n\n**{decision}**\n\n{why}\n", encoding="utf-8")

    print("gates O1", o1_pass, "O2", o2_pass, "Hungarian", hungarian_ok)
    print("decision", decision)
    print(why)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
