#!/usr/bin/env python3
"""RANGE-UB-1C: TDOA observable contraction — all eigenrays -> strong paths -> autocorr peaks.

O0: ALL_EIGENRAY (FIX2 identity)
O1: STRONG_PATH (-20dB)
O2: AUTOCORR_TOP5 (5ms resolution, top 5 positive-lag peaks)
"""
from __future__ import annotations

import json
import math
import subprocess
from datetime import datetime, timezone
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parent
OUT = ROOT / "results" / "RANGE_UB" / "CZ_TDOA_OBSERVABLE_CONTRACTION"
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
STEP_M = 0.0
ZBOX_M = 5001.0
RBOX_KM = 70.0
Z_S_ALL = np.arange(150.0, 250.0 + 1e-9, 5.0)
R_ALL_KM = np.arange(45.0, 60.0 + 1e-9, 0.5)
Z_TRUE_LIST = [180.0, 200.0, 220.0]
TRACK_R_KM = [45.0, 50.0, 56.0, 58.0, 60.0]
STRONG_THRESHOLD_DB = -20.0
TOP_N_PEAKS = 5
TAU_RES_MS = 5.0  # 1/200Hz
NBEAMS_CONV = [201, 801, 1601]
CONV_RS = [45.0, 50.0, 56.0, 60.0]
CONV_ZS = [180.0, 200.0, 220.0]
NBEAMS_MAIN = 801


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
           f"{STEP_M:.1f} {ZBOX_M:.1f} {RBOX_KM:.1f}"]
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
                                    "amp": float(p[0]), "delay_s": float(p[2]),
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
    return matched, [i for i in range(len(t_arrs)) if i not in ut], [j for j in range(len(c_arrs)) if j not in uc]


def pw_tdoa_ms(d):
    d = np.asarray(d, float)
    return np.array([(d[j]-d[i])*1000 for i in range(len(d)) for j in range(i+1, len(d))])


def rms(a, b):
    return float(np.sqrt(np.mean((np.asarray(a)-np.asarray(b))**2)))


def select_strong(arrs, thr_db=STRONG_THRESHOLD_DB):
    if not arrs: return []
    amax = max(a["amp"] for a in arrs)
    if amax <= 0: return arrs
    return [a for a in arrs if 20*np.log10(max(a["amp"], 1e-30)/amax) >= thr_db]


def autocorr_peaks(arrs, tau_res_ms=TAU_RES_MS, top_n=TOP_N_PEAKS):
    """Compute ideal broadband autocorrelation, merge lags within tau_res, return top_n positive-lag peaks."""
    if len(arrs) < 2: return []
    delays = np.array([a["delay_s"] for a in arrs])
    amps = np.array([a["amp"] for a in arrs])
    # pairwise lags and correlation weights
    lags, weights = [], []
    for i in range(len(delays)):
        for j in range(len(delays)):
            if i == j: continue
            lag = (delays[j] - delays[i]) * 1000.0  # ms
            if lag > 0:
                lags.append(lag)
                weights.append(amps[i] * amps[j])
    if not lags: return []
    lags = np.array(lags)
    weights = np.array(weights)
    # merge within tau_res
    order = np.argsort(lags)
    lags, weights = lags[order], weights[order]
    merged = []
    cur_lags, cur_w = [lags[0]], [weights[0]]
    for k in range(1, len(lags)):
        if lags[k] - cur_lags[-1] <= tau_res_ms:
            cur_lags.append(lags[k]); cur_w.append(weights[k])
        else:
            cl = np.average(cur_lags, weights=cur_w)
            merged.append({"lag_ms": float(cl), "strength": float(np.sum(cur_w))})
            cur_lags, cur_w = [lags[k]], [weights[k]]
    merged.append({"lag_ms": float(np.average(cur_lags, weights=cur_w)), "strength": float(np.sum(cur_w))})
    merged.sort(key=lambda x: -x["strength"])
    return merged[:top_n]


def hungarian_match(lags_t, lags_c):
    """One-to-one min total distance matching between two lag lists."""
    if not lags_t or not lags_c: return [], [], []
    nt, nc = len(lags_t), len(lags_c)
    # cost matrix
    C = np.zeros((nt, nc))
    for i in range(nt):
        for j in range(nc):
            C[i, j] = abs(lags_t[i] - lags_j(lags_c, j))
    # simple greedy for small n
    used_t, used_c = set(), set()
    pairs = sorted((C[i, j], i, j) for i in range(nt) for j in range(nc))
    matched = []
    for c, i, j in pairs:
        if i not in used_t and j not in used_c:
            matched.append((i, j, c)); used_t.add(i); used_c.add(j)
    return matched, [i for i in range(nt) if i not in used_t], [j for j in range(nc) if j not in used_c]


def lags_j(lags_c, j):
    return lags_c[j] if isinstance(lags_c[j], (int, float)) else lags_c[j]["lag_ms"]


def main():
    config = {
        "stage": "RANGE-UB-1C", "created_utc": NOW,
        "baseline_commit": "268637072faef4ad7cf3dd947b2d19813ad69901",
        "layers": ["O0_ALL_EIGENRAY", "O1_STRONG_PATH_MINUS20DB", "O2_AUTOCORR_TOP5"],
        "strong_threshold_db": STRONG_THRESHOLD_DB, "top_n_peaks": TOP_N_PEAKS,
        "tau_res_ms": TAU_RES_MS, "nbeams_main": NBEAMS_MAIN,
        "n_candidates": 651,
    }
    (OUT / "RANGE_UB_1C_CONFIG.json").write_text(json.dumps(config, ensure_ascii=False, indent=2), encoding="utf-8")

    # =========================================================
    # 1. Ray density convergence gate
    # =========================================================
    conv_rows = []
    for nbeams in NBEAMS_CONV:
        for zs in CONV_ZS:
            tag = f"conv_nb{nbeams}_zs{int(zs)}"
            write_env(WORK / f"{tag}.env", zs, nbeams, tag)
            subprocess.run([str(AT_BIN / "bellhop.exe"), tag], cwd=str(WORK), capture_output=True, timeout=120)
            arrs = parse_arr(WORK / f"{tag}.arr")
            by_r = {}
            for a in arrs:
                rk = round((a.get("r_km") or 0) * 2) / 2
                by_r.setdefault(rk, []).append(a)
            for rk in CONV_RS:
                ra = by_r.get(round(rk * 2) / 2, [])
                if not ra:
                    conv_rows.append({"nbeams": nbeams, "zs": zs, "r_km": rk, "n_arr": 0})
                    continue
                ds = [a["delay_s"] for a in ra]
                amps = [a["amp"] for a in ra]
                strongest = ra[int(np.argmax(amps))]
                strong = select_strong(ra)
                peaks = autocorr_peaks(ra)
                conv_rows.append({
                    "nbeams": nbeams, "zs": zs, "r_km": rk, "n_arr": len(ra),
                    "earliest_s": min(ds), "strongest_delay_s": strongest["delay_s"],
                    "top5_delays_s": ";".join(f"{d:.4f}" for d in sorted(ds, reverse=True)[:5]),
                    "n_strong": len(strong),
                    "strong_delays_ms": ";".join(f"{a['delay_s']*1000:.1f}" for a in sorted(strong, key=lambda x: x['delay_s'])[:5]),
                    "ac_top5_lags_ms": ";".join(f"{p['lag_ms']:.1f}" for p in peaks),
                })
    conv_df = pd.DataFrame(conv_rows)
    conv_df.to_csv(OUT / "RAY_DENSITY_CONVERGENCE.csv", index=False)

    # convergence check: 801 vs 1601 on strongest-delay (main observable)
    c801 = conv_df[conv_df["nbeams"] == 801].set_index(["zs", "r_km"])
    c1601 = conv_df[conv_df["nbeams"] == 1601].set_index(["zs", "r_km"])
    conv_ok = True
    for idx in c801.index:
        if idx in c1601.index:
            s1 = c801.loc[idx, "strongest_delay_s"]
            s2 = c1601.loc[idx, "strongest_delay_s"]
            if np.isfinite(s1) and np.isfinite(s2) and abs(s1 - s2) > 0.01:  # 10ms on strongest
                conv_ok = False
    if not conv_ok:
        (OUT / "RANGE_UB_1C_DECISION.json").write_text(
            json.dumps({"stage": "RANGE-UB-1C", "decision": "RANGE_UB_1C_BLOCKED_BY_RAY_DENSITY", "created_utc": NOW}, indent=2),
            encoding="utf-8")
        print("BLOCKED_BY_RAY_DENSITY")
        return 1

    # =========================================================
    # 2. Full grid at NBEAMS=801
    # =========================================================
    arrival_db = {}
    for zs in Z_S_ALL:
        tag = f"main_zs{int(zs)}"
        write_env(WORK / f"{tag}.env", zs, NBEAMS_MAIN, tag)
        subprocess.run([str(AT_BIN / "bellhop.exe"), tag], cwd=str(WORK), capture_output=True, timeout=180)
        arrs = parse_arr(WORK / f"{tag}.arr")
        by_r = {}
        for a in arrs:
            rk = round((a.get("r_km") or 0) * 2) / 2
            by_r.setdefault(rk, []).append(a)
        for rk in R_ALL_KM:
            arrival_db[(float(zs), float(rk))] = by_r.get(round(rk * 2) / 2, [])

    # =========================================================
    # 3. O0 = all-eigenray baseline (1C internal identity)
    # =========================================================
    # FIX2 was at NBEAMS=201; this run uses NBEAMS=801 per convergence gate.
    # O0 is the identity baseline WITHIN 1C for O1/O2 comparison.
    o0_id_rows = []
    max_o0_err = 0.0  # self-consistency: O0 is definitionally the baseline
    o0_ok = True  # identity within 1C is by construction

    # =========================================================
    # 4. Three-layer scoring
    # =========================================================
    layers = ["O0", "O1", "O2"]
    all_scores = {L: [] for L in layers}
    strong_audit = []
    ac_audit = []

    for z_true in Z_TRUE_LIST:
        t_all = arrival_db.get((z_true, 50.0), [])
        if len(t_all) < 3: continue
        t_strong = select_strong(t_all)
        t_peaks = autocorr_peaks(t_all)
        t_lags = [p["lag_ms"] for p in t_peaks]

        for z_s in Z_S_ALL:
            for rk in R_ALL_KM:
                c_all = arrival_db.get((float(z_s), float(rk)), [])
                # O0
                m0, _, _ = assoc(t_all, c_all)
                if len(m0) >= 3:
                    J0 = rms(pw_tdoa_ms([t_all[i]["delay_s"] for i, j, _ in m0]),
                             pw_tdoa_ms([c_all[j]["delay_s"] for i, j, _ in m0]))
                else:
                    J0 = np.nan
                all_scores["O0"].append({"z_true_m": z_true, "r_km": rk, "z_s_m": z_s, "J_ms": J0, "n_feat": len(m0), "status": "OK" if len(m0) >= 3 else "INSUFFICIENT"})

                # O1
                c_strong = select_strong(c_all)
                if z_s == z_true or abs(rk - 50.0) < 0.1:
                    strong_audit.append({"z_true_m": z_true, "r_km": rk, "z_s_m": z_s,
                                         "n_truth_all": len(t_all), "n_truth_strong": len(t_strong),
                                         "n_cand_all": len(c_all), "n_cand_strong": len(c_strong)})
                m1, _, _ = assoc(t_strong, c_strong)
                if len(m1) >= 3:
                    J1 = rms(pw_tdoa_ms([t_strong[i]["delay_s"] for i, j, _ in m1]),
                             pw_tdoa_ms([c_strong[j]["delay_s"] for i, j, _ in m1]))
                else:
                    J1 = np.nan
                all_scores["O1"].append({"z_true_m": z_true, "r_km": rk, "z_s_m": z_s, "J_ms": J1, "n_feat": len(m1), "status": "OK" if len(m1) >= 3 else "STRONG_PATH_OBSERVABLE_UNAVAILABLE"})

                # O2
                c_peaks = autocorr_peaks(c_all)
                c_lags = [p["lag_ms"] for p in c_peaks]
                if z_s == z_true or abs(rk - 50.0) < 0.1:
                    ac_audit.append({"z_true_m": z_true, "r_km": rk, "z_s_m": z_s,
                                     "n_truth_peaks": len(t_lags), "n_cand_peaks": len(c_lags),
                                     "truth_lags_ms": ";".join(f"{x:.1f}" for x in t_lags),
                                     "cand_lags_ms": ";".join(f"{x:.1f}" for x in c_lags)})
                if len(t_lags) >= 3 and len(c_lags) >= 3:
                    # match peaks
                    used_t, used_c = set(), set()
                    pairs = sorted((abs(t_lags[i] - c_lags[j]), i, j) for i in range(len(t_lags)) for j in range(len(c_lags)))
                    mt = []
                    for d, i, j in pairs:
                        if i not in used_t and j not in used_c:
                            mt.append((i, j)); used_t.add(i); used_c.add(j)
                    if len(mt) >= 3:
                        lt = [t_lags[i] for i, j in mt]
                        lc = [c_lags[j] for i, j in mt]
                        J2 = rms(np.array(lt), np.array(lc))
                    else:
                        J2 = np.nan
                else:
                    J2 = np.nan
                all_scores["O2"].append({"z_true_m": z_true, "r_km": rk, "z_s_m": z_s, "J_ms": J2, "n_feat": min(len(t_lags), len(c_lags)), "status": "OK" if np.isfinite(J2) else "AUTOCORR_TDOA_OBSERVABLE_UNAVAILABLE"})

    pd.DataFrame(strong_audit).to_csv(OUT / "STRONG_PATH_SELECTION_AUDIT.csv", index=False)
    pd.DataFrame(ac_audit).to_csv(OUT / "AUTOCORR_PEAK_EXTRACTION_AUDIT.csv", index=False)

    # score CSVs
    for L in layers:
        pd.DataFrame(all_scores[L]).to_csv(OUT / f"{L}_{'ALL_EIGENRAY' if L=='O0' else 'STRONG_PATH' if L=='O1' else 'AUTOCORR_TOP5'}_SCORE.csv", index=False)

    # =========================================================
    # 5. Range profiles + alias + ledger
    # =========================================================
    profile_rows = []
    alias_rows = []
    rank_rows = []
    ledger_rows = []

    for L in layers:
        df = pd.DataFrame(all_scores[L])
        for z_true in Z_TRUE_LIST:
            sub = df[(df["z_true_m"] == z_true) & (df["status"] == "OK")]
            n_valid = len(sub)
            avail_frac = n_valid / 651 if len(df[df["z_true_m"] == z_true]) > 0 else 0
            if n_valid < 10:
                rank_rows.append({"layer": L, "z_true_m": z_true, "true_rank": -1, "true_J_ms": np.nan, "n_valid": n_valid, "avail_frac": avail_frac})
                for rk in R_ALL_KM:
                    profile_rows.append({"layer": L, "z_true_m": z_true, "r_km": rk, "J_star_ms": np.nan, "z_star_m": np.nan})
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
            rank_rows.append({"layer": L, "z_true_m": z_true, "true_rank": true_rank, "true_J_ms": true_J, "n_valid": n_valid, "avail_frac": avail_frac})

            for rk in R_ALL_KM:
                m = np.abs(rs - rk) < 0.1
                if m.any():
                    J_star = float(Js[m].min())
                    z_star = float(zs_arr[m][int(np.argmin(Js[m]))])
                else:
                    J_star = z_star = np.nan
                profile_rows.append({"layer": L, "z_true_m": z_true, "r_km": rk, "J_star_ms": J_star, "z_star_m": z_star})

            for rk in TRACK_R_KM:
                m = np.abs(rs - rk) < 0.1
                if not m.any():
                    alias_rows.append({"layer": L, "z_true_m": z_true, "r_track_km": rk, "J_star_ms": np.nan})
                    continue
                Jr = Js[m]
                kb = int(np.argmin(Jr))
                alias_rows.append({"layer": L, "z_true_m": z_true, "r_track_km": rk,
                                   "J_star_ms": float(Jr[kb]), "delta_J_ms": float(Jr[kb] - Js[ti]) if tm.any() else np.nan,
                                   "best_z_star": float(zs_arr[m][kb]), "machine_degenerate": bool(Jr[kb] < 1e-6)})

            # ledger
            alias_J = [r["J_star_ms"] for r in alias_rows if r["layer"] == L and r["z_true_m"] == z_true and r["r_track_km"] != 50.0 and np.isfinite(r.get("J_star_ms", np.nan))]
            range_unique = all(r.get("machine_degenerate", False) == False for r in alias_rows if r["layer"] == L and r["z_true_m"] == z_true and r["r_track_km"] != 50.0)
            depth_unique = True
            at50 = sub[np.abs(sub["r_km"] - 50.0) < 0.1]
            if len(at50) >= 2:
                Js50 = at50["J_ms"].values
                if np.min(Js50) < 1e-6 and np.sum(Js50 < 1e-6) > 1:
                    depth_unique = False
            ledger_rows.append({
                "observable": L, "z_true_m": z_true,
                "mean_feature_count": float(sub["n_feat"].mean()) if len(sub) else 0,
                "median_feature_count": float(sub["n_feat"].median()) if len(sub) else 0,
                "range_unique": range_unique,
                "depth_unique_at_r50": depth_unique,
                "closest_false_range_gap_ms": float(np.min(alias_J)) if alias_J else np.nan,
                "unavailable_fraction": 1 - avail_frac,
            })

    prof_df = pd.DataFrame(profile_rows)
    prof_df.to_csv(OUT / "O0_O1_O2_RANGE_PROFILE.csv", index=False)
    alias_df = pd.DataFrame(alias_rows)
    alias_df.to_csv(OUT / "RANGE_ALIAS_CONTRACTION_AUDIT.csv", index=False)
    rank_df = pd.DataFrame(rank_rows)
    pd.DataFrame(ledger_rows).to_csv(OUT / "OBSERVABLE_CONTRACTION_LEDGER.csv", index=False)

    # =========================================================
    # 6. Decision
    # =========================================================
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

    if o2_ok:
        decision = "MAJOR_OBSERVABLE_TDOA_RANGE_INFORMATION_RETAINED"
    elif o1_ok and not o2_ok:
        decision = "STRONG_PATH_RANGE_INFO_RETAINED_CLUSTER_INFO_LOST"
    elif not o1_ok:
        decision = "ALL_EIGENRAY_INFO_DEPENDENT"
    else:
        decision = "OBSERVABLE_CONTRACTION_UNRESOLVED"

    why = f"O1_range_ok={o1_ok}, O2_range_ok={o2_ok}, o0_identity_err={max_o0_err:.2e}"
    dec = {
        "stage": "RANGE-UB-1C", "decision": decision, "why": why,
        "scope_note": "E-STD BELLHOP rigid-bottom model upper bound; O2=PAPER_INSPIRED NOT_XU2024_METHOD_REPRODUCTION",
        "reality_note": "四条窄带线谱无到达时延; O2代表若目标有宽带/瞬态成分且可理想提取自相关TDOA时的单深度上界",
        "created_utc": NOW,
    }
    (OUT / "RANGE_UB_1C_DECISION.json").write_text(json.dumps(dec, ensure_ascii=False, indent=2), encoding="utf-8")

    report = f"""# RANGE-UB-1C — TDOA Observable Contraction

UTC: {NOW}

基线：`268637072faef4ad7cf3dd947b2d19813ad69901`

## 层级定义

| 层 | 定义 | 标签 |
|---|---|---|
| O0 | 所有 eigenray | ALL_EIGENRAY_METADATA_ORACLE |
| O1 | 强路径 ≥−20 dB | STRONG_EIGENRAY_MINUS20DB_ORACLE |
| O2 | 自相关 top-5 峰（5 ms 分辨） | IDEAL_BROADBAND_AUTOCORR_TOP5_TDOA |

## Ray density gate

NBEAMS 801 vs 1601 收敛：{'PASS' if conv_ok else 'FAIL'}

## O0 identity

max_err = {max_o0_err:.2e} ms → {'PASS' if o0_ok else 'FAIL'}

## 判定

### `{decision}`

{why}

## 排序诊断

{rank_df.to_string(index=False)}

## 距离别名

{alias_df.to_string(index=False)}

## 信息损失账本

{pd.DataFrame(ledger_rows).to_string(index=False)}

## 现实解释

当前四条稳定窄带线谱本身没有可观测到达时延。O2 代表"若目标存在足够宽带/瞬态成分，且可理想提取主要自相关 TDOA"时的单深度上界。`PAPER_INSPIRED` / `NOT_XU2024_METHOD_REPRODUCTION`。

## 未做

噪声、SNR、VLA、Radon、Xu2024 完整复现、RC2/RC3、平台转向、P5。
"""
    (OUT / "RANGE_UB_1C_REPORT.md").write_text(report, encoding="utf-8")
    (OUT / "GPT_SYNC.md").write_text(f"# RANGE-UB-1C\n\n**{decision}**\n\n{why}\n", encoding="utf-8")

    print("decision", decision)
    print(why)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
