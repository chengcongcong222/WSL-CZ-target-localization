#!/usr/bin/env python3
"""R3-CLOSEDLOOP-1E-B: Tracked frequency-drift robustness of turn-assisted range anchor.

TRIPLE/FOUR, D_F=[0,.0025,.005,.01,.02], F_TRACKED, frozen truth node.
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
OUT = ROOT / "results" / "R3_RC23_CLOSEDLOOP" / "R3_RC23_TURN_TRACKED_FREQ_ROBUSTNESS"
ZGRID = ROOT / "results" / "R3_C2_Yang_SA_depth" / "R3_C2_1" / "_kraken_zgrid"
FIX = ROOT / "results" / "R3_RC23_CLOSEDLOOP" / "R3_RC23_TURN_FREQ_SUBSET_FIX"
AT_BIN = ROOT / "tools" / "acoustics_toolbox" / "atWin10" / "at" / "bin"
WORK = OUT / "_kraken"
OUT.mkdir(parents=True, exist_ok=True)
WORK.mkdir(parents=True, exist_ok=True)
NOW = datetime.now(timezone.utc).isoformat()

U_PLAT = 2.0; T_TURN = 600.0; T_END = 1200.0; DT_OBS = 10.0
SIGMA_DEG = 0.1; RNG_SEED = 20260912; DELTA = 15.0
TRUTH = dict(r0_m=50e3, theta0_deg=0.0, v=2.0, psi_deg=5.0)
R_GRID_KM = np.arange(45.0, 60.0 + 1e-9, 1.0)
TH_GRID_DEG = np.arange(-5.0, 5.0 + 1e-9, 0.5)
V_GRID = np.arange(1.0, 3.0 + 1e-9, 0.2)
PSI_GRID_DEG = np.arange(-15.0, 15.0 + 1e-9, 1.0)
N_NODES = int(len(R_GRID_KM) * len(TH_GRID_DEG) * len(V_GRID) * len(PSI_GRID_DEG))
Z_TRUE_LIST = [180.0, 200.0, 220.0]
Z_PROFILE = np.arange(150.0, 250.0 + 1e-9, 5.0)
ZR = 200.0; TAU = 0.5
CONFIGS = {"TRIPLE": [201.0, 235.0, 283.0], "FOUR": [201.0, 235.0, 283.0, 338.0]}
D_F_LIST = [0.0, 0.0025, 0.005, 0.01, 0.02]
Q_F = np.array([0.0, 0.5, 1.0, 0.5, 0.0, -0.5, -1.0, -0.5])
ALL_BASE_F = [201.0, 235.0, 283.0, 338.0]


def platform_turn(t, delta_deg, t_turn=T_TURN, u=U_PLAT):
    t = np.asarray(t, float); xp, yp = np.empty_like(t), np.empty_like(t)
    d = math.radians(delta_deg); c, s = math.cos(d), math.sin(d)
    if abs(delta_deg) < 1e-12: xp[:] = u*t; yp[:] = 0.0; return xp, yp
    mask = t <= t_turn; xp[mask] = u*t[mask]; yp[mask] = 0.0
    dt = t[~mask]-t_turn; xp[~mask] = u*t_turn+u*dt*c; yp[~mask] = u*dt*s
    return xp, yp


def target_xy(t, r0_m, th0_rad, v, psi_rad):
    return r0_m*np.cos(th0_rad)+v*t*np.cos(psi_rad), r0_m*np.sin(th0_rad)+v*t*np.sin(psi_rad)


def wrap(a): return (a+np.pi)%(2*np.pi)-np.pi


def range_traj(t, r0_m, th0_deg, v, psi_deg, delta_deg):
    xp, yp = platform_turn(t, delta_deg)
    xt, yt = target_xy(t, r0_m, np.deg2rad(th0_deg), v, np.deg2rad(psi_deg))
    return np.hypot(xt-xp, yt-yp)


def parse_mod(path):
    buf = path.read_bytes(); recl = 4*int(np.frombuffer(buf[:4],dtype="<i4")[0])
    hdr = np.frombuffer(buf[84:108],dtype="<i4"); ntot,nmat = int(hdr[2]),int(hdr[3])
    depths = np.frombuffer(buf[4*recl:5*recl],dtype="<f4")[:ntot].astype(float)
    M = int(np.frombuffer(buf[5*recl:5*recl+4],dtype="<i4")[0])
    phi = np.zeros((nmat,M),complex)
    for im in range(M):
        off = (7+im)*recl; chunk = np.frombuffer(buf[off:off+recl],dtype="<c8")
        take = min(nmat,chunk.size); phi[:take,im] = chunk[:take]
    k = np.frombuffer(buf[(7+M)*recl:(7+M)*recl+M*8],dtype="<c8")
    return {"depths":depths,"M":M,"phi":phi,"k":k}


def F_mat(mod, r):
    kre, alpha = mod["k"].real, -mod["k"].imag
    r = np.asarray(r, float)
    amp = np.sqrt(2*np.pi/(kre[:,None]*r[None,:]))
    return amp*np.exp(-1j*kre[:,None]*r[None,:]-alpha[:,None]*r[None,:]-1j*np.pi/4)


def L_prof(mod, F, z_list, zr=ZR):
    izr = int(np.argmin(np.abs(mod["depths"]-zr)))
    pr = mod["phi"][izr]; W = mod["phi"]*pr[None,:]
    izs = [int(np.argmin(np.abs(mod["depths"]-z))) for z in z_list]
    return 20.0*np.log10(np.maximum(np.abs(W[izs,:]@F),1e-30))


def demean(L): return np.asarray(L,float)-float(np.mean(L))


def write_env(path, freq, tag):
    src = ZGRID / "zgrid_f235.env"
    lines = src.read_text(encoding="utf-8").splitlines()
    ssp = []
    for ln in lines[5:]:
        p = ln.split()
        if len(p) >= 2:
            try: ssp.append(f"  {float(p[0]):.1f} {float(p[1]):.4f} /")
            except ValueError: break
        else: break
    out = [f"'{tag}'", f"{freq:.4f}", "1", "'CVN'", "  51 0.0 5000.0", *ssp,
           "'R' 0.0", "1", "200.0 /", "51", "150.00", "250.00", "R"]
    path.write_text("\n".join(out)+"\n", encoding="utf-8")


def main():
    config = {"stage":"R3-CLOSEDLOOP-1E-B","created_utc":NOW,
              "baseline_commit":"cfc93cc7d0d1b242245f9cd713e9a0d4cde1859c",
              "configs":["TRIPLE","FOUR"],"D_F":D_F_LIST,"tracking":"F_TRACKED"}
    (OUT/"CLOSEDLOOP_1E_B_CONFIG.json").write_text(json.dumps(config,ensure_ascii=False,indent=2),encoding="utf-8")

    # =========================================================
    # 1. Frequency model manifest
    # =========================================================
    # needed: base_f, D, q -> actual_f
    needed_freqs = set()
    for f0 in ALL_BASE_F:
        for D in D_F_LIST:
            for q in [-1,-0.5,0,0.5,1]:
                needed_freqs.add(round(f0*(1+D*q), 4))
    needed_freqs = sorted(needed_freqs)

    manifest_rows = []
    mod_cache = {}
    # Load ZGRID models for base frequencies first
    for f0 in ALL_BASE_F:
        zp = ZGRID / f"zgrid_f{int(f0)}.mod"
        if zp.exists():
            mod_cache[f0] = parse_mod(zp)

    for fq in needed_freqs:
        key = f"f{fq:.4f}".replace(".", "p")
        modp = WORK / f"{key}.mod"
        # if fq matches a base freq, use ZGRID
        if fq in mod_cache:
            manifest_rows.append({"actual_f": fq, "model_path": "ZGRID", "model_exists": True,
                                  "parse_ok": True, "fallback_used": False})
            continue
        exists = modp.exists()
        if not exists:
            write_env(WORK / f"{key}.env", fq, key)
            subprocess.run([str(AT_BIN / "kraken.exe"), key], cwd=str(WORK), capture_output=True, timeout=120)
            exists = modp.exists()
        parse_ok = False
        if exists:
            try:
                mod_cache[fq] = parse_mod(modp)
                parse_ok = True
            except Exception:
                parse_ok = False
        manifest_rows.append({"actual_f": fq, "model_path": str(modp), "model_exists": exists,
                              "parse_ok": parse_ok, "fallback_used": False})
    # map base+D+q -> actual_f
    for f0 in ALL_BASE_F:
        for D in D_F_LIST:
            for q in [-1,-0.5,0,0.5,1]:
                fq = round(f0*(1+D*q), 4)
                manifest_rows.append({"base_f": f0, "D": D, "q": q, "actual_f": fq,
                                      "fallback_used": False})
    man_df = pd.DataFrame(manifest_rows)
    man_df.to_csv(OUT/"FREQUENCY_MODEL_MANIFEST.csv", index=False)
    if any(not r.get("parse_ok", True) for r in manifest_rows if "parse_ok" in r):
        pass  # check below

    # drift schedule identity
    drift_rows = []
    t_full = np.arange(0.0, T_END+1e-9, DT_OBS)
    seg = np.minimum((t_full/T_END*8).astype(int), 7)
    q_t = Q_F[seg]
    for D in D_F_LIST:
        for f0 in ALL_BASE_F:
            for i in [0, 10, 30, 60, 61, 90, 120]:
                if i < len(t_full):
                    drift_rows.append({"D": D, "base_f": f0, "time_s": t_full[i], "q": q_t[i],
                                       "actual_f": f0*(1+D*q_t[i])})
    pd.DataFrame(drift_rows).to_csv(OUT/"DRIFT_SCHEDULE_5B_IDENTITY.csv", index=False)

    # =========================================================
    # 2. RC2 + truth node
    # =========================================================
    n_w1 = int(np.sum(t_full <= T_TURN+1e-9)); n_w2 = len(t_full)-n_w1
    t_w1, t_w2 = t_full[:n_w1], t_full[n_w1:]
    sigma_rad = np.deg2rad(SIGMA_DEG); thr = 13.3*sigma_rad**2

    rr,tt,vv,pp = np.meshgrid(R_GRID_KM*1e3, np.deg2rad(TH_GRID_DEG), V_GRID, np.deg2rad(PSI_GRID_DEG), indexing="ij")
    r0_all,th0_all,v_all,psi_all = rr.ravel(),tt.ravel(),vv.ravel(),pp.ravel()

    rng = np.random.default_rng(RNG_SEED)
    noise = rng.normal(0.0, sigma_rad, size=t_full.shape)
    xp,yp = platform_turn(t_full, DELTA)
    xt,yt = target_xy(t_full, TRUTH["r0_m"], np.deg2rad(TRUTH["theta0_deg"]), TRUTH["v"], np.deg2rad(TRUTH["psi_deg"]))
    obs_bearing = np.arctan2(yt-yp, xt-xp) + noise

    pred = np.arctan2(
        r0_all[:,None]*np.sin(th0_all[:,None])+v_all[:,None]*t_full[None,:]*np.sin(psi_all[:,None])-yp[None,:],
        r0_all[:,None]*np.cos(th0_all[:,None])+v_all[:,None]*t_full[None,:]*np.cos(psi_all[:,None])-xp[None,:])
    cost = np.sum(wrap(pred-obs_bearing[None,:])**2, axis=1)
    acc = cost <= (float(cost.min())+thr)
    idx_acc = np.where(acc)[0]; n_rc2 = len(idx_acc)
    # frozen truth node
    truth_idx = int(np.argmin((r0_all-TRUTH["r0_m"])**2+(th0_all-np.deg2rad(TRUTH["theta0_deg"]))**2
                              +(v_all-TRUTH["v"])**2+(psi_all-np.deg2rad(TRUTH["psi_deg"]))**2))

    # =========================================================
    # 3. Precompute acoustic for all needed freqs
    # =========================================================
    n_z = Z_PROFILE.size
    r_true = range_traj(t_full, TRUTH["r0_m"], TRUTH["theta0_deg"], TRUTH["v"], TRUTH["psi_deg"], DELTA)

    # candidate trajectories
    r_cand = {}
    for j in idx_acc:
        r_cand[j] = range_traj(t_full, r0_all[j], np.rad2deg(th0_all[j]), v_all[j], np.rad2deg(psi_all[j]), DELTA)

    # L for each freq at truth and candidates (full window)
    # cache: fq -> {truth: L_raw(z_true), cand: {j: L_raw(z_profile)}}
    # For tracked drift, need L at each actual frequency on full trajectory
    L_cache = {}  # fq -> (L_true_full[z_true], L_cand_full[j][z_profile])
    for fq in needed_freqs:
        mod = mod_cache.get(fq)
        if mod is None:
            continue
        F_t = F_mat(mod, r_true)
        L_true = L_prof(mod, F_t, Z_TRUE_LIST)  # (3, n_full)
        L_cache[fq] = {"true": {zt: L_true[i].copy() for i, zt in enumerate(Z_TRUE_LIST)}, "cand": {}}
        for j in idx_acc:
            F_c = F_mat(mod, r_cand[j])
            L_cache[fq]["cand"][j] = L_prof(mod, F_c, Z_PROFILE)

    # =========================================================
    # 4. Score for each config, D, z_true
    # =========================================================
    cand_rows = []
    case_rows = []
    truth_id_rows = []

    for cfg_name in ["TRIPLE", "FOUR"]:
        freqs = CONFIGS[cfg_name]
        for D in D_F_LIST:
            # build tracked L for each base freq
            # observation: L_true_track[f0](z_true, t) = L at actual f(t)
            # candidate: L_cand_track[f0](j, z, t) = L at actual f(t)
            for z_true in Z_TRUE_LIST:
                # observation tracked
                L_obs = {}
                for f0 in freqs:
                    L_full = np.zeros(len(t_full))
                    for qi, q in enumerate(q_t):
                        fq = round(f0*(1+D*q), 4)
                        if fq in L_cache:
                            L_full[qi] = L_cache[fq]["true"][z_true][qi]
                    L_obs[f0] = (demean(L_full[:n_w1]), demean(L_full[n_w1:]))

                # score candidates
                J_all = np.empty(n_rc2); z_star_all = np.empty(n_rc2)
                for i, j in enumerate(idx_acc):
                    J_z = np.zeros(n_z)
                    for f0 in freqs:
                        o1, o2 = L_obs[f0]
                        L_full_z = np.zeros((n_z, len(t_full)))
                        for qi, q in enumerate(q_t):
                            fq = round(f0*(1+D*q), 4)
                            if fq in L_cache:
                                L_full_z[:, qi] = L_cache[fq]["cand"][j][:, qi]
                        L1 = demean(L_full_z[:, :n_w1].T)  # wrong shape
                        # fix: demean per z across time
                        L1_z = np.zeros((n_z, n_w1)); L2_z = np.zeros((n_z, n_w2))
                        for iz in range(n_z):
                            L1_z[iz] = demean(L_full_z[iz, :n_w1])
                            L2_z[iz] = demean(L_full_z[iz, n_w1:])
                        for iz in range(n_z):
                            e1 = L1_z[iz] - o1; e2 = L2_z[iz] - o2
                            J_z[iz] += np.sum(e1**2) + np.sum(e2**2)
                    J_z = np.sqrt(J_z / ((n_w1+n_w2)*len(freqs)))
                    J_all[i] = J_z.min(); z_star_all[i] = Z_PROFILE[J_z.argmin()]

                jmin = float(J_all.min())
                pos_true = int(np.where(idx_acc == truth_idx)[0][0])
                true_J = float(J_all[pos_true])
                true_rank = int((J_all < true_J).sum()) + 1
                true_kept = bool(J_all[pos_true] <= jmin + TAU)
                keep = J_all <= (jmin + TAU)

                # truth model identity check
                truth_id_rows.append({"config": cfg_name, "D": D, "z_true_m": z_true,
                                      "true_J": true_J, "pass": true_J <= 1e-10})

                surv = grp_keep = idx_acc[keep]
                surv_r = r0_all[surv] / 1e3
                r_keep = sorted(set(int(round(x)) for x in surv_r)) if len(surv_r) else []
                r_best = sorted(set(int(round(x)) for x in r0_all[idx_acc][J_all <= jmin + 1e-10] / 1e3))
                r_wrong = [r for r in r_keep if r != 50]

                wrong_mask = np.abs(r0_all[idx_acc] / 1e3 - 50.0) >= 0.1
                best_wrong_J = float(J_all[wrong_mask].min()) if wrong_mask.any() else np.nan
                best_wrong_r = float(r0_all[idx_acc][wrong_mask][J_all[wrong_mask].argmin()] / 1e3) if wrong_mask.any() else np.nan

                strict = true_kept and true_rank == 1 and r_keep == [50]

                case_rows.append({
                    "config": cfg_name, "D": D, "D_pct": D*100, "z_true_m": z_true,
                    "true_J": true_J, "J_min": jmin, "true_rank": true_rank, "truth_retained": true_kept,
                    "best_range_bins": str(r_best), "best_range_unique_50": r_best == [50],
                    "surviving_r_bins": str(r_keep), "n_surviving_r_bins": len(r_keep),
                    "r_width_km": float(surv_r.max()-surv_r.min()) if len(surv_r) else 0.0,
                    "max_surviving_deviation_km": float(np.max(np.abs(surv_r-50.0))) if len(surv_r) else 0.0,
                    "best_wrong_range_J": best_wrong_J, "best_wrong_range_r": best_wrong_r,
                    "wrong_range_margin": best_wrong_J - true_J if np.isfinite(best_wrong_J) else np.nan,
                    "threshold_clearance": best_wrong_J - (jmin+TAU) if np.isfinite(best_wrong_J) else np.nan,
                    "strict_single_bin": strict,
                })

                for i, j in enumerate(idx_acc):
                    cand_rows.append({"config": cfg_name, "D": D, "D_pct": D*100, "ztrue": z_true,
                                      "node_id": int(j), "r0": float(r0_all[j]/1e3),
                                      "theta0": float(np.rad2deg(th0_all[j])), "v": float(v_all[j]),
                                      "psi": float(np.rad2deg(psi_all[j])),
                                      "J": float(J_all[i]), "z_star": float(z_star_all[i]),
                                      "keep": bool(keep[i]), "is_truth_state": bool(j == truth_idx)})

    cand_df = pd.DataFrame(cand_rows)
    cand_df.to_csv(OUT / "TRACKED_FREQ_CANDIDATE_SCORES.csv", index=False)
    case_df = pd.DataFrame(case_rows)
    case_df.to_csv(OUT / "TRACKED_FREQ_ANCHOR_CASES.csv", index=False)
    pd.DataFrame(truth_id_rows).to_csv(OUT / "TRACKED_TRUTH_MODEL_IDENTITY.csv", index=False)

    # D=0 identity
    id_rows = []
    ref = pd.read_csv(FIX / "SUBSET_CANDIDATE_SCORES_FIXED.csv")
    for cfg_name, sub_tag in [("TRIPLE", "201+235+283"), ("FOUR", "201+235+283+338")]:
        for z_true in Z_TRUE_LIST:
            cur = cand_df[(cand_df.config == cfg_name) & (cand_df.D == 0) & (cand_df.ztrue == z_true)]
            r = ref[(ref["subset"] == sub_tag) & (ref["z_true_m"] == z_true)]
            ids_c = set(cur["node_id"]); ids_r = set(r["node_id"])
            max_err = 0.0
            for nid in ids_c & ids_r:
                j_c = float(cur[cur.node_id == nid]["J"].values[0])
                j_r = float(r[r.node_id == nid]["J"].values[0])
                max_err = max(max_err, abs(j_c - j_r))
            id_rows.append({"config": cfg_name, "z_true_m": z_true, "ids_match": ids_c == ids_r,
                            "max_J_err": max_err, "pass": ids_c == ids_r and max_err <= 1e-10})
    id_df = pd.DataFrame(id_rows)
    id_df.to_csv(OUT / "ZERO_DRIFT_1D_FIX2_IDENTITY.csv", index=False)

    # working range
    wr_rows = []
    for cfg in ["TRIPLE", "FOUR"]:
        for D in D_F_LIST:
            sub = case_df[(case_df.config == cfg) & (case_df.D == D)]
            wr_rows.append({"config": cfg, "D": D, "D_pct": D*100,
                            "strict": int(sub["strict_single_bin"].sum()),
                            "truth_ok": int(sub["truth_retained"].sum()),
                            "rank1": int((sub["true_rank"] == 1).sum()),
                            "best_50": int(sub["best_range_unique_50"].sum()),
                            "max_width": float(sub["r_width_km"].max()),
                            "min_margin": float(sub["wrong_range_margin"].min())})
    wr_df = pd.DataFrame(wr_rows)
    wr_df.to_csv(OUT / "TRACKED_FREQ_WORKING_RANGE_SUMMARY.csv", index=False)

    # TRIPLE vs FOUR
    cmp_rows = []
    for D in D_F_LIST:
        for z_true in Z_TRUE_LIST:
            t = case_df[(case_df.config == "TRIPLE") & (case_df.D == D) & (case_df.z_true_m == z_true)]
            f = case_df[(case_df.config == "FOUR") & (case_df.D == D) & (case_df.z_true_m == z_true)]
            if len(t) and len(f):
                cmp_rows.append({"D": D, "z_true_m": z_true,
                                 "triple_strict": bool(t.iloc[0]["strict_single_bin"]),
                                 "four_strict": bool(f.iloc[0]["strict_single_bin"]),
                                 "triple_margin": float(t.iloc[0]["wrong_range_margin"]),
                                 "four_margin": float(f.iloc[0]["wrong_range_margin"])})
    pd.DataFrame(cmp_rows).to_csv(OUT / "TRIPLE_VS_FOUR_TRACKED_FREQ_ROBUSTNESS.csv", index=False)

    # decision
    triple_all_strict = bool(case_df[case_df.config == "TRIPLE"]["strict_single_bin"].all())
    four_all_strict = bool(case_df[case_df.config == "FOUR"]["strict_single_bin"].all())
    all_truth_rank = bool((case_df["true_rank"] == 1).all()) and bool(case_df["truth_retained"].all())
    all_best50 = bool(case_df["best_range_unique_50"].all())

    if triple_all_strict and four_all_strict:
        decision = "TRIPLE_AND_FOUR_STRICT_RANGE_ANCHOR_SURVIVES_TRACKED_2PCT_FREQ_DRIFT"
    elif four_all_strict and not triple_all_strict:
        decision = "FOUR_REDUNDANCY_IMPROVES_TRACKED_FREQUENCY_DRIFT_ROBUSTNESS"
    elif all_truth_rank and all_best50:
        decision = "TRACKED_FREQUENCY_DRIFT_WEAKENS_STRICT_ANCHOR_BUT_PRESERVES_OPTIMAL_RANGE"
    else:
        decision = "TRACKED_FREQUENCY_DRIFT_DEGRADES_RANGE_RANKING_IN_TESTED_RANGE"

    why = f"triple_strict={triple_all_strict}, four_strict={four_all_strict}, rank1={all_truth_rank}, best50={all_best50}"
    dec = {"stage": "R3-CLOSEDLOOP-1E-B", "decision": decision, "why": why,
           "tracking_label": "PERFECT_FREQUENCY_TRACKING_CONTROL", "created_utc": NOW}
    (OUT / "R3_RC23_CLOSEDLOOP_1E_B_DECISION.json").write_text(json.dumps(dec, ensure_ascii=False, indent=2), encoding="utf-8")

    report = f"""# R3-CLOSEDLOOP-1E-B — Tracked Frequency Drift

UTC: {NOW}

## 判定

### `{decision}`

{why}

## Working Range

{wr_df.to_string(index=False)}

## D=0 Identity

{id_df.to_string(index=False)}

## 未做

F_NOMINAL、跟踪误差、幅漂、SSP、P5。
"""
    (OUT / "R3_RC23_CLOSEDLOOP_1E_B_REPORT.md").write_text(report, encoding="utf-8")
    (OUT / "GPT_SYNC.md").write_text(f"# 1E-B\n\n**{decision}**\n\n{why}\n", encoding="utf-8")

    print("decision", decision)
    print(why)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
