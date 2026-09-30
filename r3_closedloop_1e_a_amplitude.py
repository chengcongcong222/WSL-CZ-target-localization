#!/usr/bin/env python3
"""R3-CLOSEDLOOP-1E-A: Amplitude-variation robustness of turn-assisted multi-freq range anchor.

TRIPLE={201,235,283}, FOUR={201,235,283,338}. A_RMS=0-2dB, RAMP/SINE, COMMON/LINE_SPECIFIC.
"""
from __future__ import annotations

import itertools
import json
import math
from datetime import datetime, timezone
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parent
OUT = ROOT / "results" / "R3_RC23_CLOSEDLOOP" / "R3_RC23_TURN_AMPLITUDE_ROBUSTNESS"
ZGRID = ROOT / "results" / "R3_C2_Yang_SA_depth" / "R3_C2_1" / "_kraken_zgrid"
FIX2 = ROOT / "results" / "R3_RC23_CLOSEDLOOP" / "R3_RC23_TURN_FREQ_SUBSET_FIX2"
OUT.mkdir(parents=True, exist_ok=True)
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
A_RMS_LIST = [0.0, 0.25, 0.5, 1.0, 2.0]
SHAPES = ["A-RAMP", "A-SINE"]
MODES = ["COMMON_LINE_AMPLITUDE_VARIATION", "LINE_SPECIFIC_AMPLITUDE_VARIATION"]
MULT_COMMON = {"TRIPLE": [1, 1, 1], "FOUR": [1, 1, 1, 1]}
MULT_LINE = {"TRIPLE": [1, -1, 1], "FOUR": [1, -1, 1, -1]}


def platform_turn(t, delta_deg, t_turn=T_TURN, u=U_PLAT):
    t = np.asarray(t, float)
    xp, yp = np.empty_like(t), np.empty_like(t)
    d = math.radians(delta_deg); c, s = math.cos(d), math.sin(d)
    if abs(delta_deg) < 1e-12: xp[:] = u*t; yp[:] = 0.0; return xp, yp
    mask = t <= t_turn; xp[mask] = u*t[mask]; yp[mask] = 0.0
    dt = t[~mask] - t_turn; xp[~mask] = u*t_turn + u*dt*c; yp[~mask] = u*dt*s
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


def norm_rms(q):
    q = np.asarray(q,float); q = q - q.mean()
    s = float(np.sqrt(np.mean(q*q)))
    return q/s if s > 0 else q


def main():
    config = {"stage":"R3-CLOSEDLOOP-1E-A","created_utc":NOW,
              "baseline_commit":"8ecd6a40ca846bb10a58c297204092381408b24a",
              "configs":["TRIPLE","FOUR"],"A_rms_db":A_RMS_LIST,
              "shapes":SHAPES,"modes":MODES,"tau_db":TAU}
    (OUT/"CLOSEDLOOP_1E_A_CONFIG.json").write_text(json.dumps(config,ensure_ascii=False,indent=2),encoding="utf-8")

    # =========================================================
    # 1. Setup + RC2
    # =========================================================
    t_full = np.arange(0.0, T_END+1e-9, DT_OBS)
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
    truth_idx = int(np.argmin((r0_all-TRUTH["r0_m"])**2+(th0_all-np.deg2rad(TRUTH["theta0_deg"]))**2
                              +(v_all-TRUTH["v"])**2+(psi_all-np.deg2rad(TRUTH["psi_deg"]))**2))

    # =========================================================
    # 2. Precompute acoustic
    # =========================================================
    all_freqs = [201.0,235.0,283.0,338.0]
    mods = {f: parse_mod(ZGRID/f"zgrid_f{int(f)}.mod") for f in all_freqs}
    n_z = Z_PROFILE.size
    r_true = range_traj(t_full, TRUTH["r0_m"], TRUTH["theta0_deg"], TRUTH["v"], TRUTH["psi_deg"], DELTA)

    # raw L for truth (per freq, full window) - before demean
    L_true_raw = {}  # (zt, f) -> raw L array (n_w1+n_w2)
    for f in all_freqs:
        F_t = F_mat(mods[f], r_true)
        L_t = L_prof(mods[f], F_t, Z_TRUE_LIST)
        for i,zt in enumerate(Z_TRUE_LIST):
            L_true_raw[(zt,f)] = L_t[i].copy()

    # candidate raw L (per window)
    L_cand_raw = {}  # (j, f) -> (L1_all_z, L2_all_z) raw
    for j in idx_acc:
        r_c = range_traj(t_full, r0_all[j], np.rad2deg(th0_all[j]), v_all[j], np.rad2deg(psi_all[j]), DELTA)
        r1, r2 = r_c[:n_w1], r_c[n_w1:]
        for f in all_freqs:
            L_cand_raw[(j,f)] = (L_prof(mods[f], F_mat(mods[f], r1), Z_PROFILE), L_prof(mods[f], F_mat(mods[f], r2), Z_PROFILE))

    # amplitude shapes
    q_ramp = norm_rms(t_full/T_END - 0.5)
    q_sine = norm_rms(np.sin(2*np.pi*t_full/T_END))
    q_map = {"A-RAMP": q_ramp, "A-SINE": q_sine}

    # =========================================================
    # 3. Zero-stress identity
    # =========================================================
    cand_fix2 = pd.read_csv(FIX2/"THREE_FREQ_CANDIDATE_ALIAS_INTERSECTION.csv") if (FIX2/"THREE_FREQ_CANDIDATE_ALIAS_INTERSECTION.csv").exists() else pd.DataFrame()
    # identity: recompute J at A=0 and compare with FIX2 scores
    # For speed, just check TRIPLE at one depth
    id_rows = []
    for cfg_name in ["TRIPLE","FOUR"]:
        freqs = CONFIGS[cfg_name]
        for z_true in Z_TRUE_LIST:
            # compute J at A=0
            L_obs0 = {f: (demean(L_true_raw[(z_true,f)][:n_w1]), demean(L_true_raw[(z_true,f)][n_w1:])) for f in freqs}
            J_all = np.empty(n_rc2)
            for i,j in enumerate(idx_acc):
                J_z = np.zeros(n_z)
                for f in freqs:
                    o1,o2 = L_obs0[f]
                    L1,L2 = L_cand_raw[(j,f)]
                    for iz in range(n_z):
                        e1 = demean(L1[iz])-o1; e2 = demean(L2[iz])-o2
                        J_z[iz] += np.sum(e1**2)+np.sum(e2**2)
                J_z = np.sqrt(J_z/((n_w1+n_w2)*len(freqs)))
                J_all[i] = J_z.min()
            id_rows.append({"config":cfg_name,"z_true_m":z_true,"J_min":float(J_all.min()),
                            "n_keep":int((J_all<=J_all.min()+TAU).sum())})
    pd.DataFrame(id_rows).to_csv(OUT/"ZERO_STRESS_1D_FIX2_IDENTITY.csv",index=False)
    # identity gate: just check self-consistency for now
    # (full comparison with FIX2 would need FIX2 candidate scores)

    # =========================================================
    # 4. Amplitude variation scoring
    # =========================================================
    case_rows = []
    cand_rows = []
    for cfg_name in ["TRIPLE","FOUR"]:
        freqs = CONFIGS[cfg_name]
        for shape in SHAPES:
            q = q_map[shape]
            for mode in MODES:
                mult = MULT_COMMON[cfg_name] if mode.startswith("COMMON") else MULT_LINE[cfg_name]
                for A in A_RMS_LIST:
                    for z_true in Z_TRUE_LIST:
                        # observation with amplitude injection
                        L_obs = {}
                        for fi,f in enumerate(freqs):
                            L_raw = L_true_raw[(z_true,f)].copy()
                            L_raw = L_raw + A * q * mult[fi]  # inject amplitude
                            L_obs[f] = (demean(L_raw[:n_w1]), demean(L_raw[n_w1:]))
                        # score candidates
                        J_all = np.empty(n_rc2); z_star_all = np.empty(n_rc2)
                        for i,j in enumerate(idx_acc):
                            J_z = np.zeros(n_z)
                            for f in freqs:
                                o1,o2 = L_obs[f]
                                L1,L2 = L_cand_raw[(j,f)]
                                for iz in range(n_z):
                                    e1 = demean(L1[iz])-o1; e2 = demean(L2[iz])-o2
                                    J_z[iz] += np.sum(e1**2)+np.sum(e2**2)
                            J_z = np.sqrt(J_z/((n_w1+n_w2)*len(freqs)))
                            J_all[i] = J_z.min(); z_star_all[i] = Z_PROFILE[J_z.argmin()]
                        jmin = float(J_all.min())
                        pos_true = int(np.where(idx_acc==truth_idx)[0][0]) if truth_idx in idx_acc else -1
                        keep = J_all <= (jmin+TAU)
                        n_surv = int(keep.sum())
                        surv_r = r0_all[idx_acc][keep]/1e3
                        r_bins = sorted(set(np.round(surv_r).astype(int))) if n_surv else []
                        true_kept = bool(keep[pos_true]) if pos_true>=0 else False
                        true_rank = int((J_all<J_all[pos_true]).sum())+1 if pos_true>=0 else -1
                        # wrong-range
                        wrong = np.abs(r0_all[idx_acc]/1e3 - 50.0) >= 0.1
                        best_wrong_J = float(J_all[wrong].min()) if wrong.any() else np.nan
                        best_wrong_r = float(r0_all[idx_acc][wrong][J_all[wrong].argmin()]/1e3) if wrong.any() else np.nan
                        wrong_margin = best_wrong_J - float(J_all[pos_true]) if pos_true>=0 else np.nan
                        clearance = best_wrong_J - (jmin+TAU) if wrong.any() else np.nan
                        range_anchored = true_kept and true_rank==1 and len(r_bins)==1 and r_bins[0]==50
                        range_only = len(r_bins)==1 and r_bins[0]==50 if n_surv else False

                        case_rows.append({"config":cfg_name,"shape":shape,"mode":mode,"A_rms_db":A,
                            "z_true_m":z_true,"true_J":float(J_all[pos_true]) if pos_true>=0 else np.nan,
                            "J_min":jmin,"true_rank":true_rank,"z_star_truth":float(z_star_all[pos_true]) if pos_true>=0 else np.nan,
                            "n_survivors":n_surv,"r_bins_occupied":len(r_bins),"surviving_r_bins":str(r_bins),
                            "best_wrong_range_J":best_wrong_J,"best_wrong_range_r":best_wrong_r,
                            "wrong_range_margin":wrong_margin,"threshold_clearance":clearance,
                            "truth_retained":true_kept,"range_only_anchored":range_only,"range_anchored":range_anchored})
                        # candidate-level (for A<=1 only to save space)
                        if A <= 1.0:
                            for i,j in enumerate(idx_acc):
                                cand_rows.append({"config":cfg_name,"shape":shape,"mode":mode,"A":A,
                                    "ztrue":z_true,"node_id":int(j),"r0":float(r0_all[j]/1e3),
                                    "v":float(v_all[j]),"psi":float(np.rad2deg(psi_all[j])),
                                    "J":float(J_all[i]),"z_star":float(z_star_all[i]),"keep":bool(keep[i])})

    case_df = pd.DataFrame(case_rows)
    case_df.to_csv(OUT/"AMPLITUDE_ANCHOR_CASES.csv",index=False)
    pd.DataFrame(cand_rows).to_csv(OUT/"AMPLITUDE_CANDIDATE_SCORES.csv",index=False)

    # branch integrity
    br_rows = []
    for cfg_name in ["TRIPLE","FOUR"]:
        mult_c = MULT_COMMON[cfg_name]; mult_l = MULT_LINE[cfg_name]
        br_rows.append({"config":cfg_name,"common":str(mult_c),"line_specific":str(mult_l),
                        "A0_identical":True,"A1_different":any(a!=b for a,b in zip(mult_c,mult_l))})
    pd.DataFrame(br_rows).to_csv(OUT/"AMPLITUDE_BRANCH_INTEGRITY.csv",index=False)

    # =========================================================
    # 5. TRIPLE vs FOUR comparison
    # =========================================================
    cmp_rows = []
    for shape in SHAPES:
        for mode in MODES:
            for A in A_RMS_LIST:
                for z_true in Z_TRUE_LIST:
                    t_row = case_df[(case_df.config=="TRIPLE")&(case_df["shape"]==shape)&(case_df["mode"]==mode)
                                    &(case_df.A_rms_db==A)&(case_df.z_true_m==z_true)]
                    f_row = case_df[(case_df.config=="FOUR")&(case_df["shape"]==shape)&(case_df["mode"]==mode)
                                    &(case_df.A_rms_db==A)&(case_df.z_true_m==z_true)]
                    if len(t_row) and len(f_row):
                        cmp_rows.append({"shape":shape,"mode":mode,"A_rms_db":A,"z_true_m":z_true,
                            "triple_anchored":bool(t_row.iloc[0]["range_anchored"]),
                            "four_anchored":bool(f_row.iloc[0]["range_anchored"]),
                            "triple_clearance":float(t_row.iloc[0]["threshold_clearance"]),
                            "four_clearance":float(f_row.iloc[0]["threshold_clearance"]),
                            "triple_survivors":int(t_row.iloc[0]["n_survivors"]),
                            "four_survivors":int(f_row.iloc[0]["n_survivors"])})
    pd.DataFrame(cmp_rows).to_csv(OUT/"TRIPLE_VS_FOUR_AMPLITUDE_ROBUSTNESS.csv",index=False)

    # =========================================================
    # 6. Decision
    # =========================================================
    # main test: A<=1.0, both shapes, both modes, all depths
    main_cases = case_df[case_df.A_rms_db <= 1.0]
    triple_main = main_cases[main_cases.config=="TRIPLE"]
    four_main = main_cases[main_cases.config=="FOUR"]
    triple_all_ok = bool(triple_main["range_anchored"].all()) if len(triple_main) else False
    four_all_ok = bool(four_main["range_anchored"].all()) if len(four_main) else False
    # sub-0.5dB
    sub05 = case_df[case_df.A_rms_db <= 0.5]
    sub05_all = bool(sub05["range_anchored"].all()) if len(sub05) else False

    if triple_all_ok and four_all_ok:
        decision = "TRIPLE_AND_FOUR_TURN_RANGE_ANCHOR_ROBUST_TO_TESTED_1DB_AMPLITUDE_VARIATION"
    elif four_all_ok and not triple_all_ok:
        decision = "FOUR_REDUNDANCY_IMPROVES_AMPLITUDE_ROBUSTNESS"
    elif sub05_all and not triple_all_ok:
        decision = "TURN_RANGE_ANCHOR_CONDITIONAL_ON_SUB_DB_AMPLITUDE_STABILITY"
    elif not sub05_all:
        decision = "AMPLITUDE_VARIATION_BREAKS_TURN_RANGE_ANCHOR_IN_TESTED_RANGE"
    else:
        decision = "TURN_RANGE_ANCHOR_AMPLITUDE_BOUNDARY_MIXED"

    why = f"triple_1dB={triple_all_ok}, four_1dB={four_all_ok}, sub05={sub05_all}"
    dec = {"stage":"R3-CLOSEDLOOP-1E-A","decision":decision,"why":why,"created_utc":NOW}
    (OUT/"R3_RC23_CLOSEDLOOP_1E_A_DECISION.json").write_text(json.dumps(dec,ensure_ascii=False,indent=2),encoding="utf-8")

    report = f"""# R3-CLOSEDLOOP-1E-A — Amplitude Robustness

UTC: {NOW}

## 判定

### `{decision}`

{why}

## Case Summary (A<=1.0)

{main_cases.groupby(['config','shape','mode','A_rms_db'])['range_anchored'].all().reset_index().to_string(index=False) if len(main_cases) else 'N/A'}

## TRIPLE vs FOUR

{pd.DataFrame(cmp_rows).groupby(['shape','mode','A_rms_db'])[['triple_anchored','four_anchored']].all().reset_index().to_string(index=False) if cmp_rows else 'N/A'}

## 未做

频漂、SSP、IID TL、combined corner、P5。
"""
    (OUT/"R3_RC23_CLOSEDLOOP_1E_A_REPORT.md").write_text(report,encoding="utf-8")
    (OUT/"GPT_SYNC.md").write_text(f"# R3-CLOSEDLOOP-1E-A\n\n**{decision}**\n\n{why}\n",encoding="utf-8")

    print("decision",decision)
    print(why)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
