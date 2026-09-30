#!/usr/bin/env python3
"""R3-CLOSEDLOOP-1E-C: Turn-assisted multi-freq range anchor under SSP mismatch.

E1/E2 truth vs E0 nominal template. 4G frozen SSP stress. 2 configs x 3 envs x 3 depths.
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
OUT = ROOT / "results" / "R3_RC23_CLOSEDLOOP" / "R3_RC23_TURN_SSP_MISMATCH"
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
ENVS = ["E0", "E1", "E2"]


def split_env_tail(text):
    lines = text.splitlines()
    for i, ln in enumerate(lines):
        if ln.strip().startswith("'R'") or ln.strip() == "R":
            return lines[:i], lines[i:]
    return lines, []


def read_ssp(text):
    zs, cs = [], []
    head, _ = split_env_tail(text)
    for ln in head:
        p = ln.split()
        if len(p) >= 6:
            try:
                z = float(p[0]); [float(p[i]) for i in range(6)]; c = float(p[1])
            except ValueError:
                continue
            if 0.0 <= z <= 5000.0:
                zs.append(z); cs.append(c)
    return np.array(zs, float), np.array(cs, float)


def shift_ssp(z, c, dz=0.0, dc=0.0):
    zq = np.clip(z - dz, float(z.min()), float(z.max()))
    return np.interp(zq, z, c) + dc


def write_env_pert(path, freq, z_ssp, c_ssp, env0_path):
    text0 = env0_path.read_text(encoding="utf-8")
    head, tail = split_env_tail(text0)
    z_ssp = np.asarray(z_ssp, float); c_ssp = np.asarray(c_ssp, float)
    new_head = []
    for i, ln in enumerate(head):
        if i == 0: new_head.append("'PERT'")
        elif i == 1: new_head.append(f"{freq:.3f}")
        else:
            p = ln.split()
            if len(p) >= 6:
                try:
                    z = float(p[0]); [float(p[j]) for j in range(6)]
                    c_new = float(np.interp(z, z_ssp, c_ssp))
                    new_head.append(f"{z:.3f} {c_new:.4f} 0.0 1.0 0.0 0.0")
                    continue
                except ValueError:
                    pass
            new_head.append(ln)
    path.write_text("\n".join(new_head + tail) + "\n", encoding="utf-8")
    return tail


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


def main():
    config = {"stage":"R3-CLOSEDLOOP-1E-C","created_utc":NOW,
              "baseline_commit":"8fceb304801b3027259400d93214d7816d1fd9ba",
              "ssp_stress":"FROZEN_SYNTHETIC_SSP_STRESS_FROM_REANCHOR_4G",
              "envs":ENVS,"configs":["TRIPLE","FOUR"]}
    (OUT/"CLOSEDLOOP_1E_C_CONFIG.json").write_text(json.dumps(config,ensure_ascii=False,indent=2),encoding="utf-8")

    # =========================================================
    # 1. SSP definitions + identity
    # =========================================================
    env0_text = (ZGRID / "zgrid_f235.env").read_text(encoding="utf-8")
    z0, c0 = read_ssp(env0_text)
    c_E1 = shift_ssp(z0, c0, dz=50.0, dc=0.0)
    c_E2 = shift_ssp(z0, c0, dz=150.0, dc=2.0)

    ssp_rows = []
    for i in range(len(z0)):
        ssp_rows.append({"z": z0[i], "c_E0": c0[i], "c_E1": c_E1[i], "c_E2": c_E2[i],
                         "c_E1_formula": float(np.interp(np.clip(z0[i]-50, z0.min(), z0.max()), z0, c0)),
                         "c_E2_formula": float(np.interp(np.clip(z0[i]-150, z0.min(), z0.max()), z0, c0)) + 2.0})
    ssp_df = pd.DataFrame(ssp_rows)
    ssp_df["abs_err_E1"] = np.abs(ssp_df["c_E1"] - ssp_df["c_E1_formula"])
    ssp_df["abs_err_E2"] = np.abs(ssp_df["c_E2"] - ssp_df["c_E2_formula"])
    ssp_df.to_csv(OUT/"SSP_DEFINITION_4G_IDENTITY.csv", index=False)
    ssp_id_ok = float(ssp_df[["abs_err_E1","abs_err_E2"]].max().max()) <= 1e-10

    # =========================================================
    # 2. Generate E1/E2 KRAKEN models
    # =========================================================
    ssp_variants = {"E0": (z0, c0), "E1": (z0, c_E1), "E2": (z0, c_E2)}
    mod_cache = {}  # (env, freq) -> mod
    build_rows = []
    tail_rows = []

    for env_name in ["E1", "E2"]:
        z_ssp, c_ssp = ssp_variants[env_name]
        for f in [201.0, 235.0, 283.0, 338.0]:
            key = f"{env_name}_f{int(f)}"
            envp = WORK / f"{key}.env"
            modp = WORK / f"{key}.mod"
            tail = write_env_pert(envp, f, z_ssp, c_ssp, ZGRID / f"zgrid_f{int(f)}.env")
            # tail integrity
            env0_f = ZGRID / f"zgrid_f{int(f)}.env"
            _, tail0 = split_env_tail(env0_f.read_text(encoding="utf-8"))
            tail_ok = tail == tail0
            tail_rows.append({"env": env_name, "freq": f, "tail_identical": tail_ok})

            try:
                proc = subprocess.run([str(AT_BIN/"kraken.exe"), key], cwd=str(WORK),
                                      capture_output=True, text=True, timeout=120)
                rc = proc.returncode
            except Exception:
                rc = -1
            mod_exists = modp.exists()
            parse_ok = False; M = 0; finite_k = False; finite_phi = False
            if mod_exists:
                try:
                    m = parse_mod(modp)
                    mod_cache[(env_name, f)] = m
                    M = m["M"]; parse_ok = M > 0
                    finite_k = bool(np.all(np.isfinite(m["k"])))
                    finite_phi = bool(np.all(np.isfinite(m["phi"])))
                except Exception:
                    parse_ok = False
            build_rows.append({"env": env_name, "freq": f, "return_code": rc,
                               "mod_exists": mod_exists, "mod_size": modp.stat().st_size if mod_exists else 0,
                               "parse_ok": parse_ok, "M": M, "finite_k": finite_k, "finite_phi": finite_phi,
                               "depth_min": float(m["depths"].min()) if parse_ok else np.nan,
                               "depth_max": float(m["depths"].max()) if parse_ok else np.nan,
                               "has_150_250": bool(parse_ok and m["depths"].min() <= 150 and m["depths"].max() >= 250)})

    # E0 from ZGRID
    for f in [201.0, 235.0, 283.0, 338.0]:
        mod_cache[("E0", f)] = parse_mod(ZGRID / f"zgrid_f{int(f)}.mod")
        build_rows.append({"env": "E0", "freq": f, "return_code": 0, "mod_exists": True,
                           "parse_ok": True, "M": mod_cache[("E0", f)]["M"],
                           "finite_k": True, "finite_phi": True, "has_150_250": True})

    build_df = pd.DataFrame(build_rows)
    build_df.to_csv(OUT/"SSP_MODEL_BUILD_AUDIT.csv", index=False)
    pd.DataFrame(tail_rows).to_csv(OUT/"SSP_ENV_TAIL_INTEGRITY.csv", index=False)

    ssp_build_ok = bool(build_df[build_df.env != "E0"]["parse_ok"].all()) and bool(build_df[build_df.env != "E0"]["finite_k"].all())
    if not ssp_build_ok:
        (OUT/"R3_RC23_CLOSEDLOOP_1E_C_DECISION.json").write_text(
            json.dumps({"decision":"CLOSEDLOOP_1E_C_BLOCKED_BY_SSP_MODEL_INTEGRITY","created_utc":NOW},indent=2),encoding="utf-8")
        print("BLOCKED_BY_SSP_MODEL_INTEGRITY")
        return 1

    # =========================================================
    # 3. RC2 + trajectories
    # =========================================================
    t_full = np.arange(0.0, T_END+1e-9, DT_OBS)
    n_w1 = int(np.sum(t_full <= T_TURN+1e-9)); n_w2 = len(t_full)-n_w1
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

    r_true = range_traj(t_full, TRUTH["r0_m"], TRUTH["theta0_deg"], TRUTH["v"], TRUTH["psi_deg"], DELTA)
    r_cand = {}
    for j in idx_acc:
        r_cand[j] = range_traj(t_full, r0_all[j], np.rad2deg(th0_all[j]), v_all[j], np.rad2deg(psi_all[j]), DELTA)

    # =========================================================
    # 4. Precompute L for E0 candidates and E0/E1/E2 truth
    # =========================================================
    n_z = Z_PROFILE.size
    # Candidate templates: always E0
    L_cand_e0 = {}
    for j in idx_acc:
        for f in [201.0, 235.0, 283.0, 338.0]:
            m = mod_cache[("E0", f)]
            F_c = F_mat(m, r_cand[j])
            L_cand_e0[(j, f)] = (L_prof(m, F_mat(m, r_cand[j][:n_w1]), Z_PROFILE),
                                  L_prof(m, F_mat(m, r_cand[j][n_w1:]), Z_PROFILE))

    # Truth L for each env
    L_truth = {}  # (env, f, zt) -> (w1, w2)
    for env_name in ENVS:
        for f in [201.0, 235.0, 283.0, 338.0]:
            m = mod_cache[(env_name, f)]
            F_t = F_mat(m, r_true)
            L_t = L_prof(m, F_t, Z_TRUE_LIST)
            for i, zt in enumerate(Z_TRUE_LIST):
                L_truth[(env_name, f, zt)] = (demean(L_t[i, :n_w1]), demean(L_t[i, n_w1:]))

    # =========================================================
    # 5. E0 identity check
    # =========================================================
    ref = pd.read_csv(FIX/"SUBSET_CANDIDATE_SCORES_FIXED.csv")
    e0_id_rows = []
    for cfg_name, sub_tag in [("TRIPLE","201+235+283"),("FOUR","201+235+283+338")]:
        freqs = CONFIGS[cfg_name]
        for z_true in Z_TRUE_LIST:
            L_obs = {f: L_truth[("E0", f, z_true)] for f in freqs}
            J_all = np.empty(n_rc2)
            for i,j in enumerate(idx_acc):
                J_z = np.zeros(n_z)
                for f in freqs:
                    o1,o2 = L_obs[f]
                    L1,L2 = L_cand_e0[(j,f)]
                    for iz in range(n_z):
                        e1 = demean(L1[iz])-o1; e2 = demean(L2[iz])-o2
                        J_z[iz] += np.sum(e1**2)+np.sum(e2**2)
                J_z = np.sqrt(J_z/((n_w1+n_w2)*len(freqs)))
                J_all[i] = J_z.min()
            # compare with ref
            r = ref[(ref["subset"]==sub_tag)&(ref["z_true_m"]==z_true)]
            ids_c = set(idx_acc.tolist()); ids_r = set(r["node_id"])
            max_err = 0.0
            for nid in ids_c & ids_r:
                j_c = J_all[list(idx_acc).index(nid)]
                j_r = float(r[r.node_id==nid]["J"].values[0])
                max_err = max(max_err, abs(j_c-j_r))
            e0_id_rows.append({"config":cfg_name,"z_true_m":z_true,"ids_match":ids_c==ids_r,
                               "max_J_err":max_err,"pass":ids_c==ids_r and max_err<=1e-10})
    e0_id = pd.DataFrame(e0_id_rows)
    e0_id.to_csv(OUT/"E0_1D_FIX2_IDENTITY.csv",index=False)
    if not e0_id["pass"].all():
        (OUT/"R3_RC23_CLOSEDLOOP_1E_C_DECISION.json").write_text(
            json.dumps({"decision":"CLOSEDLOOP_1E_C_BLOCKED_BY_E0_IDENTITY","created_utc":NOW},indent=2),encoding="utf-8")
        print("BLOCKED_BY_E0_IDENTITY")
        return 1

    # =========================================================
    # 6. SSP mismatch scoring
    # =========================================================
    cand_rows = []; case_rows = []; depth_rows = []

    for cfg_name in ["TRIPLE","FOUR"]:
        freqs = CONFIGS[cfg_name]
        for env_truth in ["E0","E1","E2"]:
            for z_true in Z_TRUE_LIST:
                # observation from env_truth
                L_obs = {f: L_truth[(env_truth, f, z_true)] for f in freqs}
                # candidates always E0
                J_all = np.empty(n_rc2); z_star_all = np.empty(n_rc2)
                for i,j in enumerate(idx_acc):
                    J_z = np.zeros(n_z)
                    for f in freqs:
                        o1,o2 = L_obs[f]
                        L1,L2 = L_cand_e0[(j,f)]
                        for iz in range(n_z):
                            e1 = demean(L1[iz])-o1; e2 = demean(L2[iz])-o2
                            J_z[iz] += np.sum(e1**2)+np.sum(e2**2)
                    J_z = np.sqrt(J_z/((n_w1+n_w2)*len(freqs)))
                    J_all[i] = J_z.min(); z_star_all[i] = Z_PROFILE[J_z.argmin()]

                jmin = float(J_all.min())
                pos_true = int(np.where(idx_acc==truth_idx)[0][0])
                true_J = float(J_all[pos_true])
                true_rank = int((J_all<true_J).sum())+1
                true_kept = bool(true_J<=jmin+TAU)
                keep = J_all<=(jmin+TAU)
                z_star_truth = float(z_star_all[pos_true])

                surv = idx_acc[keep]; surv_r = r0_all[surv]/1e3
                r_keep = sorted(set(int(round(x)) for x in surv_r)) if len(surv_r) else []
                r_best = sorted(set(int(round(x)) for x in r0_all[idx_acc][J_all<=jmin+1e-10]/1e3))
                wrong_mask = np.abs(r0_all[idx_acc]/1e3-50.0)>=0.1
                best_wrong_J = float(J_all[wrong_mask].min()) if wrong_mask.any() else np.nan
                strict = true_kept and true_rank==1 and r_keep==[50]

                case_rows.append({
                    "config":cfg_name,"env_truth":env_truth,"z_true_m":z_true,
                    "true_J":true_J,"J_min":jmin,"true_rank":true_rank,"truth_retained":true_kept,
                    "z_star_truth":z_star_truth,"profiled_z_shift":z_star_truth-z_true,
                    "best_range_bins":str(r_best),"best_range_unique_50":r_best==[50],
                    "surviving_r_bins":str(r_keep),"n_surviving_range_bins":len(r_keep),
                    "r_width_km":float(surv_r.max()-surv_r.min()) if len(surv_r) else 0.0,
                    "max_surviving_deviation_km":float(np.max(np.abs(surv_r-50.0))) if len(surv_r) else 0.0,
                    "best_wrong_range_J":best_wrong_J,
                    "wrong_range_margin":best_wrong_J-true_J if np.isfinite(best_wrong_J) else np.nan,
                    "threshold_clearance":best_wrong_J-(jmin+TAU) if np.isfinite(best_wrong_J) else np.nan,
                    "strict_single_bin":strict,
                })
                depth_rows.append({"config":cfg_name,"env":env_truth,"z_true":z_true,
                                   "z_star_truth":z_star_truth,"z_shift":z_star_truth-z_true,"true_J":true_J})
                for i,j in enumerate(idx_acc):
                    cand_rows.append({"config":cfg_name,"env_truth":env_truth,"z_true":z_true,
                                      "node_id":int(j),"r0":float(r0_all[j]/1e3),
                                      "theta0":float(np.rad2deg(th0_all[j])),"v":float(v_all[j]),
                                      "psi":float(np.rad2deg(psi_all[j])),
                                      "J":float(J_all[i]),"z_star":float(z_star_all[i]),
                                      "keep":bool(keep[i]),"is_truth_state":bool(j==truth_idx),
                                      "is_true_range":bool(abs(r0_all[j]/1e3-50.0)<0.1)})

    cand_df = pd.DataFrame(cand_rows)
    cand_df.to_csv(OUT/"SSP_MISMATCH_CANDIDATE_SCORES.csv",index=False)
    case_df = pd.DataFrame(case_rows)
    case_df.to_csv(OUT/"SSP_MISMATCH_ANCHOR_CASES.csv",index=False)
    pd.DataFrame(depth_rows).to_csv(OUT/"SSP_DEPTH_PROFILING_COMPENSATION.csv",index=False)

    # TRIPLE vs FOUR
    cmp_rows = []
    for env_truth in ["E0","E1","E2"]:
        for z_true in Z_TRUE_LIST:
            t = case_df[(case_df.config=="TRIPLE")&(case_df.env_truth==env_truth)&(case_df.z_true_m==z_true)]
            f = case_df[(case_df.config=="FOUR")&(case_df.env_truth==env_truth)&(case_df.z_true_m==z_true)]
            if len(t) and len(f):
                cmp_rows.append({"env_truth":env_truth,"z_true_m":z_true,
                                 "triple_strict":bool(t.iloc[0]["strict_single_bin"]),
                                 "four_strict":bool(f.iloc[0]["strict_single_bin"]),
                                 "triple_rank":int(t.iloc[0]["true_rank"]),
                                 "four_rank":int(f.iloc[0]["true_rank"])})
    pd.DataFrame(cmp_rows).to_csv(OUT/"TRIPLE_VS_FOUR_SSP_ROBUSTNESS.csv",index=False)

    # working range
    wr_rows = []
    for cfg in ["TRIPLE","FOUR"]:
        for env in ["E0","E1","E2"]:
            sub = case_df[(case_df.config==cfg)&(case_df.env_truth==env)]
            wr_rows.append({"config":cfg,"env":env,"n_cases":len(sub),
                            "strict":int(sub["strict_single_bin"].sum()),
                            "truth_ok":int(sub["truth_retained"].sum()),
                            "rank1":int((sub["true_rank"]==1).sum()),
                            "best50":int(sub["best_range_unique_50"].sum()),
                            "max_width":float(sub["r_width_km"].max()),
                            "min_margin":float(sub["wrong_range_margin"].min())})
    wr_df = pd.DataFrame(wr_rows)
    wr_df.to_csv(OUT/"SSP_WORKING_RANGE_SUMMARY.csv",index=False)

    # =========================================================
    # 7. Decision
    # =========================================================
    e12 = case_df[case_df.env_truth.isin(["E1","E2"])]
    all_strict = bool(e12["strict_single_bin"].all())
    triple_ok = bool(e12[e12.config=="TRIPLE"]["strict_single_bin"].all())
    four_ok = bool(e12[e12.config=="FOUR"]["strict_single_bin"].all())
    all_rank = bool((e12["true_rank"]==1).all()) and bool(e12["truth_retained"].all())
    all_50 = bool(e12["best_range_unique_50"].all())

    if all_strict:
        decision = "TRIPLE_AND_FOUR_TURN_RANGE_ANCHOR_SURVIVES_TESTED_SSP_STRESS"
    elif four_ok and not triple_ok:
        decision = "FOUR_REDUNDANCY_IMPROVES_SSP_MISMATCH_ROBUSTNESS"
    elif all_rank and all_50:
        decision = "SSP_MISMATCH_WEAKENS_STRICT_ANCHOR_BUT_PRESERVES_OPTIMAL_RANGE"
    else:
        decision = "SSP_MISMATCH_DEGRADES_RANGE_RANKING_IN_TESTED_STRESS"

    why = f"all_strict={all_strict}, triple={triple_ok}, four={four_ok}, rank1={all_rank}, best50={all_50}"
    dec = {"stage":"R3-CLOSEDLOOP-1E-C","decision":decision,"why":why,
           "ssp_stress":"FROZEN_SYNTHETIC_SSP_STRESS_FROM_REANCHOR_4G","created_utc":NOW}
    (OUT/"R3_RC23_CLOSEDLOOP_1E_C_DECISION.json").write_text(json.dumps(dec,ensure_ascii=False,indent=2),encoding="utf-8")

    report = f"""# R3-CLOSEDLOOP-1E-C — SSP Mismatch

UTC: {NOW}

## SSP Stress

E1: shift 50m; E2: shift 150m + 2m/s (4G frozen)

## 判定

### `{decision}`

{why}

## Working Range

{wr_df.to_string(index=False)}

## Depth Profiling Compensation

{pd.DataFrame(depth_rows).to_string(index=False) if depth_rows else 'N/A'}

## 未做

幅漂、频漂、跟踪误差、联合corner、P5。
"""
    (OUT/"R3_RC23_CLOSEDLOOP_1E_C_REPORT.md").write_text(report,encoding="utf-8")
    (OUT/"GPT_SYNC.md").write_text(f"# 1E-C\n\n**{decision}**\n\n{why}\n",encoding="utf-8")

    print("decision",decision)
    print(why)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
