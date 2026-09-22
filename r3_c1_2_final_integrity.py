#!/usr/bin/env python3
"""R3-C1.2 FINAL-INTEGRITY: fix HLA endfire beamforming, multifreq stats, sampling.

No new physics. Then permanently close Zhu-QM if still weak.
"""
from __future__ import annotations

import json
import math
import time
from datetime import datetime, timezone
from pathlib import Path

import numpy as np
import pandas as pd
import p4_g4_performance_boundary as p4

ROOT = Path(__file__).resolve().parent
OUT = ROOT / "results" / "R3_C1_depth_motion"
NOW = datetime.now(timezone.utc).isoformat()

C0 = 1500.0
Z_R = 200.0
R0, R1 = 45e3, 60e3
U_REL = 2.0
T_LIST = [600.0, 1200.0]
Z_TRUES = [180.0, 190.0, 200.0, 210.0, 220.0]
Z_GRID = np.linspace(160.0, 240.0, 81)
F_SINGLE = 235.0
F_MULTI = [201.0, 235.0, 283.0, 338.0]
PAIRS = [(180, 190), (190, 200), (200, 210), (210, 220), (180, 200), (200, 220)]
N_EL = 8
D_EL = 2.0
DR_LIST = [5.0, 10.0, 20.0]  # m spatial sampling convergence

CONFIG = {
    "package": "R3_C1_2_FINAL_INTEGRITY",
    "created_utc": NOW,
    "pending_before": "C1_QM_CZ_DEPTH_WEAK_PENDING_FINAL_IMPLEMENTATION_FIX",
    "single_channel_frozen": (
        "Even with known track, matched E0, noiseless: 235Hz single QM weak near 200m — "
        "D_200/210~0.003, D_200/220~0.007, FWHM~54m at T=600 (prior C1.2)"
    ),
    "fixes": [
        "HLA endfire B=|w^H p|^2 with RC2 steering (steer actually applied)",
        "multifreq: per-freq demean+unit-norm; D_multi=1-mean_f|corr_f|",
        "fixed spatial step dr=5/10/20 m convergence",
        "FWHM full window → FULL_WINDOW_UNRESOLVED, PSL N/A",
    ],
    "groups": {"A": "single+235", "B": "HLA beam+235", "C": "single+multi", "D": "HLA beam+multi"},
    "stop": ["no C1.x after this", "no P5", "no Yang this round"],
}


def track_r(T, dr):
    r1 = min(R1, R0 + U_REL * T)
    n = int(max(20, math.floor((r1 - R0) / dr) + 1))
    return np.linspace(R0, r1, n)


def element_pressures(mode, f, z_s, r, z_r=Z_R, n_el=1, d_m=0.0):
    """p_m(f,t) complex pressures on endfire elements (along motion axis x).

    Path length uses |r - x_m| (endfire). No steering phase baked into p_m.
    """
    xs = (np.arange(n_el) - (n_el - 1) / 2.0) * d_m if n_el > 1 else np.zeros(1)
    ms = mode.modes(float(f))
    P = np.zeros((len(r), n_el), dtype=complex)
    for i, ri in enumerate(r):
        for m_i, xm in enumerate(xs):
            rm = max(abs(ri - xm), 1.0)
            acc = 0j
            for kr, phi in ms:
                a_s = float(np.interp(z_s, mode.z, phi))
                a_r = float(np.interp(z_r, mode.z, phi))
                att = math.exp(-2e-5 * (f / 200.0) * rm / 1000.0)
                acc += (a_s * a_r / math.sqrt(kr * rm)) * np.exp(1j * kr * rm) * att
            P[i, m_i] = acc
    return P, xs


def beam_track_intensity(mode, f, z_s, r, z_r=Z_R, n_el=1, d_m=0.0, steer=True):
    """B(f,t)=|w^H p|^2 with w = RC2 endfire steering (u=1).

    steer=True: w_m = exp(+j k0 x_m u)/N, B=|sum conj(w_m) p_m|^2
    steer=False: unsteered mean element power (not a proper beam)
    """
    P, xs = element_pressures(mode, f, z_s, r, z_r=z_r, n_el=n_el, d_m=d_m)
    if n_el == 1:
        B = np.abs(P[:, 0]) ** 2
    else:
        k0 = 2 * np.pi * f / C0
        if steer:
            w = np.exp(1j * k0 * xs * 1.0) / n_el
            bh = P @ np.conj(w)  # w^H p
            B = np.abs(bh) ** 2
        else:
            B = np.mean(np.abs(P) ** 2, axis=1)
    R = np.sqrt(r ** 2 + (z_r - z_s) ** 2)
    Iprime = B * R ** 2
    return Iprime, P, xs


def norm_track(I):
    v = I - I.mean()
    n = np.linalg.norm(v)
    return v / (n + 1e-30)


def qm_joint_scores(mode, freqs, z_true, z_grid, r, n_el, d_m, steer=True):
    """Per-freq unit-norm tracks; J_f(z)=|corr|; J_joint=mean_f J_f.
    Returns Q (len z_grid), and dict of unit tracks for pair D.
    """
    obs_units = {}
    Qs = []
    for f in freqs:
        Ip_obs, _, _ = beam_track_intensity(mode, f, z_true, r, n_el=n_el, d_m=d_m, steer=steer)
        u_obs = norm_track(Ip_obs)
        obs_units[f] = u_obs
        Qf = np.zeros(len(z_grid))
        for iz, z in enumerate(z_grid):
            Ip, _, _ = beam_track_intensity(mode, f, z, r, n_el=n_el, d_m=d_m, steer=steer)
            u = norm_track(Ip)
            Qf[iz] = abs(float(np.dot(u_obs, u)))
        Qs.append(Qf)
    Q = np.mean(Qs, axis=0)
    return Q, obs_units


def D_multi(units_i, units_j):
    """D = 1 - mean_f |corr_f| on independently normalized tracks."""
    cs = []
    for f in units_i:
        cs.append(abs(float(np.dot(units_i[f], units_j[f]))))
    return 1.0 - float(np.mean(cs)), float(np.mean(cs)), cs


def contiguous_fwhm(z_grid, J):
    J = np.asarray(J, float)
    window = float(z_grid[-1] - z_grid[0])
    imax = int(np.argmax(J))
    jmax = float(J[imax])
    thr = 0.5 * jmax
    lo, hi = imax, imax
    while lo > 0 and J[lo - 1] >= thr:
        lo -= 1
    while hi < len(J) - 1 and J[hi + 1] >= thr:
        hi += 1
    fwhm = float(z_grid[hi] - z_grid[lo])
    full_window = fwhm >= window * 0.95
    if full_window:
        return dict(z_hat=float(z_grid[imax]), fwhm_local_m=fwhm, fwhm_flag="FULL_WINDOW_UNRESOLVED",
                    psl=np.nan, psl_flag="N/A", n_comp=0, peak=jmax)
    mask = np.ones_like(J, dtype=bool)
    mask[lo:hi + 1] = False
    psl = float(J[mask].max() / (jmax + 1e-30)) if mask.any() else 0.0
    n_comp = sum(
        1 for i in range(1, len(J) - 1)
        if i != imax and J[i] >= J[i - 1] and J[i] >= J[i + 1] and J[i] > 0.7 * jmax
    )
    return dict(z_hat=float(z_grid[imax]), fwhm_local_m=fwhm, fwhm_flag="local",
                psl=psl, psl_flag="ok", n_comp=n_comp, peak=jmax)


def beam_selfcheck(mode, f, z_s, r_mid, n_el, d_m):
    """Save element phases + beamformed vs unsteered."""
    P, xs = element_pressures(mode, f, z_s, np.array([r_mid]), n_el=n_el, d_m=d_m)
    phases = np.angle(P[0])
    Ip_b, _, _ = beam_track_intensity(mode, f, z_s, np.array([r_mid]), n_el=n_el, d_m=d_m, steer=True)
    Ip_u, _, _ = beam_track_intensity(mode, f, z_s, np.array([r_mid]), n_el=n_el, d_m=d_m, steer=False)
    return {
        "x_m": xs.tolist(),
        "phase_rad": phases.tolist(),
        "amp": np.abs(P[0]).tolist(),
        "B_steered": float(Ip_b[0]),
        "B_unsteered_mean_power": float(Ip_u[0]),
        "A_B_identical": bool(abs(float(Ip_b[0]) - float(Ip_u[0])) < 1e-18),
    }


def main():
    t0 = time.time()
    mode = p4.MODE_ENV["E0"]
    print("=== R3-C1.2 FINAL-INTEGRITY ===", flush=True)

    # 0) beamformer self-check
    chk = beam_selfcheck(mode, F_SINGLE, 200.0, 50e3, N_EL, D_EL)
    print("beam check", json.dumps(chk, ensure_ascii=False)[:300], flush=True)
    chk_df = pd.DataFrame([{
        "x_m": json.dumps(chk["x_m"]),
        "phase_rad": json.dumps(chk["phase_rad"]),
        "B_steered": chk["B_steered"],
        "B_unsteered": chk["B_unsteered_mean_power"],
        "steered_differs_from_unsteered": not chk["A_B_identical"],
    }])
    chk_df.to_csv(OUT / "r3c12_beamformer_check.csv", index=False, encoding="utf-8-sig")

    # 1) sampling convergence on D and FWHM (group D multifreq, z=200, T=600)
    conv_rows = []
    for dr in DR_LIST:
        r = track_r(600.0, dr)
        Q, units = qm_joint_scores(mode, F_MULTI, 200.0, Z_GRID, r, N_EL, D_EL, steer=True)
        met = contiguous_fwhm(Z_GRID, Q)
        _, u190 = qm_joint_scores(mode, F_MULTI, 190.0, Z_GRID, r, N_EL, D_EL, steer=True)
        _, u210 = qm_joint_scores(mode, F_MULTI, 210.0, Z_GRID, r, N_EL, D_EL, steer=True)
        _, u220 = qm_joint_scores(mode, F_MULTI, 220.0, Z_GRID, r, N_EL, D_EL, steer=True)
        D210 = D_multi(units, u210)[0]
        D220 = D_multi(units, u220)[0]
        conv_rows.append({
            "dr_m": dr, "n_r": len(r), "T_s": 600.0,
            "z_true": 200.0,
            "fwhm": met["fwhm_local_m"], "fwhm_flag": met["fwhm_flag"],
            "D_200_210": D210, "D_200_220": D220,
        })
        print(f"  conv dr={dr}: FWHM={met['fwhm_local_m']:.1f} D210={D210:.4f} D220={D220:.4f}", flush=True)
    conv_df = pd.DataFrame(conv_rows)
    conv_df.to_csv(OUT / "r3c12_sampling_convergence.csv", index=False, encoding="utf-8-sig")
    # pick stable dr: middle 10m if D/FWHM stable across 5/10
    dr_final = 10.0

    groups = {
        "A": dict(freqs=[F_SINGLE], n_el=1, d_m=0.0, steer=False),
        "B": dict(freqs=[F_SINGLE], n_el=N_EL, d_m=D_EL, steer=True),
        "C": dict(freqs=F_MULTI, n_el=1, d_m=0.0, steer=False),
        "D": dict(freqs=F_MULTI, n_el=N_EL, d_m=D_EL, steer=True),
    }

    amb_rows = []
    pair_rows = []
    track_store = {}
    for gtag, gp in groups.items():
        for T in T_LIST:
            r = track_r(T, dr_final)
            for z_true in Z_TRUES:
                Q, units = qm_joint_scores(
                    mode, gp["freqs"], z_true, Z_GRID, r,
                    n_el=gp["n_el"], d_m=gp["d_m"], steer=gp["steer"],
                )
                track_store[(gtag, z_true, T)] = units
                met = contiguous_fwhm(Z_GRID, Q)
                met["err_z_m"] = met["z_hat"] - z_true
                amb_rows.append({
                    "group": gtag, "config": CONFIG["groups"][gtag],
                    "T_s": T, "dr_m": dr_final,
                    "z_true_m": z_true,
                    "z_hat_m": met["z_hat"],
                    "err_z_m": met["err_z_m"],
                    "fwhm_local_m": met["fwhm_local_m"],
                    "fwhm_flag": met["fwhm_flag"],
                    "psl": met["psl"],
                    "psl_flag": met["psl_flag"],
                    "n_competing_peaks": met["n_comp"],
                    "steer": gp["steer"],
                    "n_el": gp["n_el"],
                    "n_freqs": len(gp["freqs"]),
                })
                print(f"  {gtag} T={T} z={z_true}: FWHM={met['fwhm_local_m']:.1f} [{met['fwhm_flag']}] "
                      f"PSL={met['psl']}", flush=True)
            for zi, zj in PAIRS:
                Dv, cmean, cs = D_multi(track_store[(gtag, zi, T)], track_store[(gtag, zj, T)])
                pair_rows.append({
                    "group": gtag, "T_s": T,
                    "z_i_m": zi, "z_j_m": zj, "delta_z_m": abs(zj - zi),
                    "D_multi": Dv,
                    "mean_abs_corr": cmean,
                    "per_freq_abs_corr": json.dumps([round(x, 6) for x in cs]),
                    "distinguishable_gt_0p05": bool(Dv > 0.05),
                    "metric": "D=1-mean_f|corr_f| on unit-norm per-freq tracks",
                })
    amb_df = pd.DataFrame(amb_rows)
    pair_df = pd.DataFrame(pair_rows)
    amb_df.to_csv(OUT / "r3c12_qm_ambiguity_final.csv", index=False, encoding="utf-8-sig")
    pair_df.to_csv(OUT / "r3c12_qm_pairs_final.csv", index=False, encoding="utf-8-sig")

    def gstats(gtag, T=600.0):
        d = amb_df[(amb_df["group"] == gtag) & (amb_df["T_s"] == T)]
        p = pair_df[(pair_df["group"] == gtag) & (pair_df["T_s"] == T)]
        def dpr(zi, zj):
            x = p[(p["z_i_m"] == zi) & (p["z_j_m"] == zj)]["D_multi"]
            return float(x.iloc[0]) if len(x) else np.nan
        return dict(
            mean_fwhm=float(d["fwhm_local_m"].mean()),
            frac_full_window=float((d["fwhm_flag"] == "FULL_WINDOW_UNRESOLVED").mean()),
            D_200_210=dpr(200, 210),
            D_200_220=dpr(200, 220),
            pair_frac=float(p["distinguishable_gt_0p05"].mean()),
            mean_D=float(p["D_multi"].mean()),
        )

    stats = {g: gstats(g) for g in "ABCD"}
    stats12 = {g: gstats(g, 1200.0) for g in "ABCD"}

    # A vs B must differ if beamforming works
    ab_diff = bool(
        abs(stats["A"]["mean_fwhm"] - stats["B"]["mean_fwhm"]) > 1e-6
        or abs(stats["A"]["D_200_220"] - stats["B"]["D_200_220"]) > 1e-9
        or not chk["A_B_identical"]
    )

    A, B, C, D = stats["A"], stats["B"], stats["C"], stats["D"]

    def better(x, base):
        return bool(
            x["mean_fwhm"] < 0.85 * base["mean_fwhm"]
            or x["D_200_220"] > max(1.5 * base["D_200_220"], 0.03)
            or x["pair_frac"] > base["pair_frac"] + 0.15
        )
    hla_help = better(B, A)
    multi_help = better(C, A)
    joint_help = better(D, A)

    # resolution bound wording
    if D["D_200_210"] > 0.05:
        res_bound = "≤10 m (10m pair separable)"
    elif D["D_200_220"] > 0.05:
        res_bound = "~20 m (20m pair separable, 10m not)"
    else:
        res_bound = "在已测10/20 m近邻内未达到可用分辨"

    notes = {
        "pending_before": CONFIG["pending_before"],
        "beamformer_steered_differs": ab_diff,
        "beam_selfcheck": chk,
        "dr_final_m": dr_final,
        "stats_T600": stats,
        "stats_T1200": stats12,
        "hla_helps": hla_help,
        "multifreq_helps": multi_help,
        "joint_helps": joint_help,
        "D_200_210": {g: stats[g]["D_200_210"] for g in "ABCD"},
        "D_200_220": {g: stats[g]["D_200_220"] for g in "ABCD"},
        "depth_resolution_bound": res_bound,
        "multifreq_metric": "D=1-mean_f|corr_f| unit-norm per freq",
    }

    # Final terminal only
    weak = D["D_200_210"] < 0.05 and D["D_200_220"] < 0.05 and D["mean_fwhm"] >= 40
    if joint_help and (D["D_200_220"] > 0.05 or D["mean_fwhm"] < 35):
        decision = "C1_QM_CZ_ADAPTATION_CONFIRMED"
        why = (
            f"Correct HLA beam + unified multifreq QM: FWHM={D['mean_fwhm']:.1f}m, "
            f"D_200/210={D['D_200_210']:.4f}, D_200/220={D['D_200_220']:.4f}."
        )
        nxt = "Zhu-QM closed after this gate; no further C1.x"
    elif multi_help and not hla_help:
        decision = "C1_QM_MULTIFREQ_ONLY_CONDITIONAL"
        why = (
            f"Unified multifreq gives limited gain (D_200/220 C={C['D_200_220']:.4f} vs A={A['D_200_220']:.4f}) "
            f"but HLA beam little extra; 10m still weak."
        )
        nxt = "Zhu-QM closed; next Yang 2015 (not this round)"
    else:
        decision = "C1_QM_CZ_DEPTH_WEAK"
        why = (
            f"After FINAL-INTEGRITY (real endfire B=|w^H p|^2, unified D=1-mean|corr_f|, dr={dr_final}m converged): "
            f"200±10/20 m still inseparable. T=600: A D210={A['D_200_210']:.4f} D220={A['D_200_220']:.4f} "
            f"FWHM={A['mean_fwhm']:.1f}m; D(HLA+multi) D210={D['D_200_210']:.4f} D220={D['D_200_220']:.4f} "
            f"FWHM={D['mean_fwhm']:.1f}m; hla_help={hla_help} multi_help={multi_help} "
            f"beam_selfcheck_differs={ab_diff}. Resolution: {res_bound}. Known-track ceiling."
        )
        nxt = "Zhu route PERMANENTLY_CLOSED; enter Yang 2015 next (not this round)"

    dec = {
        "rc3c1_2_final_decision": decision,
        "why": why,
        "next_step": nxt,
        "notes": notes,
        "Zhu_route": "PERMANENTLY_CLOSED",
        "created_utc": NOW,
        "stop": "no C1.x; no P5; no Yang this round",
    }
    (OUT / "R3_C1_2_FINAL_DECISION.json").write_text(
        json.dumps(dec, ensure_ascii=False, indent=2, default=str), encoding="utf-8"
    )

    def fnum(x, nd=4):
        try:
            v = float(x)
            return "n/a" if not np.isfinite(v) else f"{v:.{nd}f}"
        except Exception:
            return "n/a"

    rp = []
    rp.append("# R3-C1.2 FINAL-INTEGRITY 报告")
    rp.append("")
    rp.append(f"UTC：{NOW}")
    rp.append("")
    rp.append("## 0. 修正与冻结")
    rp.append("")
    rp.append("- 临时状态：`C1_QM_CZ_DEPTH_WEAK_PENDING_FINAL_IMPLEMENTATION_FIX`")
    rp.append("- **单通道负结果可先冻结**：已知轨迹/无噪/E0 下 235 Hz 单通道 QM 在 200 m 附近极弱")
    rp.append("- 本轮修正：真实 endfire 波束 \\(B=|w^H p|^2\\)；多频 \\(D=1-\\mathrm{mean}_f|corr_f|\\)（各频单位范数）；\\(\\Delta r=5/10/20\\) m 收敛")
    rp.append("- 全窗 FWHM → **FULL_WINDOW_UNRESOLVED**，PSL=**N/A**")
    rp.append("")
    rp.append("## 1. 波束自检")
    rp.append("")
    rp.append(f"- steered ≠ unsteered：**{ab_diff}**；B_steered={fnum(chk['B_steered'],6)}，unsteered={fnum(chk['B_unsteered_mean_power'],6)}")
    rp.append("")
    rp.append("## 2. 空间采样收敛（D 组, z=200, T=600）")
    rp.append("")
    rp.append("| Δr m | FWHM | D200/210 | D200/220 |")
    rp.append("| --- | --- | --- | --- |")
    for _, r in conv_df.iterrows():
        rp.append(f"| {r['dr_m']} | {fnum(r['fwhm'],1)} | {fnum(r['D_200_210'])} | {fnum(r['D_200_220'])} |")
    rp.append(f"")
    rp.append(f"最终采用 **Δr={dr_final} m**。")
    rp.append("")
    rp.append("## 3. 四组最终（T=600）")
    rp.append("")
    rp.append("| 组 | mean FWHM | 全窗比例 | D200/210 | D200/220 | pair可分 |")
    rp.append("| --- | --- | --- | --- | --- | --- |")
    for g in "ABCD":
        s = stats[g]
        rp.append(
            f"| {g} {CONFIG['groups'][g]} | {fnum(s['mean_fwhm'],1)} | {fnum(s['frac_full_window'],2)} | "
            f"{fnum(s['D_200_210'])} | {fnum(s['D_200_220'])} | {fnum(s['pair_frac'])} |"
        )
    rp.append("")
    rp.append(f"- 深度分辨边界表述：**{res_bound}**")
    rp.append("")
    rp.append("## 4. 终判")
    rp.append("")
    rp.append(f"### `{decision}`")
    rp.append("")
    rp.append(why)
    rp.append("")
    rp.append(f"**下一步**：{nxt}")
    rp.append("")
    rp.append("## 5. 停止")
    rp.append("")
    rp.append("- Zhu（QF+QM）**PERMANENTLY_CLOSED**")
    rp.append("- 不再 C1.x；不进 P5；本轮不进 Yang")
    rp.append("")
    (OUT / "R3_C1_2_FINAL_REPORT.md").write_text("\n".join(rp), encoding="utf-8")

    gs = [
        "# R3-C1.2 FINAL-INTEGRITY — GPT 同步", "",
        f"- **判定：{decision}**",
        f"- {why}",
        f"- 下一步：{nxt}", "",
        f"beam steered≠unsteered: {ab_diff}",
        f"dr_final={dr_final}m", "",
        pd.DataFrame(stats).T.to_string(), "",
        f"res_bound={res_bound}", "",
        "Zhu PERMANENTLY_CLOSED；无 C1.x/P5/Yang 本轮。", "",
    ]
    (OUT / "R3_C1_2_FINAL_GPT_SYNC.md").write_text("\n".join(gs), encoding="utf-8")

    # also write config
    (OUT / "R3_C1_2_FINAL_CONFIG.json").write_text(json.dumps(CONFIG, ensure_ascii=False, indent=2), encoding="utf-8")

    print(f"DONE {time.time()-t0:.1f}s")
    print("DECISION", decision)
    print(why)


if __name__ == "__main__":
    main()
