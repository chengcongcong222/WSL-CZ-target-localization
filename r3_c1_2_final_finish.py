#!/usr/bin/env python3
"""R3-C1.2 FINAL-INTEGRITY finisher: cached modes, four groups, decision."""
from __future__ import annotations

import json
import math
import time
from datetime import datetime, timezone
from pathlib import Path

import numpy as np
import pandas as pd
import p4_g4_performance_boundary as p4

OUT = Path(__file__).resolve().parent / "results" / "R3_C1_depth_motion"
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
DR_FINAL = 10.0

MODE = p4.MODE_ENV["E0"]
_MODES = {}
_PHI = {}


def get_modes(f):
    key = round(float(f), 3)
    if key not in _MODES:
        _MODES[key] = MODE.modes(float(f))
    return _MODES[key]


def phi_at(f, z):
    key = (round(float(f), 3), round(float(z), 2))
    if key not in _PHI:
        ms = get_modes(f)
        _PHI[key] = [float(np.interp(z, MODE.z, phi)) for _, phi in ms]
    return _PHI[key]


def track_r(T, dr=DR_FINAL):
    r1 = min(R1, R0 + U_REL * T)
    n = int(max(20, math.floor((r1 - R0) / dr) + 1))
    return np.linspace(R0, r1, n)


def element_p(f, z_s, ri, xm, z_r=Z_R):
    rm = max(abs(ri - xm), 1.0)
    ms = get_modes(f)
    phis = phi_at(f, z_s)
    phir = phi_at(f, z_r)
    acc = 0j
    for (kr, _), a_s, a_r in zip(ms, phis, phir):
        att = math.exp(-2e-5 * (f / 200.0) * rm / 1000.0)
        acc += (a_s * a_r / math.sqrt(kr * rm)) * np.exp(1j * kr * rm) * att
    return acc


def beam_track(f, z_s, r, n_el=1, d_m=0.0, steer=True, z_r=Z_R):
    xs = (np.arange(n_el) - (n_el - 1) / 2.0) * d_m if n_el > 1 else np.zeros(1)
    k0 = 2 * np.pi * f / C0
    B = np.zeros(len(r))
    for i, ri in enumerate(r):
        P = np.array([element_p(f, z_s, ri, xm, z_r) for xm in xs])
        if n_el == 1:
            Bi = abs(P[0]) ** 2
        elif steer:
            w = np.exp(1j * k0 * xs * 1.0) / n_el
            Bi = abs(np.dot(np.conj(w), P)) ** 2
        else:
            Bi = float(np.mean(np.abs(P) ** 2))
        R = math.sqrt(ri ** 2 + (z_r - z_s) ** 2)
        B[i] = Bi * R ** 2
    return B


def norm_track(I):
    v = I - I.mean()
    return v / (np.linalg.norm(v) + 1e-30)


def qm_scores(freqs, z_true, z_grid, r, n_el, d_m, steer):
    obs = {}
    Qs = []
    for f in freqs:
        u0 = norm_track(beam_track(f, z_true, r, n_el, d_m, steer))
        obs[f] = u0
        Qf = np.array([abs(float(np.dot(u0, norm_track(beam_track(f, z, r, n_el, d_m, steer))))) for z in z_grid])
        Qs.append(Qf)
    return np.mean(Qs, axis=0), obs


def D_multi(u1, u2):
    cs = [abs(float(np.dot(u1[f], u2[f]))) for f in u1]
    return 1.0 - float(np.mean(cs)), float(np.mean(cs)), cs


def fwhm_flags(z_grid, J):
    J = np.asarray(J, float)
    window = float(z_grid[-1] - z_grid[0])
    imax = int(np.argmax(J))
    jmax = float(J[imax])
    thr = 0.5 * jmax
    lo = hi = imax
    while lo > 0 and J[lo - 1] >= thr:
        lo -= 1
    while hi < len(J) - 1 and J[hi + 1] >= thr:
        hi += 1
    fwhm = float(z_grid[hi] - z_grid[lo])
    if fwhm >= window * 0.95:
        return dict(z_hat=float(z_grid[imax]), fwhm=fwhm, fwhm_flag="FULL_WINDOW_UNRESOLVED", psl=np.nan, psl_flag="N/A")
    mask = np.ones_like(J, bool)
    mask[lo:hi + 1] = False
    psl = float(J[mask].max() / (jmax + 1e-30)) if mask.any() else 0.0
    return dict(z_hat=float(z_grid[imax]), fwhm=fwhm, fwhm_flag="local", psl=psl, psl_flag="ok")


def main():
    t0 = time.time()
    # beam self-check
    xs = (np.arange(N_EL) - (N_EL - 1) / 2.0) * D_EL
    P = np.array([element_p(F_SINGLE, 200.0, 50e3, xm) for xm in xs])
    k0 = 2 * np.pi * F_SINGLE / C0
    w = np.exp(1j * k0 * xs) / N_EL
    B_st = abs(np.dot(np.conj(w), P)) ** 2
    B_un = float(np.mean(np.abs(P) ** 2))
    ab_diff = abs(B_st - B_un) > 1e-20
    print("beam", B_st, B_un, "differs", ab_diff, flush=True)
    pd.DataFrame([{
        "x_m": json.dumps(xs.tolist()),
        "phase_rad": json.dumps(np.angle(P).tolist()),
        "B_steered": B_st, "B_unsteered": B_un,
        "steered_differs_from_unsteered": ab_diff,
    }]).to_csv(OUT / "r3c12_beamformer_check.csv", index=False, encoding="utf-8-sig")

    groups = {
        "A": dict(freqs=[F_SINGLE], n_el=1, d_m=0.0, steer=False),
        "B": dict(freqs=[F_SINGLE], n_el=N_EL, d_m=D_EL, steer=True),
        "C": dict(freqs=F_MULTI, n_el=1, d_m=0.0, steer=False),
        "D": dict(freqs=F_MULTI, n_el=N_EL, d_m=D_EL, steer=True),
    }

    amb_rows, pair_rows = [], []
    cache = {}
    for gtag, gp in groups.items():
        for T in T_LIST:
            r = track_r(T)
            for z in Z_TRUES:
                Q, units = qm_scores(gp["freqs"], z, Z_GRID, r, gp["n_el"], gp["d_m"], gp["steer"])
                cache[(gtag, z, T)] = units
                met = fwhm_flags(Z_GRID, Q)
                amb_rows.append({
                    "group": gtag, "T_s": T, "z_true_m": z,
                    "z_hat_m": met["z_hat"], "err_z_m": met["z_hat"] - z,
                    "fwhm_local_m": met["fwhm"], "fwhm_flag": met["fwhm_flag"],
                    "psl": met["psl"], "psl_flag": met["psl_flag"],
                    "steer": gp["steer"], "n_el": gp["n_el"], "n_freqs": len(gp["freqs"]),
                    "dr_m": DR_FINAL,
                })
                print(f"  {gtag} T={T} z={z}: FWHM={met['fwhm']:.1f} {met['fwhm_flag']}", flush=True)
            for zi, zj in PAIRS:
                Dv, cm, cs = D_multi(cache[(gtag, zi, T)], cache[(gtag, zj, T)])
                pair_rows.append({
                    "group": gtag, "T_s": T, "z_i_m": zi, "z_j_m": zj,
                    "delta_z_m": abs(zj - zi), "D_multi": Dv, "mean_abs_corr": cm,
                    "per_freq_abs_corr": json.dumps([round(x, 6) for x in cs]),
                    "distinguishable_gt_0p05": bool(Dv > 0.05),
                    "metric": "D=1-mean_f|corr_f|",
                })
                print(f"  pair {gtag} T={T} {zi}/{zj}: D={Dv:.5f}", flush=True)

    amb_df = pd.DataFrame(amb_rows)
    pair_df = pd.DataFrame(pair_rows)
    amb_df.to_csv(OUT / "r3c12_qm_ambiguity_final.csv", index=False, encoding="utf-8-sig")
    pair_df.to_csv(OUT / "r3c12_qm_pairs_final.csv", index=False, encoding="utf-8-sig")

    # sampling convergence quick (D group z=200 T=600) — optional small
    conv_rows = []
    for dr in [5.0, 10.0, 20.0]:
        r = track_r(600.0, dr)
        Q, u200 = qm_scores(F_MULTI, 200.0, Z_GRID, r, N_EL, D_EL, True)
        _, u210 = qm_scores(F_MULTI, 210.0, Z_GRID, r, N_EL, D_EL, True)
        _, u220 = qm_scores(F_MULTI, 220.0, Z_GRID, r, N_EL, D_EL, True)
        met = fwhm_flags(Z_GRID, Q)
        D210 = D_multi(u200, u210)[0]
        D220 = D_multi(u200, u220)[0]
        conv_rows.append({"dr_m": dr, "n_r": len(r), "fwhm": met["fwhm"],
                          "fwhm_flag": met["fwhm_flag"], "D_200_210": D210, "D_200_220": D220})
        print("conv", dr, met["fwhm"], D210, D220, flush=True)
    pd.DataFrame(conv_rows).to_csv(OUT / "r3c12_sampling_convergence.csv", index=False, encoding="utf-8-sig")

    def gstats(gtag, T):
        d = amb_df[(amb_df.group == gtag) & (amb_df.T_s == T)]
        p = pair_df[(pair_df.group == gtag) & (pair_df.T_s == T)]

        def dpr(zi, zj):
            x = p[(p.z_i_m == zi) & (p.z_j_m == zj)].D_multi
            return float(x.iloc[0]) if len(x) else np.nan
        return dict(
            mean_fwhm=float(d.fwhm_local_m.mean()),
            frac_full=float((d.fwhm_flag == "FULL_WINDOW_UNRESOLVED").mean()),
            D_200_210=dpr(200, 210), D_200_220=dpr(200, 220),
            pair_frac=float(p.distinguishable_gt_0p05.mean()),
            mean_D=float(p.D_multi.mean()),
        )

    stats = {g: gstats(g, 600.0) for g in "ABCD"}
    stats12 = {g: gstats(g, 1200.0) for g in "ABCD"}
    A, B, C, D = stats["A"], stats["B"], stats["C"], stats["D"]

    def better(x, base):
        return bool(x["mean_fwhm"] < 0.85 * base["mean_fwhm"]
                    or x["D_200_220"] > max(1.5 * base["D_200_220"], 0.03)
                    or x["pair_frac"] > base["pair_frac"] + 0.15)

    hla_help, multi_help = better(B, A), better(C, A)
    if D["D_200_210"] > 0.05:
        res = "≤10 m"
    elif D["D_200_220"] > 0.05:
        res = "~20 m"
    else:
        res = "在已测10/20 m近邻内未达到可用分辨"

    weak = D["D_200_210"] < 0.05 and D["D_200_220"] < 0.05 and D["mean_fwhm"] >= 40
    if (not weak) and (D["D_200_220"] > 0.05 or D["mean_fwhm"] < 35):
        decision = "C1_QM_CZ_ADAPTATION_CONFIRMED"
        why = f"Joint QM usable: FWHM={D['mean_fwhm']:.1f} D210={D['D_200_210']:.4f} D220={D['D_200_220']:.4f}."
    elif multi_help and not hla_help:
        decision = "C1_QM_MULTIFREQ_ONLY_CONDITIONAL"
        why = f"Multifreq-only limited gain; HLA little extra. D220 C={C['D_200_220']:.4f} A={A['D_200_220']:.4f}."
    else:
        decision = "C1_QM_CZ_DEPTH_WEAK"
        why = (
            f"FINAL-INTEGRITY: real endfire B=|w^H p|^2 (differs={ab_diff}), unified D=1-mean|corr_f|, "
            f"dr={DR_FINAL}m; 200±10/20 m still weak. T=600 A D210={A['D_200_210']:.4f} D220={A['D_200_220']:.4f} "
            f"FWHM={A['mean_fwhm']:.1f}; D D210={D['D_200_210']:.4f} D220={D['D_200_220']:.4f} FWHM={D['mean_fwhm']:.1f}; "
            f"hla_help={hla_help} multi_help={multi_help}. {res}. Known-track ceiling."
        )
    nxt = "Zhu PERMANENTLY_CLOSED; Yang 2015 next (not this round)"

    dec = {
        "rc3c1_2_final_decision": decision, "why": why, "next_step": nxt,
        "notes": {
            "beam_steered_differs": ab_diff, "stats_T600": stats, "stats_T1200": stats12,
            "hla_helps": hla_help, "multifreq_helps": multi_help,
            "D_200_210": {g: stats[g]["D_200_210"] for g in "ABCD"},
            "D_200_220": {g: stats[g]["D_200_220"] for g in "ABCD"},
            "depth_resolution_bound": res, "dr_final_m": DR_FINAL,
            "metric": "D=1-mean_f|corr_f| unit-norm per freq",
        },
        "Zhu_route": "PERMANENTLY_CLOSED", "created_utc": NOW,
    }
    (OUT / "R3_C1_2_FINAL_DECISION.json").write_text(json.dumps(dec, ensure_ascii=False, indent=2, default=str), encoding="utf-8")

    rp = []
    rp.append("# R3-C1.2 FINAL-INTEGRITY 报告")
    rp.append("")
    rp.append(f"UTC：{NOW}")
    rp.append("")
    rp.append("修正：真实 endfire 波束；统一多频 D；固定 Δr 采样；FULL_WINDOW 标注。")
    rp.append("")
    rp.append(f"- 波束 steered≠unsteered：**{ab_diff}**")
    rp.append("")
    rp.append("## T=600 四组")
    rp.append("")
    rp.append("| 组 | FWHM | 全窗 | D200/210 | D200/220 |")
    rp.append("| --- | --- | --- | --- | --- |")
    for g in "ABCD":
        s = stats[g]
        rp.append(f"| {g} | {s['mean_fwhm']:.1f} | {s['frac_full']:.2f} | {s['D_200_210']:.5f} | {s['D_200_220']:.5f} |")
    rp.append("")
    rp.append(f"分辨边界：**{res}**")
    rp.append("")
    rp.append(f"### `{decision}`")
    rp.append("")
    rp.append(why)
    rp.append("")
    rp.append(f"下一步：{nxt}")
    rp.append("")
    rp.append("Zhu **PERMANENTLY_CLOSED**；无 C1.x / P5 / Yang 本轮。")
    (OUT / "R3_C1_2_FINAL_REPORT.md").write_text("\n".join(rp), encoding="utf-8")
    (OUT / "R3_C1_2_FINAL_GPT_SYNC.md").write_text(
        f"# R3-C1.2 FINAL\n\n**{decision}**\n\n{why}\n\nnext: {nxt}\n\nT600:\n{pd.DataFrame(stats).T}\n\nres={res}\n",
        encoding="utf-8",
    )

    print("DECISION", decision)
    print(why)
    print(f"DONE {time.time()-t0:.1f}s")


if __name__ == "__main__":
    main()
