#!/usr/bin/env python3
"""R3-C1.2: Zhu-QM limited CZ adaptation — single vs HLA beam vs multifreq.

QF direct transfer CLOSED. Known track (optimistic ceiling). Then permanently close Zhu.
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
FIG = OUT / "figures"
OUT.mkdir(parents=True, exist_ok=True)
FIG.mkdir(parents=True, exist_ok=True)
NOW = datetime.now(timezone.utc).isoformat()

C0 = 1500.0
Z_R = 200.0
R0, R1 = 45e3, 60e3
U_REL = 2.0
T_LIST = [600.0, 1200.0]
Z_TRUES = [180.0, 190.0, 200.0, 210.0, 220.0]
Z_GRID = np.linspace(160.0, 240.0, 81)  # local prior ~200±30 m search
F_SINGLE = 235.0
F_MULTI = [201.0, 235.0, 283.0, 338.0]
PAIRS = [(180, 190), (190, 200), (200, 210), (210, 220), (180, 200), (200, 220)]
N_EL = 8
D_EL = 2.0  # 14 m endfire aperture

CONFIG = {
    "package": "R3_C1_2_qm_adaptation",
    "created_utc": NOW,
    "QF_direct_transfer": "CLOSED (s=sin(arctan|zr-zs|/r)≈0 at 50km/200m — not Zhu CZ variable)",
    "QM_scope": "matched intensity structure; CZ modal I' tracks; known track = optimistic ceiling",
    "groups": {
        "A": "single + 235 Hz",
        "B": "HLA 14m endfire beamformed + 235 Hz",
        "C": "single + multi-line {201,235,283,338} Hz",
        "D": "HLA 14m + multi-line",
    },
    "metrics": {
        "D": "1-|corr|",
        "fwhm": "contiguous local FWHM around global max",
        "pairs_m": [list(p) for p in PAIRS],
    },
    "stop": ["no C1.x after this", "no Yang this round", "no P5", "no new features"],
}


def track_r(T):
    r1 = min(R1, R0 + U_REL * T)
    return np.linspace(R0, r1, 80)


def beam_track_intensity(mode, f, z_s, r, z_r=Z_R, n_el=1, d_m=0.0, steer=True):
    """Track intensity after endfire array beam toward known track bearing (RC2-aided).

    endfire: elements along motion axis x; r_m = |r - x_m|.
    steer: coherent sum with 0 phase (target on array axis endfire).
    """
    I = np.zeros(len(r))
    xs = (np.arange(n_el) - (n_el - 1) / 2.0) * d_m if n_el > 1 else np.array([0.0])
    ms = mode.modes(float(f))
    for i, ri in enumerate(r):
        acc = 0j
        for xm in xs:
            rm = max(abs(ri - xm), 1.0)
            p = 0j
            for kr, phi in ms:
                a_s = float(np.interp(z_s, mode.z, phi))
                a_r = float(np.interp(z_r, mode.z, phi))
                att = math.exp(-2e-5 * (f / 200.0) * rm / 1000.0)
                # endfire steering: plane-wave phase along x ~ kr * xm * (target direction)
                # target on +x from array: extra phase kr*xm for source at range ri from center
                p += (a_s * a_r / math.sqrt(kr * rm)) * np.exp(1j * kr * rm)
            acc += p  # steer 0 for endfire on-axis
        I[i] = abs(acc / n_el) ** 2
    # Zhu-form range compensation
    R = np.sqrt(r ** 2 + (z_r - z_s) ** 2)
    Iprime = I * R ** 2
    return r, Iprime


def qm_corr(I_obs, I_pred):
    a = I_obs - I_obs.mean()
    b = I_pred - I_pred.mean()
    return float(np.dot(a, b) / (np.linalg.norm(a) * np.linalg.norm(b) + 1e-30))


def qm_profile(mode, f, z_true, z_grid, T, n_el, d_m):
    r = track_r(T)
    _, Iobs = beam_track_intensity(mode, f, z_true, r, n_el=n_el, d_m=d_m)
    Q = np.zeros(len(z_grid))
    for iz, z in enumerate(z_grid):
        _, Ip = beam_track_intensity(mode, f, z, r, n_el=n_el, d_m=d_m)
        Q[iz] = abs(qm_corr(Iobs, Ip))
    return Q, Iobs


def contiguous_fwhm(z_grid, J):
    J = np.asarray(J, float)
    imax = int(np.argmax(J))
    jmax = float(J[imax])
    thr = 0.5 * jmax
    lo = imax
    while lo > 0 and J[lo - 1] >= thr:
        lo -= 1
    hi = imax
    while hi < len(J) - 1 and J[hi + 1] >= thr:
        hi += 1
    fwhm = float(z_grid[hi] - z_grid[lo])
    mask = np.ones_like(J, dtype=bool)
    mask[lo:hi + 1] = False
    psl = float(J[mask].max() / (jmax + 1e-30)) if mask.any() else 0.0
    n_comp = sum(
        1 for i in range(1, len(J) - 1)
        if i != imax and J[i] >= J[i - 1] and J[i] >= J[i + 1] and J[i] > 0.7 * jmax
    )
    return dict(z_hat=float(z_grid[imax]), err=None, fwhm_local_m=fwhm, psl=psl, n_comp=n_comp, peak=jmax)


def pair_D(I_a, I_b):
    return 1.0 - abs(qm_corr(I_a, I_b))


def run_group(mode, tag, z_true, T, freqs, n_el, d_m):
    """Return ambiguity metrics + track intensities per freq for pair D."""
    Q_acc = np.zeros(len(Z_GRID))
    tracks = {}
    for f in freqs:
        Qf, Iobs = qm_profile(mode, f, z_true, Z_GRID, T, n_el, d_m)
        Q_acc += Qf
        tracks[f] = Iobs
    Q_acc /= len(freqs)
    met = contiguous_fwhm(Z_GRID, Q_acc)
    met["err_z_m"] = met["z_hat"] - z_true
    return met, tracks, Q_acc


def main():
    t0 = time.time()
    mode = p4.MODE_ENV["E0"]
    print("=== R3-C1.2 QM adaptation ===", flush=True)
    (OUT / "R3_C1_2_CONFIG.json").write_text(json.dumps(CONFIG, ensure_ascii=False, indent=2), encoding="utf-8")

    groups = {
        "A": dict(freqs=[F_SINGLE], n_el=1, d_m=0.0),
        "B": dict(freqs=[F_SINGLE], n_el=N_EL, d_m=D_EL),
        "C": dict(freqs=F_MULTI, n_el=1, d_m=0.0),
        "D": dict(freqs=F_MULTI, n_el=N_EL, d_m=D_EL),
    }

    rows_amb = []
    # cache tracks for pairs at T=600
    cache_tracks = {}  # (group, z, T) -> mean-stacked track (concat freqs)

    for gtag, gp in groups.items():
        for T in T_LIST:
            for z_true in Z_TRUES:
                met, tracks, Q = run_group(mode, gtag, z_true, T, gp["freqs"], gp["n_el"], gp["d_m"])
                cache_tracks[(gtag, z_true, T)] = tracks
                rows_amb.append({
                    "group": gtag,
                    "config": CONFIG["groups"][gtag],
                    "T_s": T,
                    "z_true_m": z_true,
                    "z_hat_m": met["z_hat"],
                    "err_z_m": met["err_z_m"],
                    "fwhm_local_m": met["fwhm_local_m"],
                    "psl": met["psl"],
                    "n_competing_peaks": met["n_comp"],
                    "n_freqs": len(gp["freqs"]),
                    "n_el": gp["n_el"],
                })
                print(f"  {gtag} T={T} z={z_true}: hat={met['z_hat']:.1f} err={met['err_z_m']:.1f} "
                      f"FWHM={met['fwhm_local_m']:.1f} PSL={met['psl']:.2f}", flush=True)

    amb_df = pd.DataFrame(rows_amb)
    amb_df.to_csv(OUT / "r3c12_qm_ambiguity.csv", index=False, encoding="utf-8-sig")

    # pair D=1-|corr| on stacked multi-freq / single track signatures
    pair_rows = []
    for gtag, gp in groups.items():
        for T in T_LIST:
            for zi, zj in PAIRS:
                ti = cache_tracks[(gtag, zi, T)]
                tj = cache_tracks[(gtag, zj, T)]
                # stack freqs into one vector
                ai = np.concatenate([ti[f] for f in gp["freqs"]])
                aj = np.concatenate([tj[f] for f in gp["freqs"]])
                D = pair_D(ai, aj)
                pair_rows.append({
                    "group": gtag,
                    "T_s": T,
                    "z_i_m": zi,
                    "z_j_m": zj,
                    "delta_z_m": abs(zj - zi),
                    "D_abs": D,
                    "distinguishable_gt_0p05": bool(D > 0.05),
                    "metric": "D=1-|corr| on stacked track intensity",
                })
    pair_df = pd.DataFrame(pair_rows)
    pair_df.to_csv(OUT / "r3c12_qm_pairs.csv", index=False, encoding="utf-8-sig")

    # summary vs group A
    def group_stats(gtag, T=600.0):
        d = amb_df[(amb_df["group"] == gtag) & (amb_df["T_s"] == T)]
        p = pair_df[(pair_df["group"] == gtag) & (pair_df["T_s"] == T)]
        return dict(
            mean_fwhm=float(d["fwhm_local_m"].mean()),
            mean_psl=float(d["psl"].mean()),
            mean_abs_err=float(d["err_z_m"].abs().mean()),
            pair_D_mean=float(p["D_abs"].mean()),
            pair_frac=float(p["distinguishable_gt_0p05"].mean()),
            D_200_210=float(p[(p["z_i_m"] == 200) & (p["z_j_m"] == 210)]["D_abs"].iloc[0]) if len(p[(p["z_i_m"]==200)&(p["z_j_m"]==210)]) else np.nan,
            D_200_220=float(p[(p["z_i_m"] == 200) & (p["z_j_m"] == 220)]["D_abs"].iloc[0]) if len(p[(p["z_i_m"]==200)&(p["z_j_m"]==220)]) else np.nan,
        )

    stats = {g: group_stats(g) for g in ["A", "B", "C", "D"]}
    # also T=1200
    stats1200 = {g: group_stats(g, 1200.0) for g in ["A", "B", "C", "D"]}

    # core judgment
    A, B, C, D = stats["A"], stats["B"], stats["C"], stats["D"]
    # significant improvement vs A: FWHM shrink >=25% or D_200_220*2 or pair_frac +0.2
    def better(x, base):
        return bool(
            x["mean_fwhm"] < 0.75 * base["mean_fwhm"]
            or x["D_200_220"] > max(2 * base["D_200_220"], 0.05)
            or x["pair_frac"] > base["pair_frac"] + 0.2
        )
    hla_help = better(B, A)
    multi_help = better(C, A)
    joint_help = better(D, A)
    # 200m usable peak: FWHM < 40 and D_200_210 > 0.05 or D_200_220 > 0.08
    d200 = amb_df[(amb_df["group"] == "D") & (amb_df["T_s"] == 600) & (amb_df["z_true_m"] == 200)]
    fwhm200 = float(d200["fwhm_local_m"].iloc[0]) if len(d200) else np.nan
    usable200 = bool(np.isfinite(fwhm200) and fwhm200 < 40 and (D["D_200_210"] > 0.05 or D["D_200_220"] > 0.08))

    notes = {
        "QF": "CLOSED",
        "stats_T600": stats,
        "stats_T1200": stats1200,
        "hla_helps": hla_help,
        "multifreq_helps": multi_help,
        "joint_helps": joint_help,
        "fwhm_200m_D_T600": fwhm200,
        "D_200_210_D": D["D_200_210"],
        "D_200_220_D": D["D_200_220"],
        "depth_resolution_bound_m": (
            10.0 if D["D_200_210"] > 0.05 else (20.0 if D["D_200_220"] > 0.05 else None)
        ),
        "known_track_ceiling": True,
    }

    if usable200 or (joint_help and D["mean_fwhm"] < 45 and D["pair_frac"] >= 0.5):
        decision = "C1_QM_CZ_ADAPTATION_CONFIRMED"
        why = (
            f"Multifreq/HLA QM yields usable ~200m depth response: FWHM_200={fwhm200:.1f} m, "
            f"D_200/210={D['D_200_210']:.3f}, D_200/220={D['D_200_220']:.3f}, "
            f"pair_frac={D['pair_frac']:.2f} vs A={A['pair_frac']:.2f}."
        )
        nxt = "Zhu-QM closed after this gate either way; if confirmed only as limited adaptation, no further C1.x"
    elif multi_help and not hla_help:
        decision = "C1_QM_MULTIFREQ_ONLY_CONDITIONAL"
        why = (
            f"Multifreq joint improves somewhat (C mean FWHM={C['mean_fwhm']:.1f} vs A={A['mean_fwhm']:.1f}, "
            f"D_200/220={C['D_200_220']:.3f} vs {A['D_200_220']:.3f}) but 10m pairs still weak "
            f"(D_200/210={C['D_200_210']:.3f}). HLA beam adds little ({B})."
        )
        nxt = "permanently close Zhu route; next RC3-C = Yang 2015 (not this round)"
    else:
        decision = "C1_QM_CZ_DEPTH_WEAK"
        why = (
            f"Four QM groups fail to separate 200±10/20 m clearly on first CZ. "
            f"A: FWHM={A['mean_fwhm']:.1f}m D200/210={A['D_200_210']:.4f} D200/220={A['D_200_220']:.4f}; "
            f"D joint: FWHM={D['mean_fwhm']:.1f}m D200/210={D['D_200_210']:.4f} D200/220={D['D_200_220']:.4f}; "
            f"HLA help={hla_help}, multifreq help={multi_help}. Known-track ceiling already optimistic."
        )
        nxt = "permanently close Zhu-QM route; proceed to Yang 2015 modal synthetic aperture next (not this round)"

    dec = {
        "rc3c1_2_decision": decision,
        "why": why,
        "next_step": nxt,
        "notes": notes,
        "created_utc": NOW,
        "Zhu_route": "PERMANENTLY_CLOSED after this gate",
        "QF": "CLOSED (direct transfer)",
        "stop": "no C1.x, no Yang this round, no P5",
        "b1": "PERMANENTLY_CLOSED",
    }
    (OUT / "R3_C1_2_DECISION.json").write_text(json.dumps(dec, ensure_ascii=False, indent=2, default=str), encoding="utf-8")

    def fnum(x, nd=3):
        try:
            v = float(x)
            return "n/a" if not np.isfinite(v) else f"{v:.{nd}f}"
        except Exception:
            return "n/a"

    # report
    rp = []
    rp.append("# R3-C1.2 报告：Zhu-QM 第一 CZ 有限适配终判")
    rp.append("")
    rp.append(f"UTC：{NOW}")
    rp.append("")
    rp.append("- **QF 直接迁移：CLOSED**（CZ 中 \\(s=\\sin\\theta\\) 几何变量失效，不改造成新方法）")
    rp.append("- **QM matched intensity**：保留做本轮有限适配（HLA 波束 + 多频联合）")
    rp.append("- 已知真实轨迹 = **理论上限**（RC2 误差只会更难）")
    rp.append("")
    rp.append("## 四组（T=600 s 摘要）")
    rp.append("")
    rp.append("| 组 | mean FWHM m | mean PSL | D 200/210 | D 200/220 | pair D>0.05 |")
    rp.append("| --- | --- | --- | --- | --- | --- |")
    for g in ["A", "B", "C", "D"]:
        s = stats[g]
        rp.append(
            f"| {g} {CONFIG['groups'][g]} | {fnum(s['mean_fwhm'],1)} | {fnum(s['mean_psl'])} | "
            f"{fnum(s['D_200_210'])} | {fnum(s['D_200_220'])} | {fnum(s['pair_frac'])} |"
        )
    rp.append("")
    rp.append("## 深度对 D=1−|corr|（T=600, 组 A vs D）")
    rp.append("")
    rp.append("| pair | Δz | A | D(HLA+多频) |")
    rp.append("| --- | --- | --- | --- |")
    for zi, zj in PAIRS:
        da = pair_df[(pair_df["group"]=="A")&(pair_df["T_s"]==600)&(pair_df["z_i_m"]==zi)&(pair_df["z_j_m"]==zj)]["D_abs"]
        dd = pair_df[(pair_df["group"]=="D")&(pair_df["T_s"]==600)&(pair_df["z_i_m"]==zi)&(pair_df["z_j_m"]==zj)]["D_abs"]
        rp.append(f"| {zi}/{zj} | {zj-zi} | {fnum(da.iloc[0] if len(da) else np.nan)} | {fnum(dd.iloc[0] if len(dd) else np.nan)} |")
    rp.append("")
    rp.append(f"- 深度分辨下限：**{notes['depth_resolution_bound_m']} m**（10 m 对不过则记 >10 m；20 m 过则记 ~20 m）")
    rp.append("")
    rp.append("## 终判（Zhu 路线此后永久关闭）")
    rp.append("")
    rp.append(f"### `{decision}`")
    rp.append("")
    rp.append(why)
    rp.append("")
    rp.append(f"**下一步**：{nxt}")
    rp.append("")
    rp.append("## 停止")
    rp.append("")
    rp.append("- 不再 C1.x；不进 P5；本轮不进 Yang 2015")
    rp.append("- Zhu（QF+QM）**PERMANENTLY_CLOSED**")
    rp.append("")
    (OUT / "R3_C1_2_REPORT.md").write_text("\n".join(rp), encoding="utf-8")

    gs = [
        "# R3-C1.2 — GPT 同步", "",
        f"- **判定：{decision}**",
        f"- QF=CLOSED；QM 本轮后 Zhu 路线永久关闭",
        f"- {why}",
        f"- 下一步：{nxt}", "",
        "## T=600 stats", "",
        pd.DataFrame(stats).T.to_string(), "",
        f"D_200/210 (D组)={D['D_200_210']:.4f}  D_200/220={D['D_200_220']:.4f}  FWHM200={fwhm200:.1f}m", "",
        "停止：无 C1.x / Yang / P5。", "",
    ]
    (OUT / "R3_C1_2_GPT_SYNC.md").write_text("\n".join(gs), encoding="utf-8")

    # fig
    p = ['<svg xmlns="http://www.w3.org/2000/svg" width="720" height="360">',
         '<rect width="720" height="360" fill="#f7f4ef"/>',
         f'<text x="360" y="28" text-anchor="middle" font-family="sans-serif" font-size="14" font-weight="600">R3-C1.2 QM 组间 FWHM / D200-220</text>']
    for i, g in enumerate(["A", "B", "C", "D"]):
        x = 80 + i * 150
        h = min(stats[g]["mean_fwhm"], 120) / 120 * 180
        p.append(f'<rect x="{x}" y="{240-h:.1f}" width="40" height="{h:.1f}" fill="#b45309"/>')
        d2 = stats[g]["D_200_220"] * 400
        p.append(f'<rect x="{x+50}" y="{240-min(d2,180):.1f}" width="40" height="{min(d2,180):.1f}" fill="#0f766e"/>')
        p.append(f'<text x="{x+45}" y="260" text-anchor="middle" font-size="12" font-family="sans-serif">{g}</text>')
        p.append(f'<text x="{x+20}" y="280" text-anchor="middle" font-size="10" font-family="sans-serif">{stats[g]["mean_fwhm"]:.0f}m</text>')
        p.append(f'<text x="{x+70}" y="280" text-anchor="middle" font-size="10" font-family="sans-serif">{stats[g]["D_200_220"]:.3f}</text>')
    p.append('<text x="80" y="320" font-size="11" font-family="sans-serif">橙=mean FWHM m（≤120显示）  绿=D(200/220)×400</text></svg>')
    (FIG / "fig_r3c12_qm_groups.svg").write_text("\n".join(p), encoding="utf-8")

    print(f"DONE {time.time()-t0:.1f}s")
    print("DECISION", decision)
    print(why)


if __name__ == "__main__":
    main()
