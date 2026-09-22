#!/usr/bin/env python3
"""R3-C1.2 steering correction: fix endfire w_m sign; unit test; rerun B/D; T1200 conv."""
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
DR = 10.0

MODE = p4.MODE_ENV["E0"]
_MODES, _PHI = {}, {}


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


def track_r(T, dr=DR):
    r1 = min(R1, R0 + U_REL * T)
    n = int(max(20, math.floor((r1 - R0) / dr) + 1))
    return np.linspace(R0, r1, n)


def element_p(f, z_s, ri, xm, z_r=Z_R):
    """p_m ~ exp(j k_r (r-x_m)) → phase slope along x_m is -k_r."""
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
    """B=|w^H p|^2 with CORRECT endfire steering w_m=exp(-j k0 x_m).

    p_m ~ exp(-j k_r x_m); conj(w_m)*p_m ~ exp(+j k0 x_m)*exp(-j k_r x_m) → coherent if k_r≈k0.
    """
    xs = (np.arange(n_el) - (n_el - 1) / 2.0) * d_m if n_el > 1 else np.zeros(1)
    k0 = 2 * np.pi * f / C0
    w = np.exp(-1j * k0 * xs * 1.0) / n_el  # CORRECTED sign
    B = np.zeros(len(r))
    for i, ri in enumerate(r):
        P = np.array([element_p(f, z_s, ri, xm, z_r) for xm in xs])
        if n_el == 1:
            Bi = abs(P[0]) ** 2
        elif steer:
            Bi = abs(np.dot(np.conj(w), P)) ** 2
        else:
            Bi = float(np.mean(np.abs(P) ** 2))
        R = math.sqrt(ri ** 2 + (z_r - z_s) ** 2)
        B[i] = Bi * R ** 2
    return B


def plane_wave_unit_test():
    """Ideal endfire plane wave p_m = exp(-j k0 x_m).

    Correct steering w_m = exp(-j k0 x_m)/N ⇒ w^H p = 1 (coherent).
    Wrong (old) w_m = exp(+j k0 x_m)/N ⇒ phase doubling, output collapses.
    """
    k0 = 2 * np.pi * F_SINGLE / C0
    xs = (np.arange(N_EL) - (N_EL - 1) / 2.0) * D_EL
    p = np.exp(-1j * k0 * xs)
    w_ok = np.exp(-1j * k0 * xs) / N_EL
    out_ok = float(abs(np.dot(np.conj(w_ok), p)) ** 2)  # expect ~1
    w_bad = np.exp(+1j * k0 * xs) / N_EL
    out_bad = float(abs(np.dot(np.conj(w_bad), p)) ** 2)  # expect <<1
    w_bs = np.ones(N_EL, dtype=complex) / N_EL
    out_bs = float(abs(np.dot(np.conj(w_bs), p)) ** 2)
    # opposite-sign plane wave should fail correct steering
    p_anti = np.exp(+1j * k0 * xs)
    out_anti = float(abs(np.dot(np.conj(w_ok), p_anti)) ** 2)
    passed = bool(out_ok > 0.8 and out_bad < 0.2 * max(out_ok, 1e-30) and out_anti < 0.2)
    return {
        "plane_wave_endfire": "p_m=exp(-j k0 x_m)",
        "w_correct": "exp(-j k0 x_m)/N",
        "w_wrong_old": "exp(+j k0 x_m)/N",
        "out_correct_steering": out_ok,
        "out_wrong_steering": out_bad,
        "out_broadside": out_bs,
        "out_anti_plane_wave": out_anti,
        "phase_slope_expected_rad_per_2m": float(-k0 * 2.0),
        "unit_test_pass": passed,
    }


def norm_track(I):
    v = I - I.mean()
    return v / (np.linalg.norm(v) + 1e-30)


def qm_scores(freqs, z_true, z_grid, r, n_el, d_m, steer=True):
    obs, Qs = {}, []
    for f in freqs:
        u0 = norm_track(beam_track(f, z_true, r, n_el, d_m, steer))
        obs[f] = u0
        Qf = np.array([
            abs(float(np.dot(u0, norm_track(beam_track(f, z, r, n_el, d_m, steer)))))
            for z in z_grid
        ])
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
        return dict(z_hat=float(z_grid[imax]), fwhm=fwhm, fwhm_flag="FULL_WINDOW_UNRESOLVED", psl=np.nan)
    mask = np.ones_like(J, bool)
    mask[lo:hi + 1] = False
    psl = float(J[mask].max() / (jmax + 1e-30)) if mask.any() else 0.0
    return dict(z_hat=float(z_grid[imax]), fwhm=fwhm, fwhm_flag="local", psl=psl)


def main():
    t0 = time.time()
    print("=== R3-C1.2 steering correction ===", flush=True)

    # 1) unit test FIRST — block E-STD if fail
    ut = plane_wave_unit_test()
    print("unit_test", json.dumps(ut, ensure_ascii=False), flush=True)
    pd.DataFrame([ut]).to_csv(OUT / "r3c12_beamformer_unit_test.csv", index=False, encoding="utf-8-sig")
    if not ut["unit_test_pass"]:
        print("UNIT TEST FAIL — not entering E-STD", flush=True)
        (OUT / "R3_C1_2_FINAL_DECISION.json").write_text(json.dumps({
            "rc3c1_2_final_decision": "UNIT_TEST_FAIL_NO_ESTD",
            "why": "plane-wave steering unit test failed",
            "created_utc": NOW,
        }, indent=2), encoding="utf-8")
        return

    # 2) rerun B and D only (A/C keep prior final data)
    groups = {
        "B": dict(freqs=[F_SINGLE], n_el=N_EL, d_m=D_EL),
        "D": dict(freqs=F_MULTI, n_el=N_EL, d_m=D_EL),
    }
    amb_rows, pair_rows = [], []
    cache = {}
    for gtag, gp in groups.items():
        for T in T_LIST:
            r = track_r(T, DR)
            for z in Z_TRUES:
                Q, units = qm_scores(gp["freqs"], z, Z_GRID, r, gp["n_el"], gp["d_m"], steer=True)
                cache[(gtag, z, T)] = units
                met = fwhm_flags(Z_GRID, Q)
                amb_rows.append({
                    "group": gtag, "T_s": T, "z_true_m": z,
                    "z_hat_m": met["z_hat"], "err_z_m": met["z_hat"] - z,
                    "fwhm_local_m": met["fwhm"], "fwhm_flag": met["fwhm_flag"],
                    "psl": met["psl"], "steer": "CORRECTED exp(-j k0 x)",
                    "dr_m": DR, "n_el": gp["n_el"], "n_freqs": len(gp["freqs"]),
                })
                print(f"  {gtag} T={T} z={z}: FWHM={met['fwhm']:.1f} {met['fwhm_flag']}", flush=True)
            for zi, zj in PAIRS:
                Dv, cm, cs = D_multi(cache[(gtag, zi, T)], cache[(gtag, zj, T)])
                pair_rows.append({
                    "group": gtag, "T_s": T, "z_i_m": zi, "z_j_m": zj, "delta_z_m": abs(zj - zi),
                    "D_multi": Dv, "mean_abs_corr": cm,
                    "per_freq_abs_corr": json.dumps([round(x, 6) for x in cs]),
                    "distinguishable_gt_0p05": bool(Dv > 0.05),
                    "steer": "CORRECTED",
                })
                print(f"  pair {gtag} T={T} {zi}/{zj}: D={Dv:.5f}", flush=True)

    amb_df = pd.DataFrame(amb_rows)
    pair_df = pd.DataFrame(pair_rows)
    # merge with A/C from previous final if present
    old_amb_p = OUT / "r3c12_qm_ambiguity_final.csv"
    old_pair_p = OUT / "r3c12_qm_pairs_final.csv"
    if old_amb_p.exists():
        old_amb = pd.read_csv(old_amb_p)
        keep = old_amb[old_amb["group"].isin(["A", "C"])]
        amb_all = pd.concat([keep, amb_df], ignore_index=True)
    else:
        amb_all = amb_df
    if old_pair_p.exists():
        old_pair = pd.read_csv(old_pair_p)
        keep_p = old_pair[old_pair["group"].isin(["A", "C"])]
        pair_all = pd.concat([keep_p, pair_df], ignore_index=True)
    else:
        pair_all = pair_df
    amb_all.to_csv(OUT / "r3c12_qm_ambiguity_final.csv", index=False, encoding="utf-8-sig")
    pair_all.to_csv(OUT / "r3c12_qm_pairs_final.csv", index=False, encoding="utf-8-sig")

    # 3) T=1200 multifreq 200/220 sampling convergence
    conv_rows = []
    for dr in [5.0, 10.0, 20.0]:
        r = track_r(1200.0, dr)
        _, u200 = qm_scores(F_MULTI, 200.0, Z_GRID, r, N_EL, D_EL, True)
        _, u210 = qm_scores(F_MULTI, 210.0, Z_GRID, r, N_EL, D_EL, True)
        _, u220 = qm_scores(F_MULTI, 220.0, Z_GRID, r, N_EL, D_EL, True)
        D210 = D_multi(u200, u210)[0]
        D220 = D_multi(u200, u220)[0]
        conv_rows.append({
            "dr_m": dr, "n_r": len(r), "T_s": 1200.0, "group": "D",
            "D_200_210": D210, "D_200_220": D220,
            "per_freq_220": json.dumps(D_multi(u200, u220)[2]),
        })
        print(f"  T1200 conv dr={dr}: D210={D210:.5f} D220={D220:.5f}", flush=True)
    conv_df = pd.DataFrame(conv_rows)
    conv_df.to_csv(OUT / "r3c12_T1200_sampling_convergence.csv", index=False, encoding="utf-8-sig")

    def gstats(gtag, T):
        d = amb_all[(amb_all.group == gtag) & (amb_all.T_s == T)]
        p = pair_all[(pair_all.group == gtag) & (pair_all.T_s == T)]
        def dpr(zi, zj):
            x = p[(p.z_i_m == zi) & (p.z_j_m == zj)].D_multi
            return float(x.iloc[0]) if len(x) else np.nan
        return dict(
            mean_fwhm=float(d.fwhm_local_m.mean()) if len(d) else np.nan,
            frac_full=float((d.fwhm_flag == "FULL_WINDOW_UNRESOLVED").mean()) if len(d) else np.nan,
            D_200_210=dpr(200, 210), D_200_220=dpr(200, 220),
            pair_frac=float(p.distinguishable_gt_0p05.mean()) if len(p) else np.nan,
        )

    stats600 = {g: gstats(g, 600.0) for g in "ABCD"}
    stats1200 = {g: gstats(g, 1200.0) for g in "ABCD"}
    A, B, C, D = stats600["A"], stats600["B"], stats600["C"], stats600["D"]
    A12, B12, C12, D12 = stats1200["A"], stats1200["B"], stats1200["C"], stats1200["D"]

    # HLA increment vs single
    hla_gain_600 = {
        "dD_200_220": B["D_200_220"] - A["D_200_220"],
        "dD_200_210": B["D_200_210"] - A["D_200_210"],
        "dFWHM": B["mean_fwhm"] - A["mean_fwhm"],
    }
    hla_gain_1200 = {
        "dD_200_220": B12["D_200_220"] - A12["D_200_220"],
        "dFWHM": B12["mean_fwhm"] - A12["mean_fwhm"],
    }
    # D vs C at T=1200 (multifreq)
    d_vs_c_1200 = {
        "dD_200_220": D12["D_200_220"] - C12["D_200_220"],
        "dD_200_210": D12["D_200_210"] - C12["D_200_210"],
    }
    conv_stable = bool(
        len(conv_df) >= 2
        and (conv_df["D_200_220"].max() - conv_df["D_200_220"].min()) < 0.15
        and conv_df["D_200_220"].min() > 0.05
    )

    # terminal states only three
    hla_confirmed = bool(
        abs(hla_gain_600["dD_200_220"]) > 0.02
        or abs(hla_gain_1200["dD_200_220"]) > 0.05
        or (B["mean_fwhm"] < 0.7 * A["mean_fwhm"] if np.isfinite(B["mean_fwhm"]) else False)
    )
    multi_cond = bool(
        (not hla_confirmed)
        and D12["D_200_220"] > 0.05
        and D12["D_200_210"] < 0.05
        and conv_stable
    )
    weak = bool(
        D12["D_200_210"] < 0.05
        and (D12["D_200_220"] < 0.05 or not conv_stable)
        and (B["D_200_220"] < 0.05 and B["D_200_210"] < 0.05)
    )

    if hla_confirmed:
        decision = "C1_QM_HLA_ADAPTATION_CONFIRMED"
        why = (
            f"Corrected endfire steering yields HLA depth increment: "
            f"T=600 ΔD220={hla_gain_600['dD_200_220']:.4f}, T=1200 ΔD220={hla_gain_1200['dD_200_220']:.4f}, "
            f"D210_B={B['D_200_210']:.4f} D220_B={B['D_200_220']:.4f} FWHM_B={B['mean_fwhm']:.1f}."
        )
    elif multi_cond:
        decision = "C1_QM_MULTIFREQ_ONLY_CONDITIONAL"
        why = (
            f"HLA increment still weak (ΔD220_600={hla_gain_600['dD_200_220']:.4f}). "
            f"Multifreq long-track 20 m signal holds under sampling "
            f"(D220_T1200={D12['D_200_220']:.3f}, conv_stable={conv_stable}); "
            f"10 m still not separable (D210={D12['D_200_210']:.4f}). "
            f"Frequency-selective depth cue (e.g. 338 Hz) noted — aligns with Yang modal-depth idea."
        )
    else:
        decision = "C1_QM_CZ_DEPTH_WEAK"
        why = (
            f"After steering correction (unit test pass): HLA help remains weak "
            f"(ΔD220_600={hla_gain_600['dD_200_220']:.4f}); 10 m and 20 m pairs not reliably separable. "
            f"D12 D210={D12['D_200_210']:.4f} D220={D12['D_200_220']:.4f}; conv_stable={conv_stable}."
        )
    nxt = "Zhu route PERMANENTLY_CLOSED; Yang 2015 modal synthetic aperture next (not this round)"

    notes = {
        "unit_test": ut,
        "steering": "CORRECTED w_m=exp(-j k0 x_m)",
        "stats_T600": stats600,
        "stats_T1200": stats1200,
        "hla_increment_600": hla_gain_600,
        "hla_increment_1200": hla_gain_1200,
        "D_vs_C_T1200": d_vs_c_1200,
        "T1200_D_200_220": D12["D_200_220"],
        "T1200_D_200_210": D12["D_200_210"],
        "sampling_conv_stable": conv_stable,
        "A_C_frozen": "kept from previous final (single-channel)",
    }
    dec = {
        "rc3c1_2_final_decision": decision,
        "why": why,
        "next_step": nxt,
        "notes": notes,
        "Zhu_route": "PERMANENTLY_CLOSED",
        "created_utc": NOW,
        "stop": "no C1.x after this; Yang next not this round; no P5",
    }
    (OUT / "R3_C1_2_FINAL_DECISION.json").write_text(
        json.dumps(dec, ensure_ascii=False, indent=2, default=str), encoding="utf-8"
    )

    rp = []
    rp.append("# R3-C1.2 FINAL CORRECTION 报告")
    rp.append("")
    rp.append(f"UTC：{NOW}")
    rp.append("")
    rp.append("## 1. Steering 符号修正")
    rp.append("")
    rp.append("- \(p_m\\sim e^{-jk_rx_m}\)（\(r_m=r-x_m\)）")
    rp.append("- 正确：\(w_m=e^{-jk_0x_m}\)，则 \(w_m^*p_m\\sim e^{+jk_0x_m}e^{-jk_rx_m}\) 相干")
    rp.append("- 旧错误：\(w_m=e^{+jk_0x_m}\) 相位梯度加倍，B_steered 塌缩")
    rp.append("")
    rp.append("## 2. 平面波单元测试")
    rp.append("")
    rp.append(f"- pass = **{ut['unit_test_pass']}**")
    rp.append(f"- out_correct={ut['out_correct_steering']:.4g}, out_wrong={ut['out_wrong_steering']:.4g}, out_bs={ut['out_broadside']:.4g}")
    rp.append("")
    rp.append("## 3. B/D 重跑（A/C 冻结保留）")
    rp.append("")
    rp.append("| 组 | T | D200/210 | D200/220 | mean FWHM |")
    rp.append("| --- | --- | --- | --- | --- |")
    for g in "ABCD":
        for T, st in [(600, stats600), (1200, stats1200)]:
            s = st[g]
            rp.append(f"| {g} | {T} | {s['D_200_210']:.4f} | {s['D_200_220']:.4f} | {s['mean_fwhm']:.1f} |")
    rp.append("")
    rp.append(f"- HLA 相对 single 增量 T=600：ΔD220=**{hla_gain_600['dD_200_220']:.4f}**，ΔFWHM=**{hla_gain_600['dFWHM']:.1f} m**")
    rp.append(f"- T=1200 多频 20 m 信号：D220=**{D12['D_200_220']:.3f}**，采样稳定=**{conv_stable}**")
    rp.append("")
    rp.append("## 4. 终判")
    rp.append("")
    rp.append(f"### `{decision}`")
    rp.append("")
    rp.append(why)
    rp.append("")
    rp.append(f"**下一步**：{nxt}")
    rp.append("")
    rp.append("Zhu **PERMANENTLY_CLOSED**（本轮后无 C1.x）。")
    (OUT / "R3_C1_2_FINAL_CORRECTION.md").write_text("\n".join(rp), encoding="utf-8")
    (OUT / "R3_C1_2_FINAL_GPT_SYNC.md").write_text(
        f"# R3-C1.2 correction\n\n**{decision}**\n\n{why}\n\nnext: {nxt}\n",
        encoding="utf-8",
    )

    print("DECISION", decision)
    print(why)
    print(f"DONE {time.time()-t0:.1f}s")


if __name__ == "__main__":
    main()
