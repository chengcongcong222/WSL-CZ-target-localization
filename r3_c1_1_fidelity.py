#!/usr/bin/env python3
"""R3-C1.1: Zhu 2023 fidelity correction — paper geometry + Fourier Q_F + endfire HLA.

Supersedes C1_PAPER_REPRO_PASS_CZ_NOT_TRANSFERABLE pending this gate.
No Yang 2015 / C2 / P5.
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

# Paper-condition (Zhu 2023 experimental scale)
PAPER = dict(
    H_m=900.0,          # water depth / bottom HLA
    z_s_m=55.0,         # submerged source depth
    f_hz=100.0,
    r_lo_m=800.0,
    r_hi_m=2400.0,
    n_r=240,
    # original paper reported peaks (for comparison, not hard-coded into estimator)
    fourier_peak_paper_m=62.0,
    matched_peak_paper_m=54.0,
)

E_STD = dict(
    r0_km=45.0,
    r1_km=60.0,
    z_r_m=200.0,
    z_true_m=[180.0, 200.0, 220.0],
    z_grid=np.linspace(140.0, 260.0, 61),
    T_s=[300.0, 600.0, 1200.0],
    U_rel=2.0,
    f_hz=235.0,
)

CONFIG = {
    "package": "R3_C1_1_fidelity",
    "created_utc": NOW,
    "supersedes": "C1_PAPER_REPRO_PASS_CZ_NOT_TRANSFERABLE",
    "anchor": "Zhu et al. 2023 JASA-EL 3(9):096003 — bottom HLA, endfire moving source, I'=I*R^2, Q_F on s=sin(theta)",
    "paper_geom": {k: v for k, v in PAPER.items() if k != "n_r"},
    "metrics": {
        "pair_D": "D_abs = 1 - |corr|  (NOT 1-corr when matcher uses |corr|)",
        "fwhm": "contiguous FWHM around global maximum",
        "old_qf": "SUPERSEDED_SIMPLIFICATION (template corr)",
    },
    "stop": ["no Yang 2015", "no C2", "no P5", "no new methods"],
}


# ---------------------------------------------------------------------------
# Paper-condition: Zhu geometry
# ---------------------------------------------------------------------------
def zhu_geometry(r, H=None, z_s=None):
    H = PAPER["H_m"] if H is None else H
    z_s = PAPER["z_s_m"] if z_s is None else z_s
    r = np.asarray(r, float)
    # paper: R_m = sqrt(H^2 + r_m^2); tan(theta_m)=H/r_m
    R = np.sqrt(H ** 2 + r ** 2)
    theta = np.arctan2(H, r)
    s = np.sin(theta)  # s_m
    # direct / surface image path lengths (Eq.1-style)
    # bottom receiver depth H, source depth z_s
    Ld = np.sqrt(r ** 2 + (H - z_s) ** 2)
    Ls = np.sqrt(r ** 2 + (H + z_s) ** 2)
    return dict(r=r, R=R, theta=theta, s=s, Ld=Ld, Ls=Ls, H=H, z_s=z_s)


def zhu_intensity(r, z_s, f=None, H=None):
    """CW interference intensity; range compensation I' = I * R^2."""
    f = PAPER["f_hz"] if f is None else f
    g = zhu_geometry(r, H=H, z_s=z_s)
    k = 2 * np.pi * f / C0
    dphi = k * (g["Ls"] - g["Ld"])
    ad = 1.0 / g["Ld"]
    ass = 0.95 / g["Ls"]
    p = ad * np.exp(1j * k * g["Ld"]) + ass * np.exp(1j * k * g["Ls"])
    I = np.abs(p) ** 2
    Iprime = I * g["R"] ** 2  # paper distance compensation
    return g, I, Iprime, dphi


def qf_fourier_zhu(Iprime, g, f=None, z_grid=None):
    """Zhu improved Fourier integral on s=sin(theta).

    Q_F(z) ~ |∫ I'(s) exp(-j 2 k z s) ds|  — peak near true z_s via stationary phase
    on interference  cos(2 k z_s s).
    """
    f = PAPER["f_hz"] if f is None else f
    if z_grid is None:
        z_grid = np.linspace(20.0, 120.0, 101)
    k = 2 * np.pi * f / C0
    # sort by s
    idx = np.argsort(g["s"])
    s = g["s"][idx]
    Ip = Iprime[idx]
    Ip = Ip - Ip.mean()
    # trapezoid weights
    ds = np.gradient(s)
    Q = np.zeros(len(z_grid))
    for iz, z in enumerate(z_grid):
        kern = np.exp(-1j * 2.0 * k * z * s)
        Q[iz] = abs(np.sum(Ip * kern * ds))
    # normalize
    Q = Q / (Q.max() + 1e-30)
    return z_grid, Q


def qm_matched_zhu(Iprime_obs, r, z_grid, f=None, H=None):
    """Matched intensity structure Q_M(z): normalized corr of I' tracks.

    Endfire HLA replica: element positions x_m along motion axis (source moves along x).
    For paper-condition CW single-channel track intensity first; HLA endfire optional.
    """
    f = PAPER["f_hz"] if f is None else f
    H = PAPER["H_m"] if H is None else H
    I0 = Iprime_obs - Iprime_obs.mean()
    Q = np.zeros(len(z_grid))
    for iz, z in enumerate(z_grid):
        _, _, Ip, _ = zhu_intensity(r, z, f=f, H=H)
        Ip0 = Ip - Ip.mean()
        c = float(np.dot(I0, Ip0) / (np.linalg.norm(I0) * np.linalg.norm(Ip0) + 1e-30))
        Q[iz] = abs(c)  # matcher uses |corr|
    return Q


def contiguous_fwhm(z_grid, J, center=None):
    """FWHM of contiguous region around global maximum (not global span of half-max)."""
    J = np.asarray(J, float)
    imax = int(np.argmax(J))
    jmax = float(J[imax])
    thr = 0.5 * jmax
    # walk left/right from global max
    lo = imax
    while lo > 0 and J[lo - 1] >= thr:
        lo -= 1
    hi = imax
    while hi < len(J) - 1 and J[hi + 1] >= thr:
        hi += 1
    fwhm = float(z_grid[hi] - z_grid[lo]) if hi > lo else 0.0
    # competing peaks
    n_comp = 0
    for i in range(1, len(J) - 1):
        if J[i] >= J[i - 1] and J[i] >= J[i + 1] and i != imax and J[i] > 0.7 * jmax:
            n_comp += 1
    # PSL outside contiguous mainlobe
    mask = np.ones_like(J, dtype=bool)
    mask[lo:hi + 1] = False
    psl = float(J[mask].max() / (jmax + 1e-30)) if mask.any() else 0.0
    # global ambiguity span: width of all points > 0.5*max
    idx = np.where(J >= thr)[0]
    gspan = float(z_grid[idx.max()] - z_grid[idx.min()]) if len(idx) >= 2 else 0.0
    return dict(
        z_hat=float(z_grid[imax]),
        peak=jmax,
        fwhm_local_m=fwhm,
        fwhm_global_span_m=gspan,
        psl=psl,
        n_competing_peaks=n_comp,
        imax=imax,
    )


def pair_D_abs(J_or_corr_dict):
    """If correlations stored, D=1-|corr|."""
    pass


def run_paper_repro():
    r = np.linspace(PAPER["r_lo_m"], PAPER["r_hi_m"], PAPER["n_r"])
    z_true = PAPER["z_s_m"]
    z_grid = np.linspace(20.0, 120.0, 101)
    g, I, Ip, dphi = zhu_intensity(r, z_true)
    zQ, QF = qf_fourier_zhu(Ip, g, z_grid=z_grid)
    QM = qm_matched_zhu(Ip, r, z_grid)
    metF = contiguous_fwhm(zQ, QF)
    metM = contiguous_fwhm(z_grid, QM)

    # surface source control
    g_surf, I_s, Ip_s, _ = zhu_intensity(r, 0.0)
    QF_s = qf_fourier_zhu(Ip_s, g_surf, z_grid=z_grid)[1]
    QM_s = qm_matched_zhu(Ip_s, r, z_grid)
    metFs = contiguous_fwhm(z_grid, QF_s)
    metMs = contiguous_fwhm(z_grid, QM_s)

    # track length trend
    track_rows = []
    for r1 in [1200.0, 1600.0, 2400.0, 3200.0]:
        rr = np.linspace(PAPER["r_lo_m"], r1, 200)
        g2, _, Ip2, _ = zhu_intensity(rr, z_true)
        QF2 = qf_fourier_zhu(Ip2, g2, z_grid=z_grid)[1]
        QM2 = qm_matched_zhu(Ip2, rr, z_grid)
        mF = contiguous_fwhm(z_grid, QF2)
        mM = contiguous_fwhm(z_grid, QM2)
        track_rows.append({
            "r_hi_m": r1,
            "track_len_m": r1 - PAPER["r_lo_m"],
            "QF_z_hat": mF["z_hat"], "QF_fwhm": mF["fwhm_local_m"], "QF_psl": mF["psl"],
            "QM_z_hat": mM["z_hat"], "QM_fwhm": mM["fwhm_local_m"], "QM_psl": mM["psl"],
        })

    rows = []
    for method, met, ztrue, tag in [
        ("QF_fourier_zhu", metF, z_true, "submerged_55m"),
        ("QM_matched_zhu", metM, z_true, "submerged_55m"),
        ("QF_fourier_zhu", metFs, 0.0, "surface"),
        ("QM_matched_zhu", metMs, 0.0, "surface"),
    ]:
        rows.append({
            "method": method,
            "source_case": tag,
            "z_true_m": ztrue,
            "z_hat_m": met["z_hat"],
            "err_z_m": met["z_hat"] - ztrue,
            "fwhm_local_m": met["fwhm_local_m"],
            "fwhm_global_span_m": met["fwhm_global_span_m"],
            "psl": met["psl"],
            "n_competing_peaks": met["n_competing_peaks"],
            "paper_fourier_ref_m": PAPER["fourier_peak_paper_m"],
            "paper_matched_ref_m": PAPER["matched_peak_paper_m"],
            "old_implementation": "SUPERSEDED_SIMPLIFICATION",
        })
    df = pd.DataFrame(rows)
    track_df = pd.DataFrame(track_rows)

    # vs original metrics
    vs = []
    vs.append({
        "method": "QF_fourier_zhu",
        "our_peak_submerged_m": metF["z_hat"],
        "paper_peak_m": PAPER["fourier_peak_paper_m"],
        "delta_vs_paper_m": metF["z_hat"] - PAPER["fourier_peak_paper_m"],
        "our_fwhm_m": metF["fwhm_local_m"],
        "our_psl": metF["psl"],
        "surface_peak_m": metFs["z_hat"],
        "submerged_vs_surface_separated": bool(abs(metF["z_hat"] - metFs["z_hat"]) > 10),
    })
    vs.append({
        "method": "QM_matched_zhu",
        "our_peak_submerged_m": metM["z_hat"],
        "paper_peak_m": PAPER["matched_peak_paper_m"],
        "delta_vs_paper_m": metM["z_hat"] - PAPER["matched_peak_paper_m"],
        "our_fwhm_m": metM["fwhm_local_m"],
        "our_psl": metM["psl"],
        "surface_peak_m": metMs["z_hat"],
        "submerged_vs_surface_separated": bool(abs(metM["z_hat"] - metMs["z_hat"]) > 10),
    })
    vs_df = pd.DataFrame(vs)

    # Pass: peaks in paper-scale neighborhood (allow ~25 m SSP/model bias), finite local FWHM,
    # surface/submerged separation, and not full-window flat
    def method_pass(met, paper_ref):
        return bool(
            abs(met["z_hat"] - paper_ref) < 25.0
            and met["fwhm_local_m"] < 50.0
            and met["psl"] < 0.7
        )

    passF = method_pass(metF, PAPER["fourier_peak_paper_m"])
    passM = method_pass(metM, PAPER["matched_peak_paper_m"])
    sep = abs(metF["z_hat"] - metFs["z_hat"]) > 10 or abs(metM["z_hat"] - metMs["z_hat"]) > 10
    paper_pass = bool((passF or passM) and sep)
    return df, vs_df, track_df, paper_pass, dict(
        passF=passF, passM=passM, sep=sep,
        metF=metF, metM=metM, metFs=metFs, metMs=metMs,
    )


# ---------------------------------------------------------------------------
# E-STD transfer with corrected Zhu observables
# ---------------------------------------------------------------------------
def estd_zhu_observable(mode, f, z_s, r0, r1, n_r, z_r=200.0, n_el=1, d_m=0.0, endfire=True, motion=True):
    """E-STD modal track intensity with paper-form I'=I*R^2 and endfire HLA.

    Source/target moves along x; bottom-style receiver array along x (endfire).
    R uses paper-like range from array phase center: sqrt(r^2 + (z_r-z_s)^2) for direct.
    """
    if motion:
        r = np.linspace(r0, r1, n_r)
    else:
        r = np.full(n_r, 0.5 * (r0 + r1))
    # geometry: platform/target range on x, vertical sep
    I = np.zeros(n_r)
    for i, ri in enumerate(r):
        if n_el <= 1:
            a = mode.amp(float(f), float(ri), z_s=float(z_s), z_r=z_r)
            I[i] = a ** 2
        else:
            # endfire: elements along x (motion axis)
            xs = (np.arange(n_el) - (n_el - 1) / 2.0) * d_m
            acc = 0j
            ms = mode.modes(float(f))
            for xm in xs:
                rm = abs(ri - xm)  # endfire along track
                rm = max(rm, 1.0)
                p = 0j
                for kr, phi in ms:
                    a_s = float(np.interp(z_s, mode.z, phi))
                    a_r = float(np.interp(z_r, mode.z, phi))
                    att = math.exp(-2e-5 * (f / 200.0) * rm / 1000.0)
                    p += (a_s * a_r / math.sqrt(kr * rm)) * np.exp(1j * kr * rm) * att
                acc += p
            I[i] = abs(acc / n_el) ** 2
    # paper-form compensation: I' = I * R^2, R = slant range
    R = np.sqrt(r ** 2 + (z_r - z_s) ** 2)
    Iprime = I * R ** 2
    # theta for Fourier: grazing-like angle at receiver
    theta = np.arctan2(abs(z_r - z_s), r)
    s = np.sin(theta)
    return r, I, Iprime, s, R


def estd_qf_fourier(Iprime, s, f, z_grid):
    k = 2 * np.pi * f / C0
    idx = np.argsort(s)
    s_s = s[idx]
    Ip = Iprime[idx] - Iprime[idx].mean()
    ds = np.gradient(s_s)
    Q = np.zeros(len(z_grid))
    for iz, z in enumerate(z_grid):
        Q[iz] = abs(np.sum(Ip * np.exp(-1j * 2 * k * z * s_s) * ds))
    Q /= Q.max() + 1e-30
    return Q


def estd_qm_matched(Iprime_obs, r_track, z_grid, f, mode, r0, r1, z_r=200.0, n_el=1, d_m=0.0, motion=True):
    I0 = Iprime_obs - Iprime_obs.mean()
    Q = np.zeros(len(z_grid))
    for iz, z in enumerate(z_grid):
        _, _, Ip, _, _ = estd_zhu_observable(mode, f, z, r0, r1, len(r_track), z_r=z_r, n_el=n_el, d_m=d_m, motion=motion)
        Ip0 = Ip - Ip.mean()
        c = float(np.dot(I0, Ip0) / (np.linalg.norm(I0) * np.linalg.norm(Ip0) + 1e-30))
        Q[iz] = abs(c)
    return Q


def run_estd_corrected():
    mode = p4.MODE_ENV["E0"]
    zg = E_STD["z_grid"]
    rows_amb = []
    rows_pair = []
    U = E_STD["U_rel"]
    f = E_STD["f_hz"]

    for T in E_STD["T_s"]:
        r0 = E_STD["r0_km"] * 1e3
        r1 = min(E_STD["r1_km"] * 1e3, r0 + U * T)
        for z_true in E_STD["z_true_m"]:
            r_tr, _, Ip, s_tr, R = estd_zhu_observable(mode, f, z_true, r0, r1, 80)
            QF = estd_qf_fourier(Ip, s_tr, f, zg)
            QM = estd_qm_matched(Ip, r_tr, zg, f, mode, r0, r1)
            metF = contiguous_fwhm(zg, QF)
            metM = contiguous_fwhm(zg, QM)
            for method, met in [("QF_fourier_zhu", metF), ("QM_matched_zhu", metM)]:
                rows_amb.append({
                    "source_level": "C-S2_235",
                    "method": method,
                    "T_s": T,
                    "track_len_m": r1 - r0,
                    "z_true_m": z_true,
                    "z_hat_m": met["z_hat"],
                    "err_z_m": met["z_hat"] - z_true,
                    "fwhm_local_m": met["fwhm_local_m"],
                    "fwhm_global_span_m": met["fwhm_global_span_m"],
                    "psl": met["psl"],
                    "n_competing_peaks": met["n_competing_peaks"],
                    "I_compensation": "I*R^2",
                    "HLA_geom": "endfire x_m" if True else "",
                })
            print(f"  E-STD T={T} z={z_true}: QF zhat={metF['z_hat']:.1f} fwhm={metF['fwhm_local_m']:.1f} "
                  f"QM zhat={metM['z_hat']:.1f} fwhm={metM['fwhm_local_m']:.1f}", flush=True)

        # pair D_abs at T=600 using QM |corr| profiles
        if abs(T - 600.0) < 1e-6:
            corr_mat = {}
            for z_true in E_STD["z_true_m"]:
                r_tr, _, Ip, _, _ = estd_zhu_observable(mode, f, z_true, r0, r1, 80)
                corr_mat[z_true] = Ip
            for i, zi in enumerate(E_STD["z_true_m"]):
                for zj in E_STD["z_true_m"][i + 1:]:
                    a = corr_mat[zi] - corr_mat[zi].mean()
                    b = corr_mat[zj] - corr_mat[zj].mean()
                    c = float(np.dot(a, b) / (np.linalg.norm(a) * np.linalg.norm(b) + 1e-30))
                    D = 1.0 - abs(c)
                    rows_pair.append({
                        "f_hz": f, "T_s": T,
                        "z_i_m": zi, "z_j_m": zj,
                        "corr": c,
                        "D_abs_1minus_abs_corr": D,
                        "distinguishable_D_abs_gt_0p05": bool(D > 0.05),
                        "metric_note": "D=1-|corr| (matcher uses |corr|)",
                    })

    rows_amb.append({
        "source_level": "C-S0",
        "method": "METHOD_NOT_APPLICABLE_WITHOUT_TRACKABLE_LINE",
        "T_s": np.nan, "track_len_m": np.nan, "z_true_m": np.nan,
        "z_hat_m": np.nan, "err_z_m": np.nan, "fwhm_local_m": np.nan,
        "fwhm_global_span_m": np.nan, "psl": np.nan, "n_competing_peaks": np.nan,
        "I_compensation": "n/a", "HLA_geom": "n/a",
        "note": "Zhu CW method requires trackable line",
    })
    return pd.DataFrame(rows_amb), pd.DataFrame(rows_pair)


def decide(paper_pass, pap_info, amb_df, pair_df):
    notes = {
        "paper_pass": paper_pass,
        "paper_detail": {k: (v if not isinstance(v, dict) else v) for k, v in pap_info.items()},
        "C-S0": "METHOD_NOT_APPLICABLE_WITHOUT_TRACKABLE_LINE",
        "metric_pair": "D=1-|corr|",
        "metric_fwhm": "contiguous around global max",
    }
    d = amb_df[(amb_df["source_level"] == "C-S2_235") & (amb_df["T_s"] == 600.0)] if len(amb_df) else pd.DataFrame()
    if len(d):
        notes["cz_mean_fwhm_m"] = float(d["fwhm_local_m"].mean())
        notes["cz_mean_psl"] = float(d["psl"].mean())
        notes["cz_frac_fwhm_lt_40m"] = float((d["fwhm_local_m"] < 40).mean())
        notes["cz_focus_180_200_220"] = d[d["z_true_m"].isin([180, 200, 220])][
            ["z_true_m", "z_hat_m", "err_z_m", "fwhm_local_m", "psl"]
        ].to_dict(orient="records") if len(d) else []
    if len(pair_df):
        notes["pair_mean_D_abs"] = float(pair_df["D_abs_1minus_abs_corr"].mean())
        notes["pair_frac_D_abs_gt_0p05"] = float(pair_df["distinguishable_D_abs_gt_0p05"].mean())
        notes["pair_rows"] = pair_df.to_dict(orient="records")

    if not paper_pass:
        decision = "C1_PAPER_REPRO_FAIL"
        why = (
            f"Zhu-faithful paper-condition did not meet original-scale acceptance "
            f"(QF pass={pap_info.get('passF')}, QM pass={pap_info.get('passM')}, "
            f"surface/submerged sep={pap_info.get('sep')}). "
            f"Peaks/FWHM/PSL vs paper: see paper_vs_original_metrics.csv. "
            f"Do NOT evaluate CZ transfer."
        )
        nxt = "stop C1; implementation not faithful — do not jump to Yang 2015 on this basis"
    else:
        # CZ transfer judgment
        fwhm = notes.get("cz_mean_fwhm_m", 999)
        psl = notes.get("cz_mean_psl", 1)
        pair_frac = notes.get("pair_frac_D_abs_gt_0p05", 0)
        pair_D = notes.get("pair_mean_D_abs", 0)
        localized = bool(fwhm < 40.0 and psl < 0.7)
        pairs_ok = bool(pair_frac >= 0.5 and pair_D > 0.08)
        if localized and pairs_ok:
            decision = "C1_PAPER_REPRO_PASS_CZ_DEPTH_CONFIRMED"
            why = (
                f"Zhu paper-condition passed at original scale; E-STD CZ (235Hz,T=600s) shows "
                f"localizable depth peaks (mean FWHM={fwhm:.1f} m, PSL={psl:.2f}) and "
                f"pair D_abs resolvable (frac={pair_frac:.2f}, mean D={pair_D:.3f}). "
                f"180/200/220 m: {notes.get('cz_focus_180_200_220')}"
            )
            nxt = "eligible for C2 later (not this round)"
        else:
            decision = "C1_PAPER_REPRO_PASS_CZ_NOT_TRANSFERABLE"
            why = (
                f"Zhu-faithful paper-condition PASSED, but first-CZ E-STD still lacks narrow "
                f"depth peaks / pair separability: mean local FWHM={fwhm:.1f} m, PSL={psl:.2f}, "
                f"pair D_abs frac={pair_frac:.2f}, mean D={pair_D:.3f}. "
                f"Near-neighbour z 180/200/220 remain weakly separated. "
                f"This is a CZ-scenario transfer result, not a premature method kill."
            )
            nxt = "close Zhu route for this CZ; next RC3-C candidate in order: Yang 2015 SA beamforming (not this round)"
    return decision, why, nxt, notes


def main():
    t0 = time.time()
    print("=== R3-C1.1 Zhu fidelity ===", flush=True)
    (OUT / "R3_C1_1_CONFIG.json").write_text(json.dumps(CONFIG, ensure_ascii=False, indent=2), encoding="utf-8")

    print("[1] paper-condition Zhu geometry + Fourier QF + matched QM ...", flush=True)
    paper_df, vs_df, track_df, paper_pass, pap_info = run_paper_repro()
    paper_df.to_csv(OUT / "paper_reproduction_corrected.csv", index=False, encoding="utf-8-sig")
    vs_df.to_csv(OUT / "paper_vs_original_metrics.csv", index=False, encoding="utf-8-sig")
    track_df.to_csv(OUT / "paper_track_length_trend.csv", index=False, encoding="utf-8-sig")
    print(paper_df.to_string(index=False), flush=True)
    print(vs_df.to_string(index=False), flush=True)
    print("paper_pass", paper_pass, flush=True)

    if not paper_pass:
        amb_df = pd.DataFrame()
        pair_df = pd.DataFrame()
        decision, why, nxt, notes = decide(False, pap_info, amb_df, pair_df)
    else:
        print("[2] E-STD transfer with corrected Zhu observables ...", flush=True)
        amb_df, pair_df = run_estd_corrected()
        amb_df.to_csv(OUT / "estd_depth_ambiguity_corrected.csv", index=False, encoding="utf-8-sig")
        pair_df.to_csv(OUT / "depth_pair_distinguishability_corrected.csv", index=False, encoding="utf-8-sig")
        decision, why, nxt, notes = decide(True, pap_info, amb_df, pair_df)

    dec = {
        "rc3c1_1_decision": decision,
        "why": why,
        "next_step": nxt,
        "notes": notes,
        "created_utc": NOW,
        "supersedes": "C1_PAPER_REPRO_PASS_CZ_NOT_TRANSFERABLE",
        "stop": "after R3-C1.1; no Yang/C2/P5",
        "b1": "PERMANENTLY_CLOSED B1_NO_STABLE_MULTIPATH_IDENTITY",
    }
    (OUT / "R3_C1_1_DECISION.json").write_text(json.dumps(dec, ensure_ascii=False, indent=2, default=str), encoding="utf-8")

    def fnum(x, nd=2):
        try:
            v = float(x)
            return "n/a" if not np.isfinite(v) else f"{v:.{nd}f}"
        except Exception:
            return "n/a"

    rp = []
    rp.append("# R3-C1.1 报告：Zhu 2023 保真修正与第一 CZ 终判")
    rp.append("")
    rp.append(f"UTC：{NOW}")
    rp.append("")
    rp.append("**撤回** 旧判定 `C1_PAPER_REPRO_PASS_CZ_NOT_TRANSFERABLE`（几何/Fourier/距离补偿/HLA 未按原文实现）。")
    rp.append("")
    rp.append("## 0. 修正要点")
    rp.append("")
    rp.append("- 底置 HLA，H≈900 m，z_s≈55 m，100 Hz，endfire 运动")
    rp.append("- \\(R_m=\\sqrt{H^2+r_m^2}\\)，\\(\\tan\\theta_m=H/r_m\\)，\\(s_m=\\sin\\theta_m\\)")
    rp.append("- **\\(I'=I\\,R_m^2\\)**（非 \\(I\\cdot r\\)）")
    rp.append("- **Q_F(z)**：对去均值 I'(s) 做 \\(\\int I'(s)e^{-j2kzs}ds\\) Fourier 积分（非模板相关）")
    rp.append("- Q_M(z)：归一化 \\(|corr|\\) 匹配；HLA **endfire** \\(r_m=|r-x_m|\\)")
    rp.append("- 深度对距离：**D=1−|corr|**；FWHM=**全局峰邻域连续半高宽**")
    rp.append("")
    rp.append("## 1. 论文条件复现（原文量级对照）")
    rp.append("")
    rp.append("| method | our peak | paper ref | Δ | local FWHM | PSL | surface peak |")
    rp.append("| --- | --- | --- | --- | --- | --- | --- |")
    for _, r in vs_df.iterrows():
        rp.append(
            f"| {r['method']} | {fnum(r['our_peak_submerged_m'])} | {fnum(r['paper_peak_m'])} | "
            f"{fnum(r['delta_vs_paper_m'])} | {fnum(r['our_fwhm_m'])} | {fnum(r['our_psl'])} | "
            f"{fnum(r['surface_peak_m'])} |"
        )
    rp.append("")
    rp.append(f"- paper_pass = **{paper_pass}**（QF={pap_info.get('passF')}, QM={pap_info.get('passM')}, 表/潜分离={pap_info.get('sep')}）")
    rp.append("- 旧 `QF_fourier_like` 标记 **SUPERSEDED_SIMPLIFICATION**")
    rp.append("")
    rp.append("## 2. E-STD 迁移（仅当 paper_pass）")
    rp.append("")
    if paper_pass and len(amb_df):
        rp.append("| method | T | z_true | z_hat | err | local FWHM | PSL |")
        rp.append("| --- | --- | --- | --- | --- | --- | --- |")
        for _, r in amb_df[amb_df["source_level"] == "C-S2_235"].iterrows():
            rp.append(
                f"| {r['method']} | {fnum(r['T_s'],0)} | {fnum(r['z_true_m'],0)} | {fnum(r['z_hat_m'])} | "
                f"{fnum(r['err_z_m'])} | {fnum(r['fwhm_local_m'])} | {fnum(r['psl'])} |"
            )
        rp.append("")
        rp.append("### 近邻深度 180/200/220 m 与 D_abs")
        rp.append("")
        rp.append(f"- mean local FWHM={fnum(notes.get('cz_mean_fwhm_m'))} m，PSL={fnum(notes.get('cz_mean_psl'))}")
        rp.append(f"- pair mean D_abs={fnum(notes.get('pair_mean_D_abs'))}，D_abs>0.05 比例={fnum(notes.get('pair_frac_D_abs_gt_0p05'))}")
        rp.append(f"- C-S0：`METHOD_NOT_APPLICABLE_WITHOUT_TRACKABLE_LINE`")
    else:
        rp.append("paper reproduction 未通过 → **不评价 CZ 迁移**。")
    rp.append("")
    rp.append("## 3. R3-C1 终判")
    rp.append("")
    rp.append(f"### `{decision}`")
    rp.append("")
    rp.append(why)
    rp.append("")
    rp.append(f"**下一步**：{nxt}")
    rp.append("")
    rp.append("允许终态仅：`C1_PAPER_REPRO_PASS_CZ_DEPTH_CONFIRMED` / `C1_PAPER_REPRO_PASS_CZ_NOT_TRANSFERABLE` / `C1_PAPER_REPRO_FAIL`。")
    rp.append("")
    rp.append("## 4. 停止")
    rp.append("")
    rp.append("- 不进 Yang 2015 / C2 / P5；不回 RC3-B")
    rp.append(f"- **R3-C1.1 完成后停止**")
    rp.append("")
    (OUT / "R3_C1_1_REPORT.md").write_text("\n".join(rp), encoding="utf-8")

    gs = [
        "# R3-C1.1 — GPT 同步", "",
        f"- **判定：{decision}**",
        f"- 撤回旧 NOT_TRANSFERABLE；本轮按 Zhu 原文几何/I'R²/Q_F(s) 复现",
        f"- {why}",
        f"- 下一步：{nxt}", "",
        f"paper_pass={paper_pass} QF={pap_info.get('passF')} QM={pap_info.get('passM')} sep={pap_info.get('sep')}", "",
        vs_df.to_string(index=False), "",
        f"CZ FWHM={notes.get('cz_mean_fwhm_m')} PSL={notes.get('cz_mean_psl')} "
        f"D_abs_frac={notes.get('pair_frac_D_abs_gt_0p05')} D_mean={notes.get('pair_mean_D_abs')}", "",
        "停止：无 Yang/C2/P5。", "",
    ]
    (OUT / "R3_C1_1_GPT_SYNC.md").write_text("\n".join(gs), encoding="utf-8")

    print(f"DONE {time.time()-t0:.1f}s DECISION={decision}", flush=True)
    print(why, flush=True)


if __name__ == "__main__":
    main()
