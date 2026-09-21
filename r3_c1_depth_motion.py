#!/usr/bin/env python3
"""R3-C1: Zhu 2023 moving-source HLA depth estimation — paper repro + E-STD transfer.

No C2 / P5 / RC3-B return. Matched env, noiseless mechanism first.
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
CONFIG = {
    "package": "R3_C1_depth_motion",
    "created_utc": NOW,
    "anchor": {
        "primary": "Zhu F., Li F., Zhang Y., et al. Moving source depth estimation in deep ocean direct arrival zone with a horizontal line array. JASA Express Letters, 2023, 3(9):096003.",
        "backup": "T.C. Yang 2015 synthetic aperture beamforming for moving source depth",
        "methods": ["improved Fourier integral Q_F(z)", "matched sound intensity structure (paper internal control)"],
    },
    "source_levels": {
        "C-S2": "stable lines 201/235 Hz + multi-line joint",
        "C-S1": "in-band stable machinery lines only",
        "C-S0": "no trackable line → METHOD_NOT_APPLICABLE_WITHOUT_TRACKABLE_LINE",
    },
    "estd": {
        "r_km": [45, 60],
        "z_true_m": [150, 180, 200, 220, 250],
        "z_r_m": 200,
        "T_s": [300, 600, 1200],
        "U_rel_mps": 2.0,
        "propagation": "E-STD modal theory (frozen); algorithm stays Zhu track-intensity form",
    },
    "comparisons": ["single-channel+motion", "14m HLA+motion", "HLA static snapshot"],
    "stop": ["no C2", "no P5", "no RC3-B", "no env mismatch / MC this round"],
}


# ---------------------------------------------------------------------------
# Paper-condition: direct + surface image, endfire moving source
# ---------------------------------------------------------------------------
def paper_geometry(z_s, z_r=Z_R, r0=1000.0, r1=4000.0, n_r=200):
    """Source moves along x (endfire), array origin, constant depth z_s."""
    r = np.linspace(r0, r1, n_r)
    # direct
    Ld = r.copy()
    Ls = np.sqrt(r ** 2 + (z_s + z_r) ** 2)
    # arrival angle from horizontal (receiver→source sense)
    th_d = np.arctan2(z_s - z_r, r)
    th_s = np.arctan2(-(z_s + z_r), r)  # image below
    sin_d = np.sin(th_d)
    sin_s = np.sin(th_s)
    return dict(r=r, Ld=Ld, Ls=Ls, th_d=th_d, th_s=th_s, sin_d=sin_d, sin_s=sin_s)


def paper_intensity(z_s, f, r0=1000.0, r1=4000.0, n_r=200, A_d=1.0, A_s=0.9):
    g = paper_geometry(z_s, r0=r0, r1=r1, n_r=n_r)
    k = 2 * np.pi * f / C0
    dphi = k * (g["Ls"] - g["Ld"])
    # spherical spreading amplitudes
    ad = A_d / g["Ld"]
    ass = A_s / g["Ls"]
    p = ad * np.exp(1j * k * g["Ld"]) + ass * np.exp(1j * k * g["Ls"])
    I = np.abs(p) ** 2
    # range-compensated intensity I' as in paper-like form: I * r
    Iprime = I * g["r"]
    return g, I, Iprime, dphi


def qf_fourier(Iprime, g, f, z_grid, z_r=Z_R):
    """Improved-Fourier-like depth weight: match track structure in sin(theta) domain.

    Q_F(z) ~ |∫ I'(u) exp(-j ψ(z,u,f)) du| with ψ from two-path phase at candidate z.
    """
    u = np.sin(g["th_d"])  # use direct-arrival sin angle track
    # sort by u
    idx = np.argsort(u)
    u_s = u[idx]
    Ip = Iprime[idx]
    # remove mean
    Ip0 = Ip - Ip.mean()
    k = 2 * np.pi * f / C0
    Q = np.zeros(len(z_grid))
    for iz, z in enumerate(z_grid):
        g_z = paper_geometry(z, z_r=z_r, r0=float(g["r"][0]), r1=float(g["r"][-1]), n_r=len(g["r"]))
        u_z = np.sin(g_z["th_d"])[idx]
        dphi_z = k * (g_z["Ls"] - g_z["Ld"])[idx]
        # phase kernel: predicted interference modulation
        kern = np.cos(dphi_z)
        kern = kern - kern.mean()
        # interpolate Ip onto candidate u grid via original r order — use same r samples
        # match: correlation of I'(r) with cos(Δφ(r; z))
        num = float(np.dot(Ip0, kern))
        den = float(np.linalg.norm(Ip0) * np.linalg.norm(kern) + 1e-30)
        Q[iz] = abs(num / den)
    return Q


def matched_intensity_structure(Iprime, g_true_r, z_grid, f, z_r=Z_R, r0=1000.0, r1=4000.0):
    """Matched sound intensity structure: J(z)=corr(I'_obs(r), I'_pred(r|z))."""
    Iobs = Iprime
    Iobs0 = Iobs - Iobs.mean()
    J = np.zeros(len(z_grid))
    for iz, z in enumerate(z_grid):
        _, _, Ip, _ = paper_intensity(z, f, r0=r0, r1=r1, n_r=len(g_true_r))
        Ip0 = Ip - Ip.mean()
        den = np.linalg.norm(Iobs0) * np.linalg.norm(Ip0) + 1e-30
        J[iz] = float(np.dot(Iobs0, Ip0) / den)
    return np.abs(J)


def peak_metrics(z_grid, J, z_true):
    jmax = float(np.max(J))
    iz = int(np.argmax(J))
    z_hat = float(z_grid[iz])
    # FWHM
    thr = 0.5 * jmax
    idx = np.where(J >= thr)[0]
    fwhm = float(z_grid[idx.max()] - z_grid[idx.min()]) if len(idx) >= 2 else 0.0
    # 3dB-like
    thr3 = (10 ** (-3 / 20)) * jmax
    idx3 = np.where(J >= thr3)[0]
    w3 = float(z_grid[idx3.max()] - z_grid[idx3.min()]) if len(idx3) >= 2 else 0.0
    # sidelobe
    mask = np.ones_like(J, dtype=bool)
    mask[max(iz - 2, 0): iz + 3] = False
    side = float(J[mask].max()) if mask.any() else 0.0
    psl = side / (jmax + 1e-30)
    # multi-peak
    n_loc = 0
    for i in range(1, len(J) - 1):
        if J[i] >= J[i - 1] and J[i] >= J[i + 1] and J[i] > 0.7 * jmax:
            n_loc += 1
    return dict(
        z_hat_m=z_hat,
        err_z_m=z_hat - z_true,
        peak= jmax,
        fwhm_m=fwhm,
        width_3db_m=w3,
        psl=psl,
        n_local_peaks=n_loc,
        multimodal=bool(n_loc > 1),
    )


def run_paper_repro():
    z_grid = np.linspace(20.0, 200.0, 91)
    rows = []
    for z_true in [50.0, 55.0, 62.0, 100.0]:
        for f in [235.0]:
            for r0, r1, tag in [(1000, 4000, "track_3km"), (1000, 2000, "track_1km"), (1000, 8000, "track_7km")]:
                g, I, Ip, _ = paper_intensity(z_true, f, r0=r0, r1=r1, n_r=400)
                Q = qf_fourier(Ip, g, f, z_grid)
                metQ = peak_metrics(z_grid, Q, z_true)
                J = matched_intensity_structure(Ip, g["r"], z_grid, f, r0=r0, r1=r1)
                metJ = peak_metrics(z_grid, J, z_true)
                rows.append({
                    "method": "QF_fourier_like",
                    "z_true_m": z_true,
                    "f_hz": f,
                    "track": tag,
                    "r0_m": r0,
                    "r1_m": r1,
                    **metQ,
                })
                rows.append({
                    "method": "matched_intensity_structure",
                    "z_true_m": z_true,
                    "f_hz": f,
                    "track": tag,
                    "r0_m": r0,
                    "r1_m": r1,
                    **metJ,
                })
    df = pd.DataFrame(rows)
    df.to_csv(OUT / "paper_reproduction.csv", index=False, encoding="utf-8-sig")
    # pass criteria: peak near true z for at least one method on main track
    main = df[(df["track"] == "track_3km") & (df["z_true_m"].isin([50.0, 55.0]))]
    err_ok = main.groupby("method")["err_z_m"].apply(lambda s: float(np.mean(np.abs(s))))
    # paper-like: deeper / surface separation trends
    print(err_ok, flush=True)
    pass_repro = bool(err_ok.min() < 15.0)  # allow some bias in simplified QF
    return df, pass_repro, err_ok.to_dict()


# ---------------------------------------------------------------------------
# E-STD: modal track intensity + Zhu-form observables
# ---------------------------------------------------------------------------
def estd_track_intensity(mode, f, z_s, r0_m, r1_m, n_r, z_r=Z_R, n_el=1, d_m=0.0, motion=True):
    """I(track) using modal amp; optional HLA spatial correlation vs single channel.

    For motion=True: r spans [r0,r1] along track.
    For motion=False: single snapshot at mid-range, repeated (static HLA).
    """
    if motion:
        r = np.linspace(r0_m, r1_m, n_r)
    else:
        r = np.full(n_r, 0.5 * (r0_m + r1_m))
    I = np.zeros(n_r)
    for i, ri in enumerate(r):
        if n_el <= 1:
            a = mode.amp(float(f), float(ri), z_s=float(z_s), z_r=z_r)
            I[i] = a ** 2
        else:
            # HLA: average element powers + phase-coherent sum magnitude normalized
            ys = (np.arange(n_el) - (n_el - 1) / 2.0) * d_m
            acc = 0j
            for y in ys:
                rm = math.hypot(ri, y)  # source on x-axis at range ri
                ms = mode.modes(float(f))
                p = 0j
                for kr, phi in ms:
                    a_s = float(np.interp(z_s, mode.z, phi))
                    a_r = float(np.interp(z_r, mode.z, phi))
                    att = math.exp(-2e-5 * (f / 200.0) * rm / 1000.0)
                    p += (a_s * a_r / math.sqrt(kr * max(rm, 1))) * np.exp(1j * kr * rm) * att
                acc += p
            I[i] = abs(acc / n_el) ** 2
    Iprime = I * r
    return r, I, Iprime


def estd_depth_match(Iprime_obs, r_track, z_grid, f, mode, r0, r1, z_r=Z_R, n_el=1, d_m=0.0, motion=True):
    """Matched intensity structure on E-STD: J(z)=corr(I'_obs, I'_pred(z))."""
    I0 = Iprime_obs - Iprime_obs.mean()
    J = np.zeros(len(z_grid))
    for iz, z in enumerate(z_grid):
        _, _, Ip = estd_track_intensity(mode, f, z, r0, r1, len(r_track), z_r=z_r, n_el=n_el, d_m=d_m, motion=motion)
        Ip0 = Ip - Ip.mean()
        den = np.linalg.norm(I0) * np.linalg.norm(Ip0) + 1e-30
        J[iz] = float(np.dot(I0, Ip0) / den)
    return np.abs(J)


def structure_distance(mode, f, z_list, r0, r1, n_r=80, z_r=Z_R):
    """Pairwise D(z_i,z_j) = 1 - corr of normalized I' tracks."""
    I = {}
    for z in z_list:
        r, _, Ip = estd_track_intensity(mode, f, z, r0, r1, n_r, z_r=z_r)
        In = Ip - Ip.mean()
        I[z] = In / (np.linalg.norm(In) + 1e-30)
    D = np.zeros((len(z_list), len(z_list)))
    for i, zi in enumerate(z_list):
        for j, zj in enumerate(z_list):
            D[i, j] = 1.0 - float(np.dot(I[zi], I[zj]))
    return D, z_list


def run_estd():
    mode = p4.MODE_ENV["E0"]
    z_grid = np.linspace(140.0, 260.0, 61)
    z_trues = [150.0, 180.0, 200.0, 220.0, 250.0]
    U = 2.0
    freq_sets = {
        "C-S2_201": [201.0],
        "C-S2_235": [235.0],
        "C-S2_joint": [201.0, 235.0],
        "C-S1": [204.0, 232.0],
    }
    rows_amb = []
    rows_pair = []
    rows_track = []
    rows_cmp = []

    for T in [300.0, 600.0, 1200.0]:
        r0 = 45e3
        r1 = r0 + U * T
        if r1 > 60e3:
            r1 = 60e3
        track_len = r1 - r0
        for src_name, freqs in freq_sets.items():
            for z_true in z_trues:
                # joint matched: sum correlations over freqs
                J_sum = np.zeros(len(z_grid))
                for f in freqs:
                    r_tr, _, Ip = estd_track_intensity(mode, f, z_true, r0, r1, 60, n_el=1)
                    J = estd_depth_match(Ip, r_tr, z_grid, f, mode, r0, r1, n_el=1)
                    J_sum += J
                J_sum /= len(freqs)
                met = peak_metrics(z_grid, J_sum, z_true)
                rows_amb.append({
                    "source_level": src_name,
                    "freqs_hz": json.dumps(freqs),
                    "T_s": T,
                    "track_len_m": track_len,
                    "z_true_m": z_true,
                    "n_el": 1,
                    "motion": True,
                    **met,
                    "method": "Zhu-like matched intensity (E-STD modal I')",
                })
                print(f"  amb {src_name} T={T} z={z_true}: zhat={met['z_hat_m']:.1f} err={met['err_z_m']:.1f} fwhm={met['fwhm_m']:.1f} psl={met['psl']:.2f}", flush=True)

        # pair distinguishability at T=600, C-S2_235
        if T == 600.0:
            D, zl = structure_distance(mode, 235.0, z_trues, r0, r1)
            for i, zi in enumerate(zl):
                for j, zj in enumerate(zl):
                    if i < j:
                        rows_pair.append({
                            "f_hz": 235.0,
                            "T_s": T,
                            "z_i_m": zi,
                            "z_j_m": zj,
                            "D_1minus_corr": D[i, j],
                            "distinguishable_D_gt_0p05": bool(D[i, j] > 0.05),
                        })
            # track length boundary
            for Tt in [300.0, 600.0, 1200.0]:
                rr0 = 45e3
                rr1 = min(60e3, rr0 + U * Tt)
                for z_true in [180.0, 200.0, 220.0]:
                    r_tr, _, Ip = estd_track_intensity(mode, 235.0, z_true, rr0, rr1, 60)
                    J = estd_depth_match(Ip, r_tr, z_grid, 235.0, mode, rr0, rr1)
                    met = peak_metrics(z_grid, J, z_true)
                    rows_track.append({
                        "f_hz": 235.0,
                        "T_s": Tt,
                        "track_len_m": rr1 - rr0,
                        "z_true_m": z_true,
                        **met,
                    })
            # A/B/C: single+motion, HLA+motion, HLA static
            for z_true in [200.0]:
                for label, n_el, d_m, motion in [
                    ("A_single_motion", 1, 0.0, True),
                    ("B_hla14m_motion", 8, 2.0, True),
                    ("C_hla14m_static", 8, 2.0, False),
                ]:
                    r_tr, _, Ip = estd_track_intensity(mode, 235.0, z_true, r0, r1, 60, n_el=n_el, d_m=d_m, motion=motion)
                    J = estd_depth_match(Ip, r_tr, z_grid, 235.0, mode, r0, r1, n_el=n_el, d_m=d_m, motion=motion)
                    met = peak_metrics(z_grid, J, z_true)
                    rows_cmp.append({
                        "config": label,
                        "f_hz": 235.0,
                        "T_s": 600.0,
                        "z_true_m": z_true,
                        "n_el": n_el,
                        "aperture_m": d_m * max(n_el - 1, 0),
                        "motion": motion,
                        **met,
                    })

    # S0 not applicable note row
    rows_amb.append({
        "source_level": "C-S0",
        "freqs_hz": "[]",
        "T_s": np.nan,
        "track_len_m": np.nan,
        "z_true_m": np.nan,
        "n_el": 1,
        "motion": True,
        "z_hat_m": np.nan,
        "err_z_m": np.nan,
        "peak": np.nan,
        "fwhm_m": np.nan,
        "width_3db_m": np.nan,
        "psl": np.nan,
        "n_local_peaks": np.nan,
        "multimodal": np.nan,
        "method": "METHOD_NOT_APPLICABLE_WITHOUT_TRACKABLE_LINE",
        "note": "Zhu CW track method requires trackable stable line; S0 not forced",
    })

    return (pd.DataFrame(rows_amb), pd.DataFrame(rows_pair),
            pd.DataFrame(rows_track), pd.DataFrame(rows_cmp))


def decide(paper_pass, err_map, amb_df, pair_df):
    notes = {}
    notes["paper_repro_pass"] = bool(paper_pass)
    notes["paper_err_by_method"] = err_map
    # C-S0
    notes["C-S0"] = "METHOD_NOT_APPLICABLE_WITHOUT_TRACKABLE_LINE"
    # E-STD: use C-S2_235, T=600
    d = amb_df[(amb_df.get("source_level") == "C-S2_235") & (amb_df["T_s"] == 600.0)] if len(amb_df) else pd.DataFrame()
    if len(d):
        notes["cz_mean_abs_err_m"] = float(d["err_z_m"].abs().mean())
        notes["cz_mean_fwhm_m"] = float(d["fwhm_m"].mean())
        notes["cz_true_in_peak"] = float((d["err_z_m"].abs() < 25.0).mean())
        notes["cz_unique_peak_frac"] = float((~d["multimodal"]).mean())
        notes["cz_mean_psl"] = float(d["psl"].mean())
    if len(pair_df):
        notes["pair_mean_D"] = float(pair_df["D_1minus_corr"].mean())
        notes["pair_resolvable_frac"] = float(pair_df["distinguishable_D_gt_0p05"].mean())
    # information presence
    if len(d):
        info_present = bool(d["err_z_m"].abs().median() < 30.0 or notes.get("cz_mean_fwhm_m", 999) < 80.0)
        distinct = bool(notes.get("pair_resolvable_frac", 0) > 0.3)
    else:
        info_present = False
        distinct = False
    notes["cz_depth_info_present"] = info_present
    notes["cz_pairs_distinguishable"] = distinct

    if not paper_pass:
        decision = "C1_PAPER_REPRO_FAIL"
        why = f"Paper-condition reproduction failed (err map={err_map}). Do not evaluate CZ method."
        nxt = "fix Zhu method implementation before any CZ claim"
    elif info_present and distinct:
        decision = "C1_PAPER_REPRO_PASS_CZ_DEPTH_CONFIRMED"
        why = (
            f"Paper repro passed; E-STD C-S2_235 T=600s: mean|err_z|={notes['cz_mean_abs_err_m']:.1f} m, "
            f"mean FWHM={notes['cz_mean_fwhm_m']:.1f} m, true-in-peak={notes['cz_true_in_peak']:.2f}, "
            f"pair resolvable frac={notes.get('pair_resolvable_frac')}. "
            f"S0=NOT_APPLICABLE (no trackable line)."
        )
        nxt = "eligible for C2 (SNR/mismatch/MC) later — not this round"
    elif info_present and not distinct:
        decision = "C1_PAPER_REPRO_PASS_CZ_CONDITIONAL"
        why = (
            f"Paper repro OK; CZ shows some true-depth response "
            f"(mean|err|={notes['cz_mean_abs_err_m']:.1f} m, FWHM={notes['cz_mean_fwhm_m']:.1f} m) "
            f"but pair structure weakly distinguishable (frac={notes.get('pair_resolvable_frac')})."
        )
        nxt = "C2 only if restricted to favorable z/f/track windows"
    else:
        decision = "C1_PAPER_REPRO_PASS_CZ_NOT_TRANSFERABLE"
        why = (
            f"Paper-condition method works, but first-CZ modal track intensity does not form "
            f"clear true-depth peaks under current E-STD (mean|err|={notes.get('cz_mean_abs_err_m')}, "
            f"FWHM={notes.get('cz_mean_fwhm_m')}, pairs={notes.get('pair_resolvable_frac')})."
        )
        nxt = "next RC3-C candidate in order: Yang 2015 SA beamforming → HLA modal/k-spectrum → Emmetière 2019"
    return decision, why, nxt, notes


def _esc(s):
    return str(s).replace("&", "&amp;").replace("<", "&lt;").replace(">", "&gt;")


def write_line_svg(path, title, series, xlab, ylab, note=""):
    w, h = 760, 360
    ml, mr, mt, mb = 70, 130, 48, 50
    pw, ph = w - ml - mr, h - mt - mb
    xs, ys = [], []
    for s in series:
        xs.extend(s["x"]); ys.extend([y for y in s["y"] if y is not None and np.isfinite(y)])
    if not xs:
        xs, ys = [0, 1], [0, 1]
    x0, x1 = min(xs), max(xs)
    y0, y1 = min(ys + [0]), max(ys + [0.1])
    if x1 - x0 < 1e-12:
        x1 = x0 + 1
    cols = ["#b45309", "#0f766e", "#1d4ed8", "#8a8a8a"]
    p = [f'<svg xmlns="http://www.w3.org/2000/svg" width="{w}" height="{h}" viewBox="0 0 {w} {h}">',
         f'<rect width="{w}" height="{h}" fill="#f7f4ef"/>',
         f'<text x="{w/2}" y="26" text-anchor="middle" font-family="-apple-system,PingFang SC,Microsoft YaHei,sans-serif" font-size="14" font-weight="600">{_esc(title)}</text>',
         f'<rect x="{ml}" y="{mt}" width="{pw}" height="{ph}" fill="#fff" stroke="#ccc"/>']
    def sx(x):
        return ml + (x - x0) / (x1 - x0) * pw
    def sy(y):
        return mt + ph - (y - y0) / (y1 - y0 + 1e-15) * ph
    for i, s in enumerate(series):
        pts = [f"{sx(x):.1f},{sy(y):.1f}" for x, y in zip(s["x"], s["y"]) if y is not None and np.isfinite(y)]
        if len(pts) >= 2:
            p.append(f'<polyline fill="none" stroke="{cols[i%len(cols)]}" stroke-width="1.8" points="{" ".join(pts)}"/>')
        p.append(f'<rect x="{w-mr+8}" y="{mt+8+i*16}" width="12" height="3" fill="{cols[i%len(cols)]}"/>')
        p.append(f'<text x="{w-mr+24}" y="{mt+12+i*16}" font-size="10" font-family="sans-serif">{_esc(s["name"])}</text>')
    p.append(f'<text x="{ml+pw/2}" y="{h-16}" text-anchor="middle" font-size="11" font-family="sans-serif">{_esc(xlab)}</text>')
    p.append(f'<text x="{ml}" y="{h-4}" font-size="10" fill="#777" font-family="sans-serif">{_esc(note)} {_esc(ylab)}</text></svg>')
    path.write_text("\n".join(p), encoding="utf-8")


def main():
    t0 = time.time()
    print("=== R3-C1 depth by motion / intensity structure ===", flush=True)
    print(f"OUT={OUT}", flush=True)
    (OUT / "R3_C1_CONFIG.json").write_text(json.dumps(CONFIG, ensure_ascii=False, indent=2), encoding="utf-8")

    # Anchor method doc
    am = [
        "# ANCHOR — Zhu 2023 深海 HLA 移动源深度估计",
        "",
        f"UTC：{NOW}",
        "",
        "## 文献",
        "",
        "Zhu F., Li F., Zhang Y., et al. *Moving source depth estimation in deep ocean direct arrival zone with a horizontal line array.* JASA Express Letters, 2023, 3(9):096003.",
        "",
        "备用锚点：T.C. Yang 2015, synthetic aperture beamforming for moving source depth (PubMed 26428805)。",
        "",
        "## 方法要点（本轮实现口径）",
        "",
        "- **几何**：深海直达区；水平阵 endfire；源沿阵轴方向连续运动、深度恒定。",
        "- **观测**：沿轨迹的声强/干涉结构 I(r) 或 I'(r)=I·r（距离补偿）。",
        "- **机理**：直达 + 海面镜像路径干涉；相位差 Δφ=k(Ls−Ld) 随轨迹与深度变化 → 声强调制。",
        "- **方法 A（优先）**：improved Fourier integral 型深度权 Q_F(z)：用候选深度预测的干涉核与 I' 的匹配。",
        "- **方法 B（论文内对照）**：matched sound intensity structure：J(z)=corr(I'_obs, I'_pred(·|z))。",
        "- **变量**：可用 sin(arrival angle) 或等价轨迹变量；本轮用 r 轨迹上的 I' 相关（Zhu-form，非 P4.5 CZ 轮廓）。",
        "- **原文条件**：CW/可跟踪稳定频率；移动形成空间采样；深度搜索得到模糊峰。",
        "",
        "## paper-condition 场景",
        "",
        "- isovelocity c=1500 m/s；direct + surface image；endfire 移动；CW。",
        "- 验证：true z 峰、track 变长峰变尖、表面/水下可分趋势。",
        "",
        "## 观测条件分级",
        "",
        "| 级 | 定义 | 本轮 |",
        "| --- | --- | --- |",
        "| C-S2 | 稳定线 201/235 Hz + 联合 | 主线 |",
        "| C-S1 | 带内可跟踪机械线 | 对照 |",
        "| C-S0 | 无线谱 | **METHOD_NOT_APPLICABLE_WITHOUT_TRACKABLE_LINE** |",
        "",
        "禁止用 P4.5 简化 CZ 轮廓代替本方法观测形式。",
        "",
    ]
    (OUT / "ANCHOR_ZHU2023_METHOD.md").write_text("\n".join(am), encoding="utf-8")

    print("[1] paper-condition repro ...", flush=True)
    paper_df, paper_pass, err_map = run_paper_repro()
    print("paper pass", paper_pass, err_map, flush=True)

    print("[2] E-STD transfer ...", flush=True)
    amb_df, pair_df, track_df, cmp_df = run_estd()
    amb_df.to_csv(OUT / "estd_depth_ambiguity.csv", index=False, encoding="utf-8-sig")
    pair_df.to_csv(OUT / "depth_pair_distinguishability.csv", index=False, encoding="utf-8-sig")
    track_df.to_csv(OUT / "track_length_depth_boundary.csv", index=False, encoding="utf-8-sig")
    cmp_df.to_csv(OUT / "single_vs_hla_vs_motion.csv", index=False, encoding="utf-8-sig")

    print("[3] figures ...", flush=True)
    # fig1 paper: FWHM vs track for QF vs matched
    p1 = paper_df[(paper_df["z_true_m"] == 55.0) & (paper_df["f_hz"] == 235.0)]
    series = []
    for m in p1["method"].unique():
        d = p1[p1["method"] == m].sort_values("r1_m")
        series.append({"name": m, "x": d["track"].tolist(), "y": d["fwhm_m"].tolist()})
    write_line_svg(FIG / "fig1_paper_repro.svg", "图1  论文条件复现：FWHM vs 轨迹（z=55m）",
                   series, "track", "FWHM m", note="Zhu-like QF / matched structure")
    # fig2 estd ambiguity peak err
    d2 = amb_df[(amb_df["source_level"] == "C-S2_235") & (amb_df["T_s"] == 600.0)]
    if len(d2):
        write_line_svg(FIG / "fig2_estd_depth_ambiguity.svg",
                       "图2  E-STD 深度模糊峰误差 (C-S2 235Hz, T=600s)",
                       [{"name": "err_z", "x": d2["z_true_m"].tolist(), "y": d2["err_z_m"].tolist()},
                        {"name": "FWHM", "x": d2["z_true_m"].tolist(), "y": d2["fwhm_m"].tolist()}],
                       "z_true m", "m", note="vs P4.5 z UNRESOLVED baseline")
    # fig3 pairs
    if len(pair_df):
        write_line_svg(FIG / "fig3_depth_pair_distinguishability.svg",
                       "图3  深度对结构距离 D=1−corr (235Hz,T=600s)",
                       [{"name": "D", "x": list(range(len(pair_df))), "y": pair_df["D_1minus_corr"].tolist()}],
                       "pair idx", "D", note="D大→轨迹声强结构更可分")
    # fig4 track length
    if len(track_df):
        for zt in [200.0]:
            d4 = track_df[track_df["z_true_m"] == zt].sort_values("T_s")
            if len(d4):
                write_line_svg(FIG / "fig4_track_length_boundary.svg",
                               "图4  轨迹长度 vs 深度峰宽 (z=200m,235Hz)",
                               [{"name": "FWHM", "x": d4["track_len_m"].tolist(), "y": d4["fwhm_m"].tolist()},
                                {"name": "|err_z|", "x": d4["track_len_m"].tolist(), "y": d4["err_z_m"].abs().tolist()}],
                               "track length m", "m", note="真实相对位移 U*T，非人为孔径")
    # fig5 static vs motion
    if len(cmp_df):
        write_line_svg(FIG / "fig5_static_vs_motion.svg",
                       "图5  单通道运动 / HLA运动 / HLA静态",
                       [{"name": r["config"], "x": [0], "y": [r["fwhm_m"]]} for _, r in cmp_df.iterrows()],
                       "config", "FWHM m", note="分离物理阵列与合成孔径贡献")

    print("[4] decision ...", flush=True)
    decision, why, nxt, notes = decide(paper_pass, err_map, amb_df, pair_df)
    dec = {
        "rc3c1_decision": decision,
        "why": why,
        "next_step": nxt,
        "notes": notes,
        "created_utc": NOW,
        "b1_closed": "B1_NO_STABLE_MULTIPATH_IDENTITY / PERMANENTLY_CLOSED",
        "stop": "after R3-C1; no C2 / P5 / RC3-B",
    }
    (OUT / "R3_C1_DECISION.json").write_text(json.dumps(dec, ensure_ascii=False, indent=2, default=str), encoding="utf-8")

    def fnum(x, nd=2):
        try:
            v = float(x)
            return "n/a" if not np.isfinite(v) else f"{v:.{nd}f}"
        except Exception:
            return "n/a"

    rp = []
    rp.append("# R3-C1 报告：运动/声强结构深度信息（Zhu 2023 锚定）")
    rp.append("")
    rp.append(f"UTC：{NOW}")
    rp.append("")
    rp.append("B1 已永久关闭：`B1_NO_STABLE_MULTIPATH_IDENTITY`（当前 CZ 仅 1 条能量可观测分支，MMAC 不成立）。")
    rp.append("")
    rp.append("## 0. 锚定与禁止")
    rp.append("")
    rp.append("- 主锚点：**Zhu et al. 2023** JASA-EL 深海 HLA + 移动源深度估计")
    rp.append("- 观测形式：轨迹上的声强/干涉结构（Zhu-form）；**禁止**用 P4.5 CZ 轮廓代替")
    rp.append("- C-S0：`METHOD_NOT_APPLICABLE_WITHOUT_TRACKABLE_LINE`")
    rp.append("")
    rp.append("## 1. 论文条件复现")
    rp.append("")
    rp.append(f"- **pass = {paper_pass}**；各方法平均 \\|err_z\\|：{json.dumps(err_map, ensure_ascii=False)}")
    rp.append("")
    if len(paper_df):
        rp.append("| method | z_true | track | z_hat | err | FWHM | PSL |")
        rp.append("| --- | --- | --- | --- | --- | --- | --- |")
        for _, r in paper_df[paper_df["z_true_m"].isin([50.0, 55.0])].iterrows():
            rp.append(
                f"| {r['method']} | {r['z_true_m']} | {r['track']} | {fnum(r['z_hat_m'])} | "
                f"{fnum(r['err_z_m'])} | {fnum(r['fwhm_m'])} | {fnum(r['psl'])} |"
            )
    rp.append("")
    rp.append("## 2. E-STD 第一会聚区（模态 I' 轨迹 + Zhu 匹配）")
    rp.append("")
    rp.append("| 源 | T s | z_true | z_hat | err | FWHM | PSL | 多峰 |")
    rp.append("| --- | --- | --- | --- | --- | --- | --- | --- |")
    show = amb_df[(amb_df["source_level"] == "C-S2_235") & (amb_df["T_s"] == 600.0)] if len(amb_df) else amb_df
    for _, r in show.iterrows():
        if r.get("source_level") == "C-S0":
            rp.append(f"| C-S0 | — | — | — | — | — | — | {r.get('method')} |")
            continue
        rp.append(
            f"| {r['source_level']} | {fnum(r['T_s'],0)} | {fnum(r['z_true_m'],0)} | {fnum(r['z_hat_m'])} | "
            f"{fnum(r['err_z_m'])} | {fnum(r['fwhm_m'])} | {fnum(r['psl'])} | {r.get('multimodal')} |"
        )
    rp.append("")
    rp.append("### 与 P4.5 基线对照")
    rp.append("")
    rp.append("- P4.5 / P4：z 在简化 CZ 轮廓下总体 **UNRESOLVED / 弱约束**。")
    rp.append(f"- R3-C1 运动声强方法（C-S2 235Hz, T=600s）：mean\\|err_z\\|=**{fnum(notes.get('cz_mean_abs_err_m'))} m**，mean FWHM=**{fnum(notes.get('cz_mean_fwhm_m'))} m**，真值入峰比例=**{fnum(notes.get('cz_true_in_peak'))}**。")
    rp.append(f"- 深度对可分比例（D>0.05）=**{fnum(notes.get('pair_resolvable_frac'))}**，mean D=**{fnum(notes.get('pair_mean_D'))}**。")
    rp.append("")
    rp.append("## 3. 运动长度 / 阵列 vs 运动")
    rp.append("")
    if len(track_df):
        rp.append("| T s | track m | z_true | err_z | FWHM |")
        rp.append("| --- | --- | --- | --- | --- |")
        for _, r in track_df.iterrows():
            rp.append(f"| {fnum(r['T_s'],0)} | {fnum(r['track_len_m'],0)} | {fnum(r['z_true_m'],0)} | {fnum(r['err_z_m'])} | {fnum(r['fwhm_m'])} |")
        rp.append("")
    if len(cmp_df):
        rp.append("| config | err_z | FWHM | PSL |")
        rp.append("| --- | --- | --- | --- |")
        for _, r in cmp_df.iterrows():
            rp.append(f"| {r['config']} | {fnum(r['err_z_m'])} | {fnum(r['fwhm_m'])} | {fnum(r['psl'])} |")
    rp.append("")
    rp.append("## 4. R3-C1 判定")
    rp.append("")
    rp.append(f"### `{decision}`")
    rp.append("")
    rp.append(why)
    rp.append("")
    rp.append(f"**下一步**：{nxt}")
    rp.append("")
    rp.append("允许终态：`C1_PAPER_REPRO_PASS_CZ_DEPTH_CONFIRMED` / `C1_PAPER_REPRO_PASS_CZ_CONDITIONAL` / `C1_PAPER_REPRO_PASS_CZ_NOT_TRANSFERABLE` / `C1_PAPER_REPRO_FAIL`。")
    rp.append("")
    rp.append("## 5. 停止")
    rp.append("")
    rp.append("- 不自动进入 C2；不进 P5；不回 RC3-B")
    rp.append(f"- **R3-C1 完成后停止；判定 `{decision}`**")
    rp.append("")
    (OUT / "R3_C1_REPORT.md").write_text("\n".join(rp), encoding="utf-8")

    gs = [
        "# R3-C1 — GPT 同步稿", "",
        f"- UTC: {NOW}",
        f"- **判定：{decision}**",
        f"- {why}",
        f"- 下一步：{nxt}", "",
        f"Paper repro pass={paper_pass} err={err_map}", "",
        f"C-S0: METHOD_NOT_APPLICABLE_WITHOUT_TRACKABLE_LINE", "",
        "## E-STD C-S2_235 T=600 摘录",
        show[["z_true_m", "z_hat_m", "err_z_m", "fwhm_m", "psl", "multimodal"]].to_string(index=False) if len(show) else "", "",
        f"cz_mean_abs_err={notes.get('cz_mean_abs_err_m')} FWHM={notes.get('cz_mean_fwhm_m')} pair_frac={notes.get('pair_resolvable_frac')}", "",
        "停止：无 C2 / P5 / RC3-B。", "",
    ]
    (OUT / "R3_C1_GPT_SYNC.md").write_text("\n".join(gs), encoding="utf-8")

    print(f"DONE {time.time()-t0:.1f}s DECISION={decision}", flush=True)
    print(why, flush=True)


if __name__ == "__main__":
    main()
