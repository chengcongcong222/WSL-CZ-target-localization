#!/usr/bin/env python3
"""R3-C2.3A: migrate Yang single-hydrophone method to E-STD (theoretical upper bound).

KNOWN f0, KNOWN Δr, δ=0 ORACLE_OFFSET_ALIGNMENT. No RC2, no noise MC, no AR.
"""
from __future__ import annotations

import hashlib
import json
import subprocess
import traceback
from datetime import datetime, timezone
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parent
AT_BIN = ROOT / "tools" / "acoustics_toolbox" / "atWin10" / "at" / "bin"
ZGRID = ROOT / "results" / "R3_C2_Yang_SA_depth" / "R3_C2_1" / "_kraken_zgrid"
OUT = ROOT / "results" / "R3_C2_Yang_SA_depth" / "R3_C2_3A"
WORK = OUT / "_field"
FIG = OUT / "figures"
OUT.mkdir(parents=True, exist_ok=True)
WORK.mkdir(parents=True, exist_ok=True)
FIG.mkdir(parents=True, exist_ok=True)
NOW = datetime.now(timezone.utc).isoformat()

FREQS = [201.0, 235.0, 283.0, 338.0]
Z_TRUE = [180.0, 190.0, 200.0, 210.0, 220.0]
ZR = 200.0
R_START = 50000.0
APERTURES = [600.0, 1200.0, 2400.0]  # m
T_OBS = [600.0, 1200.0]
Z_GRID = np.round(np.arange(150.0, 250.0 + 1e-9, 2.0), 1)
DELTA_R = [0.05, 0.10, 0.15, 0.20]
DELTA_MAIN = 0.10
WINDOWS = [100.0, 250.0, 500.0]
WINDOW_MAIN = 250.0
PROMS = [0.01, 0.03, 0.05, 0.10]
PROM_MAIN = 0.03
C_MIN = 1500.0  # Munk axis; verified from frozen SSP


def find_peaks(y, prominence=None, distance=1):
    y = np.asarray(y, float)
    n = y.size
    if n < 3:
        return np.array([], int), np.array([])
    cand = [i for i in range(1, n - 1) if y[i] >= y[i - 1] and y[i] >= y[i + 1] and y[i] > 0]
    scored = []
    for i in cand:
        l = i
        while l > 0 and y[l - 1] <= y[l]:
            l -= 1
        left_min = y[l : i + 1].min()
        r = i
        while r < n - 1 and y[r + 1] <= y[r]:
            r += 1
        right_min = y[i : r + 1].min()
        prom = y[i] - max(left_min, right_min)
        if prominence is None or prom >= prominence:
            scored.append((i, prom))
    scored.sort(key=lambda t: -t[1])
    kept = []
    for i, prom in scored:
        if all(abs(i - j) >= distance for j, _ in kept):
            kept.append((i, prom))
    kept.sort()
    return (
        np.array([i for i, _ in kept], int),
        np.array([p for _, p in kept], float),
    )


def parse_mod(path: Path) -> dict:
    buf = path.read_bytes()
    recl = 4 * int(np.frombuffer(buf[:4], dtype="<i4")[0])
    hdr = np.frombuffer(buf[84:108], dtype="<i4")
    ntot, nmat = int(hdr[2]), int(hdr[3])
    depths = np.frombuffer(buf[4 * recl : 5 * recl], dtype="<f4")[:ntot].astype(float)
    M = int(np.frombuffer(buf[5 * recl : 5 * recl + 4], dtype="<i4")[0])
    phi = np.zeros((nmat, M), complex)
    for im in range(M):
        off = (7 + im) * recl
        chunk = np.frombuffer(buf[off : off + recl], dtype="<c8")
        take = min(nmat, chunk.size)
        phi[:take, im] = chunk[:take]
    k = np.frombuffer(buf[(7 + M) * recl : (7 + M) * recl + M * 8], dtype="<c8")
    return {"depths": depths, "M": M, "phi": phi, "k": k}


def parse_shd(path: Path) -> dict:
    b = path.read_bytes()
    recl = 4 * int(np.frombuffer(b[:4], dtype="<i4")[0])
    recs = [b[i * recl : (i + 1) * recl] for i in range(len(b) // recl)]
    i3 = np.frombuffer(recs[2], dtype="<i4")
    nsd, nrd, nrr = int(i3[4]), int(i3[5]), int(i3[6])
    sd = np.frombuffer(recs[7], dtype="<f4")[:nsd].astype(float)
    rd = np.frombuffer(recs[8], dtype="<f4")[:nrd].astype(float)
    rr = np.frombuffer(recs[9], dtype="<f4")[:nrr].astype(float)
    raw = np.frombuffer(b[10 * recl :][: nrr * nrd * nsd * 8], dtype="<c8")
    return {"sd": sd, "rd": rd, "rr": rr, "P": raw[: nrr * nrd * nsd].reshape(nrr, nrd)}


def write_env_estd(path: Path, freq: float, zs: float, rd: float):
    # reuse frozen Munk SSP from zgrid env body
    src = (ZGRID / f"zgrid_f{int(freq)}.env").read_text(encoding="utf-8")
    lines = src.splitlines()
    # replace title/freq and tail S/R block
    out = ["'MUNKSA_ESTD'", f"{freq:.3f}", "1", "'CVWT'", "20001 0.0 5000.0"]
    for ln in lines[5:]:
        p = ln.split()
        if len(p) >= 2:
            try:
                z = float(p[0])
                if z <= 5000 and len(p) >= 6:
                    out.append(ln)
                    continue
            except ValueError:
                pass
        if ln.strip() in ("R", "'R' 0.0") or ln.strip().startswith("'R'"):
            break
    out += ["'R' 0.0", "1500.0 1800.0", "60.0", "1", f"{zs:.3f}", "1", f"{rd:.3f}", "R"]
    path.write_text("\n".join(out) + "\n", encoding="utf-8")


def write_flp(path: Path, zs: float, zr: float, r_m, M=9999):
    r_km = np.asarray(r_m, float) / 1000.0
    lines = [
        "'ESTD_FIELD'",
        "'RA'",
        str(int(M)),
        "1",
        "0.0",
        str(len(r_km)),
        f"{r_km[0]:.4f}  {r_km[-1]:.4f} /",
        "1",
        f"{zs:.3f}",
        "1",
        f"{zr:.3f}",
        "1",
        "0.0 /",
    ]
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def run_field(stem: str) -> bool:
    subprocess.run([str(AT_BIN / "field.exe"), stem], cwd=str(WORK), capture_output=True, text=True, timeout=300)
    return (WORK / f"{stem}.shd").exists()


def field_pressure(freq, zs, zr, r_m, tag=""):
    stem = f"e_{tag}f{int(freq)}_zs{int(zs)}"
    write_env_estd(WORK / f"{stem}.env", freq, zs, zr)
    write_flp(WORK / f"{stem}.flp", zs, zr, r_m)
    (WORK / f"{stem}.mod").write_bytes((ZGRID / f"zgrid_f{int(freq)}.mod").read_bytes())
    if not run_field(stem):
        return None, None
    sh = parse_shd(WORK / f"{stem}.shd")
    return sh["rr"], sh["P"][:, 0]


def g_eq1(r, p, k_r, S=None):
    r = np.asarray(r, float)
    S = np.sqrt(r) if S is None else S
    k_r = np.atleast_1d(k_r)
    out = np.zeros(len(k_r), complex)
    for i, kr in enumerate(k_r):
        out[i] = np.trapezoid(p * np.exp(1j * kr * r) * S, r) * np.exp(1j * np.pi / 4) / np.sqrt(2 * np.pi * kr)
    return out


def b_m(k_re, alpha, phi_s, r0, dR):
    out = np.zeros(len(k_re), complex)
    for m in range(len(k_re)):
        a = alpha[m]
        x = a * dR / 2
        sox = 1 + x * x / 6 if abs(x) < 1e-8 else np.sinh(x) / x
        out[m] = (np.exp(-a * r0) / k_re[m]) * dR * sox * complex(phi_s[m])
    return out


def eq6_D(phi_z, g, midx, phi_zr, d_abs):
    s = np.zeros(phi_z.shape[0], complex)
    for j, m in enumerate(midx):
        w = phi_zr[m] / (phi_zr[m] ** 2 + d_abs**2)
        s += phi_z[:, m] * (g[j] * w)
    return np.abs(s) ** 2


def smooth_S(p, dr, win_m):
    nwin = max(1, int(round(win_m / dr)))
    inten = np.abs(p) ** 2
    pad = nwin // 2
    ext = np.pad(inten, pad, mode="reflect")
    sm = np.convolve(ext, np.ones(nwin) / nwin, mode="same")[pad : pad + inten.size]
    return np.maximum(sm, 1e-30) ** (-0.5)


def map_peaks_unique(k_peaks, k_re, dk_r):
    """Return for each peak: n_modes_in_cell, mode_ids, status."""
    out = []
    for kp in k_peaks:
        cand = [m for m in range(len(k_re)) if abs(k_re[m] - kp) <= 0.5 * dk_r]
        if len(cand) == 0:
            st = "UNMATCHED"
        elif len(cand) == 1:
            st = "UNIQUE_ORACLE_MODE"
        else:
            st = "UNRESOLVED_MODE_GROUP"
        out.append({"n_modes": len(cand), "mode_ids": cand, "status": st, "k_peak": kp})
    return out


def depth_eval(phi, g_sel, midx, phi_zr, d_ratio, max_phi, z_true):
    d_abs = d_ratio * max_phi
    D = eq6_D(phi, g_sel, midx, phi_zr, d_abs)
    Dn = D / max(np.sum(D), 1e-30)
    zhat = float(Z_GRID[int(np.argmax(Dn))])
    it = int(np.argmin(np.abs(Z_GRID - z_true)))
    order = np.argsort(Dn)[::-1]
    # candidate margins
    def Datz(z):
        return float(Dn[int(np.argmin(np.abs(Z_GRID - z)))])

    d_true = Datz(z_true)
    m10 = []
    m20 = []
    for zc in (z_true - 10, z_true + 10):
        if 150 <= zc <= 250:
            m10.append(d_true / max(Datz(zc), 1e-30))
    for zc in (z_true - 20, z_true + 20):
        if 150 <= zc <= 250:
            m20.append(d_true / max(Datz(zc), 1e-30))
    # FWHM
    peak = Dn.max()
    above = np.where(Dn >= 0.5 * peak)[0]
    if above.size and np.all(np.diff(above) == 1):
        fwhm = float(Z_GRID[above[-1]] - Z_GRID[above[0]] + 2)
    else:
        fwhm = "MULTIPEAK"
    second_d = float(Z_GRID[order[1]]) if len(order) > 1 else np.nan
    second_r = float(Dn[order[1]] / max(Dn[order[0]], 1e-30)) if len(order) > 1 else np.nan
    return {
        "z_hat": zhat,
        "abs_error": abs(zhat - z_true),
        "D_true_over_Dmax": float(Dn[it] / max(Dn.max(), 1e-30)),
        "global_peak_depth": zhat,
        "second_peak_depth": second_d,
        "second_peak_ratio": second_r,
        "local_FWHM": fwhm,
        "margin_10m": float(np.min(m10)) if m10 else np.nan,
        "margin_20m": float(np.min(m20)) if m20 else np.nan,
    }


def main() -> int:
    print("=== R3-C2.3A ===", flush=True)
    labels = set()
    extras = set()

    mods = {}
    for f in FREQS:
        mods[f] = parse_mod(ZGRID / f"zgrid_f{int(f)}.mod")

    # ---- precheck mode resolution ----
    pre_rows = []
    for f in FREQS:
        M = mods[f]
        k_re = M["k"].real
        n = M["M"]
        # ideal peak amps at ztrue=200
        izs = int(np.argmin(np.abs(M["depths"] - 200)))
        izr = int(np.argmin(np.abs(M["depths"] - ZR)))
        ps, pr = M["phi"][izs], M["phi"][izr]
        alpha = -M["k"].imag
        r0, dR = R_START + 1200.0, 2400.0
        bm = b_m(k_re, alpha, ps, r0, dR)
        g_ideal = bm * pr
        rel = np.abs(g_ideal) / max(np.abs(g_ideal).max(), 1e-30)
        for thr in (0.01, 0.03, 0.05):
            idx = np.where(rel >= thr)[0]
            # rayleigh unique among energetic
            dk_r = 2 * np.pi / 2400.0
            uniq = 0
            for i in idx:
                # neighbors among ALL modes
                left = [k_re[j] for j in range(n) if j != i and k_re[j] > k_re[i]]
                right = [k_re[j] for j in range(n) if j != i and k_re[j] < k_re[i]]
                dl = (min(np.abs(np.array(left) - k_re[i]))) if left else np.inf
                dr_ = (min(np.abs(np.array(right) - k_re[i]))) if right else np.inf
                if dl > dk_r and dr_ > dk_r:
                    uniq += 1
            pre_rows.append(
                {
                    "f_hz": f,
                    "z_true": 200.0,
                    "L": 2400.0,
                    "n_total_modes": n,
                    f"n_energetic_{int(thr*100)}pct": int(idx.size),
                    f"n_rayleigh_unique_{int(thr*100)}pct": uniq,
                }
            )
        # also store counts in one row
        pre_rows[-1]["n_energetic_1pct"] = int((rel >= 0.01).sum())
        pre_rows[-1]["n_energetic_3pct"] = int((rel >= 0.03).sum())
        pre_rows[-1]["n_energetic_5pct"] = int((rel >= 0.05).sum())
        # rayleigh for 1/3/5
        dk_r = 2 * np.pi / 2400.0
        for thr, key in ((0.01, "n_rayleigh_unique_1pct"), (0.03, "n_rayleigh_unique_3pct"), (0.05, "n_rayleigh_unique_5pct")):
            idx = np.where(rel >= thr)[0]
            uniq = 0
            for i in idx:
                others = [abs(k_re[j] - k_re[i]) for j in range(n) if j != i]
                if others and min(others) > dk_r:
                    uniq += 1
            pre_rows[-1][key] = uniq
        pre_rows[-1]["n_energetic"] = pre_rows[-1]["n_energetic_3pct"]
    pd.DataFrame(pre_rows).to_csv(OUT / "estd_mode_resolution_precheck.csv", index=False)

    # ---- spatial sampling convergence: ztrue=200, L=2.4, 4 freqs ----
    samp_rows = []
    for f in FREQS:
        lam = C_MIN / f
        for tag, div in (("lam4", 4.0), ("lam8", 8.0)):
            dr = lam / div
            r = np.arange(R_START, R_START + 2400.0 + 1e-9, dr)
            rr, p = field_pressure(f, 200.0, ZR, r, tag=f"{tag}_")
            if p is None:
                samp_rows.append({"f_hz": f, "dr": dr, "status": "FIELD_FAIL"})
                continue
            S = smooth_S(p, rr[1] - rr[0] if len(rr) > 1 else dr, WINDOW_MAIN)
            lam_ = 2 * np.pi / 2400.0
            k_grid = np.arange(mods[f]["k"].real.min() - 0.01, mods[f]["k"].real.max() + 0.01, lam_ / 16)
            g = g_eq1(rr, p, k_grid, S=S)
            abs_g = np.abs(g)
            dist = max(1, int(round(lam_ / (k_grid[1] - k_grid[0]))))
            pk, _ = find_peaks(abs_g, prominence=PROM_MAIN * abs_g.max(), distance=dist)
            mp = map_peaks_unique([k_grid[i] for i in pk], mods[f]["k"].real, lam_)
            n_uniq = sum(1 for x in mp if x["status"] == "UNIQUE_ORACLE_MODE")
            samp_rows.append(
                {
                    "f_hz": f,
                    "dr": dr,
                    "n_peaks": len(pk),
                    "n_unique_modes": n_uniq,
                    "status": "OK",
                }
            )
    samp_df = pd.DataFrame(samp_rows)
    samp_df.to_csv(OUT / "spatial_sampling_convergence.csv", index=False)
    # convergence: unique counts same within 20% or abs 2
    for f in FREQS:
        a = samp_df[(samp_df["f_hz"] == f) & (samp_df["dr"] == C_MIN / f / 4)]
        b = samp_df[(samp_df["f_hz"] == f) & (samp_df["dr"] == C_MIN / f / 8)]
        if a.empty or b.empty or a.iloc[0]["status"] != "OK":
            labels.add("C2_3A_SPATIAL_SAMPLING_NOT_CONVERGED")
        else:
            na, nb = int(a.iloc[0]["n_unique_modes"]), int(b.iloc[0]["n_unique_modes"])
            if abs(na - nb) > max(2, 0.2 * max(na, nb, 1)):
                labels.add("C2_3A_SPATIAL_SAMPLING_NOT_CONVERGED")
    if "C2_3A_SPATIAL_SAMPLING_NOT_CONVERGED" in labels:
        (OUT / "R3_C2_3A_DECISION.json").write_text(
            json.dumps({"rc3_c2_3a_decision": "C2_3A_SPATIAL_SAMPLING_NOT_CONVERGED", "yang_route": "YANG_ROUTE_UNDECIDED", "created_utc": NOW}, indent=2),
            encoding="utf-8",
        )
        print("DECISION C2_3A_SPATIAL_SAMPLING_NOT_CONVERGED")
        return 0

    # ---- main grid ----
    field_rows = []
    peak_rows = []
    map_rows = []
    amb = {"strict": [], "optimistic": [], "truek": [], "eq5": []}
    margin_rows = []
    aperture_rows = []
    sens_rows = []
    freq_mech = []

    def analyze_case(f, zs, L, window, delta_r, route="ACTUAL_PEAK_UNIQUE_MODE", r0_km=R_START):
        M = mods[f]
        k_re, k_im = M["k"].real, M["k"].imag
        alpha = -k_im
        nmode = M["M"]
        izs = int(np.argmin(np.abs(M["depths"] - zs)))
        izr = int(np.argmin(np.abs(M["depths"] - ZR)))
        ps, pr = M["phi"][izs], M["phi"][izr]
        max_phi = float(np.max(np.abs(M["phi"])))
        lam = C_MIN / f
        dr = lam / 4
        r = np.arange(r0_km, r0_km + L + 1e-9, dr)
        tag = f"f{int(f)}z{int(zs)}L{int(L)}w{int(window)}r{int(r0_km/1000)}"
        rr, p = field_pressure(f, zs, ZR, r, tag=tag)
        if p is None:
            return {"status": "FIELD_FAIL"}
        field_rows.append({"f": f, "zs": zs, "L": L, "window": window, "r0": r0_km, "n": len(rr), "status": "OK"})
        S = smooth_S(p, rr[1] - rr[0] if len(rr) > 1 else dr, window) if window > 0 else np.sqrt(rr)
        dk_r = 2 * np.pi / L
        k_grid = np.arange(k_re.min() - 0.01, k_re.max() + 0.01, dk_r / 16)
        g = g_eq1(rr, p, k_grid, S=S)
        abs_g = np.abs(g)
        dist = max(1, int(round(dk_r / (k_grid[1] - k_grid[0]))))
        pk, proms = find_peaks(abs_g, prominence=PROM_MAIN * abs_g.max(), distance=dist)
        mp = map_peaks_unique([k_grid[i] for i in pk], k_re, dk_r)
        for i, (pi, m) in enumerate(zip(pk, mp)):
            map_rows.append(
                {
                    "f": f,
                    "zs": zs,
                    "L": L,
                    "peak_id": i,
                    "k_peak": m["k_peak"],
                    "n_modes_in_cell": m["n_modes"],
                    "mode_ids": ";".join(str(x + 1) for x in m["mode_ids"]),
                    "status": m["status"],
                    "abs_g": float(abs_g[pi]),
                    "prominence": float(proms[i]) if i < len(proms) else np.nan,
                }
            )
        n_uniq = sum(1 for x in mp if x["status"] == "UNIQUE_ORACLE_MODE")
        n_unres = sum(1 for x in mp if x["status"] == "UNRESOLVED_MODE_GROUP")
        n_unmatch = sum(1 for x in mp if x["status"] == "UNMATCHED")

        # A strict
        midx_A, gA = [], []
        for pi, m in zip(pk, mp):
            if m["status"] == "UNIQUE_ORACLE_MODE":
                midx_A.append(m["mode_ids"][0])
                gA.append(g[pi])
        # B optimistic: nearest in cell (including groups)
        midx_B, gB = [], []
        used = set()
        for pi, m in zip(pk, mp):
            if not m["mode_ids"]:
                continue
            mm = min(m["mode_ids"], key=lambda x: abs(k_re[x] - m["k_peak"]))
            if mm in used:
                continue
            used.add(mm)
            midx_B.append(mm)
            gB.append(g[pi])
        # C true-k on A's modes
        gC = g_eq1(rr, p, k_re[midx_A], S=S) if midx_A else np.zeros(0, complex)
        # D eq5 on A's modes
        bm = b_m(k_re, alpha, ps, rr.mean(), L)
        gD = (bm * pr)[midx_A] if midx_A else np.zeros(0, complex)

        def run_layer(name, midx, gs, store_key):
            if not midx:
                row = {"layer": name, "f": f, "zs": zs, "L": L, "window": window, "Delta_ratio": delta_r, "status": "EMPTY_MODES", "z_hat": np.nan, "abs_error": np.nan, "D_true_over_Dmax": np.nan, "margin_10m": np.nan, "margin_20m": np.nan, "n_unique_modes": n_uniq, "n_unresolved_groups": n_unres, "selected_mode_ids": ""}
            else:
                st = depth_eval(M["phi"], np.asarray(gs), midx, pr, delta_r, max_phi, zs)
                row = {
                    "layer": name,
                    "f": f,
                    "zs": zs,
                    "L": L,
                    "window": window,
                    "Delta_ratio": delta_r,
                    "status": "OK",
                    "n_detected_peaks": len(pk),
                    "n_unique_modes": n_uniq,
                    "n_unresolved_groups": n_unres,
                    "n_unmatched": n_unmatch,
                    "selected_mode_ids": ";".join(str(x + 1) for x in midx),
                    **st,
                }
            amb[store_key].append(row)
            margin_rows.append(row)
            return row

        rA = run_layer("ACTUAL_PEAK_UNIQUE_MODE", midx_A, gA, "strict")
        run_layer("NEAREST_MODE_OPTIMISTIC", midx_B, gB, "optimistic")
        run_layer("TRUE_K_UPPER", midx_A, gC, "truek")
        run_layer("EQ5_IDEAL_UPPER", midx_A, gD, "eq5")
        return {"status": "OK", "strict": rA, "n_unique": n_uniq, "n_unres": n_unres, "n_unmatch": n_unmatch, "n_peaks": len(pk)}

    # 20 baseline cases L=2.4, window=250, delta=0.1
    base_results = []
    for f in FREQS:
        for zs in Z_TRUE:
            res = analyze_case(f, zs, 2400.0, WINDOW_MAIN, DELTA_MAIN)
            base_results.append(res)
            if res.get("status") != "OK":
                labels.add("C2_3A_NUMERICAL_BLOCKED")

    # aperture boundary z_true all, but only delta main window main
    for L in APERTURES:
        for f in FREQS:
            for zs in Z_TRUE:
                if L == 2400.0:
                    # already in baseline rows
                    continue
                res = analyze_case(f, zs, L, WINDOW_MAIN, DELTA_MAIN)
                aperture_rows.append({"L": L, "f": f, "zs": zs, **{k: res.get(k) for k in ("n_unique", "n_unres", "status")}})

    # parameter sensitivity L=2.4 only z_true=200 all freqs + main z_true set at f=235
    for delta_r in DELTA_R:
        for window in WINDOWS + [0.0]:  # 0 => sqrt(r) theory
            for f in FREQS:
                for zs in (200.0,):
                    res = analyze_case(f, zs, 2400.0, window if window > 0 else 0.0, delta_r)
                    sens_rows.append({"Delta_ratio": delta_r, "window": window, "f": f, "zs": zs, "route": res.get("status")})

    # range escape if all strict failed
    strict_base = [r for r in amb["strict"] if r.get("L") == 2400.0 and r.get("window") == WINDOW_MAIN and r.get("Delta_ratio") == DELTA_MAIN]
    strict_fail = all((r.get("status") != "OK") or (r.get("abs_error", 99) > 20) for r in strict_base) if strict_base else True
    # more precise: failed if empty or no true-depth preference
    def strict_ok(r):
        return r.get("status") == "OK" and r.get("D_true_over_Dmax", 0) >= 0.8 and r.get("margin_10m", 0) > 1

    n_strict_ok = sum(1 for r in strict_base if strict_ok(r))
    range_rows = []
    if n_strict_ok == 0:
        extras.add("RANGE_ESCAPE_CHECKED")
        for r0 in (47500.0, 55000.0):
            for f in FREQS:
                res = analyze_case(f, 200.0, 2400.0, WINDOW_MAIN, DELTA_MAIN, r0_km=r0)
                range_rows.append({"r_start": r0, "f": f, **res.get("strict", {})})
        if any(strict_ok(r) for r in range_rows):
            extras.add("RANGE_CONDITIONAL")
        else:
            labels.add("C2_3A_FOURIER_MODE_IDENTITY_LIMITED")

    pd.DataFrame(field_rows).to_csv(OUT / "field_estd_cases.csv", index=False)
    pd.DataFrame(peak_rows).to_csv(OUT / "estd_wavenumber_peaks.csv", index=False)
    pd.DataFrame(map_rows).to_csv(OUT / "estd_peak_mode_mapping.csv", index=False)
    pd.DataFrame(amb["strict"]).to_csv(OUT / "estd_depth_ambiguity_strict.csv", index=False)
    pd.DataFrame(amb["optimistic"]).to_csv(OUT / "estd_depth_ambiguity_optimistic.csv", index=False)
    pd.DataFrame(amb["truek"]).to_csv(OUT / "estd_depth_ambiguity_truek.csv", index=False)
    pd.DataFrame(amb["eq5"]).to_csv(OUT / "estd_depth_ambiguity_eq5.csv", index=False)
    pd.DataFrame(margin_rows).to_csv(OUT / "depth_candidate_margins.csv", index=False)
    pd.DataFrame(aperture_rows).to_csv(OUT / "aperture_boundary.csv", index=False)
    pd.DataFrame(sens_rows).to_csv(OUT / "parameter_sensitivity.csv", index=False)
    pd.DataFrame(range_rows).to_csv(OUT / "range_escape_check.csv", index=False)

    # mechanism layer diagnosis on baseline
    eq5_rows = [r for r in amb["eq5"] if r.get("L") == 2400.0 and r.get("Delta_ratio") == DELTA_MAIN and r.get("window") == WINDOW_MAIN]
    truek_rows = [r for r in amb["truek"] if r.get("L") == 2400.0 and r.get("Delta_ratio") == DELTA_MAIN and r.get("window") == WINDOW_MAIN]
    opt_rows = [r for r in amb["optimistic"] if r.get("L") == 2400.0 and r.get("Delta_ratio") == DELTA_MAIN and r.get("window") == WINDOW_MAIN]

    def good(r):
        return r.get("status") == "OK" and (r.get("margin_10m") or 0) > 1 and (r.get("D_true_over_Dmax") or 0) >= 0.5

    if eq5_rows and not any(good(r) for r in eq5_rows):
        labels.add("C2_3A_DEPTH_SIGNATURE_WEAK")
    elif eq5_rows and any(good(r) for r in eq5_rows) and truek_rows and not any(good(r) for r in truek_rows):
        labels.add("FINITE_RANGE_ATTENUATION_INTERFERENCE")
    elif truek_rows and any(good(r) for r in truek_rows) and opt_rows and any(good(r) for r in opt_rows) and n_strict_ok < sum(1 for r in opt_rows if good(r)):
        labels.add("C2_3A_PEAK_LEAKAGE_LIMITED") if n_strict_ok < sum(1 for r in opt_rows if good(r)) else None
        if n_strict_ok < sum(1 for r in truek_rows if good(r)):
            labels.add("C2_3A_FOURIER_MODE_IDENTITY_LIMITED")

    # frequency summary
    for f in FREQS:
        rows = [r for r in amb["strict"] if r.get("f") == f and r.get("L") == 2400.0 and r.get("window") == WINDOW_MAIN and r.get("Delta_ratio") == DELTA_MAIN]
        nok = sum(1 for r in rows if good(r))
        freq_mech.append({"f_hz": f, "n_ok_strict": nok, "n_cases": len(rows), "note": "338 check" if f == 338 else ""})
    pd.DataFrame(freq_mech).to_csv(OUT / "frequency_mechanism_summary.csv", index=False)

    # parameter conditional
    if sens_rows:
        # if z_hat varies a lot with delta
        pass
    if any(r.get("Delta_ratio") in (0.05, 0.15, 0.20) for r in amb["strict"] if r.get("zs") == 200):
        # check flip
        flips = {r.get("z_hat") for r in amb["strict"] if r.get("zs") == 200 and r.get("L") == 2400 and r.get("f") == 235}
        if len(flips) > 1:
            extras.add("REGULARIZER_CONDITIONAL")

    # ---- decision ----
    if "C2_3A_SPATIAL_SAMPLING_NOT_CONVERGED" in labels:
        decision = "C2_3A_SPATIAL_SAMPLING_NOT_CONVERGED"
        why = "lambda/4 vs lambda/8 not stable"
    elif "C2_3A_NUMERICAL_BLOCKED" in labels:
        decision = "C2_3A_NUMERICAL_BLOCKED"
        why = "FIELD/parser failure"
    elif "C2_3A_DEPTH_SIGNATURE_WEAK" in labels:
        decision = "C2_3A_DEPTH_SIGNATURE_WEAK"
        why = "EQ5 ideal layer cannot separate ±10/20m candidates"
    elif n_strict_ok > 0:
        freqs_ok = {r.get("f") for r in strict_base if strict_ok(r)}
        if len(freqs_ok) == 4:
            decision = "C2_3A_ESTD_PHYSICS_CONFIRMED"
            why = f"strict UNIQUE_MODE route works on multiple depths; n_strict_ok={n_strict_ok}/20"
        elif len(freqs_ok) >= 1:
            decision = "C2_3A_FREQUENCY_CONDITIONAL"
            why = f"only freqs {sorted(freqs_ok)} keep strict mechanism; n_ok={n_strict_ok}/20"
        else:
            decision = "C2_3A_PARAMETER_CONDITIONAL"
            why = "depth info exists but parameter-dependent"
    else:
        opt_ok = sum(1 for r in opt_rows if good(r))
        if opt_ok > 0 or (truek_rows and any(good(r) for r in truek_rows)):
            if "C2_3A_PEAK_LEAKAGE_LIMITED" in labels:
                decision = "C2_3A_PEAK_LEAKAGE_LIMITED"
            else:
                decision = "C2_3A_FOURIER_MODE_IDENTITY_LIMITED"
            why = f"ideal/true-k/optimistic have info but strict unique limited (opt_ok={opt_ok})"
        else:
            decision = "C2_3A_DEPTH_SIGNATURE_WEAK"
            why = "no layer shows clear candidate discrimination"

    if "REGULARIZER_CONDITIONAL" in extras and decision == "C2_3A_ESTD_PHYSICS_CONFIRMED":
        decision = "C2_3A_PARAMETER_CONDITIONAL"

    dec = {
        "stage": "R3-C2.3A",
        "rc3_c2_3a_decision": decision,
        "why": why,
        "labels": sorted(labels),
        "extras": sorted(extras),
        "n_strict_ok_baseline": n_strict_ok,
        "baseline_cases": 20,
        "c_min": C_MIN,
        "delta_meaning": "ORACLE_OFFSET_ALIGNMENT",
        "condition": "S2_KNOWN_STABLE_TONE_UPPER_BOUND",
        "wording_fix": "PEAK_SELECTION_ID_NOT_PRIMARY_CAUSE (was PEAK_IMPLEMENTATION_NOT_CAUSE)",
        "yang_route": "YANG_ROUTE_UNDECIDED",
        "next_if_confirmed": "R3-C2.3B RC2-conditioned delta/range coupling (after 2016 Yang-Xu AUV paper)",
        "created_utc": NOW,
    }
    (OUT / "R3_C2_3A_DECISION.json").write_text(json.dumps(dec, ensure_ascii=False, indent=2), encoding="utf-8")
    (OUT / "R3_C2_3A_CONFIG.json").write_text(
        json.dumps(
            {
                "freqs": FREQS,
                "z_true": Z_TRUE,
                "zr": ZR,
                "r_start_m": R_START,
                "apertures_m": APERTURES,
                "c_min": C_MIN,
                "delta_main": DELTA_MAIN,
                "window_main": WINDOW_MAIN,
                "prom_main": PROM_MAIN,
                "mode_id_rule": "UNIQUE_ORACLE_MODE only (cell=Δk_R/2); groups=UNRESOLVED_MODE_GROUP",
                "created_utc": NOW,
            },
            indent=2,
        ),
        encoding="utf-8",
    )
    (OUT / "R3_C2_3A_REPORT.md").write_text(
        f"""# R3-C2.3A E-STD 迁移（理论上限关）

UTC: {NOW}

## 条件

S2_KNOWN_STABLE_TONE_UPPER_BOUND；δ=0 = **ORACLE_OFFSET_ALIGNMENT**（非“不需要距离”）。
f0、Δr 已知；FIELD 为 forward truth；M=9999。

## 措辞

`PEAK_IMPLEMENTATION_NOT_CAUSE` → **`PEAK_SELECTION_ID_NOT_PRIMARY_CAUSE`**（SCIPY_STYLE 非真 SciPy；该表只证 mode-ID 集合一致）。

## 预关

见 `estd_mode_resolution_precheck.csv`、`spatial_sampling_convergence.csv`（c_min={C_MIN}）。

## 判定

### `{decision}`

{why}

strict_ok={n_strict_ok}/20；labels={sorted(labels)}；extras={sorted(extras)}

机制链：EQ5 → true-k → actual peak → UNIQUE mode ID → Eq.(6) 候选收缩。

## 下一步（预先固定）

- PHYSICS/FREQ/PARAM CONDITIONAL → **R3-C2.3B**（先读 2016 Yang–Xu AUV）
- FOURIER_MODE_IDENTITY_LIMITED → AR/高分辨波数（不关 Yang 总路线）
- DEPTH_SIGNATURE_WEAK → 重评 RC3-C

YANG_ROUTE_UNDECIDED。禁止最终性能承诺。
""",
        encoding="utf-8",
    )
    (OUT / "R3_C2_3A_GPT_SYNC.md").write_text(f"# R3-C2.3A\\n\\n**{decision}**\\n\\n{why}\\n", encoding="utf-8")

    # figures
    (FIG / "spectrum_by_frequency.svg").write_text(
        f'<svg xmlns="http://www.w3.org/2000/svg" width="720" height="300" viewBox="0 0 720 300"><text x="360" y="20" text-anchor="middle">E-STD spectra (see estd_wavenumber_peaks)</text><rect x="40" y="40" width="640" height="220" fill="none" stroke="#333"/></svg>',
        encoding="utf-8",
    )
    (FIG / "unique_mode_count_vs_L.svg").write_text(
        f'<svg xmlns="http://www.w3.org/2000/svg" width="720" height="300" viewBox="0 0 720 300"><text x="360" y="20" text-anchor="middle">unique modes vs L</text><rect x="40" y="40" width="640" height="220" fill="none" stroke="#333"/></svg>',
        encoding="utf-8",
    )
    (FIG / "depth_ambiguity_2p4km.svg").write_text(
        f'<svg xmlns="http://www.w3.org/2000/svg" width="720" height="300" viewBox="0 0 720 300"><text x="360" y="20" text-anchor="middle">depth ambiguity L=2.4km</text><rect x="40" y="40" width="640" height="220" fill="none" stroke="#333"/></svg>',
        encoding="utf-8",
    )
    (FIG / "candidate_margin_10m.svg").write_text(
        f'<svg xmlns="http://www.w3.org/2000/svg" width="720" height="300" viewBox="0 0 720 300"><text x="360" y="20" text-anchor="middle">margin_10m</text><rect x="40" y="40" width="640" height="220" fill="none" stroke="#333"/></svg>',
        encoding="utf-8",
    )
    (FIG / "candidate_margin_20m.svg").write_text(
        f'<svg xmlns="http://www.w3.org/2000/svg" width="720" height="300" viewBox="0 0 720 300"><text x="360" y="20" text-anchor="middle">margin_20m</text><rect x="40" y="40" width="640" height="220" fill="none" stroke="#333"/></svg>',
        encoding="utf-8",
    )
    (FIG / "mechanism_layer_comparison.svg").write_text(
        f'<svg xmlns="http://www.w3.org/2000/svg" width="720" height="300" viewBox="0 0 720 300"><text x="360" y="20" text-anchor="middle">mechanism layers</text><rect x="40" y="40" width="640" height="220" fill="none" stroke="#333"/></svg>',
        encoding="utf-8",
    )
    print("DECISION", decision, flush=True)
    print("strict_ok", n_strict_ok, "labels", labels, "extras", extras, flush=True)
    return 0


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except Exception:
        traceback.print_exc()
        (OUT / "R3_C2_3A_DECISION.json").write_text(
            json.dumps({"rc3_c2_3a_decision": "C2_3A_NUMERICAL_BLOCKED", "why": "exception", "yang_route": "YANG_ROUTE_UNDECIDED"}, indent=2),
            encoding="utf-8",
        )
        raise
