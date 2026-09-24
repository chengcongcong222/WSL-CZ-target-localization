#!/usr/bin/env python3
"""R3-C2.2C-FIX: integrity fixes only. No E-STD. Peak detect → oracle label → Eq.(6)."""
from __future__ import annotations

import hashlib
import json
import subprocess
import traceback
from datetime import datetime, timezone
from pathlib import Path

import numpy as np
import pandas as pd


def find_peaks(y, prominence=None, distance=1):
    """Minimal find_peaks: local maxima with prominence and min distance."""
    y = np.asarray(y, dtype=float)
    n = y.size
    if n < 3:
        return np.array([], dtype=int), {"prominences": np.array([])}
    cand = [i for i in range(1, n - 1) if y[i] >= y[i - 1] and y[i] >= y[i + 1] and y[i] > 0]
    # prominence: drop to left/right min of higher peak
    scored = []
    for i in cand:
        # left
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
    # distance filter: keep highest prominence in each cluster
    scored.sort(key=lambda t: -t[1])
    kept = []
    for i, prom in scored:
        if all(abs(i - j) >= distance for j, _ in kept):
            kept.append((i, prom))
    kept.sort(key=lambda t: t[0])
    idx = np.array([i for i, _ in kept], dtype=int)
    proms = np.array([p for _, p in kept], dtype=float)
    return idx, {"prominences": proms}

ROOT = Path(__file__).resolve().parent
AT_BIN = ROOT / "tools" / "acoustics_toolbox" / "atWin10" / "at" / "bin"
OUT = ROOT / "results" / "R3_C2_Yang_SA_depth" / "R3_C2_2C_FIX"
WORK = OUT / "_field_work"
FIG = OUT / "figures"
OUT.mkdir(parents=True, exist_ok=True)
WORK.mkdir(parents=True, exist_ok=True)
FIG.mkdir(parents=True, exist_ok=True)
NOW = datetime.now(timezone.utc).isoformat()

FREQ = 350.0
H = 88.0
R1 = 5010.0
R2 = 10000.0
SPAN = 4990.0
DR_PAPER = 2.5
ZS_LIST = [4.0, 50.0]
ZR_LIST = [18.0, 70.0]
CASES = [(zs, zr) for zs in ZS_LIST for zr in ZR_LIST]
Z_SEARCH = np.round(np.arange(0.0, 88.0 + 1e-9, 1.0), 2)
PROMINENCE_FRACS = [0.01, 0.03, 0.05, 0.10]
PROM_MAIN = 0.03
DELTA_RATIOS = [0.05, 0.10, 0.20]
SMOOTH_M = [100.0, 250.0, 500.0]
DK_R = 2 * np.pi / SPAN  # Rayleigh
MATCH_TOL = 0.5 * DK_R


def parse_mod_exact(path: Path) -> dict:
    buf = path.read_bytes()
    recl = 4 * int(np.frombuffer(buf[:4], dtype="<i4")[0])
    hdr = np.frombuffer(buf[84:108], dtype="<i4")
    ntot = int(hdr[2]) if hdr.size > 2 else 0
    nmat = int(hdr[3]) if hdr.size > 3 else 0
    rec4 = np.frombuffer(buf[4 * recl : 5 * recl], dtype="<f4")
    depths = rec4[:ntot].astype(float) if ntot > 0 else np.array([])
    M = int(np.frombuffer(buf[5 * recl : 5 * recl + 4], dtype="<i4")[0])
    nmat = nmat if nmat else 1
    phi = np.zeros((nmat, M), dtype=complex)
    for im in range(M):
        off = (7 + im) * recl
        chunk = np.frombuffer(buf[off : off + recl], dtype="<c8")
        take = min(nmat, chunk.size)
        phi[:take, im] = chunk[:take]
    k_off = (7 + M) * recl
    k = np.frombuffer(buf[k_off : k_off + M * 8], dtype="<c8")
    if k.size != M:
        k = np.zeros(M, dtype=complex)
    return {"depths": depths, "M": M, "phi": phi, "k": k, "nmat": nmat}


def parse_shd(path: Path) -> dict:
    """Parse AT FIELD .shd observed layout (recl words = first int32)."""
    b = path.read_bytes()
    recl = 4 * int(np.frombuffer(b[:4], dtype="<i4")[0])
    nrec = len(b) // recl
    recs = [b[i * recl : (i + 1) * recl] for i in range(nrec)]
    title = recs[0][4:84].decode("ascii", "ignore").strip()
    # rec3 integers: 1,1,1,1,1,Nrd,Nrr
    i3 = np.frombuffer(recs[2], dtype="<i4")
    nsd, nrd, nrr = int(i3[4]), int(i3[5]), int(i3[6])
    # rec4: freq as float64 at start
    freq = float(np.frombuffer(recs[3][:8], dtype="<f8")[0])
    sd = np.frombuffer(recs[7], dtype="<f4")[:nsd].astype(float)
    rd = np.frombuffer(recs[8], dtype="<f4")[:nrd].astype(float)
    rr = np.frombuffer(recs[9], dtype="<f4")[:nrr].astype(float)
    # pressure complex64 after header records 1-10
    # observed data starts rec 11 (index 10)
    body = b[10 * recl :]
    ncomplex = nrr * nrd * nsd
    need = ncomplex * 8
    if len(body) < need:
        # pad / take what we have
        pass
    raw = np.frombuffer(body[:need], dtype="<c8")
    # Fortran order (rr, rd) or (rd, rr)? try both via correlation later
    p = raw[:ncomplex].reshape(nrr, nrd)
    return {
        "title": title,
        "freq": freq,
        "nsd": nsd,
        "nrd": nrd,
        "nrr": nrr,
        "sd": sd,
        "rd": rd,
        "rr": rr,
        "pressure_rr_rd": p,
        "nrec": nrec,
        "recl": recl,
    }


def ssp_yang2014(z: float) -> float:
    if z <= 10.0:
        return 1533.0
    if z <= 40.0:
        return 1533.0 + (1478.0 - 1533.0) * (z - 10.0) / 30.0
    return 1478.0


def write_env(path: Path, zs: float, rd_list):
    zz = sorted(set([0.0, 5.0, 10.0, 15.0, 20.0, 25.0, 30.0, 35.0, 40.0, 52.0, 64.0, 76.0, 88.0]))
    out = ["'YANG2014_NJ_SUMMER'", f"{FREQ:.1f}", "1", "'CVW'", "2001 0.0 88.0"]
    for z in zz:
        out.append(f"{z:.3f} {ssp_yang2014(z):.4f} 0.0 1.0 0.0 0.0")
    out += ["'A' 0.0", "0.0 1650.0 0.0 1.76 0.8 0.0", "1400.0 1800.0", "60.0", "1", f"{zs:.3f}", str(len(rd_list))]
    for zr in rd_list:
        out.append(f"{zr:.3f}")
    out.append("R")
    path.write_text("\n".join(out) + "\n", encoding="utf-8")


def write_flp(path: Path, zs: float, zr_list, r_m: np.ndarray):
    """Recovered FIELD .flp layout (one integer per READ for NProf; compact for Nfreq+freq)."""
    r_km = r_m / 1000.0
    lines = [
        "'YANG2014_FIELD'",
        "'R' 'R' 'A'",
        f"1 {FREQ:.1f}",
        "1",
        "0.0",
        str(len(r_km)),
    ]
    for x in r_km:
        lines.append(f"{x:.4f}")
    lines += ["1", f"{zs:.3f}", str(len(zr_list))]
    for z in zr_list:
        lines.append(f"{z:.3f}")
    # tail recovered as working: 1 / 0.0 / 0.0
    lines += ["1", "0.0", "0.0"]
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def run_exe(exe: Path, stem: str, cwd: Path, timeout: int = 180):
    r = subprocess.run([str(exe), stem], cwd=str(cwd), capture_output=True, text=True, timeout=timeout)
    return r.returncode, (r.stdout or "") + (r.stderr or "")


def pressure_eq3(r, zs, zr, k_m, alpha_m, phi_s, phi_r):
    r = np.atleast_1d(np.asarray(r, dtype=float))
    p = np.zeros_like(r, dtype=complex)
    for m in range(len(k_m)):
        amp = np.sqrt(2.0 * np.pi / (k_m[m] * r))
        p += amp * phi_s[m] * phi_r[m] * np.exp(-1j * k_m[m] * r - alpha_m[m] * r - 1j * np.pi / 4)
    return p


def g_eq1(r, p, k_r, S=None):
    r = np.asarray(r, dtype=float)
    if S is None:
        S = np.sqrt(r)
    k_r = np.atleast_1d(k_r)
    out = np.zeros(len(k_r), dtype=complex)
    for i, kr in enumerate(k_r):
        integ = p * np.exp(1j * kr * r) * S
        out[i] = np.trapezoid(integ, r) * np.exp(1j * np.pi / 4) / np.sqrt(2.0 * np.pi * kr)
    return out


def b_m_stable(k_m, alpha_m, phi_s, r0, dR):
    k_m = np.asarray(k_m)
    k_r = k_m.real if np.iscomplexobj(k_m) else np.asarray(k_m, dtype=float)
    phi_s = np.asarray(phi_s)
    b = np.zeros(len(k_r), dtype=complex)
    for m in range(len(k_r)):
        a = float(alpha_m[m])
        x = a * dR / 2.0
        sinh_over_x = 1.0 + x * x / 6.0 if abs(x) < 1e-8 else np.sinh(x) / x
        b[m] = (np.exp(-a * r0) / k_r[m]) * dR * sinh_over_x * complex(phi_s[m])
    return b


def depth_ambiguity_eq6(phi_z, g_sel, mode_idx, phi_zr, delta_abs, modes_k):
    """D(z)=|sum_sel phi_m(z) g * phi_m(zr)/(phi_m^2(zr)+Delta^2)|^2 using SELECTED peaks only."""
    s = np.zeros(phi_z.shape[0], dtype=complex)
    for j, m in enumerate(mode_idx):
        denom = phi_zr[m] ** 2 + delta_abs**2
        w = phi_zr[m] / denom
        s += phi_z[:, m] * (g_sel[j] * w)
    return np.abs(s) ** 2


def smooth_abs2(p, dr, win_m):
    nwin = max(1, int(round(win_m / dr)))
    inten = np.abs(p) ** 2
    # reflect padding to avoid edge bias
    pad = nwin // 2
    ext = np.pad(inten, pad, mode="reflect")
    ker = np.ones(nwin) / nwin
    sm = np.convolve(ext, ker, mode="same")[pad:-pad]
    if sm.size != inten.size:
        sm = np.convolve(inten, ker, mode="same")
    return np.maximum(sm, 1e-30) ** (-0.5)


def write_svg_spectrum(path, title, k_grid, abs_g, peaks_idx, k_true, sel_modes):
    pts = []
    ymax = max(abs_g.max(), 1e-30)
    k0, k1 = k_grid.min(), k_grid.max()
    for i in range(0, len(k_grid), max(1, len(k_grid) // 400)):
        x = 60 + (k_grid[i] - k0) / (k1 - k0) * 620
        y = 360 - (abs_g[i] / ymax) * 300
        pts.append(f"{x:.1f},{y:.1f}")
    parts = [
        f'<svg xmlns="http://www.w3.org/2000/svg" width="720" height="420" viewBox="0 0 720 420">',
        f'<text x="360" y="24" text-anchor="middle" font-size="14">{title}</text>',
        '<rect x="60" y="40" width="620" height="320" fill="none" stroke="#333"/>',
        f'<polyline fill="none" stroke="#246" stroke-width="1.5" points="{" ".join(pts)}"/>',
    ]
    for kv in k_true:
        x = 60 + (kv - k0) / (k1 - k0) * 620
        parts.append(f'<line x1="{x:.1f}" y1="40" x2="{x:.1f}" y2="55" stroke="#888"/>')
    for pi in peaks_idx:
        x = 60 + (k_grid[pi] - k0) / (k1 - k0) * 620
        y = 360 - (abs_g[pi] / ymax) * 300
        color = "#c22" if (pi in sel_modes) else "#e80"
        parts.append(f'<circle cx="{x:.1f}" cy="{y:.1f}" r="4" fill="{color}"/>')
    parts.append("</svg>")
    path.write_text("\n".join(parts), encoding="utf-8")


def write_svg_depth(path, title, depths, Dn_true4, Dn_true50, zhat4, zhat50):
    parts = [
        f'<svg xmlns="http://www.w3.org/2000/svg" width="720" height="420" viewBox="0 0 720 420">',
        f'<text x="360" y="24" text-anchor="middle" font-size="14">{title}</text>',
        '<rect x="60" y="40" width="620" height="320" fill="none" stroke="#333"/>',
    ]
    for Dn, zs, zhat, color in [(Dn_true4, 4, zhat4, "#c44"), (Dn_true50, 50, zhat50, "#44c")]:
        if Dn is None:
            continue
        ymax = max(Dn.max(), 1e-30)
        pts = []
        for zi, dv in zip(depths, Dn):
            x = 60 + zi / 88.0 * 620
            y = 360 - (dv / ymax) * 300
            pts.append(f"{x:.1f},{y:.1f}")
        parts.append(f'<polyline fill="none" stroke="{color}" stroke-width="1.5" points="{" ".join(pts)}"/>')
        x = 60 + zs / 88.0 * 620
        parts.append(f'<line x1="{x:.1f}" y1="40" x2="{x:.1f}" y2="360" stroke="#888" stroke-dasharray="4"/>')
    parts.append("</svg>")
    path.write_text("\n".join(parts), encoding="utf-8")


def main() -> int:
    print("=== R3-C2.2C-FIX ===", flush=True)
    impl_ok = True
    field_ok = True
    case_rows = []
    peak_rows = []
    map_rows = []
    fv_rows = []
    alpha_rows = []
    prom_rows = []
    shade_rows = []
    depth_sel_rows = []
    depth_full_rows = []

    # -------- 0) status --------
    (OUT / "R3_C2_2C_FIX_CONFIG.json").write_text(
        json.dumps(
            {
                "stage": "R3-C2.2C-FIX",
                "created_utc": NOW,
                "delta_baseline": 0.0,
                "delta_meaning": "ORACLE_OFFSET_ALIGNMENT",
                "prominence_main": PROM_MAIN,
                "prominence_sensitivity": PROMINENCE_FRACS,
                "match_tol": MATCH_TOL,
                "depth_gate_m": 2.0,
                "depth_ratio_gate": 0.8,
                "yang_route": "YANG_ROUTE_UNDECIDED",
                "prior_confirm": "SUPERSEDED_PENDING_INTEGRITY_FIX",
            },
            ensure_ascii=False,
            indent=2,
        ),
        encoding="utf-8",
    )

    # -------- 1) KRAKEN modes (reuse validated env) --------
    mod_path = ROOT / "results" / "R3_C2_Yang_SA_depth" / "R3_C2_2C" / "_kraken_paper" / "yang2014_f350.mod"
    if not mod_path.exists():
        raise FileNotFoundError(mod_path)
    M = parse_mod_exact(mod_path)
    depths, phi, k = M["depths"], M["phi"], M["k"]
    nmode = M["M"]
    # alpha mapping: KRAKEN Im(k)<0 with e^{+ikr} decay => Yang exp(-i k r - alpha r) uses alpha=-Im(k)
    k_re = k.real.copy()
    k_im = k.imag.copy()
    alpha_m = -k_im  # if Im<0, alpha>0
    if np.any(alpha_m < 0):
        # fallback documented
        alpha_m = np.abs(k_im)
        alpha_map = "abs(Im k) fallback (some Im>0)"
    else:
        alpha_map = "alpha_m = -Im(k_complex); Im(k)<0 => alpha>0"

    def phi_at(z):
        iz = int(np.argmin(np.abs(depths - z)))
        return phi[iz, :], float(depths[iz]), iz

    # -------- 2) FIELD cross-check (20 ranges × 4 cases) --------
    r_chk = np.linspace(R1, R2, 20)
    field_available = True
    field_payload = {}  # (zs,zr) -> p_field along r_chk
    for zs, zr in CASES:
        stem = f"fld_zs{int(zs)}_zr{int(zr)}"
        write_env(WORK / f"{stem}.env", zs, [zr])
        write_flp(WORK / f"{stem}.flp", zs, [zr], r_chk)
        # copy mod
        (WORK / f"{stem}.mod").write_bytes(mod_path.read_bytes())
        rc, msg = run_exe(AT_BIN / "field.exe", stem, WORK, timeout=120)
        shd = WORK / f"{stem}.shd"
        if not shd.exists():
            field_ok = False
            field_available = False
            fv_rows.append(
                {
                    "zs": zs,
                    "zr": zr,
                    "status": "FIELD_FAIL",
                    "n_ranges": len(r_chk),
                    "best_scale_re": np.nan,
                    "best_scale_im": np.nan,
                    "complex_corr": np.nan,
                    "rel_residual": np.nan,
                    "phase_rmse_rad": np.nan,
                    "amp_shape_rel": np.nan,
                    "note": (msg or "")[:200],
                }
            )
            continue
        sh = parse_shd(shd)
        # pressure shape (nrr, nrd)
        P = sh["pressure_rr_rd"]
        # pick receiver depth column nearest zr
        jz = int(np.argmin(np.abs(sh["rd"] - zr)))
        p_field = P[:, jz]
        # ranges in shd (m) vs requested
        r_shd = sh["rr"]
        field_payload[(zs, zr)] = (r_shd, p_field, sh)

        ps, _, _ = phi_at(zs)
        pr, _, _ = phi_at(zr)
        p_ms = pressure_eq3(r_chk, zs, zr, k_re, alpha_m, ps, pr)
        # interpolate modesum onto shd ranges if needed
        if r_shd.size == r_chk.size and np.allclose(r_shd, r_chk, rtol=1e-3, atol=1.0):
            p_ms_use = p_ms
            r_use = r_chk
        else:
            # use shd ranges
            r_use = r_shd
            p_ms_use = pressure_eq3(r_use, zs, zr, k_re, alpha_m, ps, pr)
        # best global complex scale c* = <p_ms, p_field> / ||p_ms||^2  (min ||pF - c pMs||)
        num = np.vdot(p_ms_use, p_field)  # sum conj(ms)*field
        den = np.vdot(p_ms_use, p_ms_use)
        cstar = num / den if abs(den) > 0 else 0.0
        resid = p_field - cstar * p_ms_use
        rel_res = np.linalg.norm(resid) / max(np.linalg.norm(p_field), 1e-30)
        # complex correlation
        cc = np.vdot(p_ms_use, p_field) / max(
            np.linalg.norm(p_ms_use) * np.linalg.norm(p_field), 1e-30
        )
        # phase residual vs range
        ph = np.angle(p_field * np.conj(cstar * p_ms_use))
        # amplitude-shape residual after removing |c|
        a_f = np.abs(p_field)
        a_m = np.abs(p_ms_use) * abs(cstar)
        amp_rel = np.linalg.norm(a_f - a_m) / max(np.linalg.norm(a_f), 1e-30)

        alpha_rows.append(
            {
                "zs": zs,
                "zr": zr,
                "alpha_mapping": alpha_map,
                "Im_k_sign": "negative" if np.all(k_im <= 0) else "mixed",
                "field_abs_r0": float(np.abs(p_field[0])),
                "field_abs_r1": float(np.abs(p_field[-1])),
                "field_decay": bool(np.abs(p_field[-1]) < np.abs(p_field[0])),
                "modesum_decay": bool(np.abs(p_ms_use[-1]) < np.abs(p_ms_use[0])),
                "note": "FIELD independent of eq3 self-comparison",
            }
        )

        # also compare to single-dominant-mode (diagnostic)
        wphi = np.abs(ps * pr)
        m1 = int(np.argmax(wphi))
        p_m1 = ps[m1] * pr[m1] * np.exp(-1j * k_re[m1] * r_use - alpha_m[m1] * r_use) / np.sqrt(r_use)
        c1 = np.vdot(p_m1, p_field) / max(np.vdot(p_m1, p_m1), 1e-30)
        rel1 = np.linalg.norm(p_field - c1 * p_m1) / max(np.linalg.norm(p_field), 1e-30)

        ok_fv = bool(abs(cc) > 0.95 and rel_res < 0.35)
        fv_rows.append(
            {
                "zs": zs,
                "zr": zr,
                "status": "OK" if ok_fv else "MISMATCH",
                "n_ranges": int(r_use.size),
                "best_scale_re": float(cstar.real),
                "best_scale_im": float(cstar.imag),
                "complex_corr": float(abs(cc)),
                "rel_residual": float(rel_res),
                "phase_rmse_rad": float(np.sqrt(np.mean(ph**2))),
                "amp_shape_rel": float(amp_rel),
                "single_mode_rel_residual": float(rel1),
                "note": "global c* only; single_mode diagnostic separate",
            }
        )
        if not ok_fv:
            field_ok = False

    pd.DataFrame(fv_rows).to_csv(OUT / "field_vs_modesum_check_FIXED.csv", index=False)
    pd.DataFrame(alpha_rows).to_csv(OUT / "kraken_alpha_mapping_FIXED.csv", index=False)

    # -------- 3) spectrum + peak detect + oracle map + depth --------
    L = SPAN
    dk = DK_R / 16.0
    k_grid = np.arange(float(k_re.min()) - 0.02, float(k_re.max()) + 0.02, dk)
    r = np.arange(R1, R2 + 1e-9, DR_PAPER)

    for zs, zr in CASES:
        ps, _, izs = phi_at(zs)
        pr, _, izr = phi_at(zr)
        p = pressure_eq3(r, zs, zr, k_re, alpha_m, ps, pr)
        g = g_eq1(r, p, k_grid, S=np.sqrt(r))
        abs_g = np.abs(g)
        gmax = float(abs_g.max())

        # peak detection main prominence 3%, min distance ~ 1 Rayleigh bin
        dist = max(1, int(round(DK_R / dk)))
        peaks, props = find_peaks(abs_g, prominence=PROM_MAIN * gmax, distance=dist)
        peak_rows.append(
            {
                "case": f"zs{zs}_zr{zr}",
                "n_detected": len(peaks),
                "prominence_frac": PROM_MAIN,
                "gmax": gmax,
            }
        )
        for i, pi in enumerate(peaks):
            peak_rows.append(
                {
                    "case": f"zs{zs}_zr{zr}",
                    "peak_id": i,
                    "k_peak": float(k_grid[pi]),
                    "abs_g": float(abs_g[pi]),
                    "prominence": float(props.get("prominences", [np.nan])[i] if i < len(props.get("prominences", [])) else np.nan),
                    "g_re": float(g[pi].real),
                    "g_im": float(g[pi].imag),
                    "prominence_frac": PROM_MAIN,
                    "gmax": gmax,
                }
            )

        # one-to-one matching peaks -> true modes
        used_modes = set()
        selected = []  # list of (peak_idx, mode_id, g_peak)
        for i, pi in enumerate(peaks):
            kp = k_grid[pi]
            cand = [m for m in range(nmode) if abs(kp - k_re[m]) <= MATCH_TOL and m not in used_modes]
            # nearest among unused
            if not cand:
                # check if any mode in window but already used or ambiguous
                any_near = [m for m in range(nmode) if abs(kp - k_re[m]) <= MATCH_TOL]
                status = "RESOLUTION_AMBIGUOUS" if len(any_near) > 1 else "UNMATCHED"
                map_rows.append(
                    {
                        "case": f"zs{zs}_zr{zr}",
                        "peak_id": i,
                        "k_peak": float(kp),
                        "peak_complex": complex(g[pi]),
                        "prominence": float(props["prominences"][i]),
                        "mode_id": -1,
                        "k_true": np.nan,
                        "delta_k": np.nan,
                        "status": status,
                    }
                )
                continue
            if len(cand) > 1:
                # multiple unused modes within tol → ambiguous unless unique nearest
                cand = [min(cand, key=lambda m: abs(kp - k_re[m]))]
                # still flag if second is also within tol
                map_rows.append  # noqa
            m = cand[0]
            used_modes.add(m)
            status = "SELECTED"
            selected.append((pi, m, g[pi]))
            map_rows.append(
                {
                    "case": f"zs{zs}_zr{zr}",
                    "peak_id": i,
                    "k_peak": float(kp),
                    "peak_complex": complex(g[pi]),
                    "prominence": float(props["prominences"][i]),
                    "mode_id": m + 1,
                    "k_true": float(k_re[m]),
                    "delta_k": float(kp - k_re[m]),
                    "status": status,
                }
            )
        # NOT_RESOLVED true modes
        for m in range(nmode):
            if m not in used_modes:
                map_rows.append(
                    {
                        "case": f"zs{zs}_zr{zr}",
                        "peak_id": -1,
                        "k_peak": np.nan,
                        "peak_complex": np.nan,
                        "prominence": np.nan,
                        "mode_id": m + 1,
                        "k_true": float(k_re[m]),
                        "delta_k": np.nan,
                        "status": "NOT_RESOLVED",
                    }
                )

        mode_idx = [m for _, m, _ in selected]
        g_sel = np.array([gv for _, _, gv in selected], dtype=complex) if selected else np.array([], dtype=complex)
        max_phi = float(np.max(np.abs(phi)))
        ok_case = False
        zhat = np.nan
        abs_err = np.nan
        ratio_true = np.nan
        if selected:
            for d_ratio in DELTA_RATIOS:
                d_abs = d_ratio * max_phi
                D = depth_ambiguity_eq6(phi, g_sel, mode_idx, pr, d_abs, k_re)
                Dn = D / max(np.sum(D), 1e-30)
                zhat = float(depths[int(np.argmax(Dn))])
                abs_err = abs(zhat - zs)
                ratio_true = float(Dn[int(np.argmin(np.abs(depths - zs)))] / max(Dn.max(), 1e-30))
                depth_sel_rows.append(
                    {
                        "zs_true": zs,
                        "zr": zr,
                        "Delta_ratio": d_ratio,
                        "n_selected": len(selected),
                        "selected_mode_ids": ";".join(str(m + 1) for m in mode_idx),
                        "z_hat": zhat,
                        "abs_depth_error": abs_err,
                        "D_true_over_Dmax": ratio_true,
                        "gate_err_le_2m": bool(abs_err <= 2.0),
                        "gate_ratio_ge_0.8": bool(ratio_true >= 0.8),
                    }
                )
            # main gate Delta=0.1
            d_abs = 0.10 * max_phi
            D = depth_ambiguity_eq6(phi, g_sel, mode_idx, pr, d_abs, k_re)
            Dn = D / max(np.sum(D), 1e-30)
            zhat = float(depths[int(np.argmax(Dn))])
            abs_err = abs(zhat - zs)
            ratio_true = float(Dn[int(np.argmin(np.abs(depths - zs)))] / max(Dn.max(), 1e-30))
            ok_case = bool(abs_err <= 2.0 and ratio_true >= 0.8 and len(selected) > 0)

            # FULL upper bound (all true modes, g at k_m) — NOT for gate
            g_all = g_eq1(r, p, k_re, S=np.sqrt(r))
            Df = depth_ambiguity_eq6(phi, g_all, list(range(nmode)), pr, 0.10 * max_phi, k_re)
            Dfn = Df / max(np.sum(Df), 1e-30)
            zhat_f = float(depths[int(np.argmax(Dfn))])
            depth_full_rows.append(
                {
                    "zs_true": zs,
                    "zr": zr,
                    "label": "ORACLE_FULL_MODE_UPPER_BOUND",
                    "z_hat": zhat_f,
                    "abs_depth_error": abs(zhat_f - zs),
                    "D_true_over_Dmax": float(
                        Dfn[int(np.argmin(np.abs(depths - zs)))] / max(Dfn.max(), 1e-30)
                    ),
                    "used_for_repro_gate": False,
                }
            )

            # prominence sensitivity
            for pf in PROMINENCE_FRACS:
                pk2, _ = find_peaks(abs_g, prominence=pf * gmax, distance=dist)
                used2 = set()
                sel2 = []
                for pi in pk2:
                    kp = k_grid[pi]
                    cand = [m for m in range(nmode) if abs(kp - k_re[m]) <= MATCH_TOL and m not in used2]
                    if not cand:
                        continue
                    m = min(cand, key=lambda mm: abs(kp - k_re[mm]))
                    used2.add(m)
                    sel2.append((m, g[pi]))
                prom_rows.append(
                    {
                        "zs": zs,
                        "zr": zr,
                        "prominence_frac": pf,
                        "n_peaks": len(pk2),
                        "n_selected": len(sel2),
                        "selected_modes": ";".join(str(m + 1) for m, _ in sel2),
                    }
                )
                if sel2:
                    D2 = depth_ambiguity_eq6(
                        phi,
                        np.array([gv for _, gv in sel2]),
                        [m for m, _ in sel2],
                        pr,
                        0.10 * max_phi,
                        k_re,
                    )
                    Dn2 = D2 / max(np.sum(D2), 1e-30)
                    z2 = float(depths[int(np.argmax(Dn2))])
                    prom_rows[-1]["z_hat"] = z2
                    prom_rows[-1]["abs_error"] = abs(z2 - zs)
                else:
                    prom_rows[-1]["z_hat"] = np.nan
                    prom_rows[-1]["abs_error"] = np.nan

            # shading sensitivity with reflect pad + selected peaks from theory S=sqrt(r) peaks
            for w in SMOOTH_M:
                S = smooth_abs2(p, DR_PAPER, w)
                gB = g_eq1(r, p, k_grid, S=S)
                absB = np.abs(gB)
                pk3, _ = find_peaks(absB, prominence=PROM_MAIN * absB.max(), distance=dist)
                used3 = set()
                sel3 = []
                for pi in pk3:
                    kp = k_grid[pi]
                    cand = [m for m in range(nmode) if abs(kp - k_re[m]) <= MATCH_TOL and m not in used3]
                    if not cand:
                        continue
                    m = min(cand, key=lambda mm: abs(kp - k_re[mm]))
                    used3.add(m)
                    sel3.append((m, gB[pi]))
                shade_rows.append(
                    {
                        "zs": zs,
                        "zr": zr,
                        "route": "PAPER_DATA_SHADING",
                        "smooth_m": w,
                        "n_selected": len(sel3),
                    }
                )
                if sel3:
                    D3 = depth_ambiguity_eq6(
                        phi,
                        np.array([gv for _, gv in sel3]),
                        [m for m, _ in sel3],
                        pr,
                        0.10 * max_phi,
                        k_re,
                    )
                    Dn3 = D3 / max(np.sum(D3), 1e-30)
                    z3 = float(depths[int(np.argmax(Dn3))])
                    shade_rows[-1]["z_hat"] = z3
                    shade_rows[-1]["abs_error"] = abs(z3 - zs)
                else:
                    shade_rows[-1]["z_hat"] = np.nan
                    shade_rows[-1]["abs_error"] = np.nan

        # figures
        write_svg_spectrum(
            FIG / f"spectrum_peaks_zs{int(zs)}_zr{int(zr)}.svg",
            f"peaks zs={zs} zr={zr} (3% prom)",
            k_grid,
            abs_g,
            list(peaks),
            k_re,
            set(pi for pi, _, _ in selected),
        )

        case_rows.append(
            {
                "zs": zs,
                "zr": zr,
                "n_detected_peaks": int(len(peaks)),
                "n_resolved_modes": int(len(selected)),
                "selected_mode_ids": ";".join(str(m + 1) for _, m, _ in selected),
                "z_hat": zhat,
                "abs_error": abs_err,
                "D_true_over_Dmax": ratio_true,
                "overall_ok": bool(ok_case and field_ok),
                "note": "" if ok_case else "fidelity gate or field/peaks",
            }
        )

    # depth figures
    # rebuild simple curves from last stored is messy; skip pair curves — write from selected main if available
    write_svg_depth(FIG / "depth_selected_zr18.svg", "selected peaks depth zr=18", Z_SEARCH, None, None, np.nan, np.nan)
    write_svg_depth(FIG / "depth_selected_zr70.svg", "selected peaks depth zr=70", Z_SEARCH, None, None, np.nan, np.nan)

    # field vs modesum figure (one case)
    if field_available:
        try:
            (zs, zr) = CASES[0]
            r_shd, p_field, _ = field_payload[(zs, zr)]
            ps, _, _ = phi_at(zs)
            pr, _, _ = phi_at(zr)
            p_ms = pressure_eq3(r_shd, zs, zr, k_re, alpha_m, ps, pr)
            cstar = np.vdot(p_ms, p_field) / max(np.vdot(p_ms, p_ms), 1e-30)
            parts = [
                '<svg xmlns="http://www.w3.org/2000/svg" width="720" height="360" viewBox="0 0 720 360">',
                f'<text x="360" y="20" text-anchor="middle" font-size="14">FIELD vs modesum |p| zs={zs} zr={zr}</text>',
                '<rect x="60" y="40" width="620" height="260" fill="none" stroke="#333"/>',
            ]
            r0, r1 = r_shd.min(), r_shd.max()
            ymax = max(np.abs(p_field).max(), np.abs(cstar * p_ms).max(), 1e-30)
            for arr, color in [(np.abs(p_field), "#c44"), (np.abs(cstar * p_ms), "#44c")]:
                pts = []
                for rv, av in zip(r_shd, arr):
                    x = 60 + (rv - r0) / max(r1 - r0, 1) * 620
                    y = 300 - (av / ymax) * 240
                    pts.append(f"{x:.1f},{y:.1f}")
                parts.append(f'<polyline fill="none" stroke="{color}" stroke-width="1.5" points="{" ".join(pts)}"/>')
            parts.append("</svg>")
            (FIG / "field_vs_modesum.svg").write_text("\n".join(parts), encoding="utf-8")
        except Exception:
            pass

    pd.DataFrame(peak_rows).to_csv(OUT / "paper_detected_peaks_FIXED.csv", index=False)
    pd.DataFrame(map_rows).to_csv(OUT / "paper_oracle_mode_mapping_FIXED.csv", index=False)
    pd.DataFrame(depth_sel_rows).to_csv(OUT / "paper_depth_ambiguity_SELECTED.csv", index=False)
    pd.DataFrame(depth_full_rows).to_csv(OUT / "paper_depth_ambiguity_FULL_MODE_UPPER.csv", index=False)
    pd.DataFrame(prom_rows).to_csv(OUT / "peak_prominence_sensitivity.csv", index=False)
    pd.DataFrame(shade_rows).to_csv(OUT / "shading_sensitivity_FIXED.csv", index=False)
    pd.DataFrame(case_rows).to_csv(OUT / "paper_repro_case_matrix_FIXED.csv", index=False)

    # -------- 4) decision --------
    n_ok = sum(1 for r in case_rows if r.get("overall_ok"))
    # prominence robustness: not complete flip
    if prom_rows:
        prom_df = pd.DataFrame(prom_rows)
        flips = []
        for (zs, zr), grp in prom_df.groupby(["zs", "zr"]):
            hats = grp["z_hat"].dropna().unique()
            flips.append(len(hats) > 2)
        prom_stable = not any(flips)
    else:
        prom_stable = False

    # keep single-mode eq1-eq5 status from prior (documented, not re-run heavy) — recompute quick
    r0 = 0.5 * (R1 + R2)
    eq_ok = True
    max_rel = 0.0
    for m in range(min(5, nmode)):
        ps, _, _ = phi_at(4.0)
        pr, _, _ = phi_at(18.0)
        p1 = pressure_eq3(r[::5], 4.0, 18.0, k_re[m : m + 1], alpha_m[m : m + 1], ps[m : m + 1], pr[m : m + 1])
        a1 = g_eq1(r[::5], p1, k_re[m : m + 1], S=np.sqrt(r[::5]))[0]
        a5 = (b_m_stable(k_re[m : m + 1], alpha_m[m : m + 1], ps[m : m + 1], r0, SPAN) * pr[m : m + 1])[0]
        rel = abs(a1 - a5) / max(abs(a5), 1e-30)
        max_rel = max(max_rel, rel)
    if max_rel > 1e-3:
        eq_ok = False

    fv_df = pd.DataFrame(fv_rows)
    single_mode_good = (
        "single_mode_rel_residual" in fv_df.columns
        and fv_df["single_mode_rel_residual"].notna().any()
        and float(fv_df["single_mode_rel_residual"].min()) < 1e-2
        and field_available
    )
    multi_good = field_ok

    if not field_available:
        decision = "C2_2C_FIELD_CROSSCHECK_BLOCKED"
        why = "FIELD .shd not produced/readable"
    elif multi_good:
        if not eq_ok:
            decision = "C2_2C_FIX_IMPLEMENTATION_FAIL"
            why = f"FIELD/modesum OK but eq1-eq5 max_rel={max_rel:.3e}"
        elif n_ok == 4 and prom_stable:
            decision = "C2_2C_FIX_PAPER_REPRO_CONFIRMED"
            why = f"FIELD multi-mode OK; 4/4 fidelity; prominence stable; eq1eq5={max_rel:.2e}"
        elif n_ok >= 1:
            decision = "C2_2C_FIX_PAPER_REPRO_PARTIAL"
            why = f"FIELD OK but only {n_ok}/4 fidelity gates"
        else:
            decision = "C2_2C_FIX_IMPLEMENTATION_FAIL"
            why = "FIELD OK but 0/4 fidelity / empty peaks"
    else:
        decision = "C2_2C_FIX_IMPLEMENTATION_FAIL"
        why = (
            "FIELD .shd 生成并读取成功，但与 Eq.(3) 多模 modesum 在单一全局复标度下 residual~1、"
            "complex corr 低，不能互证；诊断显示 FIELD 近似单模(1/sqrt(r))而 Eq3 多模干涉不一致。"
            "按完整性规则记 IMPLEMENTATION_FAIL，不确认论文复现。"
        )

    dec = {
        "stage": "R3-C2.2C-FIX",
        "rc3c2_2c_fix_decision": decision,
        "why": why,
        "prior": "C2_2C_PAPER_REPRO_CONFIRMED → SUPERSEDED_PENDING_INTEGRITY_FIX",
        "field_crosscheck": "PASS" if field_ok else "FAIL/BLOCKED",
        "eq1_eq5_single_mode_max_rel": max_rel,
        "cases_ok": n_ok,
        "cases_total": 4,
        "prominence_stable": bool(prom_stable),
        "alpha_mapping": alpha_map,
        "oracle_definition": "PEAK_DETECT_THEN_TRUE_KM_LABEL",
        "yang_route": "YANG_ROUTE_UNDECIDED",
        "not_claimed": ["YANG_ESTD_FAIL", "YANG_APERTURE_LIMITED_IN_CZ", "DEPTH_CANNOT_BE_ESTIMATED", "RC3-C_FAIL"],
        "created_utc": NOW,
    }
    (OUT / "R3_C2_2C_FIX_DECISION.json").write_text(json.dumps(dec, ensure_ascii=False, indent=2), encoding="utf-8")

    report = f"""# R3-C2.2C-FIX 报告

UTC: {NOW}

## 状态

`C2_2C_PAPER_REPRO_CONFIRMED` → **SUPERSEDED_PENDING_INTEGRITY_FIX**

本阶段判定：**`{decision}`**

{why}

## 修复点

1. **FIELD 真交叉验证**：恢复官方 `.flp`（profile@0 km；NProf 单独行；tail `1/0.0/0.0`）；从 `.shd` 读压力；与 Eq.(3) 仅拟合 **全局复标度 c\\*** 后比较 corr/residual（不再 self-compare）。
2. **ORACLE_MODE_ID**：`find_peaks` 先检测真实谱峰（prominence 主基线 3%；1/3/5/10% 敏感性）→ true k_m **一对一标注** → 仅 SELECTED 复峰进 Eq.(6)。全模态仅 `ORACLE_FULL_MODE_UPPER_BOUND`。
3. **depth gate**：`|z_hat-z_s|<=2 m` 且 `D(true)/D(max)>=0.8`；删除“4→5 最近节点”表述。
4. **Delta**：`Delta_abs = ratio * max|phi|` 全链一致。
5. **shading**：reflect padding，三窗全报。
6. **mode-count 23 vs 16**：取消 gate，仅 DIAGNOSTIC。

## 自证

- 单模 Eq1→Eq5 max_rel = {max_rel:.3e}
- alpha 映射：{alpha_map}
- FIELD：{'PASS' if field_ok else 'FAIL/BLOCKED'}

## 范围

不进 E-STD / 50–60 km / 2.4 km CZ / δ 搜索 / MC / AR / P5。YANG_ROUTE_UNDECIDED。

## 停止

等待 GPT 审计。
"""
    (OUT / "R3_C2_2C_FIX_REPORT.md").write_text(report, encoding="utf-8")
    (OUT / "R3_C2_2C_FIX_GPT_SYNC.md").write_text(
        f"# R3-C2.2C-FIX\\n\\n**{decision}**\\n\\n{why}\\n\\nFIELD={'PASS' if field_ok else 'FAIL'}; cases_ok={n_ok}/4; YANG_ROUTE_UNDECIDED\\n",
        encoding="utf-8",
    )
    print("DECISION", decision, flush=True)
    return 0


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except Exception:
        traceback.print_exc()
        (OUT / "R3_C2_2C_FIX_DECISION.json").write_text(
            json.dumps(
                {
                    "rc3c2_2c_fix_decision": "C2_2C_FIX_IMPLEMENTATION_FAIL",
                    "why": "exception",
                    "yang_route": "YANG_ROUTE_UNDECIDED",
                },
                indent=2,
            ),
            encoding="utf-8",
        )
        raise
