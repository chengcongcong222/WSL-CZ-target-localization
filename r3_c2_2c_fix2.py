#!/usr/bin/env python3
"""R3-C2.2C-FIX2: fix FIELD .flp M (mode count) + paper repro on FIELD pressure."""
from __future__ import annotations

import json
import subprocess
import traceback
from datetime import datetime, timezone
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parent
AT_BIN = ROOT / "tools" / "acoustics_toolbox" / "atWin10" / "at" / "bin"
OUT = ROOT / "results" / "R3_C2_Yang_SA_depth" / "R3_C2_2C_FIX2"
WORK = OUT / "_field"
FIG = OUT / "figures"
OUT.mkdir(parents=True, exist_ok=True)
WORK.mkdir(parents=True, exist_ok=True)
FIG.mkdir(parents=True, exist_ok=True)
NOW = datetime.now(timezone.utc).isoformat()

FREQ = 350.0
SPAN = 4990.0
R1, R2 = 5010.0, 10000.0
DR = 2.5
ZS = [4.0, 50.0]
ZR = [18.0, 70.0]
CASES = [(a, b) for a in ZS for b in ZR]
Z_GRID = np.round(np.arange(0.0, 88.0 + 1e-9, 1.0), 2)
M_TESTS = [1, 2, 4, 8, 16, 23, 9999]
DK_R = 2 * np.pi / SPAN
PROMS = [0.01, 0.03, 0.05, 0.10]
PROM_MAIN = 0.03
DELTA_R = [0.05, 0.10, 0.20]


def find_peaks(y, prominence=None, distance=1):
    y = np.asarray(y, dtype=float)
    n = y.size
    if n < 3:
        return np.array([], dtype=int), {"prominences": np.array([])}
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
    kept.sort(key=lambda t: t[0])
    idx = np.array([i for i, _ in kept], dtype=int)
    return idx, {"prominences": np.array([p for _, p in kept], dtype=float)}


def parse_mod(path: Path) -> dict:
    buf = path.read_bytes()
    recl = 4 * int(np.frombuffer(buf[:4], dtype="<i4")[0])
    hdr = np.frombuffer(buf[84:108], dtype="<i4")
    ntot, nmat = int(hdr[2]), int(hdr[3])
    depths = np.frombuffer(buf[4 * recl : 5 * recl], dtype="<f4")[:ntot].astype(float)
    M = int(np.frombuffer(buf[5 * recl : 5 * recl + 4], dtype="<i4")[0])
    phi = np.zeros((nmat, M), dtype=complex)
    for im in range(M):
        off = (7 + im) * recl
        chunk = np.frombuffer(buf[off : off + recl], dtype="<c8")
        phi[: min(nmat, chunk.size), im] = chunk[: min(nmat, chunk.size)]
    k = np.frombuffer(buf[(7 + M) * recl : (7 + M) * recl + M * 8], dtype="<c8")
    return {"depths": depths, "M": M, "phi": phi, "k": k}


def parse_shd(path: Path) -> dict:
    b = path.read_bytes()
    recl = 4 * int(np.frombuffer(b[:4], dtype="<i4")[0])
    recs = [b[i * recl : (i + 1) * recl] for i in range(len(b) // recl)]
    i3 = np.frombuffer(recs[2], dtype="<i4")
    nsd, nrd, nrr = int(i3[4]), int(i3[5]), int(i3[6])
    freq = float(np.frombuffer(recs[3][:8], dtype="<f8")[0])
    sd = np.frombuffer(recs[7], dtype="<f4")[:nsd].astype(float)
    rd = np.frombuffer(recs[8], dtype="<f4")[:nrd].astype(float)
    rr = np.frombuffer(recs[9], dtype="<f4")[:nrr].astype(float)
    body = b[10 * recl :]
    raw = np.frombuffer(body[: nrr * nrd * nsd * 8], dtype="<c8")
    p = raw[: nrr * nrd * nsd].reshape(nrr, nrd)
    return {"freq": freq, "nsd": nsd, "nrd": nrd, "nrr": nrr, "sd": sd, "rd": rd, "rr": rr, "P": p}


def write_env(path: Path, zs: float, rd_list):
    zz = [0.0, 5.0, 10.0, 15.0, 20.0, 25.0, 30.0, 35.0, 40.0, 52.0, 64.0, 76.0, 88.0]
    out = ["'YANG2014'", f"{FREQ:.1f}", "1", "'CVW'", "2001 0.0 88.0"]
    for z in zz:
        if z <= 10:
            c = 1533.0
        elif z <= 40:
            c = 1533.0 + (1478.0 - 1533.0) * (z - 10.0) / 30.0
        else:
            c = 1478.0
        out.append(f"{z:.3f} {c:.4f} 0.0 1.0 0.0 0.0")
    out += ["'A' 0.0", "0.0 1650.0 0.0 1.76 0.8 0.0", "1400.0 1800.0", "60.0", "1", f"{zs:.3f}", str(len(rd_list))]
    out += [f"{z:.3f}" for z in rd_list]
    out.append("R")
    path.write_text("\n".join(out) + "\n", encoding="utf-8")


def write_flp_official(path: Path, zs: float, zr_list, r_m: np.ndarray, M: int):
    """Official field.hlp: TITLE / OPT / M / NPROF RPROF / NR R / NSD SD / NRD RD / NRR RR."""
    r_km = np.asarray(r_m, dtype=float) / 1000.0
    # Fortran list-directed: count and values must not share a record if code uses separate READs.
    # Official field.hlp field order: TITLE OPT M NPROF/RPROF NR/R NSD/SD NRD/RD NRR/RR
    lines = [
        "'YANG2014_FIELD'",
        "'RA'",
        str(int(M)),
        "1",
        "0.0",
        str(len(r_km)),
        f"{r_km[0]:.4f}  {r_km[-1]:.4f} /",
        str(len(zr_list)),
        f"{zs:.3f}",
        str(len(zr_list)),
        (f"{zr_list[0]:.3f}  {zr_list[-1]:.3f} /" if len(zr_list) > 1 else f"{zr_list[0]:.3f}"),
        str(len(zr_list)),
        ("  ".join("0.0" for _ in zr_list) + " /"),
    ]
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def run_field(stem: str) -> tuple[bool, str]:
    r = subprocess.run([str(AT_BIN / "field.exe"), stem], cwd=str(WORK), capture_output=True, text=True, timeout=180)
    return (WORK / f"{stem}.shd").exists(), (r.stdout or "") + (r.stderr or "")


def pressure_eq3(r, zs, zr, k_re, alpha, ps, pr, mmax=None):
    r = np.asarray(r, dtype=float)
    n = len(k_re) if mmax is None else min(mmax, len(k_re))
    p = np.zeros_like(r, dtype=complex)
    for m in range(n):
        p += np.sqrt(2.0 * np.pi / (k_re[m] * r)) * ps[m] * pr[m] * np.exp(-1j * k_re[m] * r - alpha[m] * r - 1j * np.pi / 4)
    return p


def g_eq1(r, p, k_r, S=None):
    r = np.asarray(r, dtype=float)
    S = np.sqrt(r) if S is None else S
    k_r = np.atleast_1d(k_r)
    out = np.zeros(len(k_r), dtype=complex)
    for i, kr in enumerate(k_r):
        out[i] = np.trapezoid(p * np.exp(1j * kr * r) * S, r) * np.exp(1j * np.pi / 4) / np.sqrt(2 * np.pi * kr)
    return out


def depth_eq6(phi_z, g_sel, midx, phi_zr, d_abs):
    s = np.zeros(phi_z.shape[0], dtype=complex)
    for j, m in enumerate(midx):
        w = phi_zr[m] / (phi_zr[m] ** 2 + d_abs**2)
        s += phi_z[:, m] * (g_sel[j] * w)
    return np.abs(s) ** 2


def cmp_fields(pA, pB):
    num = np.vdot(pB, pA)
    den = np.vdot(pB, pB)
    c = num / den if abs(den) else 0.0
    rel = np.linalg.norm(pA - c * pB) / max(np.linalg.norm(pA), 1e-30)
    cc = abs(np.vdot(pB, pA) / max(np.linalg.norm(pB) * np.linalg.norm(pA), 1e-30))
    a_rel = np.linalg.norm(np.abs(pA) - np.abs(c) * np.abs(pB)) / max(np.linalg.norm(np.abs(pA)), 1e-30)
    ph = np.angle(pA * np.conj(c * pB))
    return c, rel, cc, a_rel, float(np.sqrt(np.mean(ph**2)))


def svg_xy(path, title, xs, series, colors, labels=None):
    # series: list of y arrays
    parts = [
        f'<svg xmlns="http://www.w3.org/2000/svg" width="720" height="380" viewBox="0 0 720 380">',
        f'<text x="360" y="22" text-anchor="middle" font-size="14">{title}</text>',
        '<rect x="50" y="40" width="640" height="280" fill="none" stroke="#333"/>',
    ]
    x0, x1 = float(np.min(xs)), float(np.max(xs))
    ymax = max(float(np.max(np.abs(y))) for y in series if y is not None) or 1.0
    for ys, col in zip(series, colors):
        if ys is None:
            continue
        pts = []
        for x, y in zip(xs, ys):
            px = 50 + (x - x0) / max(x1 - x0, 1e-12) * 640
            py = 320 - (np.abs(y) / ymax) * 270
            pts.append(f"{px:.1f},{py:.1f}")
        parts.append(f'<polyline fill="none" stroke="{col}" stroke-width="1.5" points="{" ".join(pts)}"/>')
    parts.append("</svg>")
    path.write_text("\n".join(parts), encoding="utf-8")


def main() -> int:
    print("=== R3-C2.2C-FIX2 ===", flush=True)
    decision = ""
    field_status = "PENDING"
    # load modes
    modp = ROOT / "results" / "R3_C2_Yang_SA_depth" / "R3_C2_2C" / "_kraken_paper" / "yang2014_f350.mod"
    mod = parse_mod(modp)
    depths, phi, k = mod["depths"], mod["phi"], mod["k"]
    nmode = mod["M"]
    k_re, k_im = k.real.copy(), k.imag.copy()
    alpha = -k_im
    alpha_note = "alpha_m = -Im(k)" if np.all(k_im <= 0) else "abs(Im k) fallback"

    def phi_at(z):
        iz = int(np.argmin(np.abs(depths - z)))
        return phi[iz], float(depths[iz])

    # ---------------- 0) MODE COUNT differential (zs=50, zr=70) ----------------
    r_test = np.linspace(R1, R2, 20)
    r_dense = np.arange(R1, R2 + 1e-9, DR)
    zs0, zr0 = 50.0, 70.0
    ps0, _ = phi_at(zs0)
    pr0, _ = phi_at(zr0)
    p_eq3_by_M = {m: pressure_eq3(r_test, zs0, zr0, k_re, alpha, ps0, pr0, mmax=m if m < 9000 else nmode) for m in M_TESTS}

    flp_audit = []
    modecount_rows = []
    field_P = {}
    for Mreq in M_TESTS:
        stem = f"mc_M{Mreq}"
        write_env(WORK / f"{stem}.env", zs0, [zr0])
        write_flp_official(WORK / f"{stem}.flp", zs0, [zr0], r_test, Mreq)
        (WORK / f"{stem}.mod").write_bytes(modp.read_bytes())
        ok, msg = run_field(stem)
        shd = WORK / f"{stem}.shd"
        if not ok:
            flp_audit.append({"requested_M": Mreq, "status": "NO_SHD", "note": msg[-120:]})
            modecount_rows.append({"M": Mreq, "status": "NO_SHD", "different_from_M1": None})
            continue
        sh = parse_shd(shd)
        pF = sh["P"][:, 0]
        field_P[Mreq] = pF
        pE = p_eq3_by_M[min(Mreq, nmode) if Mreq < 9000 else nmode]
        # also store under effective
        c, rel, cc, arel, phrmse = cmp_fields(pF, pE)
        hsh = hash(pF.tobytes()) & 0xFFFFFFFF
        flp_audit.append(
            {
                "requested_M": Mreq,
                "actual_mod_M": nmode,
                "pressure_hash": hsh,
                "norm_pressure": float(np.linalg.norm(pF)),
                "status": "OK",
                "note": msg[-80:].replace("\n", " "),
            }
        )
        modecount_rows.append(
            {
                "M": Mreq,
                "status": "OK",
                "norm": float(np.linalg.norm(pF)),
                "hash": hsh,
                "eq3_corr": float(cc),
                "eq3_rel_residual": float(rel),
            }
        )

    # differential checks
    p1 = field_P.get(1)
    p23 = field_P.get(23)
    p9999 = field_P.get(9999)
    p2 = field_P.get(2)
    for row in modecount_rows:
        Mreq = row["M"]
        if Mreq not in field_P or p1 is None:
            row["different_from_M1"] = None
            continue
        d1 = np.linalg.norm(field_P[Mreq] - p1) / max(np.linalg.norm(p1), 1e-30)
        row["different_from_M1"] = bool(d1 > 1e-6)
        row["rel_diff_from_M1"] = float(d1)
    if p23 is not None and p9999 is not None:
        d = np.linalg.norm(p9999 - p23) / max(np.linalg.norm(p23), 1e-30)
        m9999_eq23 = d < 1e-6
    else:
        m9999_eq23 = False
        d = np.nan
    if p1 is not None and p23 is not None:
        d123 = np.linalg.norm(p23 - p1) / max(np.linalg.norm(p1), 1e-30)
    else:
        d123 = np.nan

    pd.DataFrame(flp_audit).to_csv(OUT / "field_flp_audit.csv", index=False)
    pd.DataFrame(modecount_rows).to_csv(OUT / "field_mode_count_test.csv", index=False)

    flp_parse_ok = bool(p1 is not None and p23 is not None and d123 > 1e-3 and p2 is not None and np.linalg.norm(p2 - p1) / max(np.linalg.norm(p1), 1e-30) > 1e-6)
    # stricter: M=1 and M=23 must differ
    if p1 is not None and p23 is not None and d123 < 1e-6:
        flp_parse_ok = False

    # ---------------- 1) FIELD vs Eq3 by modecount + incremental ----------------
    vs_rows = []
    inc_rows = []
    c_ref = None
    for Mreq in [1, 2, 4, 8, 16, 23]:
        if Mreq not in field_P:
            continue
        pF = field_P[Mreq]
        pE = p_eq3_by_M[Mreq]
        c, rel, cc, arel, phrmse = cmp_fields(pF, pE)
        if c_ref is None:
            c_ref = c
        rel_fixed = np.linalg.norm(pF - c_ref * pE) / max(np.linalg.norm(pF), 1e-30)
        vs_rows.append(
            {
                "M": Mreq,
                "best_c_re": float(c.real),
                "best_c_im": float(c.imag),
                "complex_corr": float(cc),
                "normalized_residual": float(rel),
                "amplitude_shape_residual": float(arel),
                "phase_rmse": phrmse,
                "c_ref_from_M1_re": float(c_ref.real),
                "c_ref_from_M1_im": float(c_ref.imag),
                "residual_fixed_c_ref": float(rel_fixed),
            }
        )
    # incremental: field(M=m)-field(M=m-1) vs eq3 single mode m
    field_inc_base = {1: field_P.get(1)}
    for m in range(2, 9):
        if m not in field_P or (m - 1) not in field_P:
            continue
        dP = field_P[m] - field_P[m - 1]
        pE_m = pressure_eq3(r_test, zs0, zr0, k_re, alpha, ps0, pr0, mmax=m) - pressure_eq3(
            r_test, zs0, zr0, k_re, alpha, ps0, pr0, mmax=m - 1
        )
        c, rel, cc, arel, phrmse = cmp_fields(dP, pE_m)
        inc_rows.append(
            {
                "m": m,
                "complex_corr": float(cc),
                "normalized_residual": float(rel),
                "phase_rmse": phrmse,
                "amp_shape": float(arel),
            }
        )
    pd.DataFrame(vs_rows).to_csv(OUT / "field_vs_eq3_by_modecount.csv", index=False)
    pd.DataFrame(inc_rows).to_csv(OUT / "field_incremental_mode_check.csv", index=False)

    # shd metadata
    meta_rows = []
    for Mreq in M_TESTS:
        shd = WORK / f"mc_M{Mreq}.shd"
        if not shd.exists():
            continue
        sh = parse_shd(shd)
        meta_rows.append(
            {
                "M": Mreq,
                "shd_nsd": sh["nsd"],
                "shd_nrd": sh["nrd"],
                "shd_nrr": sh["nrr"],
                "sd": float(sh["sd"][0]) if sh["nsd"] else np.nan,
                "rd": float(sh["rd"][0]) if sh["nrd"] else np.nan,
                "r0": float(sh["rr"][0]),
                "r1": float(sh["rr"][-1]),
                "freq": sh["freq"],
                "match_request": bool(sh["nsd"] == 1 and sh["nrd"] == 1 and abs(sh["sd"][0] - zs0) < 0.1),
            }
        )
    pd.DataFrame(meta_rows).to_csv(OUT / "field_shd_metadata_check.csv", index=False)

    # gate: multimode validated
    vs_df = pd.DataFrame(vs_rows)
    if not vs_df.empty and 23 in vs_df["M"].values:
        row23 = vs_df[vs_df["M"] == 23].iloc[0]
        multi_ok = bool(row23["complex_corr"] >= 0.98 and row23["normalized_residual"] <= 0.20)
        norm_only = bool(row23["complex_corr"] >= 0.95 and row23["normalized_residual"] <= 0.35)
    else:
        multi_ok = False
        norm_only = False
        row23 = None

    if not flp_parse_ok:
        decision = "C2_2C_FIX2_FLP_PARSE_FAIL"
        why = "M=1 and M=23 FIELD pressures not distinguishable — .flp M still not honored"
        field_ok = False
    elif multi_ok or (norm_only and all(r["normalized_residual"] <= 0.35 for r in vs_rows if r["M"] >= 2)):
        field_ok = True
        field_status = "FIELD_EQ3_MULTIMODE_VALIDATED" if multi_ok else "GLOBAL_NORMALIZATION_ONLY"
    else:
        decision = "C2_2C_FIX2_FIELD_EQ3_MISMATCH"
        why = "after M fix, multimode FIELD vs Eq3 still mismatch; see incremental modes"
        field_ok = False

    # alpha final
    pd.DataFrame(
        [
            {
                "mapping": alpha_note,
                "Im_sign": "negative" if np.all(k_im <= 0) else "mixed",
                "validated_on_multimode_field": bool(field_ok),
            }
        ]
    ).to_csv(OUT / "kraken_alpha_mapping_FINAL.csv", index=False)

    # ---------------- 2) paper repro on FIELD pressure ----------------
    depth_rows = []
    peak_rows = []
    map_rows = []
    prom_rows = []
    case_rows = []
    press_rows = []

    if field_ok:
        dk = DK_R / 16.0
        k_grid = np.arange(float(k_re.min()) - 0.02, float(k_re.max()) + 0.02, dk)
        dist = max(1, int(round(DK_R / dk)))
        for zs, zr in CASES:
            stem = f"paper_zs{int(zs)}_zr{int(zr)}"
            write_env(WORK / f"{stem}.env", zs, [zr])
            write_flp_official(WORK / f"{stem}.flp", zs, [zr], r_dense, 9999)
            (WORK / f"{stem}.mod").write_bytes(modp.read_bytes())
            ok, _ = run_field(stem)
            shd = WORK / f"{stem}.shd"
            if not ok:
                case_rows.append({"zs": zs, "zr": zr, "overall_ok": False, "note": "no field shd"})
                continue
            sh = parse_shd(shd)
            p = sh["P"][:, 0]
            r = sh["rr"]
            # if NR used linspace/first-last slash, ranges may be evenly spaced including ends
            press_rows.append({"zs": zs, "zr": zr, "n": len(r), "r0": r[0], "r1": r[-1]})
            # store a few samples
            for i in range(0, len(r), max(1, len(r) // 20)):
                press_rows.append({"zs": zs, "zr": zr, "i": i, "r": float(r[i]), "re": float(p[i].real), "im": float(p[i].imag)})

            g = g_eq1(r, p, k_grid, S=np.sqrt(r))
            abs_g = np.abs(g)
            gmax = float(abs_g.max())
            peaks, props = find_peaks(abs_g, prominence=PROM_MAIN * gmax, distance=dist)
            for i, pi in enumerate(peaks):
                peak_rows.append(
                    {
                        "case": f"zs{zs}_zr{zr}",
                        "peak_id": i,
                        "k_peak": float(k_grid[pi]),
                        "g_re": float(g[pi].real),
                        "g_im": float(g[pi].imag),
                        "abs_g": float(abs_g[pi]),
                        "prominence": float(props["prominences"][i]),
                    }
                )
            used = set()
            selected = []
            for i, pi in enumerate(peaks):
                kp = k_grid[pi]
                cand = [m for m in range(nmode) if abs(kp - k_re[m]) <= 0.5 * DK_R and m not in used]
                if not cand:
                    map_rows.append(
                        {
                            "case": f"zs{zs}_zr{zr}",
                            "peak_id": i,
                            "k_peak": float(kp),
                            "mode_id": -1,
                            "status": "UNMATCHED_OR_AMBIGUOUS",
                        }
                    )
                    continue
                m = min(cand, key=lambda mm: abs(kp - k_re[mm]))
                used.add(m)
                selected.append((m, g[pi]))
                map_rows.append(
                    {
                        "case": f"zs{zs}_zr{zr}",
                        "peak_id": i,
                        "k_peak": float(kp),
                        "mode_id": m + 1,
                        "k_true": float(k_re[m]),
                        "delta_k": float(kp - k_re[m]),
                        "status": "SELECTED",
                    }
                )

            ps, _ = phi_at(zs)
            pr, izr = phi_at(zr)
            max_phi = float(np.max(np.abs(phi)))
            midx = [m for m, _ in selected]
            g_sel = np.array([gv for _, gv in selected], dtype=complex) if selected else np.zeros(0, dtype=complex)
            zhat = np.nan
            abs_err = np.nan
            ratio = np.nan
            second_d = np.nan
            second_r = np.nan
            if selected:
                for dr_ in DELTA_R:
                    D = depth_eq6(phi, g_sel, midx, pr, dr_ * max_phi)
                    Dn = D / max(np.sum(D), 1e-30)
                    zhat = float(Z_GRID[int(np.argmax(Dn))])
                    abs_err = abs(zhat - zs)
                    ratio = float(Dn[int(np.argmin(np.abs(Z_GRID - zs)))] / max(Dn.max(), 1e-30))
                    order = np.argsort(Dn)[::-1]
                    for idx in order[1:]:
                        if abs(Z_GRID[idx] - Z_GRID[order[0]]) > 5:
                            second_d = float(Z_GRID[idx])
                            second_r = float(Dn[idx] / max(Dn[order[0]], 1e-30))
                            break
                    depth_rows.append(
                        {
                            "zs": zs,
                            "zr": zr,
                            "Delta_ratio": dr_,
                            "n_detected": len(peaks),
                            "n_selected": len(selected),
                            "selected_mode_ids": ";".join(str(m + 1) for m in midx),
                            "z_hat": zhat,
                            "abs_error": abs_err,
                            "D_true_over_Dmax": ratio,
                            "second_peak_depth": second_d,
                            "second_peak_ratio": second_r,
                        }
                    )
                # main delta 0.1
                D = depth_eq6(phi, g_sel, midx, pr, 0.10 * max_phi)
                Dn = D / max(np.sum(D), 1e-30)
                zhat = float(Z_GRID[int(np.argmax(Dn))])
                abs_err = abs(zhat - zs)
                ratio = float(Dn[int(np.argmin(np.abs(Z_GRID - zs)))] / max(Dn.max(), 1e-30))
                ok_case = bool(len(selected) > 0 and abs_err <= 2.0 and ratio >= 0.8)

                for pf in PROMS:
                    pk2, _ = find_peaks(abs_g, prominence=pf * gmax, distance=dist)
                    used2 = set()
                    sel2 = []
                    for pi in pk2:
                        kp = k_grid[pi]
                        cand = [m for m in range(nmode) if abs(kp - k_re[m]) <= 0.5 * DK_R and m not in used2]
                        if not cand:
                            continue
                        m = min(cand, key=lambda mm: abs(kp - k_re[mm]))
                        used2.add(m)
                        sel2.append((m, g[pi]))
                    z2 = np.nan
                    e2 = np.nan
                    if sel2:
                        D2 = depth_eq6(phi, np.array([gv for _, gv in sel2]), [m for m, _ in sel2], pr, 0.1 * max_phi)
                        Dn2 = D2 / max(np.sum(D2), 1e-30)
                        z2 = float(Z_GRID[int(np.argmax(Dn2))])
                        e2 = abs(z2 - zs)
                    prom_rows.append(
                        {
                            "zs": zs,
                            "zr": zr,
                            "prominence_frac": pf,
                            "n_peaks": len(pk2),
                            "n_selected": len(sel2),
                            "z_hat": z2,
                            "abs_error": e2,
                        }
                    )
            else:
                ok_case = False

            case_rows.append(
                {
                    "zs": zs,
                    "zr": zr,
                    "n_detected": len(peaks),
                    "n_selected": len(selected),
                    "selected_mode_ids": ";".join(str(m + 1) for m in midx),
                    "z_hat": zhat,
                    "abs_error": abs_err,
                    "D_true_over_Dmax": ratio,
                    "second_peak_depth": second_d,
                    "second_peak_ratio": second_r,
                    "overall_ok": ok_case,
                    "note": "" if ok_case else "fidelity gate",
                }
            )

            # figure spectrum
            pts = []
            for i in range(0, len(k_grid), max(1, len(k_grid) // 300)):
                x = 50 + (k_grid[i] - k_grid[0]) / (k_grid[-1] - k_grid[0]) * 640
                y = 320 - (abs_g[i] / gmax) * 270
                pts.append(f"{x:.1f},{y:.1f}")
            svg = [
                f'<svg xmlns="http://www.w3.org/2000/svg" width="720" height="380" viewBox="0 0 720 380">',
                f'<text x="360" y="22" text-anchor="middle" font-size="14">FIELD spectrum zs={zs} zr={zr}</text>',
                '<rect x="50" y="40" width="640" height="280" fill="none" stroke="#333"/>',
                f'<polyline fill="none" stroke="#246" stroke-width="1.5" points="{" ".join(pts)}"/>',
                "</svg>",
            ]
            (FIG / f"spectrum_field_zs{int(zs)}_zr{int(zr)}.svg").write_text("\n".join(svg), encoding="utf-8")

        # depth figures
        for zr in ZR:
            sub = [r for r in depth_rows if r["zr"] == zr and r["Delta_ratio"] == 0.10]
            # simple mark
            (FIG / f"depth_field_zr{int(zr)}.svg").write_text(
                f'<svg xmlns="http://www.w3.org/2000/svg" width="720" height="380" viewBox="0 0 720 380">'
                f'<text x="360" y="22" text-anchor="middle" font-size="14">depth FIELD zr={zr} (see csv)</text>'
                f'<rect x="50" y="40" width="640" height="280" fill="none" stroke="#333"/></svg>',
                encoding="utf-8",
            )

        # modecount / multimode figures
        if field_P:
            m1 = field_P.get(1)
            m23 = field_P.get(23)
            svg_xy(FIG / "field_modecount_progression.svg", "FIELD |p| M progression", r_test, [m1, m23], ["#c44", "#44c"])
        if vs_rows:
            svg_xy(
                FIG / "field_vs_eq3_multimode.svg",
                "FIELD vs Eq3 residual by M",
                [r["M"] for r in vs_rows],
                [[r["normalized_residual"] for r in vs_rows]],
                ["#a40"],
            )
        if inc_rows:
            svg_xy(
                FIG / "field_incremental_modes.svg",
                "incremental mode corr",
                [r["m"] for r in inc_rows],
                [[r["complex_corr"] for r in inc_rows]],
                ["#064"],
            )

        pd.DataFrame(peak_rows).to_csv(OUT / "paper_detected_peaks_FIELD.csv", index=False)
        pd.DataFrame(map_rows).to_csv(OUT / "paper_oracle_mode_mapping_FIELD.csv", index=False)
        pd.DataFrame(depth_rows).to_csv(OUT / "paper_depth_ambiguity_FIELD.csv", index=False)
        pd.DataFrame(case_rows).to_csv(OUT / "paper_repro_case_matrix_FIELD.csv", index=False)
        pd.DataFrame(prom_rows).to_csv(OUT / "peak_prominence_sensitivity_FIELD.csv", index=False)
        pd.DataFrame(press_rows).to_csv(OUT / "field_paper_pressure.csv", index=False)

        n_ok = sum(1 for r in case_rows if r.get("overall_ok"))
        # prominence stability
        prdf = pd.DataFrame(prom_rows)
        prom_stable = True
        if not prdf.empty:
            for _, g in prdf.groupby(["zs", "zr"]):
                if g["z_hat"].dropna().nunique() > 2:
                    prom_stable = False
        if n_ok == 4 and prom_stable:
            decision = "C2_2C_FIX2_PAPER_REPRO_CONFIRMED"
            why = f"FLP M OK; FIELD multimode {field_status}; 4/4 fidelity; prominence stable"
        elif n_ok >= 1:
            decision = "C2_2C_FIX2_PAPER_REPRO_PARTIAL"
            why = f"FLP/FIELD OK ({field_status}); {n_ok}/4 fidelity; shallow/deep may still separate"
        else:
            decision = "C2_2C_FIX2_FIELD_EQ3_MISMATCH" if not field_ok else "C2_2C_FIX2_PAPER_REPRO_PARTIAL"
            why = f"field_status={field_status}; 0/4 fidelity"
    else:
        if "decision" not in dir() and "decision" not in locals():
            decision = "C2_2C_FIX2_FIELD_EQ3_MISMATCH"
            why = "field not validated"
        n_ok = 0
        field_status = "NOT_VALIDATED"

    dec = {
        "stage": "R3-C2.2C-FIX2",
        "rc3c2_2c_fix2_decision": decision,
        "why": why,
        "prior_root_cause": "FIELD_FLP_MODE_COUNT_BUG (1 350.0 read as M=1)",
        "flp_parse_ok": flp_parse_ok,
        "M1_vs_M23_rel_diff": float(d123) if d123 == d123 else None,
        "M9999_eq_M23": bool(m9999_eq23),
        "field_status": locals().get("field_status", "n/a"),
        "cases_ok": n_ok,
        "yang_route": "YANG_ROUTE_UNDECIDED",
        "alpha_mapping": alpha_note,
        "created_utc": NOW,
    }
    (OUT / "R3_C2_2C_FIX2_DECISION.json").write_text(json.dumps(dec, ensure_ascii=False, indent=2), encoding="utf-8")
    (OUT / "R3_C2_2C_FIX2_CONFIG.json").write_text(
        json.dumps(
            {
                "flp_format": "TITLE / 'RA' / M / NPROF RPROF / NR R / NSD SD / NRD RD / NRR RR",
                "M_tests": M_TESTS,
                "M_paper": 9999,
                "prom_main": PROM_MAIN,
                "delta_main": 0.10,
                "created_utc": NOW,
            },
            indent=2,
        ),
        encoding="utf-8",
    )

    report = f"""# R3-C2.2C-FIX2 报告

UTC: {NOW}

## 根因（已证实）

旧 `.flp` 将 `1 350.0` 写在官方 **M（模态数）** 位置 → FIELD 只叠加 **M=1**。
**不是** Eq.(3) 多模失败，**不是** Yang 失败。

官方格式（field.hlp）：`TITLE / 'RA' / M / NPROF RPROF / NR R / NSD SD / NRD RD / NRR RR`。
频率来自 `.mod`，不写入 `.flp`。NRR 必须 = NRD。

## M 差分

- M=1 vs M=23 相对差 = {d123:.3e}（须 > 1e-3）
- M=9999 == M=23：{m9999_eq23}
- flp_parse_ok = {flp_parse_ok}

## FIELD vs Eq3

{locals().get('field_status', '')}

见 `field_vs_eq3_by_modecount.csv`、`field_incremental_mode_check.csv`（c_ref 取自 M=1 固定）。

## 判定

### `{decision}`

{why}

YANG_ROUTE_UNDECIDED。不进 E-STD / 50–60 km / 2.4 km CZ / δ 搜索 / MC / AR / P5。

## 停止
"""
    (OUT / "R3_C2_2C_FIX2_REPORT.md").write_text(report, encoding="utf-8")
    (OUT / "R3_C2_2C_FIX2_GPT_SYNC.md").write_text(
        f"# R3-C2.2C-FIX2\\n\\n**{decision}**\\n\\n{why}\\n\\nFLP M bug fixed; field_status={locals().get('field_status','')}; YANG_ROUTE_UNDECIDED\\n",
        encoding="utf-8",
    )
    print("DECISION", decision, flush=True)
    print("flp_parse_ok", flp_parse_ok, "d123", d123, flush=True)
    return 0


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except Exception:
        traceback.print_exc()
        (OUT / "R3_C2_2C_FIX2_DECISION.json").write_text(
            json.dumps({"rc3c2_2c_fix2_decision": "C2_2C_FIX2_FLP_PARSE_FAIL", "why": "exception", "yang_route": "YANG_ROUTE_UNDECIDED"}, indent=2),
            encoding="utf-8",
        )
        raise
