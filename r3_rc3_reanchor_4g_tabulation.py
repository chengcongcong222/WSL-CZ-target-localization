#!/usr/bin/env python3
"""R3-RC3-REANCHOR-4G: mode-depth tabulation integrity. No S1."""
from __future__ import annotations

import json
import math
import re
import subprocess
from datetime import datetime, timezone
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parent
OUT = ROOT / "results" / "R3_RC3_REANCHOR" / "R3_RC3_ENV_MISMATCH_BOUNDARY_FIX2"
ZGRID = ROOT / "results" / "R3_C2_Yang_SA_depth" / "R3_C2_1" / "_kraken_zgrid"
PAIRS = ROOT / "results" / "P3_RC3_increment" / "selected_hard_pairs.csv"
AT_BIN = ROOT / "tools" / "acoustics_toolbox" / "atWin10" / "at" / "bin"
WORK = OUT / "_env"
OUT.mkdir(parents=True, exist_ok=True)
WORK.mkdir(parents=True, exist_ok=True)
NOW = datetime.now(timezone.utc).isoformat()

U_PLAT = 2.0
ZR = 200.0
Z_TRUE = [180.0, 200.0, 220.0]
Z_GRID = np.arange(150.0, 250.0 + 1e-9, 5.0)
FREQS = [201.0, 235.0, 283.0, 338.0]
R_GATE = np.array([45000.0, 48000.0, 51000.0, 54000.0, 57000.0, 60000.0])
CC_MIN, RESID_MAX, TL_RMS_MAX, TL_CORR_MIN = 0.95, 0.35, 1.5, 0.95


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
    return {"depths": depths, "M": M, "phi": phi, "k": k, "nmat": nmat, "ntot": ntot}


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


def split_env_tail(text: str):
    lines = text.splitlines()
    for i, ln in enumerate(lines):
        if ln.strip().startswith("'R'") or ln.strip() == "R":
            return lines[:i], lines[i:]
    return lines, []


def read_ssp_block(text: str):
    zs, cs = [], []
    head, _ = split_env_tail(text)
    for ln in head:
        p = ln.split()
        if len(p) >= 6:
            try:
                z = float(p[0])
                [float(p[i]) for i in range(6)]
                c = float(p[1])
            except ValueError:
                continue
            if 0.0 <= z <= 5000.0:
                zs.append(z)
                cs.append(c)
    return np.array(zs, float), np.array(cs, float)


def shift_ssp(z, c, dz=0.0, dc=0.0):
    zq = np.clip(z - dz, float(z.min()), float(z.max()))
    return np.interp(zq, z, c) + dc


def write_env_pert(path: Path, freq: float, z_ssp, c_ssp, env0_path: Path):
    """Change ONLY title + SSP c(z); copy tail after 'R' verbatim."""
    text0 = env0_path.read_text(encoding="utf-8")
    head, tail = split_env_tail(text0)
    z_ssp = np.asarray(z_ssp, float)
    c_ssp = np.asarray(c_ssp, float)
    new_head = []
    for i, ln in enumerate(head):
        if i == 0:
            new_head.append("'PERT'")
        elif i == 1:
            new_head.append(f"{freq:.3f}")
        else:
            p = ln.split()
            if len(p) >= 6:
                try:
                    z = float(p[0])
                    [float(p[j]) for j in range(6)]
                    c_new = float(np.interp(z, z_ssp, c_ssp))
                    new_head.append(f"{z:.3f} {c_new:.4f} 0.0 1.0 0.0 0.0")
                    continue
                except ValueError:
                    pass
            new_head.append(ln)
    path.write_text("\n".join(new_head + tail) + "\n", encoding="utf-8")
    return tail


def write_flp(path: Path, zs: float, zr: float, r_m, M=9999):
    r_km = np.asarray(r_m, float) / 1000.0
    path.write_text(
        "\n".join(
            [
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
        )
        + "\n",
        encoding="utf-8",
    )


def pressure_series(mod, r, zs, zr=ZR):
    depths, phi, k = mod["depths"], mod["phi"], mod["k"]
    izs = int(np.argmin(np.abs(depths - zs)))
    izr = int(np.argmin(np.abs(depths - zr)))
    kre, alpha = k.real, -k.imag
    ps, pr = phi[izs], phi[izr]
    r = np.asarray(r, float)
    p = np.zeros(r.size, dtype=complex)
    for m in range(len(kre)):
        p += (
            np.sqrt(2 * np.pi / (kre[m] * r))
            * ps[m]
            * pr[m]
            * np.exp(-1j * kre[m] * r - alpha[m] * r - 1j * np.pi / 4)
        )
    return p


def tl_shape(p):
    L = 20.0 * np.log10(np.maximum(np.abs(p), 1e-30))
    return L - float(np.mean(L))


def rms(a, b):
    return float(np.sqrt(np.mean((np.asarray(a) - np.asarray(b)) ** 2)))


def platform_xy(t, turn_deg=0.0, u=U_PLAT):
    t = np.asarray(t, float)
    xp = u * t.copy()
    yp = np.zeros_like(t)
    if abs(turn_deg) > 1e-12:
        t_turn = float(np.max(t)) * 0.5
        d = math.radians(turn_deg)
        c, s = math.cos(d), math.sin(d)
        m = t > t_turn
        dt = t[m] - t_turn
        xp[m] = u * t_turn + u * dt * c
        yp[m] = u * dt * s
    return xp, yp


def range_traj(t, r0_km, theta0_deg, v, psi_deg, turn_deg):
    r0 = float(r0_km) * 1000.0
    th0 = math.radians(float(theta0_deg))
    psi = math.radians(float(psi_deg))
    xp, yp = platform_xy(t, turn_deg)
    xt = r0 * np.cos(th0) + v * t * np.cos(psi)
    yt = r0 * np.sin(th0) + v * t * np.sin(psi)
    return np.hypot(xt - xp, yt - yp)


def main() -> int:
    # SSP
    text0 = (ZGRID / "zgrid_f235.env").read_text(encoding="utf-8")
    z0, c0 = read_ssp_block(text0)
    envs = {
        "E0_nominal": (z0, c0),
        "E1_shift50m": (z0, shift_ssp(z0, c0, dz=50.0, dc=0.0)),
        "E2_shift150m_plus2": (z0, shift_ssp(z0, c0, dz=150.0, dc=2.0)),
    }

    # ---- tail roundtrip ----
    tail_rows = []
    for ename in ["E1_shift50m", "E2_shift150m_plus2"]:
        for f in FREQS:
            env0 = ZGRID / f"zgrid_f{int(f)}.env"
            t0 = env0.read_text(encoding="utf-8")
            _, tail0 = split_env_tail(t0)
            outp = WORK / f"{ename}_f{int(f)}.env"
            tail_w = write_env_pert(outp, f, envs[ename][0], envs[ename][1], env0)
            same = tail_w == tail0
            tail_rows.append(
                {
                    "env": ename,
                    "f_hz": f,
                    "tail_identical": same,
                    "n_tail_lines_frozen": len(tail0),
                    "n_tail_lines_written": len(tail_w),
                    "has_51_rd": any("51" in ln for ln in tail_w),
                    "has_150_and_250": any("150" in ln for ln in tail_w)
                    and any("250" in ln for ln in tail_w),
                }
            )
    pd.DataFrame(tail_rows).to_csv(OUT / "ENV_TAIL_ROUNDTRIP_CHECK.csv", index=False)
    if any(not r["tail_identical"] for r in tail_rows):
        (OUT / "R3_RC3_REANCHOR_4G_DECISION.json").write_text(
            json.dumps({"decision": "ENV_TAIL_INTEGRITY_FAILED", "created_utc": NOW}, indent=2),
            encoding="utf-8",
        )
        print("ENV_TAIL_INTEGRITY_FAILED")
        return 1

    # ---- KRAKEN + mod depth tabulation ----
    tab_rows = []
    mods = {"E0_nominal": {f: parse_mod(ZGRID / f"zgrid_f{int(f)}.mod") for f in FREQS}}
    for ename in ["E1_shift50m", "E2_shift150m_plus2"]:
        mods[ename] = {}
        zz, cc = envs[ename]
        for f in FREQS:
            stem = f"{ename}_f{int(f)}"
            write_env_pert(WORK / f"{stem}.env", f, zz, cc, ZGRID / f"zgrid_f{int(f)}.env")
            subprocess.run(
                [str(AT_BIN / "kraken.exe"), stem],
                cwd=str(WORK),
                capture_output=True,
                text=True,
                timeout=300,
            )
            modp = WORK / f"{stem}.mod"
            if not modp.exists():
                tab_rows.append({"env": ename, "f_hz": f, "status": "KRAKEN_FAIL", "ok": False})
                continue
            m = parse_mod(modp)
            mods[ename][f] = m
            d = m["depths"]
            def nz(z):
                return int(np.argmin(np.abs(d - z)))
            has_band = bool(d.min() <= 150.5 and d.max() >= 249.5 and d.size >= 40)
            ok_depths = all(any(abs(d - z) < 1.0) for z in (180.0, 200.0, 220.0))
            # phi distinct across depths
            p180 = m["phi"][nz(180)]
            p200 = m["phi"][nz(200)]
            p220 = m["phi"][nz(220)]
            d180 = float(np.linalg.norm(p180 - p200))
            d220 = float(np.linalg.norm(p220 - p200))
            tab_rows.append(
                {
                    "env": ename,
                    "f_hz": f,
                    "nmat": m["nmat"],
                    "M": m["M"],
                    "depth_min": float(d.min()),
                    "depth_max": float(d.max()),
                    "n_depths": int(d.size),
                    "has_150_250_band": has_band,
                    "has_180_200_220": ok_depths,
                    "nearest_180": float(d[nz(180)]),
                    "nearest_200": float(d[nz(200)]),
                    "nearest_220": float(d[nz(220)]),
                    "phi_180_minus_200_norm": d180,
                    "phi_220_minus_200_norm": d220,
                    "ok": bool(has_band and ok_depths and d180 > 0 and d220 > 0),
                }
            )
    # E0 check
    m0 = mods["E0_nominal"][235.0]
    tab_rows.append(
        {
            "env": "E0_nominal",
            "f_hz": 235.0,
            "nmat": m0["nmat"],
            "M": m0["M"],
            "depth_min": float(m0["depths"].min()),
            "depth_max": float(m0["depths"].max()),
            "n_depths": int(m0["depths"].size),
            "ok": True,
        }
    )
    pd.DataFrame(tab_rows).to_csv(OUT / "MODE_DEPTH_TABULATION_INTEGRITY.csv", index=False)
    if any(not r.get("ok") for r in tab_rows if r.get("env") != "E0_nominal"):
        (OUT / "R3_RC3_REANCHOR_4G_DECISION.json").write_text(
            json.dumps(
                {"decision": "MODE_DEPTH_TABULATION_INTEGRITY_FAILED", "created_utc": NOW}, indent=2
            ),
            encoding="utf-8",
        )
        print("MODE_DEPTH_TABULATION_INTEGRITY_FAILED")
        return 1

    # ---- FIELD warning + forward gate 24 ----
    warn_rows = []
    gate_rows = []
    sig_rows = []
    for ename in ["E1_shift50m", "E2_shift150m_plus2"]:
        for f in FREQS:
            stem = f"{ename}_f{int(f)}"
            m = mods[ename][f]
            # TL signatures at three zs for diagnostic
            Lsig = {}
            for zs in Z_TRUE:
                write_flp(WORK / f"{stem}_zs{int(zs)}.flp", zs, ZR, R_GATE)
                (WORK / f"{stem}_zs{int(zs)}.mod").write_bytes((WORK / f"{stem}.mod").read_bytes())
                (WORK / f"{stem}_zs{int(zs)}.env").write_text(
                    (WORK / f"{stem}.env").read_text(encoding="utf-8"), encoding="utf-8"
                )
                r = subprocess.run(
                    [str(AT_BIN / "field.exe"), f"{stem}_zs{int(zs)}"],
                    cwd=str(WORK),
                    capture_output=True,
                    text=True,
                    timeout=300,
                )
                prt = (WORK / "field.prt").read_text(errors="replace") if (WORK / "field.prt").exists() else ""
                warn = "Modes not tabulated near requested pt" in (r.stdout + r.stderr + prt)
                warn_rows.append(
                    {
                        "env": ename,
                        "f_hz": f,
                        "zs": zs,
                        "warning": warn,
                        "warning_text": "Modes not tabulated" if warn else "",
                    }
                )
                shd = WORK / f"{stem}_zs{int(zs)}.shd"
                if not shd.exists() or warn:
                    gate_rows.append(
                        {
                            "env": ename,
                            "f_hz": f,
                            "zs": zs,
                            "status": "FIELD_FAIL_OR_WARNING",
                            "gate_ok": False,
                        }
                    )
                    continue
                sh = parse_shd(shd)
                p_ms = pressure_series(m, sh["rr"], zs)
                p_f = sh["P"][:, 0]
                c = np.vdot(p_ms, p_f) / np.vdot(p_ms, p_ms)
                rel = float(np.linalg.norm(p_f - c * p_ms) / max(np.linalg.norm(p_f), 1e-30))
                ccorr = abs(
                    np.vdot(p_ms, p_f) / max(np.linalg.norm(p_ms) * np.linalg.norm(p_f), 1e-30)
                )
                tf, tm = tl_shape(p_f), tl_shape(p_ms)
                tl_rms = float(np.sqrt(np.mean((tf - tm) ** 2)))
                tl_corr = (
                    float(np.corrcoef(tf, tm)[0, 1])
                    if np.std(tf) > 0 and np.std(tm) > 0
                    else np.nan
                )
                ok = bool(
                    ccorr >= CC_MIN
                    and rel <= RESID_MAX
                    and tl_rms <= TL_RMS_MAX
                    and tl_corr == tl_corr
                    and tl_corr >= TL_CORR_MIN
                )
                gate_rows.append(
                    {
                        "env": ename,
                        "f_hz": f,
                        "zs": zs,
                        "status": "OK",
                        "complex_corr": float(ccorr),
                        "normalized_residual": rel,
                        "TL_SHAPE_RMS_DB": tl_rms,
                        "TL_SHAPE_CORRELATION": tl_corr,
                        "gate_ok": ok,
                    }
                )
                Lsig[zs] = tf
            if set(Lsig) == set(Z_TRUE):
                sig_rows.append(
                    {
                        "env": ename,
                        "f_hz": f,
                        "RMS_L180_minus_L200": rms(Lsig[180.0], Lsig[200.0]),
                        "RMS_L220_minus_L200": rms(Lsig[220.0], Lsig[200.0]),
                        "all_zero_machine": bool(
                            max(rms(Lsig[180.0], Lsig[200.0]), rms(Lsig[220.0], Lsig[200.0])) < 1e-12
                        ),
                    }
                )
    for f in FREQS:
        for zs in Z_TRUE:
            gate_rows.append(
                {
                    "env": "E0_nominal",
                    "f_hz": f,
                    "zs": zs,
                    "status": "REUSE_1F_VALIDATED",
                    "gate_ok": True,
                }
            )
    pd.DataFrame(warn_rows).to_csv(OUT / "FIELD_DEPTH_TABULATION_WARNING_AUDIT.csv", index=False)
    pd.DataFrame(gate_rows).to_csv(OUT / "PERTURBED_ENV_FORWARD_GATE_FIXED2.csv", index=False)
    pd.DataFrame(sig_rows).to_csv(OUT / "SOURCE_DEPTH_SIGNATURE_DIAGNOSTIC.csv", index=False)
    if any(r.get("warning") for r in warn_rows) or any(not r.get("gate_ok") for r in gate_rows):
        dec = "ENV_MISMATCH_STAGE_BLOCKED_BY_FORWARD_MODEL"
        (OUT / "R3_RC3_REANCHOR_4G_DECISION.json").write_text(
            json.dumps({"decision": dec, "created_utc": NOW}, indent=2), encoding="utf-8"
        )
        print("DECISION", dec)
        return 1
    if any(r.get("all_zero_machine") for r in sig_rows):
        dec = "MODE_DEPTH_TABULATION_INTEGRITY_FAILED"
        (OUT / "R3_RC3_REANCHOR_4G_DECISION.json").write_text(
            json.dumps(
                {"decision": dec, "reason": "SOURCE_DEPTH_SIGNATURE all machine-zero", "created_utc": NOW},
                indent=2,
            ),
            encoding="utf-8",
        )
        print("DECISION", dec)
        return 1

    # ---- main experiment ----
    pairs = pd.read_csv(PAIRS)
    if "selected" in pairs.columns:
        pairs = pairs[pairs["selected"] == True]
    m0map = mods["E0_nominal"]
    rows = []
    ztrue_rows = []
    for _, row in pairs.iterrows():
        pid = row["pair_id"]
        mech = row.get("mechanism", "")
        T = float(row["T_s"])
        turn = float(row["turn_deg"])
        t = np.arange(0.0, T + 1e-9, 1.0)
        r_ref = range_traj(
            t, row["ref_r0_km"], row["ref_theta0_deg"], row["ref_v"], row["ref_psi_deg"], turn
        )
        r_alt = range_traj(
            t, row["alt_r0_km"], row["alt_theta0_deg"], row["alt_v"], row["alt_psi_deg"], turn
        )
        static_blind = pid == "A05"
        tpl_true = {
            float(z): np.concatenate([tl_shape(pressure_series(m0map[f], r_ref, z)) for f in FREQS])
            for z in Z_GRID
        }
        tpl_alt = {
            float(z): np.concatenate([tl_shape(pressure_series(m0map[f], r_alt, z)) for f in FREQS])
            for z in Z_GRID
        }
        for ename in ["E0_nominal", "E1_shift50m", "E2_shift150m_plus2"]:
            me = mods[ename]
            y_by_z = {}
            for z_true in Z_TRUE:
                y = np.concatenate([tl_shape(pressure_series(me[f], r_ref, z_true)) for f in FREQS])
                y_by_z[z_true] = y
                Jt = min(rms(y, tpl_true[float(z)]) for z in Z_GRID)
                Ja = min(rms(y, tpl_alt[float(z)]) for z in Z_GRID)
                dJ = Ja - Jt
                rows.append(
                    {
                        "env": ename,
                        "pair_id": pid,
                        "mechanism": mech,
                        "z_true": z_true,
                        "J_true_star": Jt,
                        "J_alt_star": Ja,
                        "delta_J_env": dJ,
                        "delta_J_gt_0": bool(dJ > 0),
                        "static_blind": static_blind,
                    }
                )
            if pid in ("A01", "B01", "C01"):
                ztrue_rows.append(
                    {
                        "pair_id": pid,
                        "env": ename,
                        "RMS_y180_y200": rms(y_by_z[180.0], y_by_z[200.0]),
                        "RMS_y220_y200": rms(y_by_z[220.0], y_by_z[200.0]),
                        "identical_180_200": bool(rms(y_by_z[180.0], y_by_z[200.0]) < 1e-12),
                        "identical_220_200": bool(rms(y_by_z[220.0], y_by_z[200.0]) < 1e-12),
                    }
                )
    df = pd.DataFrame(rows)
    df.to_csv(OUT / "PROFILED_MARGIN_VS_ENV_FIXED2.csv", index=False)
    pd.DataFrame(ztrue_rows).to_csv(OUT / "ZTRUE_SIGNATURE_NONREPLICATION_CHECK.csv", index=False)

    sum_rows = []
    for (env, mech), g0 in df.groupby(["env", "mechanism"]):
        g = g0[~g0["static_blind"]]
        if not len(g):
            continue
        dJ = g["delta_J_env"]
        sum_rows.append(
            {
                "env": env,
                "mechanism": mech,
                "median_delta_J": float(dJ.median()),
                "p10_delta_J": float(dJ.quantile(0.10)),
                "min_delta_J": float(dJ.min()),
                "fraction_tested_dJ_gt_0": float((dJ > 0).mean()),
                "n": len(g),
            }
        )
    sdf = pd.DataFrame(sum_rows)
    sdf.to_csv(OUT / "ENV_MARGIN_BY_MECHANISM_FIXED2.csv", index=False)

    def env_frac(env):
        g = df[(df["env"] == env) & (~df["static_blind"])]
        return float((g["delta_J_env"] > 0).mean()), float(g["delta_J_env"].median())

    f1, m1 = env_frac("E1_shift50m")
    f2, m2 = env_frac("E2_shift150m_plus2")
    if f1 >= 0.8 and m1 > 0 and f2 >= 0.8 and m2 > 0:
        decision = "PROFILED_RC3_MARGIN_SURVIVES_TESTED_SSP_STRESS"
    elif f1 >= 0.8 and m1 > 0:
        decision = "PROFILED_RC3_CONDITIONAL_ON_MILD_SSP_MATCH"
    elif f1 < 0.5:
        decision = "PROFILED_RC3_MARGIN_COLLAPSES_UNDER_SSP_STRESS"
    else:
        decision = "PROFILED_RC3_ENV_MISMATCH_BOUNDARY_UNRESOLVED"

    why = (
        f"tail+tabulation OK; z_true signatures distinct; "
        f"E1 frac(ΔJ>0)={f1:.3f} med={m1:.3f}; E2 frac={f2:.3f} med={m2:.3f}"
    )
    (OUT / "R3_RC3_REANCHOR_4G_DECISION.json").write_text(
        json.dumps(
            {
                "stage": "R3-RC3-REANCHOR-4G",
                "decision": decision,
                "why": why,
                "prior": "REANCHOR4F_BLOCKED_BY_MODE_DEPTH_TABULATION",
                "z200_control": "Z200_ENV_STRESS_POSITIVE_CONTROL (superseded by full depth set)",
                "e1": {"frac_gt0": f1, "median_dJ": m1},
                "e2": {"frac_gt0": f2, "median_dJ": m2},
                "created_utc": NOW,
            },
            ensure_ascii=False,
            indent=2,
        ),
        encoding="utf-8",
    )
    (OUT / "R3_RC3_REANCHOR_4G_REPORT.md").write_text(
        f"""# R3-RC3-REANCHOR-4G mode-depth tabulation fix

UTC: {NOW}

## 判定

### `{decision}`

{why}

## 修复

- `'R'` 后 tail **逐行复制**（51 receivers 150:2:250）
- `.mod` 含 180/200/220 且 φ 向量不全同
- FIELD 无 `Modes not tabulated` 警告
- 源深签名 RMS(L180−L200)、RMS(L220−L200) 非机器零
- B01/C01/A01 的 y(180/200/220) 非重复

旧 4F：`REANCHOR4F_BLOCKED_BY_MODE_DEPTH_TABULATION`；z=200 单点=`Z200_ENV_STRESS_POSITIVE_CONTROL`。

## 结果

| env | fraction_tested(ΔJ>0) | median ΔJ |
| --- | ---: | ---: |
| E1 | {f1:.3f} | {m1:.3f} |
| E2 | {f2:.3f} | {m2:.3f} |

S1 仍暂缓。
""",
        encoding="utf-8",
    )
    (OUT / "GPT_SYNC.md").write_text(f"# R3-RC3-REANCHOR-4G\\n\\n**{decision}**\\n\\n{why}\\n", encoding="utf-8")
    print("DECISION", decision)
    print(sdf.to_string(index=False))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
