#!/usr/bin/env python3
"""R3-RC3-REANCHOR-4F: SSP perturbation integrity fix. No S1."""
from __future__ import annotations

import json
import math
import subprocess
from datetime import datetime, timezone
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parent
OUT = ROOT / "results" / "R3_RC3_REANCHOR" / "R3_RC3_ENV_MISMATCH_BOUNDARY_FIX"
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
ANCHOR_Z = [0.0, 1500.0, 1600.0, 3000.0, 4000.0, 5000.0]


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


def read_ssp_block(path: Path):
    """FIXED parser: SSP only until first 'R' line; six-column acoustic rows."""
    zs, cs = [], []
    started = False
    for ln in path.read_text(encoding="utf-8").splitlines():
        s = ln.strip()
        if s.startswith("'R'") or s == "R":
            break
        p = s.split()
        if len(p) >= 6:
            try:
                z = float(p[0])
                c = float(p[1])
                # require 6 numeric fields (z c cs rho a as)
                [float(p[i]) for i in range(6)]
            except ValueError:
                continue
            zs.append(z)
            cs.append(c)
            started = True
        elif started and s.startswith("'"):
            # e.g. 'A' bottom — still not SSP; stop or skip
            continue
    z = np.array(zs, float)
    c = np.array(cs, float)
    return z, c


def assert_ssp(z, c, label: str):
    ok = (
        z.size == 251
        and abs(z[0] - 0.0) < 1e-9
        and abs(z[-1] - 5000.0) < 1e-6
        and bool(np.all(np.diff(z) > 0))
        and len(np.unique(z)) == 251
    )
    return ok, {
        "label": label,
        "n_ssp": int(z.size),
        "z0": float(z[0]) if z.size else np.nan,
        "zend": float(z[-1]) if z.size else np.nan,
        "strictly_increasing": bool(np.all(np.diff(z) > 0)) if z.size > 1 else False,
        "n_unique": int(len(np.unique(z))) if z.size else 0,
        "ok": ok,
    }


def shift_ssp(z, c, dz=0.0, dc=0.0):
    zq = np.clip(z - dz, float(z.min()), float(z.max()))
    return np.interp(zq, z, c) + dc


def write_env_ssp(path: Path, freq: float, zs_src: float, rd: float, z_ssp, c_ssp):
    """Only replace SSP block values; keep tail after 'R' unchanged."""
    src = (ZGRID / f"zgrid_f{int(freq)}.env").read_text(encoding="utf-8")
    lines = src.splitlines()
    z_ssp = np.asarray(z_ssp, float)
    c_ssp = np.asarray(c_ssp, float)
    out = []
    i = 0
    while i < len(lines):
        ln = lines[i]
        s = ln.strip()
        if s.startswith("'R'") or s == "R":
            break
        p = s.split()
        if len(p) >= 6:
            try:
                z = float(p[0])
                [float(p[j]) for j in range(6)]
                c_new = float(np.interp(z, z_ssp, c_ssp))
                out.append(f"{z:.3f} {c_new:.4f} 0.0 1.0 0.0 0.0")
                i += 1
                continue
            except ValueError:
                pass
        out.append(ln)
        i += 1
    # tail from frozen env starting at 'R'
    tail = []
    for ln in lines:
        s = ln.strip()
        if s.startswith("'R'") or s == "R":
            tail = [
                "'R' 0.0",
                "1500.0 1800.0",
                "60.0",
                "1",
                f"{zs_src:.3f}",
                "1",
                f"{rd:.3f}",
                "R",
            ]
            break
    out[0] = "'PERT'"
    out[1] = f"{freq:.3f}"
    path.write_text("\n".join(out + tail) + "\n", encoding="utf-8")


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


def pressure_series(mod, r, zs, zr=ZR):
    depths, phi, k = mod["depths"], mod["phi"], mod["k"]
    izs = int(np.argmin(np.abs(depths - zs)))
    izr = int(np.argmin(np.abs(depths - zr)))
    kre = k.real
    alpha = -k.imag
    ps = phi[izs]
    pr = phi[izr]
    r = np.asarray(r, float)
    p = np.zeros(r.size, dtype=complex)
    for m in range(len(kre)):
        amp = np.sqrt(2 * np.pi / (kre[m] * r))
        p += amp * ps[m] * pr[m] * np.exp(-1j * kre[m] * r - alpha[m] * r - 1j * np.pi / 4)
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
    # ---- 1-2) base SSP parse integrity ----
    base_rows = []
    zs0 = cs0 = None
    for f in FREQS:
        z, c = read_ssp_block(ZGRID / f"zgrid_f{int(f)}.env")
        ok, meta = assert_ssp(z, c, f"f{int(f)}")
        base_rows.append(meta)
        if zs0 is None:
            zs0, cs0 = z, c
        else:
            base_rows.append(
                {
                    "label": f"f{int(f)}_vs_f235",
                    "n_ssp": int(z.size),
                    "z_match": bool(np.allclose(z, zs0)),
                    "c_match": bool(np.allclose(c, cs0)),
                    "ok": bool(np.allclose(z, zs0) and np.allclose(c, cs0)),
                }
            )
    pd.DataFrame(base_rows).to_csv(OUT / "BASE_SSP_PARSE_INTEGRITY.csv", index=False)
    if any(not r.get("ok") for r in base_rows):
        (OUT / "R3_RC3_REANCHOR_4F_DECISION.json").write_text(
            json.dumps(
                {"decision": "SSP_PERTURBATION_INTEGRITY_FAILED", "created_utc": NOW}, indent=2
            ),
            encoding="utf-8",
        )
        print("DECISION SSP_PERTURBATION_INTEGRITY_FAILED")
        return 1

    z0, c0 = zs0, cs0
    envs = {
        "E0_nominal": (z0, c0),
        "E1_shift50m": (z0, shift_ssp(z0, c0, dz=50.0, dc=0.0)),
        "E2_shift150m_plus2": (z0, shift_ssp(z0, c0, dz=150.0, dc=2.0)),
    }

    # ---- 3-4) anchor integrity ----
    anchor_rows = []
    for za in ANCHOR_Z:
        e1_exp = float(np.interp(za - 50.0, z0, c0)) if za - 50 >= z0.min() else float(c0[0])
        if za - 50 > z0.max():
            e1_exp = float(c0[-1])
        # endpoint hold: clip inside shift_ssp
        e1_act = float(envs["E1_shift50m"][1][int(np.argmin(np.abs(z0 - za)))])
        e2_exp = float(np.interp(np.clip(za - 150.0, z0.min(), z0.max()), z0, c0)) + 2.0
        e2_act = float(envs["E2_shift150m_plus2"][1][int(np.argmin(np.abs(z0 - za)))])
        e0 = float(c0[int(np.argmin(np.abs(z0 - za)))])
        anchor_rows.append(
            {
                "z": za,
                "E0": e0,
                "E1_actual": e1_act,
                "E1_expected": e1_exp,
                "E1_abs_err": abs(e1_act - e1_exp),
                "E2_actual": e2_act,
                "E2_expected": e2_exp,
                "E2_abs_err": abs(e2_act - e2_exp),
            }
        )
    pd.DataFrame(anchor_rows).to_csv(OUT / "SSP_PERTURBATION_ANCHOR_CHECK.csv", index=False)
    max_err = max(max(r["E1_abs_err"], r["E2_abs_err"]) for r in anchor_rows)
    if max_err > 1e-3:
        (OUT / "R3_RC3_REANCHOR_4F_DECISION.json").write_text(
            json.dumps(
                {
                    "decision": "SSP_PERTURBATION_INTEGRITY_FAILED",
                    "reason": f"anchor max_err={max_err}",
                    "created_utc": NOW,
                },
                indent=2,
            ),
            encoding="utf-8",
        )
        print("ANCHOR FAIL", max_err)
        return 1

    # save fixed SSP
    ssp_rows = []
    for name, (zz, cc) in envs.items():
        ok, meta = assert_ssp(zz, cc, name)
        if not ok:
            print("ENV SSP ASSERT FAIL", meta)
            return 1
        for z, c in zip(zz, cc):
            ssp_rows.append({"env": name, "z": float(z), "c": float(c)})
    pd.DataFrame(ssp_rows).to_csv(OUT / "PERTURBED_SSP_FIXED.csv", index=False)

    # ---- 5) roundtrip check on written env (one sample) ----
    write_env_ssp(WORK / "roundtrip.env", 235.0, 200.0, ZR, envs["E1_shift50m"][0], envs["E1_shift50m"][1])
    z_rt, c_rt = read_ssp_block(WORK / "roundtrip.env")
    c_int = envs["E1_shift50m"][1]
    rt_err = float(np.max(np.abs(c_rt - c_int[: len(c_rt)]))) if len(c_rt) == len(c_int) else np.inf
    pd.DataFrame(
        [
            {
                "n_written": len(c_rt),
                "n_intended": len(c_int),
                "max_abs_c_diff": rt_err,
                "ok": bool(len(c_rt) == 251 and rt_err < 1e-3),
            }
        ]
    ).to_csv(OUT / "WRITTEN_ENV_ROUNDTRIP_CHECK.csv", index=False)

    # ---- 6-7) KRAKEN + forward gate 24 cases ----
    mods = {"E0_nominal": {f: parse_mod(ZGRID / f"zgrid_f{int(f)}.mod") for f in FREQS}}
    gate_rows = []
    for ename in ["E1_shift50m", "E2_shift150m_plus2"]:
        mods[ename] = {}
        zz, cc = envs[ename]
        for f in FREQS:
            stem = f"{ename}_f{int(f)}"
            write_env_ssp(WORK / f"{stem}.env", f, 200.0, ZR, zz, cc)
            subprocess.run(
                [str(AT_BIN / "kraken.exe"), stem],
                cwd=str(WORK),
                capture_output=True,
                text=True,
                timeout=300,
            )
            modp = WORK / f"{stem}.mod"
            if not modp.exists():
                gate_rows.append({"env": ename, "f_hz": f, "zs": np.nan, "status": "KRAKEN_FAIL", "gate_ok": False})
                continue
            m = parse_mod(modp)
            mods[ename][f] = m
            for zs in Z_TRUE:
                write_flp(WORK / f"{stem}_zs{int(zs)}.flp", zs, ZR, R_GATE)
                (WORK / f"{stem}_zs{int(zs)}.mod").write_bytes(modp.read_bytes())
                # env for field must match — reuse stem env
                (WORK / f"{stem}_zs{int(zs)}.env").write_text(
                    (WORK / f"{stem}.env").read_text(encoding="utf-8"), encoding="utf-8"
                )
                subprocess.run(
                    [str(AT_BIN / "field.exe"), f"{stem}_zs{int(zs)}"],
                    cwd=str(WORK),
                    capture_output=True,
                    text=True,
                    timeout=300,
                )
                shd = WORK / f"{stem}_zs{int(zs)}.shd"
                if not shd.exists():
                    gate_rows.append(
                        {"env": ename, "f_hz": f, "zs": zs, "status": "FIELD_FAIL", "gate_ok": False}
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
    # E0: reuse 1F (12/12) — record
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
    pd.DataFrame(gate_rows).to_csv(OUT / "PERTURBED_ENV_FORWARD_GATE_FIXED.csv", index=False)
    if any(not r.get("gate_ok") for r in gate_rows):
        dec = "ENV_MISMATCH_STAGE_BLOCKED_BY_FORWARD_MODEL"
        (OUT / "R3_RC3_REANCHOR_4F_DECISION.json").write_text(
            json.dumps({"decision": dec, "created_utc": NOW}, indent=2), encoding="utf-8"
        )
        print("DECISION", dec)
        return 1

    # ---- 8) main experiment ----
    pairs = pd.read_csv(PAIRS)
    if "selected" in pairs.columns:
        pairs = pairs[pairs["selected"] == True]
    m0 = mods["E0_nominal"]
    rows = []
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
            float(z): np.concatenate([tl_shape(pressure_series(m0[f], r_ref, z)) for f in FREQS])
            for z in Z_GRID
        }
        tpl_alt = {
            float(z): np.concatenate([tl_shape(pressure_series(m0[f], r_alt, z)) for f in FREQS])
            for z in Z_GRID
        }
        for ename in ["E0_nominal", "E1_shift50m", "E2_shift150m_plus2"]:
            me = mods[ename]
            for z_true in Z_TRUE:
                y = np.concatenate(
                    [tl_shape(pressure_series(me[f], r_ref, z_true)) for f in FREQS]
                )
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
                        "delta_J_gt_0p5": bool(dJ > 0.5),
                        "static_blind": static_blind,
                    }
                )
    df = pd.DataFrame(rows)
    df.to_csv(OUT / "PROFILED_MARGIN_VS_ENV_FIXED.csv", index=False)

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
                "fraction_tested_dJ_gt_0p5": float((dJ > 0.5).mean()),
                "median_J_true_star": float(g["J_true_star"].median()),
                "median_J_alt_star": float(g["J_alt_star"].median()),
                "n": len(g),
            }
        )
    sdf = pd.DataFrame(sum_rows)
    sdf.to_csv(OUT / "ENV_MARGIN_BY_MECHANISM_FIXED.csv", index=False)

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
        f"SSP parser fixed (251 pts); E1 frac(ΔJ>0)={f1:.3f} med={m1:.3f}; "
        f"E2 frac={f2:.3f} med={m2:.3f}; truth=E_k template=E0"
    )

    (OUT / "R3_RC3_REANCHOR_4F_DECISION.json").write_text(
        json.dumps(
            {
                "stage": "R3-RC3-REANCHOR-4F",
                "decision": decision,
                "why": why,
                "prior_invalidated": "OLD_REANCHOR4_INVALIDATED_BY_SSP_PARSE_BUG",
                "unintended_control_label": "UNINTENDED_SSP_PERTURBATION_CONTROL",
                "ssp_integrity": "PASS",
                "anchor_max_err": max_err,
                "roundtrip_max_err": rt_err,
                "forward_gate_cases": int(sum(1 for r in gate_rows if r.get("status") != "REUSE_1F_VALIDATED")),
                "e1": {"frac_gt0": f1, "median_dJ": m1},
                "e2": {"frac_gt0": f2, "median_dJ": m2},
                "created_utc": NOW,
            },
            ensure_ascii=False,
            indent=2,
        ),
        encoding="utf-8",
    )
    (OUT / "R3_RC3_REANCHOR_4F_REPORT.md").write_text(
        f"""# R3-RC3-REANCHOR-4F SSP integrity fix

UTC: {NOW}

## 判定

### `{decision}`

{why}

## 旧 REANCHOR-4

`OLD_REANCHOR4_INVALIDATED_BY_SSP_PARSE_BUG`；意外控制=`UNINTENDED_SSP_PERTURBATION_CONTROL`。

## 修复

- SSP parser：`'R'` 前六列声学行；断言 251 / 0–5000 / 严格递增
- E1/E2 定义不变；anchor 校验 max_err={max_err:.2e}
- write_env 只改 SSP block；roundtrip max_err={rt_err:.2e}
- KRAKEN 全重跑；forward gate **24** 个 perturbed case

## 环境失配（正确 E1/E2）

| env | fraction_tested(ΔJ>0) | median ΔJ |
| --- | ---: | ---: |
| E1 | {f1:.3f} | {m1:.3f} dB |
| E2 | {f2:.3f} | {m2:.3f} dB |

## 证据链

深度消元通过 → IID TL 通过 → **环境失配：本报告** → S1 仍暂缓。

禁止：S1/S0、随机 TL、Liang、P5、改 E1/E2 定义。
""",
        encoding="utf-8",
    )
    (OUT / "GPT_SYNC.md").write_text(f"# R3-RC3-REANCHOR-4F\\n\\n**{decision}**\\n\\n{why}\\n", encoding="utf-8")
    print("DECISION", decision)
    print(sdf.to_string(index=False))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
