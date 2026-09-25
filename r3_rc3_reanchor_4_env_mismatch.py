#!/usr/bin/env python3
"""R3-RC3-REANCHOR-4: env-mismatch boundary for depth-profiled RC3. No S1, no random TL error."""
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
OUT = ROOT / "results" / "R3_RC3_REANCHOR" / "R3_RC3_ENV_MISMATCH_BOUNDARY"
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


def read_ssp_from_env(path: Path):
    zs, cs = [], []
    for ln in path.read_text(encoding="utf-8").splitlines():
        p = ln.split()
        if len(p) >= 2:
            try:
                z = float(p[0])
                c = float(p[1])
            except ValueError:
                continue
            if 0.0 <= z <= 5000.0:
                zs.append(z)
                cs.append(c)
    return np.array(zs), np.array(cs)


def shift_ssp(z, c, dz=0.0, dc=0.0):
    """c_new(z) = c_old(z-dz) + dc; endpoint extrapolation."""
    zq = np.clip(z - dz, float(z.min()), float(z.max()))
    cnew = np.interp(zq, z, c) + dc
    return cnew


def write_env_ssp(path: Path, freq: float, zs_src: float, rd: float, z_ssp, c_ssp):
    """Clone frozen zgrid env skeleton; only replace SSP c(z) values."""
    src = (ZGRID / f"zgrid_f{int(freq)}.env").read_text(encoding="utf-8")
    lines = src.splitlines()
    z_ssp = np.asarray(z_ssp, float)
    c_ssp = np.asarray(c_ssp, float)
    out = []
    for ln in lines:
        p = ln.split()
        # SSP lines: z c ...
        if len(p) >= 2:
            try:
                z = float(p[0])
                c_old = float(p[1])
            except ValueError:
                out.append(ln)
                continue
            if 0.0 <= z <= 5000.0 and len(p) >= 6:
                c_new = float(np.interp(z, z_ssp, c_ssp))
                out.append(f"{z:.3f} {c_new:.4f} 0.0 1.0 0.0 0.0")
                continue
        out.append(ln)
    # fix title and tail source/receiver block
    out[0] = "'PERT'"
    out[1] = f"{freq:.3f}"
    # replace last R-section: find first line starting with 'R'
    for i, ln in enumerate(out):
        if ln.strip().startswith("'R'"):
            out = out[:i] + [
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


def run_kraken_field(env_stem: str, freq, zs, zr, r_m, z_ssp, c_ssp, mod_src: Path):
    write_env_ssp(WORK / f"{env_stem}.env", freq, zs, zr, z_ssp, c_ssp)
    # kraken
    # reuse zgrid env structure: actually run kraken on our env
    r = subprocess.run(
        [str(AT_BIN / "kraken.exe"), env_stem],
        cwd=str(WORK),
        capture_output=True,
        text=True,
        timeout=180,
    )
    modp = WORK / f"{env_stem}.mod"
    if not modp.exists():
        return None, None, f"KRAKEN_FAIL {r.stderr[-80:] if r.stderr else ''}"
    write_flp(WORK / f"{env_stem}.flp", zs, zr, r_m)
    (WORK / f"{env_stem}.mod").write_bytes(modp.read_bytes())
    r2 = subprocess.run(
        [str(AT_BIN / "field.exe"), env_stem],
        cwd=str(WORK),
        capture_output=True,
        text=True,
        timeout=300,
    )
    shd = WORK / f"{env_stem}.shd"
    if not shd.exists():
        return parse_mod(modp), None, f"FIELD_FAIL {r2.stderr[-80:] if r2.stderr else ''}"
    sh = parse_shd(shd)
    return parse_mod(modp), sh, "OK"


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
    # SSP from frozen E-STD env (E0)
    z0, c0 = read_ssp_from_env(ZGRID / "zgrid_f235.env")
    envs = {
        "E0_nominal": (z0, c0),
        "E1_shift50m": (z0, shift_ssp(z0, c0, dz=50.0, dc=0.0)),
        "E2_shift150m_plus2": (z0, shift_ssp(z0, c0, dz=150.0, dc=2.0)),
    }
    # save SSPs
    ssp_rows = []
    for name, (zz, cc) in envs.items():
        for z, c in zip(zz, cc):
            ssp_rows.append({"env": name, "z": float(z), "c": float(c)})
    pd.DataFrame(ssp_rows).to_csv(OUT / "PERTURBED_SSP.csv", index=False)

    (OUT / "ENV_MISMATCH_AXIS_LOCK.md").write_text(
        f"""# ENV_MISMATCH_AXIS_LOCK

UTC: {NOW}

`PROJECT_PREDEFINED_SSP_MISMATCH_STRESS_AXIS`（工程敏感性轴，非海况概率）

| 档 | 定义 |
| --- | --- |
| E0 | nominal E-STD（冻结 zgrid SSP） |
| E1 | `SSP_VERTICAL_SHIFT_50M`：c(z)→c(z−50 m)，端点外推 |
| E2 | `SSP_VERTICAL_SHIFT_150M_PLUS_2MPS`：c(z)→c(z−150 m)+2 m/s |

信号固定 `S2_LIKE_TRACKABLE_FOUR_TONE_UPPER_BOUND`；ε=0（无随机 TL 误差）。

**真值用 E_k；估计器模板只用 E0**（否则不是失配）。
""",
        encoding="utf-8",
    )

    # build env-specific modes + forward gate + store mods
    mods = {"E0_nominal": {}}
    gate_rows = []
    for ename in ["E1_shift50m", "E2_shift150m_plus2"]:
        mods[ename] = {}
    for ename, (zz, cc) in envs.items():
        for f in FREQS:
            stem = f"{ename}_f{int(f)}"
            if ename == "E0_nominal":
                # use frozen zgrid mods
                m = parse_mod(ZGRID / f"zgrid_f{int(f)}.mod")
                mods[ename][f] = m
                # gate with FIELD vs frozen mode sum (already validated) — still record
                ok_field = True
                cc_max, res_max = 1.0, 0.0
                # E0: reuse 1F-validated FIELD vs mode-sum; no re-field required
                gate_rows.append(
                    {
                        "env": ename,
                        "f_hz": f,
                        "status": "REUSE_1F_VALIDATED",
                        "complex_corr": np.nan,
                        "normalized_residual": np.nan,
                        "TL_SHAPE_RMS_DB": np.nan,
                        "TL_SHAPE_CORRELATION": np.nan,
                        "gate_ok": True,
                    }
                )
                continue
            m, sh, st = run_kraken_field(
                stem, f, 200.0, ZR, R_GATE, zz, cc, ZGRID / f"zgrid_f{int(f)}.mod"
            )
            if sh is None or m is None:
                gate_rows.append(
                    {
                        "env": ename,
                        "f_hz": f,
                        "status": st,
                        "gate_ok": False,
                    }
                )
                continue
            mods[ename][f] = m
            # compare field vs mode-sum in THIS env (self-consistency of new KRAKEN)
            p_ms = pressure_series(m, sh["rr"], 200.0)
            p_f = sh["P"][:, 0]
            c = np.vdot(p_ms, p_f) / np.vdot(p_ms, p_ms)
            rel = float(np.linalg.norm(p_f - c * p_ms) / max(np.linalg.norm(p_f), 1e-30))
            ccorr = abs(np.vdot(p_ms, p_f) / max(np.linalg.norm(p_ms) * np.linalg.norm(p_f), 1e-30))
            tf, tm = tl_shape(p_f), tl_shape(p_ms)
            tl_rms = float(np.sqrt(np.mean((tf - tm) ** 2)))
            tl_corr = (
                float(np.corrcoef(tf, tm)[0, 1]) if np.std(tf) > 0 and np.std(tm) > 0 else np.nan
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
                    "status": st,
                    "complex_corr": float(ccorr),
                    "normalized_residual": rel,
                    "TL_SHAPE_RMS_DB": tl_rms,
                    "TL_SHAPE_CORRELATION": tl_corr,
                    "gate_ok": ok,
                }
            )
    pd.DataFrame(gate_rows).to_csv(OUT / "PERTURBED_ENV_FORWARD_GATE.csv", index=False)

    if any(not r.get("gate_ok") for r in gate_rows):
        dec = "ENV_MISMATCH_STAGE_BLOCKED_BY_FORWARD_MODEL"
        (OUT / "R3_RC3_REANCHOR_4_DECISION.json").write_text(
            json.dumps({"decision": dec, "created_utc": NOW}, indent=2), encoding="utf-8"
        )
        print("DECISION", dec)
        return 1

    # load E0 truth/matching mods (templates)
    m0 = mods["E0_nominal"]

    pairs = pd.read_csv(PAIRS)
    if "selected" in pairs.columns:
        pairs = pairs[pairs["selected"] == True]

    rows = []
    for _, row in pairs.iterrows():
        pid = row["pair_id"]
        mech = row.get("mechanism", "")
        T = float(row["T_s"])
        turn = float(row["turn_deg"])
        t = np.arange(0.0, T + 1e-9, 1.0)
        r_ref = range_traj(t, row["ref_r0_km"], row["ref_theta0_deg"], row["ref_v"], row["ref_psi_deg"], turn)
        r_alt = range_traj(t, row["alt_r0_km"], row["alt_theta0_deg"], row["alt_v"], row["alt_psi_deg"], turn)
        static_blind = pid == "A05"

        # E0 templates for matching (all z)
        tpl_true = {}
        tpl_alt = {}
        for z in Z_GRID:
            tpl_true[float(z)] = np.concatenate(
                [tl_shape(pressure_series(m0[f], r_ref, z)) for f in FREQS]
            )
            tpl_alt[float(z)] = np.concatenate(
                [tl_shape(pressure_series(m0[f], r_alt, z)) for f in FREQS]
            )

        for ename in ["E0_nominal", "E1_shift50m", "E2_shift150m_plus2"]:
            me = mods[ename]
            for z_true in Z_TRUE:
                # obs from true env E_k
                y = np.concatenate([tl_shape(pressure_series(me[f], r_ref, z_true)) for f in FREQS])
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
                        "blind_label": "RELATIVE_TL_STATIC_RANGE_BLIND_CASE" if static_blind else "",
                    }
                )

    df = pd.DataFrame(rows)
    df.to_csv(OUT / "PROFILED_MARGIN_VS_ENV.csv", index=False)

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
    sdf.to_csv(OUT / "ENV_MARGIN_BY_MECHANISM.csv", index=False)

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
        f"E1 fraction(ΔJ>0)={f1:.3f} median={m1:.3f} dB; "
        f"E2 fraction={f2:.3f} median={m2:.3f} dB; "
        f"truth=E_k, template=E0; S2 four-tone; no random TL error"
    )

    (OUT / "R3_RC3_REANCHOR_4_DECISION.json").write_text(
        json.dumps(
            {
                "stage": "R3-RC3-REANCHOR-4",
                "decision": decision,
                "why": why,
                "axis": "PROJECT_PREDEFINED_SSP_MISMATCH_STRESS_AXIS",
                "signal": "S2_LIKE_TRACKABLE_FOUR_TONE_UPPER_BOUND",
                "e1": {"frac_gt0": f1, "median_dJ": m1},
                "e2": {"frac_gt0": f2, "median_dJ": m2},
                "roles_frozen": {
                    "science": "DEPTH_ROLE_MIXED_OR_UNRESOLVED",
                    "system": "DEPTH_PROFILED_NUISANCE_VARIABLE",
                },
                "reanchor3_note": "TESTED_ERROR=IID_ZERO_MEAN_TL_PERTURBATION; fraction_tested not probability",
                "not_done": ["S0/S1", "random TL error", "Liang", "P5", "depth RMSE"],
                "created_utc": NOW,
            },
            ensure_ascii=False,
            indent=2,
        ),
        encoding="utf-8",
    )
    (OUT / "R3_RC3_REANCHOR_4_REPORT.md").write_text(
        f"""# R3-RC3-REANCHOR-4 env-mismatch boundary

UTC: {NOW}

## 判定

### `{decision}`

{why}

## 轴

E0 nominal · E1 SSP 垂向平移 50 m · E2 平移 150 m + 2 m/s  
`PROJECT_PREDEFINED_SSP_MISMATCH_STRESS_AXIS`（工程敏感性，非海况概率）

**真值 E_k，模板恒为 E0**（ε=0，无随机 TL 误差；S2 四稳定线谱）

## 结果

| env | fraction_tested(ΔJ>0) | median ΔJ |
| --- | ---: | ---: |
| E1 | {f1:.3f} | {m1:.3f} dB |
| E2 | {f2:.3f} | {m2:.3f} dB |

A05 继续 `RELATIVE_TL_STATIC_RANGE_BLIND_CASE`（不进主统计）。

## 口径（相对 REANCHOR-3）

- `fraction_tested` 非概率（5 seeds 非独立 MC）
- 仅 IID 零均值 TL 扰动已测；慢变/环境误差另论（本阶段即环境）
- 2 dB 最差 ΔJ≈0.11 dB（A10）仍为正，但不宽裕

## 禁止

S0/S1、随机 TL、Liang、P5、重挑样本、改 E1/E2。
""",
        encoding="utf-8",
    )
    (OUT / "GPT_SYNC.md").write_text(f"# R3-RC3-REANCHOR-4\\n\\n**{decision}**\\n\\n{why}\\n", encoding="utf-8")
    print("DECISION", decision)
    print(sdf.to_string(index=False))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
