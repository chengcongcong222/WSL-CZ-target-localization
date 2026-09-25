#!/usr/bin/env python3
"""R3-RC3-REANCHOR-1F: E-STD forward integrity gate only. No rerun of 90 depth cases."""
from __future__ import annotations

import json
import subprocess
from datetime import datetime, timezone
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parent
OUT = ROOT / "results" / "R3_RC3_REANCHOR" / "R3_RC3_DEPTH_ROLE_AUDIT_FIX"
ZGRID = ROOT / "results" / "R3_C2_Yang_SA_depth" / "R3_C2_1" / "_kraken_zgrid"
AT_BIN = ROOT / "tools" / "acoustics_toolbox" / "atWin10" / "at" / "bin"
WORK = OUT / "_field"
OUT.mkdir(parents=True, exist_ok=True)
WORK.mkdir(parents=True, exist_ok=True)
NOW = datetime.now(timezone.utc).isoformat()

ZR = 200.0
FREQS = [201.0, 235.0, 283.0, 338.0]
ZS_LIST = [180.0, 200.0, 220.0]
R_POINTS = np.array([45000.0, 48000.0, 51000.0, 54000.0, 57000.0, 60000.0])
# gate thresholds: reuse spirit of FIELD_EQ3_MULTIMODE_VALIDATED
CC_MIN = 0.95
RESID_MAX = 0.35
TL_RMS_MAX_DB = 1.5
TL_CORR_MIN = 0.95


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
    """Verified parser from r3_c2_3a_estd.py — do not rewrite offsets."""
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
    src = (ZGRID / f"zgrid_f{int(freq)}.env").read_text(encoding="utf-8")
    lines = src.splitlines()
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
        if ln.strip().startswith("'R'") or ln.strip() == "R":
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
    subprocess.run(
        [str(AT_BIN / "field.exe"), stem],
        cwd=str(WORK),
        capture_output=True,
        text=True,
        timeout=300,
    )
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


def main() -> int:
    mods = {f: parse_mod(ZGRID / f"zgrid_f{int(f)}.mod") for f in FREQS}
    rows = []
    tl_rows = []
    gate_pass = True

    for f in FREQS:
        for zs in ZS_LIST:
            rr, p_fld = field_pressure(f, zs, ZR, R_POINTS, tag="1f_")
            p_ms = pressure_series(mods[f], R_POINTS, zs)
            if p_fld is None:
                rows.append({"f": f, "zs": zs, "status": "FIELD_FAIL"})
                gate_pass = False
                continue
            # one global complex scale over whole range vector
            c = np.vdot(p_ms, p_fld) / np.vdot(p_ms, p_ms)
            resid = p_fld - c * p_ms
            rel_res = float(np.linalg.norm(resid) / max(np.linalg.norm(p_fld), 1e-30))
            cc = abs(np.vdot(p_ms, p_fld) / max(np.linalg.norm(p_ms) * np.linalg.norm(p_fld), 1e-30))
            a_rel = float(
                np.linalg.norm(np.abs(p_fld) - np.abs(c) * np.abs(p_ms))
                / max(np.linalg.norm(np.abs(p_fld)), 1e-30)
            )
            ph = np.angle(p_fld * np.conj(c * p_ms))
            ph_rmse = float(np.sqrt(np.mean(ph**2)))
            # TL shape (observable actually used)
            tf = tl_shape(p_fld)
            tm = tl_shape(p_ms)
            tl_rms = float(np.sqrt(np.mean((tf - tm) ** 2)))
            tl_corr = float(np.corrcoef(tf, tm)[0, 1]) if np.std(tf) > 0 and np.std(tm) > 0 else np.nan
            ok = bool(
                cc >= CC_MIN
                and rel_res <= RESID_MAX
                and tl_rms <= TL_RMS_MAX_DB
                and (tl_corr == tl_corr and tl_corr >= TL_CORR_MIN)
            )
            if not ok:
                gate_pass = False
            rows.append(
                {
                    "f_hz": f,
                    "zs": zs,
                    "n_ranges": int(len(R_POINTS)),
                    "c_star_re": float(c.real),
                    "c_star_im": float(c.imag),
                    "complex_corr": float(cc),
                    "normalized_complex_residual": rel_res,
                    "amplitude_shape_residual": a_rel,
                    "phase_rmse_rad": ph_rmse,
                    "TL_SHAPE_RMS_DB": tl_rms,
                    "TL_SHAPE_CORRELATION": tl_corr,
                    "gate_ok": ok,
                    "status": "OK" if ok else "MISMATCH",
                }
            )
            tl_rows.append(
                {
                    "f_hz": f,
                    "zs": zs,
                    "TL_SHAPE_RMS_DB": tl_rms,
                    "TL_SHAPE_CORRELATION": tl_corr,
                    "observable": "SOURCE_LEVEL_FREE_RELATIVE_TL_SHAPE",
                }
            )

    pd.DataFrame(rows).to_csv(OUT / "FIELD_VS_MODESUM_MULTIRANGE.csv", index=False)
    pd.DataFrame(tl_rows).to_csv(OUT / "RELATIVE_TL_SHAPE_VALIDATION.csv", index=False)

    n_ok = sum(1 for r in rows if r.get("gate_ok"))
    n_tot = len(rows)
    decision = (
        "ESTD_RELATIVE_TL_FORWARD_MODEL_VALIDATED"
        if gate_pass and n_ok == n_tot
        else "DEPTH_ROLE_BLOCKED_BY_FORWARD_MODEL"
    )
    why = (
        f"multi-range FIELD vs mode-sum: {n_ok}/{n_tot} cases pass "
        f"(cc>={CC_MIN}, resid<={RESID_MAX}, TL_rms<={TL_RMS_MAX_DB} dB, TL_corr>={TL_CORR_MIN}); "
        f"one global c* per (f,zs); no per-point scale; REANCHOR-1 90 cases not recomputed"
        if gate_pass
        else f"forward gate failed {n_ok}/{n_tot}"
    )

    (OUT / "FORWARD_GATE_AUDIT.md").write_text(
        f"""# FORWARD_GATE_AUDIT

UTC: {NOW}

## 旧 Gate 无效

单点 `c=p_fld/p_ms` ⇒ residual 恒为 0，`field_ok=True` 无信息量。

## 修正

- 复用 `r3_c2_3a_estd.py` 的 `parse_shd` / `write_env_estd` / `write_flp` / `field_pressure`
- **整条 range vector** 一个 `c* = p_MS^H p_FIELD / p_MS^H p_MS`
- 禁止逐点 scale
- 报告 complex corr / normalized residual / amp-shape / phase + **TL_SHAPE**（真正使用的 observable）
- Gate 失败 **必须停止**（不再 `if fail: pass`）

## Gate 阈值（沿用 FIELD_EQ3_MULTIMODE_VALIDATED 精神）

- complex_corr ≥ {CC_MIN}
- normalized residual ≤ {RESID_MAX}
- TL_SHAPE_RMS_DB ≤ {TL_RMS_MAX_DB}
- TL_SHAPE_CORRELATION ≥ {TL_CORR_MIN}

## 结果

{decision}

{why}

REANCHOR-1 的 90 case **不重算**；仅解除 forward Gate。
""",
        encoding="utf-8",
    )

    (OUT / "R3_RC3_REANCHOR_1F_DECISION.json").write_text(
        json.dumps(
            {
                "stage": "R3-RC3-REANCHOR-1F",
                "decision": decision,
                "why": why,
                "cases_ok": n_ok,
                "cases_total": n_tot,
                "gate": {
                    "complex_corr_min": CC_MIN,
                    "resid_max": RESID_MAX,
                    "tl_rms_db_max": TL_RMS_MAX_DB,
                    "tl_corr_min": TL_CORR_MIN,
                },
                "reanchor1_cases_recomputed": False,
                "temp_science_label": "NO_STRONG_DEPTH_COMPENSATION_OBSERVED",
                "created_utc": NOW,
            },
            ensure_ascii=False,
            indent=2,
        ),
        encoding="utf-8",
    )
    (OUT / "GPT_SYNC.md").write_text(
        f"# R3-RC3-REANCHOR-1F\\n\\n**{decision}**\\n\\n{why}\\n",
        encoding="utf-8",
    )
    print("DECISION", decision, f"{n_ok}/{n_tot}")
    for r in rows:
        print(r.get("f_hz"), r.get("zs"), r.get("complex_corr"), r.get("normalized_complex_residual"), r.get("TL_SHAPE_RMS_DB"), r.get("gate_ok"))
    return 0 if gate_pass else 1


if __name__ == "__main__":
    raise SystemExit(main())
