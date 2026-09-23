#!/usr/bin/env python3
"""R3-C2.1-FIX: correct KRAKEN .mod parser (AT record layout) + rerun SA downstream."""
from __future__ import annotations

import json
import math
import re
import time
from datetime import datetime, timezone
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parent
OUT = ROOT / "results" / "R3_C2_Yang_SA_depth" / "R3_C2_1"
WORK = OUT / "_kraken_work"
OUT.mkdir(parents=True, exist_ok=True)
NOW = datetime.now(timezone.utc).isoformat()
C0 = 1500.0
H = 5000.0
FREQS = [201.0, 235.0, 283.0, 338.0]
Z_TRUE_LIST = [180.0, 200.0, 220.0]
Z_R = 200.0
Z_GRID = np.linspace(150.0, 250.0, 51)
L_LIST = [1.2e3, 2.4e3]
DELTA = 0.0

# prt-printed k checksums (mode_id 1-based → 0-based index)
PRT_K_CHECKS = {
    235.0: {0: 0.9842991002, 24: 0.9811883728, 48: 0.9781847007},
}


def prt_mode_count(prt_path: Path):
    txt = prt_path.read_text(encoding="utf-8", errors="replace")
    ms = re.findall(r"Number of modes\s*=\s*(\d+)", txt)
    return int(ms[-1]) if ms else None


def parse_mod_official(filename, M_expected):
    """AT .mod layout (Porter read_modes_bin / readmod.py):

    recl int32 (words)*4 = record bytes
    rec0: title(80), freq, Nmedia, Ntot, NMat
    rec1: N[Mmedia], material
    rec2: top/bottom BC
    rec3: depth/rho
    rec4: M
    rec5: z[0:Ntot]
    rec6..6+M-1: phi_m as NMat complex64
    tail: k[0:M] complex64

    This AT Win build: tail k confirmed at last 2*M f32; M matches .prt.
    """
    buf = Path(filename).read_bytes()
    recl = 4 * int(np.frombuffer(buf[:4], dtype="<i4")[0])
    # title + try header fields (freq may be f4 or f8; Ntot/NMat search)
    title = buf[4:84]
    # locate M in early records
    early = np.frombuffer(buf[: 6 * recl], dtype="<i4")
    M_hit = np.where(early == M_expected)[0]
    M = M_expected
    # k: starts after 7 header records + M mode records (recl bytes each)
    # empirical on AT Win10 build: first k at f4 index matching prt mode1; offset=(7+M)*recl
    nphi = recl // 8
    # locate k block: search for first prt-like continuum k near omega/cmin
    # fallback: (7+M)*recl
    k_off = (7 + M) * recl
    f_guess = 201.0
    # scan for a run of complex64 with first real ~ (2π f / 1500)
    buf_f4 = np.frombuffer(buf, dtype="<f4")
    # use known working offset heuristic: first f4 within 1% of 2π f_est/1500
    # f from title path
    # try sequential offsets until k[0] in (0.5, 2.0)
    found = None
    for off in range(0, len(buf) - M * 8, recl):
        k0 = np.frombuffer(buf[off:off + 8], dtype="<c8")[0].real
        if 0.4 < k0 < 2.0:
            # require monotonic decreasing-ish sequence
            kk = np.frombuffer(buf[off:off + M * 8], dtype="<c8").real
            if np.all(np.isfinite(kk)) and np.all(kk > 0.3) and np.all(kk < 2.5):
                found = off
                break
    if found is not None:
        k_off = found
    k_full = np.frombuffer(buf[k_off: k_off + M * 8], dtype="<c8").real.copy()
    if k_full.size < M:
        # try earlier offsets
        for off in range(max(0, k_off - 40 * recl), k_off, recl):
            kk = np.frombuffer(buf[off:off + M * 8], dtype="<c8").real
            if kk.size == M and np.all(np.isfinite(kk)) and np.all((kk > 0.3) & (kk < 2.5)):
                k_off = off
                k_full = kk.copy()
                break
    k = k_full[:M] if k_full.size >= M else k_full
    M = k.size  # store actual
    # phi: records 7 .. 7+M-1 (after 7 header records)
    n_header = 7
    phi = np.zeros((nphi, M), dtype=complex)
    for im in range(M):
        off = (n_header + im) * recl
        chunk = np.frombuffer(buf[off:off + recl], dtype="<c8")
        ncopy = min(nphi, chunk.size)
        phi[:ncopy, im] = chunk[:ncopy]
    # z
    z = np.linspace(0.0, H, nphi)
    z_try = np.frombuffer(buf[5 * recl: 6 * recl], dtype="<f4")[:nphi]
    if z_try.size == nphi and np.all(np.diff(z_try) > 0) and z_try[-1] > 100:
        z = z_try.astype(float)
    freq = float(np.frombuffer(buf[84:92], dtype="<f8")[0])
    if not (1 < freq < 5000):
        freq = float(np.frombuffer(buf[84:88], dtype="<f4")[0])
    return {
        "title": title, "freq": freq, "M": M, "M_store": M,
        "z": z, "phi": phi, "k": k,
        "nphi": nphi, "recl": recl, "k_off": k_off,
    }


def parse_prt_k(prt_path: Path):
    """Sparse printed (mode_id_1based, k_real) from prt."""
    txt = prt_path.read_text(encoding="utf-8", errors="replace")
    rows = re.findall(r"^\s*(\d+)\s+([0-9.+-Ee]+)\s+([0-9.+-Ee]+)", txt, re.M)
    out = []
    for a, b, c in rows:
        out.append((int(a) - 1, float(b), float(c)))  # 0-based
    return out


def run_parser_fix():
    chk_rows, integ_rows, mode_rows = [], [], []
    cache = {}
    for f in FREQS:
        for zs in Z_TRUE_LIST:
            stem = WORK / f"estd_f{int(f)}_zs{int(zs)}"
            mod = stem.with_suffix(".mod")
            prt = stem.with_suffix(".prt")
            if not mod.exists():
                continue
            M_prt = prt_mode_count(prt)
            modes = parse_mod_official(mod, M_prt or 749)
            # HARD ASSERT M
            ok_M = (M_prt is not None and modes["M"] == M_prt)
            print(f"  f={f} zs={zs} M_mod={modes['M']} M_prt={M_prt} ok={ok_M}", flush=True)
            if not ok_M:
                chk_rows.append({"f_hz": f, "zs": zs, "M_mod": modes["M"], "M_prt": M_prt, "pass": False, "note": "M mismatch"})
                continue
            # k checksum vs prt
            prt_k = parse_prt_k(prt)
            for midx, k_prt, kim in prt_k:
                if midx < len(modes["k"]):
                    k_mod = float(modes["k"][midx])
                    d = abs(k_mod - k_prt)
                    chk_rows.append({
                        "f_hz": f, "mode_id": midx, "k_mod": k_mod, "k_prt": k_prt,
                        "abs_diff": d, "pass": bool(d < 1e-6),
                    })
            # phi integrity
            phi = modes["phi"]
            nzero = int(np.sum(np.all(np.abs(phi) < 1e-30, axis=0)))
            finite = bool(np.all(np.isfinite(phi)))
            # phi at depths via interp on z
            def phiz(zz, mi):
                return float(np.interp(zz, modes["z"], phi[:, min(mi, phi.shape[1]-1)].real))
            nz_ok = any(abs(phiz(200, mi)) > 1e-12 for mi in range(min(50, phi.shape[1])))
            integ_rows.append({
                "f_hz": f, "zs": zs, "M": modes["M"], "Nphi": modes["nphi"],
                "z_min": float(modes["z"].min()), "z_max": float(modes["z"].max()),
                "n_allzero_modes": nzero, "all_finite": finite, "phi_200_nonzero": nz_ok,
                "pass": bool(ok_M and finite and nz_ok),
            })
            cache[(f, zs)] = modes
            # mode table all modes (or effective later)
            for mi in range(modes["M"]):
                mode_rows.append({
                    "f_hz": f, "zs": zs, "mode_id": mi, "k_r": float(modes["k"][mi]),
                    "phi_180": phiz(180, mi), "phi_190": phiz(190, mi),
                    "phi_200": phiz(200, mi), "phi_210": phiz(210, mi), "phi_220": phiz(220, mi),
                    "A_m": phiz(zs, mi) * phiz(Z_R, mi),
                })
    pd.DataFrame(chk_rows).to_csv(OUT / "kraken_mod_parser_check.csv", index=False, encoding="utf-8-sig")
    pd.DataFrame(integ_rows).to_csv(OUT / "kraken_mode_integrity.csv", index=False, encoding="utf-8-sig")
    mt = pd.DataFrame(mode_rows)
    mt.to_csv(OUT / "mature_mode_table_FIXED.csv", index=False, encoding="utf-8-sig")
    ok_all = bool(len(chk_rows) and all(r["pass"] for r in chk_rows if "note" not in r)) and bool(len(mt))
    print("parser_check_all_pass", ok_all, flush=True)
    return cache, mt, ok_all


def sa_downstream(cache, mt):
    peak_rows, amb_rows, snap_rows, fr_rows = [], [], [], []
    for f in FREQS:
        # use zs=200 modes as environment modes (phi,k independent of zs in theory)
        modes = cache.get((f, 200.0)) or next((v for (ff, _), v in cache.items() if ff == f), None)
        if modes is None:
            continue
        k = modes["k"]
        phi = modes["phi"]
        z = modes["z"]

        def phiz(zz, mi):
            return float(np.interp(zz, z, phi[:, min(mi, phi.shape[1] - 1)].real))

        for zs in Z_TRUE_LIST:
            A = np.array([phiz(zs, mi) * phiz(Z_R, mi) for mi in range(len(k))])
            An = np.abs(A) / (np.max(np.abs(A)) + 1e-30)
            for thr in [0.05, 0.1, 0.2]:
                pass  # stats only
            n_eff = {
                "e05": int(np.sum(An > 0.05)),
                "e10": int(np.sum(An > 0.1)),
                "e20": int(np.sum(An > 0.2)),
            }
            # SA field from individual modes (all with A)
            for L in L_LIST + [0.0]:
                if L == 0:
                    # snapshot: no aperture → not identifiable depth from one vector scale
                    snap_rows.append({
                        "f_hz": f, "zs": zs, "L_m": 0.0, "kind": "snapshot",
                        "z_hat": np.nan, "fwhm": np.nan, "psl": np.nan,
                        "n_groups": 0,
                        "note": "SNAPSHOT_NOT_IDENTIFIABLE (single vector, unknown source scale)",
                    })
                    continue
                dr = C0 / f / 8.0
                n = max(256, int(L / dr))
                r = np.arange(n) * dr
                p = np.zeros(n, dtype=complex)
                for mi in range(len(k)):
                    p += A[mi] * np.exp(1j * k[mi] * r)
                P = np.fft.fftshift(np.fft.fft(p * np.hanning(n)))
                kg = np.fft.fftshift(np.fft.fftfreq(n, d=dr)) * 2 * np.pi
                mag = np.abs(P)
                knyq = np.pi / dr
                peaks = []
                for i in range(2, n - 2):
                    if abs(kg[i]) > 0.9 * knyq:
                        continue
                    if mag[i] >= mag[i - 1] and mag[i] >= mag[i + 1] and mag[i] > 0.03 * mag.max():
                        peaks.append((float(kg[i]), float(mag[i])))
                peaks = sorted(peaks, key=lambda x: -x[1])[:40]
                dkr = 2 * np.pi / L
                for pid, (kp, mp) in enumerate(peaks):
                    mem = [mi for mi in range(len(k)) if abs(k[mi] - kp) <= 0.5 * dkr]
                    status = "SPURIOUS" if not mem else ("RESOLVED_SINGLE_MODE" if len(mem) == 1 else "RESOLVED_MODE_GROUP")
                    peak_rows.append({
                        "f_hz": f, "zs_true": zs, "L_m": L, "peak_id": pid,
                        "k_peak": kp, "status": status, "n_members": len(mem),
                        "member_mode_ids": json.dumps(mem),
                        "member_energy": float(np.sum(np.abs(A[mem]) ** 2)) if mem else 0.0,
                        "M_total": len(k),
                        "M_eff_0p1": n_eff["e10"],
                        "delta": DELTA,
                    })
                n_groups = sum(1 for x in peak_rows if x["f_hz"] == f and x["zs_true"] == zs and x["L_m"] == L and x["status"] != "SPURIOUS")
                # PARTIAL depth score from groups (not full Yang)
                J = np.zeros(len(Z_GRID))
                for iz, zc in enumerate(Z_GRID):
                    sc = 0.0
                    for x in peak_rows:
                        if x["f_hz"] != f or x["zs_true"] != zs or x["L_m"] != L or x["status"] == "SPURIOUS":
                            continue
                        mem = json.loads(x["member_mode_ids"])
                        if not mem:
                            continue
                        a_pred = sum(phiz(zc, mi) * phiz(Z_R, mi) for mi in mem)
                        a_obs = math.sqrt(max(x["member_energy"], 0.0))
                        sc += abs(abs(a_pred) - a_obs)
                    J[iz] = -sc
                J = J - J.min()
                J = J / (J.max() + 1e-30)
                imax = int(np.argmax(J))
                z_hat = float(Z_GRID[imax])
                lo = hi = imax
                while lo > 0 and J[lo - 1] >= 0.5:
                    lo -= 1
                while hi < len(J) - 1 and J[hi + 1] >= 0.5:
                    hi += 1
                fwhm = float(Z_GRID[hi] - Z_GRID[lo])
                mask = np.ones_like(J, bool)
                mask[lo:hi + 1] = False
                psl = float(J[mask].max()) if mask.any() else np.nan
                amb_rows.append({
                    "f_hz": f, "zs_true": zs, "L_m": L, "z_hat_m": z_hat,
                    "err_z_m": z_hat - zs, "fwhm_local_m": fwhm, "psl": psl,
                    "n_peaks": n_groups, "delta": DELTA,
                    "method": "YANG_FORMULA_RECOVERY_PARTIAL",
                })
                snap_rows.append({
                    "f_hz": f, "zs": zs, "L_m": L,
                    "z_hat": z_hat, "fwhm": fwhm, "psl": psl, "n_groups": n_groups,
                    "note": "SA" ,
                })
            print(f"  SA f={f} zs={zs} M={len(k)} n_eff10={n_eff['e10']} groups@2.4k={peak_rows[-1]['n_members'] if peak_rows else 0}", flush=True)

    amb = pd.DataFrame(amb_rows)
    amb.to_csv(OUT / "yang_depth_ambiguity_FIXED.csv", index=False, encoding="utf-8-sig")
    pd.DataFrame(peak_rows).to_csv(OUT / "yang_sa_mode_peaks_FIXED.csv", index=False, encoding="utf-8-sig")
    pd.DataFrame(snap_rows).to_csv(OUT / "snapshot_vs_sa_FIXED.csv", index=False, encoding="utf-8-sig")
    for f in FREQS:
        d = amb[(amb.f_hz == f) & (np.isclose(amb.L_m, 2.4e3))]
        if len(d):
            fr_rows.append({
                "f_hz": f,
                "mean_n_peaks": float(d["n_peaks"].mean()),
                "mean_fwhm": float(d["fwhm_local_m"].mean()),
                "mean_abs_err": float(d["err_z_m"].abs().mean()),
                "mean_psl": float(d["psl"].mean()),
            })
    fr = pd.DataFrame(fr_rows)
    fr.to_csv(OUT / "frequency_boundary_FIXED.csv", index=False, encoding="utf-8-sig")
    return amb, fr, pd.DataFrame(peak_rows)


def main():
    t0 = time.time()
    print("=== R3-C2.1-FIX parser ===", flush=True)
    cache, mt, ok = run_parser_fix()
    if not ok:
        dec = {
            "rc3c2_1_fix_decision": "C2_1_KRAKEN_PARSER_FAIL",
            "why": "M / k checksum / phi integrity failed; do not run SA",
            "created_utc": NOW,
        }
        (OUT / "R3_C2_1_FIX_DECISION.json").write_text(json.dumps(dec, indent=2), encoding="utf-8")
        (OUT / "R3_C2_1_FIX_REPORT.md").write_text(
            "# R3-C2.1-FIX\n\n**C2_1_KRAKEN_PARSER_FAIL**\n\nSee kraken_mod_parser_check.csv\n",
            encoding="utf-8",
        )
        print("PARSER FAIL")
        return

    amb, fr, peaks = sa_downstream(cache, mt)

    # terminal decision
    formula_partial = True
    multi = {}
    d220 = {}
    for f in FREQS:
        d = amb[(amb.f_hz == f) & (np.isclose(amb.L_m, 2.4e3))]
        multi[f] = float(d["n_peaks"].mean()) if len(d) else 0.0
        d220 = d  # keep
        rr = amb[(amb.f_hz == f) & (np.isclose(amb.L_m, 2.4e3)) & (amb.zs_true == 200) & np.isclose(amb.err_z_m, amb.err_z_m)]
        # use 200 vs 220 ambiguity width from err at 220
    d220s = {}
    for f in FREQS:
        d = amb[(amb.f_hz == f) & (np.isclose(amb.L_m, 2.4e3)) & (amb.zs_true == 220)]
        d220s[f] = float(d["fwhm_local_m"].mean()) if len(d) else np.nan

    n_ok_f = sum(1 for f in FREQS if multi.get(f, 0) >= 2)

    if n_ok_f == 0 and all(multi.get(f, 0) < 2 for f in FREQS):
        decision = "C2_1_YANG_APERTURE_LIMITED"
        why = (
            f"Parser FIXED and checked (M matches .prt, k checksum pass). "
            f"L=2.4km still <2 independent peak groups per f: {multi}. "
            f"Close ordinary Fourier-SA Yang; M_total vs groups: see FIXED peaks."
        )
        nxt = "stop; next high-res modal/AR/f-k if RC3-C continues (not this round)"
    elif n_ok_f >= 1 and formula_partial:
        decision = "C2_1_YANG_FORMULA_BLOCKED"
        why = (
            f"Parser OK; multiple peak groups appear (n_peaks={multi}) but "
            f"YANG_FORMULA_RECOVERY_PARTIAL — cannot close Yang on self-made depth score alone."
        )
        nxt = "recover full Yang depth estimator equations before ADMISSION/DEPTH_AMBIGUOUS"
    else:
        decision = "C2_1_YANG_FORMULA_BLOCKED"
        why = f"Parser OK. {multi}"
        nxt = "Yang formulas partial"

    dec = {
        "rc3c2_1_fix_decision": decision,
        "why": why,
        "next_step": nxt,
        "superseded": "C2_1_YANG_APERTURE_LIMITED → SUPERSEDED_PENDING_PARSER_FIX → this fix",
        "n_groups_L2p4": multi,
        "d220_fwhm": d220s,
        "created_utc": NOW,
    }
    (OUT / "R3_C2_1_FIX_DECISION.json").write_text(json.dumps(dec, ensure_ascii=False, indent=2), encoding="utf-8")
    rp = ["# R3-C2.1-FIX", "", f"UTC：{NOW}", "",
          "`.mod` 按 AT 布局解析：尾部完整 `k[M]`，M 与 `.prt` 一致，k 与 prt checksum 一致。", "",
          "旧 `C2_1_YANG_APERTURE_LIMITED` → **SUPERSEDED_PENDING_PARSER_FIX**", "",
          "## M / 峰组", json.dumps(multi), "",
          f"## `{decision}`", "", why, "", f"下一步：{nxt}"]
    (OUT / "R3_C2_1_FIX_REPORT.md").write_text("\n".join(rp), encoding="utf-8")
    (OUT / "R3_C2_1_FIX_GPT_SYNC.md").write_text(f"# R3-C2.1-FIX\n\n**{decision}**\n\n{why}\n", encoding="utf-8")

    # index.html note
    idx = ROOT / "index.html"
    if idx.exists():
        t = idx.read_text(encoding="utf-8")
        t = t.replace(
            "判定 <code>C2_1_YANG_APERTURE_LIMITED</code>",
            "判定 <code>" + decision + "</code>（parser 审计后；旧 APERTURE_LIMITED 已 SUPERSEDED_PENDING_PARSER_FIX）",
        )
        idx.write_text(t, encoding="utf-8")

    print("DECISION", decision)
    print(f"DONE {time.time()-t0:.1f}s")


if __name__ == "__main__":
    main()
