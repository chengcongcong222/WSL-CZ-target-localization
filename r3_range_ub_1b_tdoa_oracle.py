#!/usr/bin/env python3
"""RANGE-UB-1B: E-STD single-depth oracle multipath-delay range information upper bound.

BELLHOP eigenray arrivals at z_r=200m, source depths 150:5:250m, range 45:0.5:60km.
Source-time-free delay signature. Zero error. r-z identifiability.
NOT reproducing Xu & Guo 2024.
"""
from __future__ import annotations

import json
import math
import subprocess
from datetime import datetime, timezone
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parent
OUT = ROOT / "results" / "RANGE_UB" / "CZ_TDOA_SINGLE_DEPTH_ORACLE"
ZGRID = ROOT / "results" / "R3_C2_Yang_SA_depth" / "R3_C2_1" / "_kraken_zgrid"
AT_BIN = ROOT / "tools" / "acoustics_toolbox" / "atWin10" / "at" / "bin"
WORK = OUT / "_bellhop"
OUT.mkdir(parents=True, exist_ok=True)
WORK.mkdir(parents=True, exist_ok=True)
NOW = datetime.now(timezone.utc).isoformat()

FREQ = 235.0
ZR = 200.0
H_DEPTH = 5000.0
Z_S_LIST = np.arange(150.0, 250.0 + 1e-9, 5.0)  # 21
R_LIST_KM = np.arange(45.0, 60.0 + 1e-9, 0.5)  # 31
Z_TRUE_LIST = [180.0, 200.0, 220.0]
TRACK_R_KM = [45.0, 50.0, 56.0, 58.0, 60.0]


def write_bellhop_env(path: Path, freq: float, zs: float, zr: float, r_list_km, tag: str):
    """Write BELLHOP env for arrival computation using E-STD SSP."""
    src_env = ZGRID / f"zgrid_f{int(freq)}.env"
    lines = src_env.read_text(encoding="utf-8").splitlines()
    ssp_lines = []
    for ln in lines[5:]:
        parts = ln.split()
        if len(parts) >= 2:
            try:
                z = float(parts[0])
                c = float(parts[1])
                ssp_lines.append(f"  {z:.1f} {c:.4f} /")
            except ValueError:
                break
        else:
            break

    n_ssp = len(ssp_lines)
    r_km = float(r_list_km[0]) if len(r_list_km) == 1 else float(max(r_list_km))
    r_max_m = r_km * 1000.0 + 5000.0  # extra margin

    out = []
    out.append(f"'{tag}'")
    out.append(f"{freq:.1f}")
    out.append("1")
    out.append("'CVN'")
    out.append(f"  51 0.0 {H_DEPTH:.1f}")
    out.extend(ssp_lines)
    out.append("'A' 0.0")
    out.append(f"  {H_DEPTH:.1f} 1600.0 0.0 1.0 0.0 0.0 /")
    out.append("1")
    out.append(f"{zs:.1f} /")
    out.append("1")
    out.append(f"{zr:.1f} /")
    out.append(f"{len(r_list_km)}")
    if len(r_list_km) == 1:
        out.append(f"{r_km:.1f} /")
    else:
        out.append(f"{min(r_list_km):.1f} {r_km:.1f} /")
    out.append("'A'")  # arrivals ASCII
    out.append("201")
    out.append("-85.0 85.0 /")
    out.append(f"0.0 {r_max_m:.1f} {H_DEPTH:.1f}")
    path.write_text("\n".join(out) + "\n", encoding="utf-8")


def parse_arr_ascii(arr_path: Path) -> list[dict]:
    """Parse BELLHOP .arr ASCII file. Returns list of arrivals.

    Format after header counts:
    amp  delay  angle_src  x  angle_rcv  y  ntop  nbot
    """
    arrivals = []
    if not arr_path.exists():
        return arrivals
    text = arr_path.read_text(encoding="utf-8", errors="replace")
    lines = text.splitlines()
    for ln in lines:
        parts = ln.split()
        if len(parts) >= 6:
            try:
                amp = float(parts[0])
                delay = float(parts[1])
                a_src = float(parts[2])
                a_rcv = float(parts[4]) if len(parts) > 4 else 0.0
                ntop = int(float(parts[6])) if len(parts) > 6 else 0
                nbot = int(float(parts[7])) if len(parts) > 7 else 0
                arrivals.append(
                    {
                        "amp": abs(amp),
                        "delay": delay,
                        "angle_src_deg": a_src,
                        "angle_rcv_deg": a_rcv,
                        "n_top": ntop,
                        "n_bot": nbot,
                    }
                )
            except ValueError:
                pass
    return arrivals


def delay_signature(arrivals: list[dict], min_amp_ratio: float = 1e-6) -> np.ndarray | None:
    """Source-time-free delay signature: sorted delay differences from earliest."""
    if len(arrivals) < 2:
        return None
    amps = np.array([a["amp"] for a in arrivals])
    delays = np.array([a["delay"] for a in arrivals])
    # keep arrivals above threshold
    mask = amps >= min_amp_ratio * amps.max()
    delays = np.sort(delays[mask])
    if len(delays) < 2:
        return None
    return delays - delays[0]  # source-time-free


def rms(a, b):
    return float(np.sqrt(np.mean((np.asarray(a) - np.asarray(b)) ** 2)))


def main() -> int:
    config = {
        "stage": "RANGE-UB-1B",
        "created_utc": NOW,
        "baseline_commit": "dda967d77c1870217be8973d279c6e0d4533dad1",
        "question": "Does ideal multipath delay at z_r=200m alone break r-z degeneracy?",
        "not_reproducing": "Xu & Guo 2024 paper formulas",
        "freq_hz": FREQ,
        "zr_m": ZR,
        "z_s_list_m": Z_S_LIST.tolist(),
        "r_list_km": R_LIST_KM.tolist(),
        "oracle_assumptions": [
            "PERFECT_MULTIPATH_DETECTION",
            "PERFECT_ARRIVAL_ASSOCIATION",
            "ZERO_DELAY_MEASUREMENT_ERROR",
        ],
        "observable": "SOURCE_TIME_FREE_MULTIPATH_DELAY_SIGNATURE",
        "forbidden": [
            "2024 paper formula reproduction",
            "vertical array",
            "noise sweep",
            "RC2/RC3 fusion",
            "platform turn",
            "P5",
        ],
    }
    (OUT / "RANGE_UB_1B_CONFIG.json").write_text(
        json.dumps(config, ensure_ascii=False, indent=2), encoding="utf-8"
    )

    # =========================================================
    # 1. BELLHOP arrivals for each z_s
    # =========================================================
    # Store: (z_s, r_km) -> list of arrivals
    arrival_db: dict[tuple[float, float], list[dict]] = {}
    integrity_rows = []

    for zs in Z_S_LIST:
        tag = f"oracle_zs{int(zs)}"
        env_path = WORK / f"{tag}.env"
        write_bellhop_env(env_path, FREQ, float(zs), ZR, R_LIST_KM, tag)

        # Run BELLHOP with arrivals
        rc = subprocess.run(
            [str(AT_BIN / "bellhop.exe"), tag],
            cwd=str(WORK),
            capture_output=True,
            text=True,
            timeout=120,
        )
        arr_path = WORK / f"{tag}.arr"
        if not arr_path.exists():
            # try .arr from different naming
            arr_candidates = list(WORK.glob(f"{tag}*.arr"))
            if arr_candidates:
                arr_path = arr_candidates[0]

        arrivals_all = parse_arr_ascii(arr_path)

        # Distribute arrivals to ranges (BELLHOP .arr may have all in one stream)
        # For simplicity, use eigenray mode per-range if needed
        # First try: parse all arrivals as belonging to the target range
        # If multiple ranges, we need to map them
        # For now, re-run per key range for arrival integrity
        pass

    # Better approach: run BELLHOP per (z_s, r) for the key points first
    # to check arrival integrity, then expand

    # Arrival integrity check on subset
    key_zs = [180.0, 200.0, 220.0]
    key_rs = list(np.arange(45.0, 60.5, 1.0))  # 16 ranges

    for zs in key_zs:
        for r_km in key_rs:
            tag = f"chk_zs{int(zs)}_r{r_km:.0f}"
            env_path = WORK / f"{tag}.env"
            write_bellhop_env(env_path, FREQ, float(zs), ZR, [r_km], tag)
            rc = subprocess.run(
                [str(AT_BIN / "bellhop.exe"), tag],
                cwd=str(WORK),
                capture_output=True,
                text=True,
                timeout=60,
            )
            arr_path = WORK / f"{tag}.arr"
            if not arr_path.exists():
                cands = list(WORK.glob(f"{tag}*.arr"))
                arr_path = cands[0] if cands else arr_path

            arrivals = parse_arr_ascii(arr_path)
            n_arr = len(arrivals)
            if n_arr > 0:
                delays = sorted([a["delay"] for a in arrivals])
                d_spread = delays[-1] - delays[0]
                earliest, latest = delays[0], delays[-1]
            else:
                d_spread = earliest = latest = np.nan

            integrity_rows.append(
                {
                    "z_s_m": zs,
                    "r_km": r_km,
                    "n_arrivals": n_arr,
                    "earliest_s": earliest,
                    "latest_s": latest,
                    "delay_spread_s": d_spread,
                    "has_2plus": n_arr >= 2,
                    "bellhop_rc": rc.returncode,
                }
            )
            arrival_db[(zs, r_km)] = arrivals

    integ_df = pd.DataFrame(integrity_rows)
    integ_df.to_csv(OUT / "ARRIVAL_STRUCTURE_MAP.csv", index=False)

    n_with_2 = int(integ_df["has_2plus"].sum())
    n_total = len(integ_df)
    frac_2plus = n_with_2 / n_total if n_total > 0 else 0.0

    # If most points have <2 arrivals, TDOA not available
    if frac_2plus < 0.3:
        decision = "ORACLE_TDOA_NOT_AVAILABLE_IN_ESTD_PROPAGATION"
        (OUT / "RANGE_UB_1B_DECISION.json").write_text(
            json.dumps(
                {
                    "stage": "RANGE-UB-1B",
                    "decision": decision,
                    "fraction_with_2plus_arrivals": frac_2plus,
                    "n_checked": n_total,
                    "created_utc": NOW,
                },
                indent=2,
            ),
            encoding="utf-8",
        )
        (OUT / "ARRIVAL_INTEGRITY_REPORT.md").write_text(
            f"# ARRIVAL_INTEGRITY_REPORT\n\nUTC: {NOW}\n\n"
            f"检查点数: {n_total}, ≥2 arrivals: {n_with_2} ({frac_2plus*100:.1f}%)\n\n"
            f"**{decision}**\n\n第一CZ单深度不足以产生稳定多途时延签名。\n",
            encoding="utf-8",
        )
        print("decision", decision, "frac_2plus", frac_2plus)
        return 0

    (OUT / "ARRIVAL_INTEGRITY_REPORT.md").write_text(
        f"# ARRIVAL_INTEGRITY_REPORT\n\nUTC: {NOW}\n\n"
        f"检查点数: {n_total}, ≥2 arrivals: {n_with_2} ({frac_2plus*100:.1f}%)\n\n"
        f"Gate: PASS (≥30% 有多途)\n\n{integ_df.to_string(index=False)}\n",
        encoding="utf-8",
    )

    # =========================================================
    # 2. Full r-z scoring with delay signatures
    # =========================================================
    # Need arrivals for all (z_s, r) pairs
    # Use key_zs key_rs already computed; extend to full grid
    full_zs = [180.0, 200.0, 220.0]  # truth depths
    full_rs = list(R_LIST_KM)  # 31 ranges

    sig_db: dict[tuple[float, float], np.ndarray | None] = {}
    for zs in full_zs:
        for r_km in full_rs:
            key = (zs, r_km)
            if key in arrival_db:
                sig_db[key] = delay_signature(arrival_db[key])
            else:
                # compute
                tag = f"sc_zs{int(zs)}_r{r_km:.1f}"
                env_path = WORK / f"{tag}.env"
                write_bellhop_env(env_path, FREQ, float(zs), ZR, [r_km], tag)
                subprocess.run(
                    [str(AT_BIN / "bellhop.exe"), tag],
                    cwd=str(WORK),
                    capture_output=True,
                    text=True,
                    timeout=60,
                )
                arr_path = WORK / f"{tag}.arr"
                if not arr_path.exists():
                    cands = list(WORK.glob(f"{tag}*.arr"))
                    arr_path = cands[0] if cands else arr_path
                arrivals = parse_arr_ascii(arr_path)
                arrival_db[key] = arrivals
                sig_db[key] = delay_signature(arrivals)

    # Also need candidate z values (profile) — but this is source depth in arrivals
    # For r-z identifiability: candidates are (r, z_s) pairs
    # Truth: (r=50km, z_s in {180,200,220})

    score_rows = []
    alias_rows = []
    profile_rows = []
    rank_rows = []

    for z_true in full_zs:
        sig_true = sig_db.get((z_true, 50.0))
        if sig_true is None:
            continue

        for r_km in full_rs:
            for z_s in Z_S_LIST:
                sig_cand = sig_db.get((z_s, r_km))
                if sig_cand is None:
                    # try key grid
                    key_r = round(r_km * 2) / 2  # snap to 0.5
                    sig_cand = sig_db.get((z_s, key_r))
                if sig_cand is None or sig_true is None:
                    score_rows.append(
                        {
                            "z_true_m": z_true,
                            "r_km": r_km,
                            "z_s_m": z_s,
                            "J_tau": np.nan,
                            "n_common": 0,
                            "status": "ARRIVAL_STRUCTURE_MISMATCH",
                        }
                    )
                    continue

                # compare delay signatures (oracle: same number of arrivals assumed)
                n_common = min(len(sig_true), len(sig_cand))
                if n_common < 2:
                    score_rows.append(
                        {
                            "z_true_m": z_true,
                            "r_km": r_km,
                            "z_s_m": z_s,
                            "J_tau": np.nan,
                            "n_common": n_common,
                            "status": "ARRIVAL_STRUCTURE_MISMATCH",
                        }
                    )
                    continue

                J = rms(sig_true[:n_common], sig_cand[:n_common])
                score_rows.append(
                    {
                        "z_true_m": z_true,
                        "r_km": r_km,
                        "z_s_m": z_s,
                        "J_tau": J,
                        "n_common": n_common,
                        "status": "OK",
                    }
                )

        # rank and profile for this z_true
        sub = [r for r in score_rows if r["z_true_m"] == z_true and r["status"] == "OK"]
        if not sub:
            continue
        Js = np.array([r["J_tau"] for r in sub])
        rs = np.array([r["r_km"] for r in sub])
        zs_arr = np.array([r["z_s_m"] for r in sub])
        order = np.argsort(Js)
        # find true
        true_mask = (np.abs(rs - 50.0) < 0.1) & (np.abs(zs_arr - z_true) < 0.1)
        if true_mask.any():
            true_idx = int(np.where(true_mask)[0][int(np.argmin(Js[true_mask]))])
            true_rank = int(np.where(order == true_idx)[0][0]) + 1
            true_J = float(Js[true_idx])
        else:
            true_rank = -1
            true_J = np.nan

        best_false = None
        for k in order:
            if not (np.abs(rs[k] - 50.0) < 0.1 and np.abs(zs_arr[k] - z_true) < 0.1):
                best_false = {"J": float(Js[k]), "r_km": float(rs[k]), "z_s_m": float(zs_arr[k])}
                break

        rank_rows.append(
            {
                "z_true_m": z_true,
                "true_rank": true_rank,
                "true_J": true_J,
                "best_false_J": best_false["J"] if best_false else np.nan,
                "best_false_r_km": best_false["r_km"] if best_false else np.nan,
                "best_false_z_s_m": best_false["z_s_m"] if best_false else np.nan,
                "n_scored": len(sub),
            }
        )

        # J_tau*(r) = min_z J(r,z)
        for r_km in full_rs:
            mask_r = np.abs(rs - r_km) < 0.1
            if mask_r.any():
                J_star = float(Js[mask_r].min())
                z_star = float(zs_arr[mask_r][int(np.argmin(Js[mask_r]))])
            else:
                J_star = np.nan
                z_star = np.nan
            profile_rows.append(
                {"z_true_m": z_true, "r_km": r_km, "J_tau_star": J_star, "z_star_m": z_star}
            )

        # range alias audit
        for r_km in TRACK_R_KM:
            mask_r = np.abs(rs - r_km) < 0.1
            if not mask_r.any():
                alias_rows.append(
                    {
                        "z_true_m": z_true,
                        "r_track_km": r_km,
                        "present": False,
                        "J_tau_star": np.nan,
                        "best_z_star": np.nan,
                    }
                )
                continue
            J_r = Js[mask_r]
            k_best = int(np.argmin(J_r))
            alias_rows.append(
                {
                    "z_true_m": z_true,
                    "r_track_km": r_km,
                    "present": True,
                    "J_tau_star": float(J_r[k_best]),
                    "delta_J": float(J_r[k_best] - Js[true_idx]) if true_mask.any() else np.nan,
                    "best_z_star": float(zs_arr[mask_r][k_best]),
                    "n_at_r": int(mask_r.sum()),
                }
            )

    score_df = pd.DataFrame(score_rows)
    score_df.to_csv(OUT / "ORACLE_TDOA_RZ_SCORE.csv", index=False)
    pd.DataFrame(profile_rows).to_csv(OUT / "ORACLE_TDOA_RANGE_PROFILE.csv", index=False)
    alias_df = pd.DataFrame(alias_rows)
    alias_df.to_csv(OUT / "TDOA_RANGE_ALIAS_AUDIT.csv", index=False)
    rank_df = pd.DataFrame(rank_rows)
    rank_df.to_csv(OUT / "ORACLE_TDOA_SIGNATURES.csv", index=False)  # rank diagnostics

    # =========================================================
    # 3. Decision
    # =========================================================
    all_rank1 = bool((rank_df["true_rank"] == 1).all()) if len(rank_df) else False

    # Check if J_tau*(r) has unique minimum at r=50 for all z_true
    unique_min = True
    for z_true in full_zs:
        prof = [p for p in profile_rows if p["z_true_m"] == z_true and np.isfinite(p["J_tau_star"])]
        if not prof:
            unique_min = False
            continue
        Js_p = [p["J_tau_star"] for p in prof]
        rs_p = [p["r_km"] for p in prof]
        k_min = int(np.argmin(Js_p))
        if abs(rs_p[k_min] - 50.0) > 0.1:
            unique_min = False
        # check no near-degenerate at aliases
        jmin = Js_p[k_min]
        for r_alias in [45.0, 56.0, 58.0, 60.0]:
            for p in prof:
                if abs(p["r_km"] - r_alias) < 0.1 and np.isfinite(p["J_tau_star"]):
                    if p["J_tau_star"] < jmin + 1e-6:
                        unique_min = False

    # Check r-z alias persistence
    alias_persists = False
    for _, a in alias_df.iterrows():
        if a.get("present") and np.isfinite(a.get("delta_J", np.nan)):
            if a["delta_J"] < 1e-6 and abs(a["r_track_km"] - 50.0) > 0.1:
                alias_persists = True

    if all_rank1 and unique_min and not alias_persists:
        decision = "SINGLE_DEPTH_ORACLE_TDOA_RANGE_INFORMATION_CONFIRMED"
    elif alias_persists or not unique_min:
        if all_rank1:
            decision = "SINGLE_DEPTH_ORACLE_TDOA_RANGE_INFORMATION_PARTIAL"
        else:
            decision = "SINGLE_DEPTH_ORACLE_TDOA_RZ_DEGENERACY_PERSISTS"
    else:
        decision = "SINGLE_DEPTH_ORACLE_TDOA_RANGE_INFORMATION_PARTIAL"

    why = (
        f"rank1_all={all_rank1}, unique_min_at_50={unique_min}, alias_persists={alias_persists}; "
        f"ranks={list(rank_df['true_rank'].values) if len(rank_df) else 'NA'}; "
        f"frac_2plus_arrivals={frac_2plus:.2f}"
    )

    dec = {
        "stage": "RANGE-UB-1B",
        "decision": decision,
        "why": why,
        "fraction_2plus_arrivals": frac_2plus,
        "all_true_rank1": all_rank1,
        "unique_min_at_true_r": unique_min,
        "range_alias_persists": alias_persists,
        "oracle_label": "PERFECT_MULTIPATH_DETECTION_PERFECT_ASSOCIATION_ZERO_ERROR",
        "created_utc": NOW,
    }
    (OUT / "RANGE_UB_1B_DECISION.json").write_text(
        json.dumps(dec, ensure_ascii=False, indent=2), encoding="utf-8"
    )

    # report
    report = f"""# RANGE-UB-1B — E-STD 单深度 Oracle 多途时延距离上界

UTC: {NOW}

基线：`dda967d77c1870217be8973d279c6e0d4533dad1`

**不是复现徐嘉璘2024论文。** 不使用论文未恢复的公式。

## 问题

> 在 E-STD、单接收深度 z_r=200 m，若理想获得多途到达时延，时延差本身能否打破 r-z 简并？

## 判定

### `{decision}`

{why}

## Oracle 假设

- PERFECT_MULTIPATH_DETECTION
- PERFECT_ARRIVAL_ASSOCIATION
- ZERO_DELAY_MEASUREMENT_ERROR

## Arrival 完整性

检查点数 {n_total}，≥2 arrivals: {n_with_2} ({frac_2plus*100:.1f}%)

## 排序诊断

{rank_df.to_string(index=False)}

## 距离别名审计

{alias_df.to_string(index=False)}

## 未做

噪声扫描、RC2/RC3 融合、垂直阵、论文公式复现、平台转向、P5。
"""
    (OUT / "RANGE_UB_1B_REPORT.md").write_text(report, encoding="utf-8")
    (OUT / "GPT_SYNC.md").write_text(
        f"# RANGE-UB-1B\n\n**{decision}**\n\n{why}\n", encoding="utf-8"
    )

    print("decision", decision)
    print(why)
    print("frac_2plus", frac_2plus)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
