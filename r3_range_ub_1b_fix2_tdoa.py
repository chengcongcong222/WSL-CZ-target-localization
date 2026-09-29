#!/usr/bin/env python3
"""RANGE-UB-1B-FIX2: BELLHOP BOX correction + oracle revalidation.

Fix: BOX line must be STEP(m) ZBOX(m) RBOX(km), not STEP ZBOX_m RBOX_m.
Correct: 0.0 5001.0 70.0
Also: .prt integrity gate, old-vs-new arrival comparison, eigenray complexity.
"""
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
OUT = ROOT / "results" / "RANGE_UB" / "CZ_TDOA_SINGLE_DEPTH_ORACLE_FIX2"
OLD = ROOT / "results" / "RANGE_UB" / "CZ_TDOA_SINGLE_DEPTH_ORACLE_FIX"
ZGRID = ROOT / "results" / "R3_C2_Yang_SA_depth" / "R3_C2_1" / "_kraken_zgrid"
AT_BIN = ROOT / "tools" / "acoustics_toolbox" / "atWin10" / "at" / "bin"
WORK = OUT / "_bellhop"
OUT.mkdir(parents=True, exist_ok=True)
WORK.mkdir(parents=True, exist_ok=True)
NOW = datetime.now(timezone.utc).isoformat()

FREQ = 235.0
ZR = 200.0
H_DEPTH = 5000.0
STEP_M = 0.0
ZBOX_M = 5001.0
RBOX_KM = 70.0
Z_S_ALL = np.arange(150.0, 250.0 + 1e-9, 5.0)
R_ALL_KM = np.arange(45.0, 60.0 + 1e-9, 0.5)
Z_TRUE_LIST = [180.0, 200.0, 220.0]
TRACK_R_KM = [45.0, 50.0, 56.0, 58.0, 60.0]


def write_bellhop_env(path: Path, freq: float, zs: float, zr: float, r_list_km, tag: str):
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

    r_km = float(r_list_km[0]) if len(r_list_km) == 1 else float(max(r_list_km))
    out = []
    out.append(f"'{tag}'")
    out.append(f"{freq:.1f}")
    out.append("1")
    out.append("'CVN'")
    out.append(f"  51 0.0 {H_DEPTH:.1f}")
    out.extend(ssp_lines)
    out.append("'R' 0.0")
    out.append("1")
    out.append(f"{zs:.1f} /")
    out.append("1")
    out.append(f"{zr:.1f} /")
    out.append(f"{len(r_list_km)}")
    if len(r_list_km) == 1:
        out.append(f"{r_km:.1f} /")
    else:
        out.append(f"{min(r_list_km):.1f} {r_km:.1f} /")
    out.append("'A'")
    out.append("201")
    out.append("-85.0 85.0 /")
    # CORRECT: STEP(m) ZBOX(m) RBOX(km)
    out.append(f"{STEP_M:.1f} {ZBOX_M:.1f} {RBOX_KM:.1f}")
    path.write_text("\n".join(out) + "\n", encoding="utf-8")


def parse_arr_structured(arr_path: Path) -> dict:
    if not arr_path.exists():
        return {"arrivals": [], "error": "file_not_found"}
    text = arr_path.read_text(encoding="utf-8", errors="replace")
    lines = [ln.strip() for ln in text.splitlines() if ln.strip()]
    if not lines:
        return {"arrivals": [], "error": "empty"}
    hdr = lines[0].split()
    Nsd = int(float(hdr[1]))
    Nrd = int(float(hdr[2]))
    Nrr = int(float(hdr[3]))
    idx = 1
    s_depth = [float(x) for x in lines[idx].split()[:Nsd]]
    idx += 1
    r_depth = [float(x) for x in lines[idx].split()[:Nrd]]
    idx += 1
    r_range = [float(x) for x in lines[idx].split()[:Nrr]]
    idx += 1
    all_arrivals = []
    for isd in range(Nsd):
        if idx >= len(lines):
            break
        idx += 1  # Narrmx
        for ird in range(Nrd):
            for ir in range(Nrr):
                if idx >= len(lines):
                    break
                narr = int(float(lines[idx]))
                idx += 1
                for k in range(narr):
                    if idx >= len(lines):
                        break
                    parts = lines[idx].split()
                    idx += 1
                    if len(parts) >= 8:
                        all_arrivals.append(
                            {
                                "isd": isd, "ird": ird, "ir": ir,
                                "zs": s_depth[isd] if isd < len(s_depth) else None,
                                "zr": r_depth[ird] if ird < len(r_depth) else None,
                                "r_km": r_range[ir] / 1000.0 if ir < len(r_range) else None,
                                "amp": float(parts[0]),
                                "phase_deg": float(parts[1]),
                                "delay_s": float(parts[2]),
                                "delay_imag_s": float(parts[3]),
                                "src_angle_deg": float(parts[4]),
                                "rcv_angle_deg": float(parts[5]),
                                "n_top": int(float(parts[6])),
                                "n_bot": int(float(parts[7])),
                            }
                        )
    return {"arrivals": all_arrivals, "n_arr": len(all_arrivals)}


def parse_prt(prt_path: Path) -> dict:
    if not prt_path.exists():
        return {"parsed_ZBOX_m": None, "parsed_RBOX_km": None, "step_m": None, "n_insufficient": -1, "n_other_warn": -1}
    text = prt_path.read_text(encoding="utf-8", errors="replace")
    zbox = None
    rbox = None
    step = None
    m = re.search(r"Maximum ray depth.*?=\s*([\d.]+)", text)
    if m:
        zbox = float(m.group(1))
    m = re.search(r"Maximum ray range.*?=\s*([\d.]+)", text)
    if m:
        rbox_m = float(m.group(1))
        rbox = rbox_m / 1000.0 if rbox_m > 100 else rbox_m
    m = re.search(r"Step length.*?=\s*([\d.]+)", text)
    if m:
        step = float(m.group(1))
    n_insuf = text.count("Insufficient storage")
    n_other = text.count("WARNING") - n_insuf if "WARNING" in text else 0
    return {"parsed_ZBOX_m": zbox, "parsed_RBOX_km": rbox, "step_m": step, "n_insufficient": n_insuf, "n_other_warn": max(0, n_other)}


def path_key(a):
    return (a["n_top"], a["n_bot"], int(np.sign(a["src_angle_deg"])), int(np.sign(a["rcv_angle_deg"])))


def associate_paths(truth_arrs, cand_arrs):
    from collections import defaultdict
    t_by, c_by = defaultdict(list), defaultdict(list)
    for i, a in enumerate(truth_arrs):
        t_by[path_key(a)].append(i)
    for j, a in enumerate(cand_arrs):
        c_by[path_key(a)].append(j)
    matched, used_t, used_c = [], set(), set()
    for key in set(t_by) | set(c_by):
        tl = [i for i in t_by.get(key, []) if i not in used_t]
        cl = [j for j in c_by.get(key, []) if j not in used_c]
        if not tl or not cl:
            continue
        pairs = sorted((abs(truth_arrs[i]["src_angle_deg"] - cand_arrs[j]["src_angle_deg"]), i, j) for i in tl for j in cl)
        for da, i, j in pairs:
            if i not in used_t and j not in used_c:
                matched.append((i, j, da))
                used_t.add(i)
                used_c.add(j)
    return matched, [i for i in range(len(truth_arrs)) if i not in used_t], [j for j in range(len(cand_arrs)) if j not in used_c]


def pairwise_tdoa_ms(delays):
    d = np.asarray(delays, float)
    return np.array([(d[j] - d[i]) * 1000.0 for i in range(len(d)) for j in range(i + 1, len(d))])


def rms(a, b):
    return float(np.sqrt(np.mean((np.asarray(a) - np.asarray(b)) ** 2)))


def main() -> int:
    config = {
        "stage": "RANGE-UB-1B-FIX2",
        "created_utc": NOW,
        "baseline_commit": "6b203f06b1b12381187491c6a314585c6640834c",
        "fix": "BOX: STEP(m) ZBOX(m) RBOX(km) = 0.0 5001.0 70.0",
        "oracle_label": "ALL_EIGENRAY_ZERO_ERROR_SINGLE_DEPTH_TDOA",
        "n_candidates": 651,
    }
    (OUT / "RANGE_UB_1B_FIX2_CONFIG.json").write_text(json.dumps(config, ensure_ascii=False, indent=2), encoding="utf-8")

    # =========================================================
    # 1. BELLHOP runs with corrected BOX
    # =========================================================
    arrival_db = {}
    box_rows = []
    for zs in Z_S_ALL:
        tag = f"f2_zs{int(zs)}"
        env_path = WORK / f"{tag}.env"
        write_bellhop_env(env_path, FREQ, float(zs), ZR, list(R_ALL_KM), tag)
        rc = subprocess.run([str(AT_BIN / "bellhop.exe"), tag], cwd=str(WORK), capture_output=True, text=True, timeout=180)
        arr_path = WORK / f"{tag}.arr"
        prt_path = WORK / f"{tag}.prt"
        parsed = parse_arr_structured(arr_path)
        prt = parse_prt(prt_path)
        by_r = {}
        for a in parsed.get("arrivals", []):
            rk = round(a.get("r_km", 0) * 2) / 2
            by_r.setdefault(rk, []).append(a)
        for rk in R_ALL_KM:
            arrival_db[(float(zs), float(rk))] = by_r.get(round(rk * 2) / 2, [])
        box_rows.append({
            "zs": zs, "return_code": rc.returncode,
            "parsed_ZBOX_m": prt["parsed_ZBOX_m"], "parsed_RBOX_km": prt["parsed_RBOX_km"],
            "step_m": prt["step_m"], "n_insufficient_storage": prt["n_insufficient"],
            "n_other_warning": prt["n_other_warn"], "n_arrivals_total": parsed.get("n_arr", 0),
        })

    box_df = pd.DataFrame(box_rows)
    box_df.to_csv(OUT / "BELLHOP_BOX_AND_WARNING_AUDIT.csv", index=False)
    ok_box = bool((box_df["n_insufficient_storage"] == 0).all() and (box_df["return_code"] == 0).all())
    ok_zbox = bool(box_df["parsed_ZBOX_m"].apply(lambda x: x is not None and abs(x - 5001.0) < 10).all()) if box_df["parsed_ZBOX_m"].notna().any() else False

    (OUT / "BELLHOP_NUMERICAL_INTEGRITY.md").write_text(
        f"""# BELLHOP_NUMERICAL_INTEGRITY

UTC: {NOW}

## BOX correction

Before: `0.0 r_max_m H_DEPTH` → Box%z=70000m, Box%r=5000km (WRONG)
After: `0.0 {ZBOX_M} {RBOX_KM}` → STEP=auto, ZBOX={ZBOX_M}m, RBOX={RBOX_KM}km

## .prt audit (21 runs)

- All return_code=0: {(box_df['return_code']==0).all()}
- Zero insufficient_storage: {(box_df['n_insufficient_storage']==0).all()}
- ZBOX≈5001m: {ok_zbox}

Status: {'PASS' if ok_box and ok_zbox else 'FAIL'}
""",
        encoding="utf-8",
    )
    if not (ok_box and ok_zbox):
        (OUT / "RANGE_UB_1B_FIX2_DECISION.json").write_text(
            json.dumps({"stage": "RANGE-UB-1B-FIX2", "decision": "RANGE_UB_1B_FIX2_BLOCKED_BY_BELLHOP_NUMERICAL_INTEGRITY", "created_utc": NOW}, indent=2),
            encoding="utf-8",
        )
        print("BLOCKED_BY_NUMERICAL_INTEGRITY")
        return 1

    # =========================================================
    # 2. Parser sample + old-vs-new comparison
    # =========================================================
    sample = arrival_db.get((200.0, 50.0), [])
    if sample:
        pd.DataFrame(sample).to_csv(OUT / "PARSED_ARRIVAL_SAMPLE_50KM_ZS200_FIX2.csv", index=False)

    # old vs new
    old_map = {}
    old_arr_dir = OLD / "_bellhop"
    for zs in Z_S_ALL:
        op = old_arr_dir / f"fix_zs{int(zs)}.arr"
        if op.exists():
            old_parsed = parse_arr_structured(op)
            for a in old_parsed.get("arrivals", []):
                rk = round(a.get("r_km", 0) * 2) / 2
                old_map.setdefault((float(zs), rk), []).append(a)

    cmp_rows = []
    for r_km in TRACK_R_KM:
        for zs in Z_TRUE_LIST:
            old_a = old_map.get((zs, r_km), [])
            new_a = arrival_db.get((zs, r_km), [])
            def stats(arrs):
                if not arrs:
                    return {"n": 0, "earliest": np.nan, "latest": np.nan, "spread": np.nan, "max_top": np.nan, "max_bot": np.nan}
                ds = [a["delay_s"] for a in arrs]
                return {"n": len(arrs), "earliest": min(ds), "latest": max(ds), "spread": max(ds)-min(ds),
                        "max_top": max(a["n_top"] for a in arrs), "max_bot": max(a["n_bot"] for a in arrs)}
            so, sn = stats(old_a), stats(new_a)
            cmp_rows.append({"r_km": r_km, "zs_m": zs, **{f"old_{k}": v for k, v in so.items()}, **{f"new_{k}": v for k, v in sn.items()}})
    cmp_df = pd.DataFrame(cmp_rows)
    cmp_df.to_csv(OUT / "OLD_VS_FIXED2_ARRIVAL_COMPARISON.csv", index=False)
    affected = bool((cmp_df["old_n"] != cmp_df["new_n"]).any())

    # =========================================================
    # 3. Eigenray complexity diagnostic (651 points)
    # =========================================================
    cx_rows = []
    for zs in Z_S_ALL:
        for rk in R_ALL_KM:
            arrs = arrival_db.get((float(zs), float(rk)), [])
            if arrs:
                amps = [a["amp"] for a in arrs]
                ds = [a["delay_s"] for a in arrs]
                cx_rows.append({
                    "zs_m": zs, "r_km": rk, "n_arrivals": len(arrs),
                    "max_n_top": max(a["n_top"] for a in arrs), "max_n_bot": max(a["n_bot"] for a in arrs),
                    "delay_spread_s": max(ds) - min(ds), "amp_max": max(amps), "amp_min": min(amps),
                    "amp_ratio": max(amps) / min(amps) if min(amps) > 0 else np.inf,
                })
            else:
                cx_rows.append({"zs_m": zs, "r_km": rk, "n_arrivals": 0})
    pd.DataFrame(cx_rows).to_csv(OUT / "EIGENRAY_COMPLEXITY_DIAGNOSTIC.csv", index=False)

    # =========================================================
    # 4. TDOA scoring
    # =========================================================
    score_rows, alias_rows, profile_rows, rank_rows = [], [], [], []
    for z_true in Z_TRUE_LIST:
        truth_arrs = arrival_db.get((z_true, 50.0), [])
        if len(truth_arrs) < 3:
            continue
        for z_s in Z_S_ALL:
            for r_km in R_ALL_KM:
                cand = arrival_db.get((float(z_s), float(r_km)), [])
                matched, _, _ = associate_paths(truth_arrs, cand)
                n_m = len(matched)
                if n_m >= 3:
                    t_d = [truth_arrs[i]["delay_s"] for i, j, _ in matched]
                    c_d = [cand[j]["delay_s"] for i, j, _ in matched]
                    J = rms(pairwise_tdoa_ms(t_d), pairwise_tdoa_ms(c_d))
                    status = "OK"
                else:
                    J, status = np.nan, "INSUFFICIENT_COMMON_PATHS"
                score_rows.append({"z_true_m": z_true, "r_km": r_km, "z_s_m": z_s, "J_tau_ms": J, "n_matched": n_m, "status": status})

        sub = [r for r in score_rows if r["z_true_m"] == z_true and r["status"] == "OK"]
        if not sub:
            continue
        Js = np.array([r["J_tau_ms"] for r in sub])
        rs = np.array([r["r_km"] for r in sub])
        zs_arr = np.array([r["z_s_m"] for r in sub])
        order = np.argsort(Js)
        true_mask = (np.abs(rs - 50.0) < 0.1) & (np.abs(zs_arr - z_true) < 0.1)
        if true_mask.any():
            ti = int(np.where(true_mask)[0][int(np.argmin(Js[true_mask]))])
            true_rank = int(np.where(order == ti)[0][0]) + 1
            true_J = float(Js[ti])
        else:
            true_rank, true_J = -1, np.nan
        for rk in R_ALL_KM:
            m = np.abs(rs - rk) < 0.1
            if m.any():
                J_star = float(Js[m].min())
                z_star = float(zs_arr[m][int(np.argmin(Js[m]))])
            else:
                J_star = z_star = np.nan
            profile_rows.append({"z_true_m": z_true, "r_km": rk, "J_tau_star_ms": J_star, "z_star_m": z_star})
        for rk in TRACK_R_KM:
            m = np.abs(rs - rk) < 0.1
            if not m.any():
                alias_rows.append({"z_true_m": z_true, "r_track_km": rk, "present": False, "J_tau_star_ms": np.nan})
                continue
            Jr = Js[m]
            kb = int(np.argmin(Jr))
            alias_rows.append({"z_true_m": z_true, "r_track_km": rk, "present": True,
                "J_tau_star_ms": float(Jr[kb]), "delta_J_ms": float(Jr[kb] - Js[ti]) if true_mask.any() else np.nan,
                "best_z_star": float(zs_arr[m][kb]), "machine_degenerate": bool(Jr[kb] < 1e-6)})
        rank_rows.append({"z_true_m": z_true, "true_rank": true_rank, "true_J_ms": true_J, "n_valid_scored": len(sub)})

    score_df = pd.DataFrame(score_rows)
    score_df.to_csv(OUT / "ORACLE_TDOA_RZ_SCORE_FIX2.csv", index=False)
    pd.DataFrame(profile_rows).to_csv(OUT / "ORACLE_TDOA_RANGE_PROFILE_FIX2.csv", index=False)
    alias_df = pd.DataFrame(alias_rows)
    alias_df.to_csv(OUT / "TDOA_RANGE_ALIAS_AUDIT_FIX2.csv", index=False)
    rank_df = pd.DataFrame(rank_rows)
    rank_df.to_csv(OUT / "ORACLE_TDOA_SIGNATURES_FIX2.csv", index=False)

    # =========================================================
    # 5. Decision (narrowed label)
    # =========================================================
    range_ok = True
    for z_true in Z_TRUE_LIST:
        prof = [p for p in profile_rows if p["z_true_m"] == z_true and np.isfinite(p.get("J_tau_star_ms", np.nan))]
        if not prof:
            range_ok = False
            continue
        k = int(np.argmin([p["J_tau_star_ms"] for p in prof]))
        if abs(prof[k]["r_km"] - 50.0) > 0.1:
            range_ok = False
        for p in prof:
            if abs(p["r_km"] - 50.0) > 0.1 and p["J_tau_star_ms"] < 1e-6:
                range_ok = False

    all_rank1 = bool((rank_df["true_rank"] == 1).all()) if len(rank_df) else False
    n_valid_ok = bool((rank_df["n_valid_scored"] == 651).all()) if len(rank_df) else False

    if range_ok and all_rank1 and n_valid_ok:
        decision = "ALL_EIGENRAY_ZERO_ERROR_SINGLE_DEPTH_TDOA_RANGE_INFORMATION_CONFIRMED"
    else:
        decision = "ALL_EIGENRAY_ZERO_ERROR_TDOA_RANGE_INFORMATION_NOT_CONFIRMED"

    why = f"range_ok={range_ok}, all_rank1={all_rank1}, n_valid=651x3={n_valid_ok}, box_affected_old={affected}"

    dec = {
        "stage": "RANGE-UB-1B-FIX2", "decision": decision, "why": why,
        "oracle_layer": "ALL_EIGENRAY_PERFECT_DETECTION_ORACLE (physical upper bound)",
        "not_equivalent_to_paper_TDOA": True,
        "old_arrival_set_affected": "OLD_1B_FIX_ARRIVAL_SET_AFFECTED_BY_BOX_CONFIG" if affected else "UNCHANGED",
        "created_utc": NOW,
    }
    (OUT / "RANGE_UB_1B_FIX2_DECISION.json").write_text(json.dumps(dec, ensure_ascii=False, indent=2), encoding="utf-8")

    report = f"""# RANGE-UB-1B-FIX2 — BELLHOP BOX correction + oracle revalidation

UTC: {NOW}

基线：`6b203f06b1b12381187491c6a314585c6640834c`

## BOX 修复

| | 旧 | 新 |
|---|---|---|
| 写法 | `0.0 r_max_m H_DEPTH` | `0.0 {ZBOX_M} {RBOX_KM}` |
| Box%z | 70000 m（错） | {ZBOX_M} m |
| Box%r | 5000 km（错） | {RBOX_KM} km |
| insufficient storage | 大量 | **0** |

## .prt 完整性

21 runs 全部 return_code=0，零 insufficient storage。

## 判定

### `{decision}`

{why}

**这是所有模型 eigenray 完美可知的物理信息上限**，不等价于徐嘉璘 2024 论文的可提取 TDOA。

## 排序诊断

{rank_df.to_string(index=False)}

## 距离别名

{alias_df.to_string(index=False)}

## 新旧 arrival 对比

{'OLD_1B_FIX_ARRIVAL_SET_AFFECTED_BY_BOX_CONFIG' if affected else 'arrival set unchanged'}

{cmp_df.head(10).to_string(index=False)}

## 未做

噪声、VLA、RC2/RC3、paper-like cluster 简化、amplitude cutoff、P5。
"""
    (OUT / "RANGE_UB_1B_FIX2_REPORT.md").write_text(report, encoding="utf-8")
    (OUT / "GPT_SYNC.md").write_text(f"# RANGE-UB-1B-FIX2\n\n**{decision}**\n\n{why}\n", encoding="utf-8")

    print("box_ok", ok_box, "zbox_ok", ok_zbox)
    print("decision", decision)
    print(why)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
