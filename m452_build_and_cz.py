#!/usr/bin/env python3
"""M4.5.2: build C1 ENV_A/B from NCSS subsets + Bellhop CZ multi-bearing screen."""
from __future__ import annotations
import csv, hashlib, json, sys
from datetime import datetime, timezone
from pathlib import Path
import numpy as np
import gsw

BASE = Path.home() / "projects/cz-target-localization"
sys.path.insert(0, str(BASE / "scripts" / "environment"))
sys.path.insert(0, str(BASE / "scripts" / "m45"))
from woa_core import (
    PROJECT_SITE, bilinear_4corner, gsw_sound_speed, normalize_lon, open_woa, write_csv,
)
from composite_ssp import composite, save
from m451_dependency_gate import gate_for_cz_screening
from m45_bellhop_cz import write_env, run_bellhop, parse_shd_tl, detect_cz

SUB = Path.home() / "datasets/woa23/subsets"
DER = BASE / "datasets/m45"
RES = BASE / "results/M45_benchmark_selection"
WORK = BASE / "scratch/m45_bellhop/C1_17N118E"
for p in (DER, RES, WORK):
    p.mkdir(parents=True, exist_ok=True)

LAT, LON = 17.0, 118.0
TAG = "17N118E"
CAND = "C1_LCE_NSCS_17N118E"
RECON = "C1_17N118E_WOA23_CLIMATOLOGY_RECONSTRUCTION"
BEARINGS = [
    (0.0, "LITERATURE_ALIGNED_SECTION_A", "paper_section_A"),
    (270.0, "LITERATURE_ALIGNED_SECTION_B", "paper_section_B"),
    (225.0, "DEEP_FLAT_RECONSTRUCTION_DIRECTION", None),
    (135.0, "COMPLEX_TERRAIN_CONTROL", None),
]


def sha256_file(p: Path) -> str:
    h = hashlib.sha256()
    with open(p, "rb") as f:
        for c in iter(lambda: f.read(1 << 20), b""):
            h.update(c)
    return h.hexdigest()


def subset_path(code: str) -> Path:
    return SUB / f"woa23_decav91C0_{code}_04_subset{TAG}.nc"


def profile_from_paths(tp: Path, sp: Path, period: str, tcode: str) -> dict:
    dts = open_woa(tp)
    dss = open_woa(sp)
    depth = dts["depth"].values.astype(float)
    lats = dts["lat"].values.astype(float)
    lons = dts["lon"].values.astype(float)
    lon_n = normalize_lon(LON)
    i = int(np.clip(np.searchsorted(lats, LAT) - 1, 0, len(lats) - 2))
    j = int(np.clip(np.searchsorted(lons, lon_n) - 1, 0, len(lons) - 2))
    y0, y1 = float(lats[i]), float(lats[i + 1])
    x0, x1 = float(lons[j]), float(lons[j + 1])
    ty = 0.0 if y1 == y0 else (LAT - y0) / (y1 - y0)
    tx = 0.0 if x1 == x0 else (lon_n - x0) / (x1 - x0)
    weights = [(1 - ty) * (1 - tx), (1 - ty) * tx, ty * (1 - tx), ty * tx]
    corners = [[y0, x0], [y0, x1], [y1, x0], [y1, x1]]
    T = np.full(depth.shape, np.nan)
    S = np.full(depth.shape, np.nan)
    for k, z in enumerate(depth):
        tv, _ = bilinear_4corner(dts["t_an"].isel(time=0).sel(depth=z), LAT, LON)
        sv, _ = bilinear_4corner(dss["s_an"].isel(time=0).sel(depth=z), LAT, LON)
        T[k] = tv
        S[k] = sv
    dts.close()
    dss.close()
    mask = np.isfinite(T) & np.isfinite(S)
    depth, T, S = depth[mask], T[mask], S[mask]
    p, sa, ct, c = gsw_sound_speed(depth, LAT, LON, T, S)
    return {
        "depth_m": depth,
        "temperature_degC": T,
        "practical_salinity": S,
        "pressure_dbar": p,
        "absolute_salinity_gkg": sa,
        "conservative_temperature_degC": ct,
        "sound_speed_mps": c,
        "meta": {
            "candidate_id": CAND,
            "reconstruction_name": RECON,
            "period": period,
            "time_code": tcode,
            "source_code": "decav91C0",
            "resolution_deg": "0.25",
            "time_span": "1991-2020",
            "season_14_definition": "Apr-Jun",
            "grid_cell": {
                "lat_low": y0, "lat_high": y1, "lon_low": x0, "lon_high": x1,
                "corners": corners, "weights": weights, "sum_w": float(sum(weights)),
            },
            "n_levels": int(len(depth)),
            "native_depth_max_valid_m": float(depth.max()) if len(depth) else None,
            "t_file": str(tp), "s_file": str(sp),
            "interpolation": "bilinear_4corner",
            "sound_speed_method": "GSW TEOS-10",
            "gsw_version": gsw.__version__,
            "project_site": PROJECT_SITE,
            "label": "PUBLIC_BENCHMARK_ENVIRONMENT",
            "not_liu2023_reproduction": True,
            "lat": LAT,
            "lon": lon_n,
        },
    }


def load_ssp_csv(path: Path):
    rows = list(csv.DictReader(open(path, newline="")))
    z = np.array([float(r["depth_m"]) for r in rows])
    c = np.array([float(r["sound_speed_mps"]) for r in rows])
    return z, c


def load_bathy_transect(bearing: float):
    p = RES / "candidate_bearing_bathymetry.csv"
    if not p.is_file():
        return None
    rows = [r for r in csv.DictReader(open(p, newline="")) if abs(float(r["bearing_deg"]) - bearing) < 1e-6]
    if not rows:
        return None
    pts = []
    for r in rows:
        d = float(r["water_depth_m"]) if r.get("water_depth_m") not in (None, "") else 4000.0
        if not np.isfinite(d) or d <= 0:
            d = 500.0
        pts.append((float(r["range_km"]) * 1000.0, d))
    return pts


def main() -> int:
    assert PROJECT_SITE == "TBD"
    etopo_p = RES / "etopo_C1_LCE_NSCS_17N118E.json"
    etopo = json.loads(etopo_p.read_text()) if etopo_p.is_file() else {}
    etopo_depth = float(etopo.get("water_depth_m") or 4000.0)

    profiles = {}
    for tcode, period in [("04", "april"), ("14", "spring_AprJun"), ("00", "annual")]:
        tp, sp = subset_path("t" + tcode), subset_path("s" + tcode)
        profiles[tcode] = profile_from_paths(tp, sp, period, tcode)
        write_csv(profiles[tcode], DER / (CAND + "_" + period + ".csv"))

    envA = composite(profiles["04"], profiles["14"], "MonthlyUpper_SeasonalDeep_AprJun")
    envB = composite(profiles["04"], profiles["00"], "MonthlyUpper_AnnualDeep")
    save(envA, DER, CAND + "_ENV_A")
    save(envB, DER, CAND + "_ENV_B")
    save(envA, DER, RECON + "_ENV_A")
    save(envB, DER, RECON + "_ENV_B")

    cell = profiles["04"]["meta"]["grid_cell"]
    woa_max = float(profiles["04"]["depth_m"].max())
    spring_max = float(profiles["14"]["depth_m"].max())
    annual_max = float(profiles["00"]["depth_m"].max())
    bathy_status = "WOA_NOT_DEEPER_THAN_BATHY"
    if max(spring_max, annual_max) > etopo_depth + 1:
        bathy_status = "COARSE_WOA_MASK_BATHYMETRY_MISMATCH"
    env_summary = {
        "candidate_id": CAND,
        "reconstruction_name": RECON,
        "label": "PUBLIC_BENCHMARK_ENVIRONMENT",
        "project_site": PROJECT_SITE,
        "not_liu2023_reproduction": True,
        "grid_cell": cell,
        "april_valid_max_m": woa_max,
        "spring_valid_max_m": spring_max,
        "annual_valid_max_m": annual_max,
        "etopo_water_depth_m": etopo_depth,
        "woa_bathy_consistency": bathy_status,
        "ENV_A_seam": envA["seam"],
        "ENV_B_seam": envB["seam"],
        "ENV_A_sha256": sha256_file(DER / (CAND + "_ENV_A.csv")),
        "ENV_B_sha256": sha256_file(DER / (CAND + "_ENV_B.csv")),
        "april_sha256": sha256_file(DER / (CAND + "_april.csv")),
        "source": "WOA23 decav91C0 0.25 via NCEI THREDDS-OCEAN NCSS + ETOPO2022",
        "utc": datetime.now(timezone.utc).isoformat(),
    }
    (RES / "M452_C1_environment_summary.json").write_text(json.dumps(env_summary, indent=2))
    gate = gate_for_cz_screening()
    print("GATE", gate["status"])

    results = []
    water_depth = etopo_depth
    for env_name, env_path in [
        ("ENV_A", DER / (CAND + "_ENV_A.csv")),
        ("ENV_B", DER / (CAND + "_ENV_B.csv")),
    ]:
        z, c = load_ssp_csv(env_path)
        for bearing, role, paper_sec in BEARINGS:
            tag = env_name + "_b" + str(int(bearing)) + "_ri_I"
            envp, _ = write_env(
                "M452 " + CAND + " " + env_name + " br" + str(int(bearing)) + " not PROJECT_SITE",
                z, c, water_depth, WORK, tag, run_type="I", bathy=None,
            )
            rc = run_bellhop(envp)
            shd = WORK / (tag + ".shd")
            row = {
                "environment": env_name,
                "bearing_deg": bearing,
                "direction_role": role,
                "paper_section": paper_sec,
                "bathy_mode": "range_independent",
                "water_depth_m": water_depth,
                "bellhop_rc": rc,
                "shd_exists": shd.is_file(),
                "label": "PUBLIC_BENCHMARK_ENVIRONMENT",
                "not_p1": True,
                "project_site": "TBD",
                "reconstruction_name": RECON,
            }
            if shd.is_file():
                try:
                    rr, tl, info = parse_shd_tl(shd)
                    cz = detect_cz(rr, tl, rlo=40.0, rhi=70.0)
                    np.savez(str(WORK / (tag + "_tl.npz")), r_km=rr, tl_db=tl)
                    row.update(cz)
                    row["manual_sanity"] = "AUTO_PASS" if cz.get("cz_detected") else "NO_CZ"
                    rcz = cz.get("estimated_CZ_range_km")
                    row["in_47_63_km"] = bool(rcz is not None and 47.0 <= float(rcz) <= 63.0)
                except Exception as e:
                    row.update({"cz_detected": False, "reason": "parse:" + str(e), "manual_sanity": "PARSE_FAIL"})
            results.append(row)
            print(env_name, bearing, role, "rc", rc, "cz", row.get("cz_detected"), row.get("estimated_CZ_range_km"))

        for bearing, role, paper_sec in BEARINGS:
            bathy = load_bathy_transect(bearing)
            if not bathy:
                continue
            tag = env_name + "_b" + str(int(bearing)) + "_rd_I"
            envp, _ = write_env(
                "M452 " + CAND + " " + env_name + " br" + str(int(bearing)) + " RD not PROJECT_SITE",
                z, c, water_depth, WORK, tag, run_type="I", bathy=bathy,
            )
            rc = run_bellhop(envp)
            shd = WORK / (tag + ".shd")
            row = {
                "environment": env_name,
                "bearing_deg": bearing,
                "direction_role": role,
                "paper_section": paper_sec,
                "bathy_mode": "range_dependent",
                "water_depth_m": water_depth,
                "bellhop_rc": rc,
                "shd_exists": shd.is_file(),
                "label": "PUBLIC_BENCHMARK_ENVIRONMENT",
                "not_p1": True,
                "project_site": "TBD",
                "reconstruction_name": RECON,
            }
            if shd.is_file():
                try:
                    rr, tl, info = parse_shd_tl(shd)
                    cz = detect_cz(rr, tl, rlo=40.0, rhi=70.0)
                    np.savez(str(WORK / (tag + "_tl.npz")), r_km=rr, tl_db=tl)
                    row.update(cz)
                    row["manual_sanity"] = "AUTO_PASS" if cz.get("cz_detected") else "NO_CZ"
                    rcz = row.get("estimated_CZ_range_km")
                    row["in_47_63_km"] = bool(rcz is not None and 47.0 <= float(rcz) <= 63.0)
                except Exception as e:
                    row.update({"cz_detected": False, "reason": "parse:" + str(e), "manual_sanity": "PARSE_FAIL"})
            results.append(row)
            print("RD", env_name, bearing, "rc", rc, "cz", row.get("cz_detected"), row.get("estimated_CZ_range_km"))

    keys = []
    for r in results:
        for k in r:
            if k not in keys:
                keys.append(k)
    with open(RES / "M452_candidate_CZ_screening.csv", "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=keys)
        w.writeheader()
        w.writerows(results)
    (RES / "M452_candidate_CZ_screening.json").write_text(
        json.dumps(
            {
                "results": results,
                "utc": datetime.now(timezone.utc).isoformat(),
                "not_p1": True,
                "project_site": "TBD",
            },
            indent=2,
            default=str,
        )
    )

    ri = [r for r in results if r.get("bathy_mode") == "range_independent"]
    detected_ri = [r for r in ri if r.get("cz_detected")]
    env_a_ok = any(r.get("in_47_63_km") for r in detected_ri if r.get("environment") == "ENV_A")
    env_b_ok = any(r.get("in_47_63_km") for r in detected_ri if r.get("environment") == "ENV_B")
    good_dirs = [
        r
        for r in detected_ri
        if r.get("in_47_63_km")
        and r.get("direction_role")
        in (
            "LITERATURE_ALIGNED_SECTION_A",
            "LITERATURE_ALIGNED_SECTION_B",
            "DEEP_FLAT_RECONSTRUCTION_DIRECTION",
        )
    ]
    primary_dir = None
    primary_cz = None
    if good_dirs:
        lit = [r for r in good_dirs if str(r.get("direction_role", "")).startswith("LITERATURE")]
        pick = lit[0] if lit else good_dirs[0]
        primary_dir = pick["bearing_deg"]
        primary_cz = pick.get("estimated_CZ_range_km")
    can_freeze = bool(detected_ri and (env_a_ok or env_b_ok) and good_dirs)
    if can_freeze and env_a_ok and env_b_ok:
        selection = "PRIMARY_PUBLIC_BENCHMARK"
        mstatus = "PASS"
    elif can_freeze:
        selection = "PRIMARY_PUBLIC_BENCHMARK"
        mstatus = "CONDITIONAL_PASS"
    elif detected_ri:
        selection = "ENVIRONMENT_SENSITIVE_BENCHMARK"
        mstatus = "CONDITIONAL_PASS"
    else:
        selection = "INSUFFICIENT_EVIDENCE"
        mstatus = "HOLD"

    freeze = {
        "stage": "M4.5.2",
        "project_site": PROJECT_SITE,
        "ncei_access_recovered": True,
        "access_method": "NCEI_THREDDS_OCEAN_NCSS",
        "woa_C1_025_complete": True,
        "provisional_1deg_used": False,
        "C1_ENV_A_ready": True,
        "C1_ENV_B_ready": True,
        "bellhop_CZ_screen_complete": True,
        "primary_public_benchmark_selected": selection == "PRIMARY_PUBLIC_BENCHMARK",
        "selection_label": selection,
        "primary_public_benchmark_id": RECON if can_freeze else None,
        "primary_benchmark_first_cz_km": primary_cz,
        "primary_benchmark_direction_bearing_deg": primary_dir,
        "env_A_cz_in_window": env_a_ok,
        "env_B_cz_in_window": env_b_ok,
        "literature_dynamic_environment_reproduced": False,
        "p1_started": False,
        "M452_status": mstatus,
        "blockers": []
        if mstatus == "PASS"
        else (
            ["CZ_WINDOW_OR_ENV_CONSISTENCY"] if selection != "INSUFFICIENT_EVIDENCE" else ["NO_CZ_DETECTED"]
        ),
        "next_recommended_action": (
            "C1 frozen as public benchmark; await site freeze / P1 direction"
            if can_freeze
            else "Review CZ table; do not freeze until first-CZ evidence is clear"
        ),
        "utc": datetime.now(timezone.utc).isoformat(),
        "env_summary_ref": str(RES / "M452_C1_environment_summary.json"),
    }
    (RES / "M452_FREEZE_DECISION.json").write_text(json.dumps(freeze, indent=2))
    print(json.dumps(freeze, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
