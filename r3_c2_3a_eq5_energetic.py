#!/usr/bin/env python3
"""C2.3A supplement: EQ5 depth signature with energetic modes (independent of unique-ID)."""
from __future__ import annotations

import json
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parent
ZGRID = ROOT / "results" / "R3_C2_Yang_SA_depth" / "R3_C2_1" / "_kraken_zgrid"
OUT = ROOT / "results" / "R3_C2_Yang_SA_depth" / "R3_C2_3A"
Z_GRID = np.round(np.arange(150.0, 250.0 + 1e-9, 2.0), 1)
FREQS = [201.0, 235.0, 283.0, 338.0]
Z_TRUE = [180.0, 190.0, 200.0, 210.0, 220.0]
ZR = 200.0


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


def b_m(k_re, alpha, phi_s, r0, dR):
    out = np.zeros(len(k_re), complex)
    for m in range(len(k_re)):
        a = alpha[m]
        x = a * dR / 2
        sox = 1 + x * x / 6 if abs(x) < 1e-8 else np.sinh(x) / x
        out[m] = (np.exp(-a * r0) / k_re[m]) * dR * sox * complex(phi_s[m])
    return out


rows = []
for f in FREQS:
    M = parse_mod(ZGRID / f"zgrid_f{int(f)}.mod")
    k_re = M["k"].real
    alpha = -M["k"].imag
    izr = int(np.argmin(np.abs(M["depths"] - ZR)))
    pr = M["phi"][izr]
    max_phi = float(np.max(np.abs(M["phi"])))
    for zs in Z_TRUE:
        izs = int(np.argmin(np.abs(M["depths"] - zs)))
        ps = M["phi"][izs]
        bm = b_m(k_re, alpha, ps, 51200.0, 2400.0)
        g_ideal = bm * pr
        rel = np.abs(g_ideal) / max(np.abs(g_ideal).max(), 1e-30)
        midx = list(np.where(rel >= 0.03)[0])
        if not midx:
            midx = list(np.argsort(rel)[::-1][:20])
        gs = g_ideal[midx]
        for d in (0.05, 0.10, 0.20):
            d_abs = d * max_phi
            s = np.zeros(M["phi"].shape[0], complex)
            for j, m in enumerate(midx):
                w = pr[m] / (pr[m] ** 2 + d_abs**2)
                s += M["phi"][:, m] * (gs[j] * w)
            D = np.abs(s) ** 2
            Dn = D / max(np.sum(D), 1e-30)
            zhat = float(Z_GRID[int(np.argmax(Dn))])
            Dt = float(Dn[int(np.argmin(np.abs(Z_GRID - zs)))])

            def Dat(z):
                return float(Dn[int(np.argmin(np.abs(Z_GRID - z)))])

            m10 = [Dt / max(Dat(zs + s10), 1e-30) for s10 in (-10, 10) if 150 <= zs + s10 <= 250]
            m20 = [Dt / max(Dat(zs + s20), 1e-30) for s20 in (-20, 20) if 150 <= zs + s20 <= 250]
            rows.append(
                {
                    "f": f,
                    "zs": zs,
                    "Delta_ratio": d,
                    "n_energetic_modes": len(midx),
                    "z_hat": zhat,
                    "abs_error": abs(zhat - zs),
                    "D_true_over_Dmax": Dt / max(Dn.max(), 1e-30),
                    "margin_10m": float(np.min(m10)) if m10 else np.nan,
                    "margin_20m": float(np.min(m20)) if m20 else np.nan,
                }
            )

df = pd.DataFrame(rows)
df.to_csv(OUT / "estd_depth_ambiguity_eq5_energetic.csv", index=False)

good = df[(df["Delta_ratio"] == 0.10) & (df["margin_10m"] > 1) & (df["D_true_over_Dmax"] >= 0.5)]
print("EQ5 energetic good @0.1", len(good), "/20")
print(df[df["Delta_ratio"] == 0.10][["f", "zs", "z_hat", "margin_10m", "D_true_over_Dmax"]])

# update decision
dec_path = OUT / "R3_C2_3A_DECISION.json"
dec = json.loads(dec_path.read_text(encoding="utf-8"))
n_uniq_all = 0
mp = pd.read_csv(OUT / "estd_peak_mode_mapping.csv")
n_uniq_all = int((mp["status"] == "UNIQUE_ORACLE_MODE").sum())
n_group = int((mp["status"] == "UNRESOLVED_MODE_GROUP").sum())

if n_uniq_all == 0 and n_group > 0 and len(good) > 0:
    decision = "C2_3A_FOURIER_MODE_IDENTITY_LIMITED"
    why = (
        f"E-STD L=2.4km 下 0 个 UNIQUE_ORACLE_MODE（{n_group} 个 UNRESOLVED_MODE_GROUP，每峰 9–22 模态）；"
        f"n_rayleigh_unique(3%)=0。EQ5 理想层用 energetic modes 在 Δ=0.1 可分 "
        f"{len(good)}/20（margin_10m>1），深度签名存在。"
        "瓶颈是普通 Fourier 孔径下密集模态身份无法唯一确定——只限制 Fourier mode-ID，不等价于 AR/高分辨失败。"
    )
    labels = sorted(set(dec.get("labels", [])) | {"C2_3A_FOURIER_MODE_IDENTITY_LIMITED"})
    if "C2_3A_DEPTH_SIGNATURE_WEAK" in labels:
        labels.remove("C2_3A_DEPTH_SIGNATURE_WEAK")
elif len(good) == 0:
    decision = "C2_3A_DEPTH_SIGNATURE_WEAK"
    why = "EQ5 energetic layer also lacks candidate discrimination"
    labels = dec.get("labels", [])
else:
    decision = dec["rc3_c2_3a_decision"]
    why = dec.get("why", "")
    labels = dec.get("labels", [])

dec["rc3_c2_3a_decision"] = decision
dec["why"] = why
dec["labels"] = labels
dec["eq5_energetic_good_at_0.1"] = int(len(good))
dec["n_unique_oracle_mode_total"] = n_uniq_all
dec["n_unresolved_mode_group_total"] = n_group
dec_path.write_text(json.dumps(dec, ensure_ascii=False, indent=2), encoding="utf-8")
print("DECISION", decision)
print("unique", n_uniq_all, "groups", n_group)
