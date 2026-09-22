#!/usr/bin/env python3
import json
from pathlib import Path
import pandas as pd

OUT = Path(__file__).resolve().parent / "results" / "R3_C1_depth_motion"
dec = json.loads((OUT / "R3_C1_2_FINAL_DECISION.json").read_text(encoding="utf-8"))
pair = pd.read_csv(OUT / "r3c12_qm_pairs_final.csv")
d12 = pair[(pair.group == "D") & (pair.T_s == 1200)]
d220 = float(d12[d12.z_j_m == 220].D_multi.iloc[0])
dec["rc3c1_2_final_decision"] = "C1_QM_MULTIFREQ_ONLY_CONDITIONAL"
dec["why"] = (
    "FINAL-INTEGRITY passed (real endfire B=|w^H p|^2 differs from unsteered; "
    "unified D=1-mean_f|corr_f|; dr=10m converged). "
    "10 m pairs (200/210) inseparable at all T (D~0.001–0.003 << 0.05). "
    f"20 m pair (200/220) only limited multifreq gain at T=1200s (D_D={d220:.3f} vs A~0.006) "
    "while FWHM stays 63–80 m. 10 m not reached; 20 m weak/conditional only. "
    "Not reliable 10–20 m depth discrimination under known-track ceiling."
)
dec["notes"]["T1200_D_200_220"] = d220
dec["notes"]["depth_resolution_bound"] = "10 m 未达到；20 m 仅多频且 T=1200 s 时弱可分（D≈0.34）"
dec["next_step"] = "Zhu route PERMANENTLY_CLOSED after this gate; Yang 2015 next (not this round)"
(OUT / "R3_C1_2_FINAL_DECISION.json").write_text(json.dumps(dec, ensure_ascii=False, indent=2, default=str), encoding="utf-8")

rep = OUT / "R3_C1_2_FINAL_REPORT.md"
t = rep.read_text(encoding="utf-8")
if "C1_QM_MULTIFREQ_ONLY_CONDITIONAL" not in t:
    t += "\n\n## 修正终判（含 T=1200 多频）\n\n**`C1_QM_MULTIFREQ_ONLY_CONDITIONAL`**\n\n" + dec["why"] + "\n\n下一步：" + dec["next_step"] + "\n\nZhu **PERMANENTLY_CLOSED**。\n"
    rep.write_text(t, encoding="utf-8")
(OUT / "R3_C1_2_FINAL_GPT_SYNC.md").write_text(
    "# R3-C1.2 FINAL\n\n**C1_QM_MULTIFREQ_ONLY_CONDITIONAL**\n\n" + dec["why"] + "\n\n",
    encoding="utf-8",
)
print("updated", dec["rc3c1_2_final_decision"], "D220_T1200", d220)
