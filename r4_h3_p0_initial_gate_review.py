"""Read-only supplemental review of the frozen pre-propagation stop. No new solver."""
from pathlib import Path
import json
import numpy as np
import r4_h3_p0 as E
import r4_h3_g0 as G
import r4_h3_g0_audit as C
import r4_e2_nine_pair_pilot_audit as B

def main():
    o=E.O
    reg=G.rows(o/"NOMINAL_REGRESSION_AND_CHAIN_CONTROLS.csv")
    idx={r["check"]:r for r in reg}
    replay=[r for r in reg if r["check"].startswith(("nominal_","prior_information_"))]
    assert len(replay)==126 and all(r["PASS"]=="True" for r in replay)
    checks=[];derivatives=[];zeros=[]
    oldpath=G.O/"NUMERICAL_STABILITY.csv"
    old=[r for r in G.rows(oldpath) if r["calibration"]=="C1b" and r["scene_id"].endswith("_M")]
    oldidx={(r["scene_id"],r["resource"],r["noise_model"]):r for r in old}
    for grid in E.GRIDS:
        mods=[B.mod(E.modpath(f,grid,0)) for f in E.F]
        for sc in E.scenes():
            for res in E.RES:
                p,j,h=C.cold(mods,sc["state"],res)
                j=j/p[...,None];h=h/p[...,None]
                # Actual free modal source coefficient absorption, reused nominal model only.
                pp,jp,hp,_,err=G.field(mods,sc["state"],res)
                zeros.append({"scene":sc["scene_id"],"resource":res,"mesh":grid,"offset_mps":0,
                    "control":"C1g_FREE_NOMINAL_SOURCE_DEPTH_ABSORPTION","residual":err,
                    "limit":1e-10,"PASS":err<=1e-10,"depth_information":0,"uncertainty":"INF",
                    "scope":"NOMINAL_ONLY;OFFSET_NOT_EXECUTED"})
                for noise in E.NOISE:
                    y,n=G.system(p,j,noise,"C1b");yh,nh=G.system(p,h,noise,"C1b")
                    z=C.qrremove(y,n);zh=C.qrremove(yh,nh)
                    error=np.linalg.norm(z[:,2]-zh[:,2])/np.linalg.norm(zh[:,2])
                    key=f"source_depth_steps_{sc['scene_id']}_{res}_{grid}_{noise}"
                    row=idx[key]
                    mismatch=abs(error-float(row["value"]))
                    checks.append({"check":key,"independent_QR_relative":float(error),
                        "saved_SVD_relative":float(row["value"]),"reconstruction_difference":mismatch,
                        "reconstruction_PASS":mismatch<1e-9,"frozen_threshold":float(row["limit"]),
                        "frozen_science_gate_PASS":bool(error<=float(row["limit"]))})
                    prior=oldidx[sc["scene_id"],res,noise]
                    derivatives.append({"scene":sc["scene_id"],"resource":res,"mesh":grid,"noise_model":noise,
                       "P0_depth_column_step_relative":float(error),"P0_frozen_depth_column_limit":.02,
                       "P0_depth_column_gate_PASS":bool(error<=.02),
                       "G0_whole_J_step_relative":float(prior["depth_step_relative"]),
                       "G0_joint_grid_step_effective_Iz_relative":float(prior["depth_information_relative"]),
                       "G0_registered_numeric_status":prior["status"],
                       "comparison":"DIFFERENT_NORMALIZATION;NO_RETROACTIVE_GATE_CHANGE"})
                    # C2 arbitrary real/imag per-element/frequency/time tangent spans entire data space.
                    raw,_=G.system(p,j,noise,"C0",True)
                    saturated=raw-np.eye(raw.shape[0])@raw
                    ce=float(np.linalg.norm(saturated))
                    zeros.append({"scene":sc["scene_id"],"resource":res,"mesh":grid,"offset_mps":0,
                       "control":"C2_NOMINAL_FREE_ALL_RESPONSES_"+noise,"residual":ce,"limit":0,
                       "PASS":ce==0,"depth_information":0,"uncertainty":"INF",
                       "scope":"NOMINAL_ONLY;OFFSET_NOT_EXECUTED"})
    assert len(checks)==36 and all(r["reconstruction_PASS"] for r in checks)
    assert sum(not r["frozen_science_gate_PASS"] for r in checks)==20
    assert all(r["PASS"] for r in zeros)
    G.write(o/"INITIAL_GATE_INDEPENDENT_QR_REVIEW.csv",checks)
    G.write(o/"G0_VS_P0_DEPTH_STEP_METRIC.csv",derivatives)
    G.write(o/"UNKNOWN_GAIN_ZERO_CONTROLS.csv",zeros)
    failed=[r for r in reg if r["PASS"]!="True"]
    d=E.read(o/"H3_P0_DECISION.json")
    d.update(nominal_replay_original_values_PASS=True,
        nominal_replay_checks=len(replay),initial_gate_checks=len(reg),
        initial_gate_failed_checks=[r["check"] for r in failed],
        initial_gate_failure_count=len(failed),
        source_depth_only_gate_checks=36,source_depth_only_gate_failures=20,
        max_source_depth_only_step_relative=max(float(r["value"]) for r in reg if r["check"].startswith("source_depth_steps_")),
        actual_stop_reason="NEW_PREREGISTERED_DEPTH_COLUMN_STEP_2PCT_GATE",
        SSP_provider_attempted=False,SSP_information="NOT_EVALUATED",
        original_G0_gate="UNCHANGED;WHOLE_JACOBIAN_AND_EFFECTIVE_INFORMATION_DIAGNOSTICS_DIFFERENT",
        independent_initial_gate_QR_review="PASS_RECONSTRUCTION;20_SCIENCE_GATE_FAILURES_RETAINED")
    E.dump(o/"H3_P0_DECISION.json",d)
    E.dump(o/"INITIAL_GATE_REVIEW_VALIDATION.json",{
        "PASS":True,"independent_depth_column_reconstruction_checks":36,
        "max_QR_SVD_difference":max(r["reconstruction_difference"] for r in checks),
        "numerical_science_gate_FAIL":20,"historical_numeric_artifact_path":str(oldpath).replace("\\","/"),
        "historical_numeric_artifact_sha256":E.sha(oldpath),
        "historical_G0_status_changed":False,"thresholds_changed":False,
        "new_KRAKEN":0,"new_random_draws":0,"original_nominal_replay_checks_PASS":126})
    print(json.dumps({"reconstruction_PASS":True,"depth_gate_FAILED":20,
                     "max_depth_difference":d["max_source_depth_only_step_relative"],"new_calls":0}))
if __name__=="__main__":main()
