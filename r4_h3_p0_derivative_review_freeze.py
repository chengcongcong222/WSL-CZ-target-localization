"""Create design at parent only; no saved-field analysis."""
import r4_h3_p0_derivative_review as R
E=R.E;G=R.G;O=R.O
def main():
    assert G.git("rev-parse","HEAD")==R.PARENT
    assert not O.exists()
    O.mkdir(parents=True)
    files=[E.ROOT/n for n in ["r4_h3_p0_derivative_review.py","r4_h3_p0_derivative_review_audit.py","r4_h3_p0_derivative_review_freeze.py","r4_h3_p0_derivative_review_report.py","r4_h3_p0.py","r4_h3_g0.py","r4_h3_g0_audit.py","r4_e2_nine_pair_pilot_audit.py","r4_e2_g0.py","r4_e2_paired_recovery.py"]]
    files += [E.O/n for n in ["DESIGN_FREEZE.json","H3_P0_DECISION.json","NOMINAL_REGRESSION_AND_CHAIN_CONTROLS.csv"]]
    files += [G.O/n for n in ["DESIGN_FREEZE.json","VERTICAL_INFORMATION_BY_SCENE.csv","NUMERICAL_STABILITY.csv"]]
    for g in E.GRIDS:
        for f in E.F:
            files += [E.modpath(f,g,0),E.modpath(f,g,0).with_suffix(".env")]
            files += [E.modpath(f,g,d).with_suffix(".env") for d in [-1,1]]
        files += [G.O/f"INPUT_{sc['scene_id']}_{res}_n{g}.npz" for sc in E.scenes() for res in E.RES]
    bindings=[{"path":str(p.resolve().relative_to(E.ROOT)).replace("\\","/"),"sha256":E.sha(p)} for p in files]
    E.dump(O/"DESIGN_FREEZE.json",{
        "stage":"R4_H3_P0_DEPTH_DERIVATIVE_AND_MODAL_WINDOW_REVIEW","parent_SHA":R.PARENT,
        "research_lead_parent_acceptance":"H3_P0_NUMERICAL_OR_FORMULA_INCOMPLETE",
        "accepted_failure_scope":"NEW_PREREGISTERED_DEPTH_COLUMN_STEP_2PCT_GATE",
        "scope":"READ_ONLY_SAVED_DATA_REVIEW;NO_P0_RERUN;NO_NEW_PROPAGATION",
        "source_depths_m":R.DEPTH,"methods":R.METHODS,"coefficients_per_m":R.COEF.tolist(),
        "method_count_note":"THREE_SECOND_ORDER_PLUS_TWO_RICHARDSON;ALL_FIVE_RETAINED",
        "frequencies":E.F,"grids":E.GRIDS,"scenes":E.scenes(),"resources":E.RES,
        "snapshots_s":[0,600,1200],"calibration":"C1b:unknown complex source per frequency/time;fixed complex gain per element/frequency across snapshots",
        "noise":E.NOISE,"sigma":E.SIGMAS,"absolute_floor_reference_amplitude":G.REF,
        "noise_scope":"FROZEN_DESIGN_ASSUMPTIONS;NOT_MEASURED_SNR",
        "source_samples":"EXACT_SAVED_DEPTHS;NO_INTERPOLATION_OR_NEW_DEPTHS",
        "projection":"real whitened tangent;C1b nuisance SVD removal;then four horizontal columns;independent pivoted QR",
        "rank_policy":"normalized nuisance columns cutoff1e-12;relative SVD cutoff1e-10;depth null loading>1e-8 => zero information and INF scale;no pseudoinverse zero variance",
        "comparisons":{"stencil":"ALL_TEN_METHOD_PAIRS_AT_BOTH_PROJECTION_STAGES;relative norm(a-b)/max(norm(a),norm(b),1e-30)",
            "truncation_law":"D2-D1 vs 4*(D1-D0.5);fitted ratio and vector residual;conditional smooth h2 expansion only",
            "grid":"same method/noise/scene/resource between80001/160001;both directions and effective information",
            "information":"squared horizontal-profiled depth norm divided by sigma squared;report relative to oldD1;rank and scale retained",
            "cancellation":"sum absolute modal terms / absolute sum;source-stencil coefficient-weighted pressure sum / absolute derivative",
            "sample_dependency":"R1,R2 share D1 and some samples;agreement is NOT independent derivative accuracy certification"},
        "internal_consistency_diagnostic_only":{
            "nominal_pressure_D1_Dhalf_replay_tolerance":1e-9,"nominal_D1_information_replay_tolerance":1e-6,
            "max_R1_R2_relative_direction":.01,"max_R1_R2_effective_information_relative":.05,
            "max_R1_R2_same_method_grid_direction":.02,"max_R1_R2_same_method_grid_information":.1,
            "policy":"INTERNAL_REVIEW_CLASSIFICATION_ONLY;NOT_RELAXED_ORIGINAL_2PCT_GATE;NO_PROJECT_PERFORMANCE_PASS"},
        "classification_rules":{
            "CONSISTENCY_SUPPORTED":"all internal diagnostics above satisfied",
            "ACCURACY_UNCERTIFIED":"no independent reference derivative or rigorous truncation/serialization remainder bound exists;may coexist with consistency",
            "MODAL_WINDOW_REDESIGN_REQUIRED":"CLOW excludes some physically admissible speed range under negative SSP and no excluded-mode contribution certificate",
            "modal_count":"NO_PREDICTION_FROM_NOMINAL_COUNTS;NO_OBSERVED_LOSS_WITHOUT_SOLVES",
            "future_window":"all SSP shifts AND nominal must share any redesigned window;new nominal baseline mandatory"},
        "budget":{"primary_analysis_wall_s":900,"independent_review_wall_s":900,"output_bytes":64000000,
                  "new_KRAKEN":0,"new_FIELD":0,"new_MC":0,"new_audio":0,"once_only":True},
        "stop":"after independent reconstruction,CommitB,push,remote verify;noSSP;noexpandeddepths;noautomaticnextstage",
        "original_P0":"NUMERICAL_OR_FORMULA_INCOMPLETE","original_H3_G0_G1_G2E":"FROZEN_UNCHANGED",
        "SSP_depth_information":"NOT_EVALUATED","R4_percent":0,
        "manual_source":"https://oalib-acoustics.org/website_resources/AcousticsToolbox/manual/node47.html",
        "bindings":bindings})
    (O/"DESIGN_NOTES.md").write_text(
        "# H3-P0 saved-depth and modal-window review\n\n"
        "The research lead authorized this stage only. The old 2% source-column Gate remains failed; original126/126 G0 replay is retained.\n\n"
        "All five stencils are frozen before execution. Shared samples cannot certify true derivative accuracy. Source-depth location, C-linear interpolation and serialized complex-float32 mode shapes will be audited as independent limitations.\n\n"
        "No new KRAKEN/FIELD, MC, receive recording, SSP fields or environment-information evaluation. Both consistency and accuracy-uncertified labels may apply. Unknown accuracy stops the current local Fisher environment route.\n",
        encoding="utf-8")
    print("FROZEN",len(bindings),"bindings",flush=True)
if __name__=="__main__":main()
