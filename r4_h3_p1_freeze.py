"""Generate and hash-lock all60 environments before solver execution."""
import numpy as np
import r4_h3_p1 as P
E=P.E;G=P.G;O=P.O
def main():
    assert G.git("rev-parse","HEAD")==P.PARENT;assert not O.exists()
    (O/"modal").mkdir(parents=True)
    envrows=[];paths=[]
    for c,g,f,d in P.configs():
        base=E.modpath(f,g,0).with_suffix(".env")
        lines=base.read_text(encoding="utf-8").splitlines();out=[];changed=0
        for i,line in enumerate(lines):
            a=line.split()
            try:
                water=len(a)==6 and 0<=float(a[0])<=5000 and float(a[1])>1000
            except ValueError:water=False
            if water:
                a[1]=f"{float(a[1])+d:.7f}";out.append(" ".join(a));changed+=1
            elif i>0 and lines[i-1].strip()=="'R' 0.0":out.append(f"{c:.1f} 1800.0")
            else:out.append(line)
        assert changed==251
        p=P.path(c,g,f,d).with_suffix(".env");p.write_text("\n".join(out)+"\n",encoding="utf-8")
        s,cl,ch,opt=P.R.environment(p);sb,_,_,_=P.R.environment(base)
        assert np.array_equal(s[:,0],sb[:,0]) and np.max(abs(s[:,1]-sb[:,1]-d))<1e-10
        assert cl==c and ch==1800
        assert float(out[1])==f and int(out[4].split()[0])==g
        # Only sound-speed rows and the common lower modal window may change.
        for x,y in zip(lines,out):
            if x.strip()==y.strip():continue
            xx=x.split();yy=y.split()
            if len(xx)==len(yy)==6:
                assert all(float(xx[i])==float(yy[i]) for i in [0,2,3,4,5])
            else:assert len(xx)==len(yy)==2 and float(yy[0])==c and float(yy[1])==1800
        assert len(lines)==len(out)
        envrows.append({"CLOW":c,"CHIGH":ch,"mesh":g,"frequency":f,"delta_c":d,
            "water_rows":changed,"min_sound_speed":float(s[:,1].min()),"max_sound_speed":float(s[:,1].max()),
            "only_registered_changes":True,"all_13_exact_output_depths":True,"source_depth_m":200.,
            "env_sha256":E.sha(p),"PASS":True})
        paths.extend([p,base])
    files=[E.ROOT/n for n in ["r4_h3_p1.py","r4_h3_p1_freeze.py","r4_h3_p1_audit.py","r4_h3_p1_report.py",
        "r4_h3_p0.py","r4_h3_g0.py","r4_h3_g0_audit.py","r4_e2_g0.py","r4_e2_paired_recovery.py",
        "r4_e2_nine_pair_pilot_audit.py","r4_h3_p0_derivative_review.py","r4_h3_p0_derivative_review_audit.py"]]
    files += [E.M.AT/"kraken.exe",P.R.O/"DESIGN_FREEZE.json",P.R.O/"REVIEW_DECISION.json"]
    files += [P.R.O/f"STENCILS_{sc['scene_id']}_{res}_n{g}.npz" for sc in E.scenes() for res in P.RES for g in P.GRIDS]
    files += [G.O/"DESIGN_FREEZE.json"]
    files += [E.modpath(f,g,0) for f in P.F for g in P.GRIDS]
    unique=sorted(set(p.resolve() for p in files+paths))
    bindings=[{"path":str(p.relative_to(E.ROOT)).replace("\\","/"),"sha256":E.sha(p)} for p in unique]
    E.write(O/"UNIFIED_ENVIRONMENT_AUDIT.csv",envrows)
    E.dump(O/"DESIGN_FREEZE.json",{
        "stage":"R4_H3_P1_UNIFIED_MODAL_WINDOW_AND_LOCAL_SSP_RESPONSE","parent_SHA":P.PARENT,
        "research_lead_accepted_previous":["H3_P0_DERIVATIVE_STENCIL_CONSISTENCY_SUPPORTED","H3_P0_DEPTH_DERIVATIVE_ACCURACY_UNCERTIFIED","H3_P0_MODAL_WINDOW_REDESIGN_REQUIRED"],
        "frequencies_hz":P.F,"grids":P.GRIDS,"CLOW":P.LOWS,"CHIGH":1800,
        "delta_c_mps":P.DELTA,"exact_depths_m":E.M.DEPTHS.tolist(),"source_depth_m":200.,
        "scenes":E.scenes(),"resources":P.RES,"snapshots_s":[0,600,1200],
        "scenario_scope":"INHERIT_H3_P0_H01_H06_H12_M;NO_NEW_SCENARIOS",
        "environment":"E-STD;uniform water sound-speed offsets only;all other profiles,boundaries,attenuation,depths,CVWT unchanged",
        "historical_1500":"HASHED_HISTORY_ONLY;NOT_USED_AS_NEW_NOMINAL",
        "nominal_provider":"NEW_DELTA_ZERO_FOR_EACH_WINDOW_GRID_FREQUENCY",
        "field":"full sum phi_source*phi_receiver sqrt(2pi/(k*r)) exp(-ikr-i*pi/4);complex k retained;no mode deletion",
        "noise":"inherited RELATIVE and ABSOLUTE_FLOOR;REF="+str(G.REF)+";noise1pct=.01;design assumptions only",
        "unknowns":"C0 source per frequency/snapshot;C1b additionally fixed unknown complex gain per element/frequency shared across all3snapshots",
        "projection":"real whitened tangent SVD;independent QR/chart;no formal source-depth/SSP Fisher profiling",
        "priority":["provider_and_independent_formula","complete_mode_sets_and_grid_window_response","secant_locality"],
        "pre_registered_gates":{
            "provider":"exact depths,finite complex fields,zero failed solver/binary contracts;independent field rel<=1e-9",
            "grid_and_window_raw_relative":.01,"grid_and_window_source_projector_to_noise1":.5,
            "grid_and_window_C1b_whitened_noise1_RMS":.5,
            "grid_and_window_C1b_SSP_direction_relative":.02,
            "window_unmatched_field_fraction":.005,"lower_and_upper_phase_speed_margin_mps":1.,
            "C1b_01_vs_02_centered_secant_direction_relative":.05,
            "C1b_forward_backward_secant_relative":.2,
            "source_invariant_principal_phase_max_rad":float(np.pi/2),
            "cold_QR_derived_diagnostic_reconstruction":1e-7,
            "scope":"CONJUNCTIVE_NUMERICAL_SCREEN;NOT_NEW_PROJECT_PERFORMANCE_PASS;CHIGH_NOT_INDEPENDENTLY_CERTIFIED"},
        "comparisons":{
            "relative_norm":"norm(a-b)/max(norm(a),norm(b),1e-30)",
            "source_projector":"normalize each t/f spatial complex vector;average squared Frobenius projector difference;noise inherited exact G0 formula",
            "fixed_gain_response":"C1b projection of q/p-1, divided by .01*sqrt(2*original complex sample count)",
            "secants":"(P(+h)-P(-h))/(2h), forward(P(+h)-P0)/h,backward(P0-P(-h))/h;divide by same new nominal field then joint nuisance projection",
            "step_locality":"centered .01 vs .02 relative direction AND positive/negative agreement at BOTHsteps;all retained",
            "phase":"principal phase of source-invariant element/reference response vs new nominal;not unwrapped derivative certificate",
            "mode_matching":"Hungarian assignment1-correlation+0.05*clipped k-distance/local gap;13sampled-depth shape correlation;corr<.95 flagged;not exact identity certificate",
            "unmatched":"all apparent unmatched full modal terms retained and summed;window omitted contribution norm/fullfield norm<=.005;count equality not required",
            "ambiguous":"low-correlation modal fields retained/report;mode assignment not used to truncate or construct secants"},
        "classification":{
            "A":"all provider,cold,window/grid and locality checks pass;only three frequencies/tested window steps",
            "B":"any window/grid/full-response/margin/unmatched or grid/window secant check fails",
            "C":"any step/forward-backward/principal-phase locality check fails",
            "D":"provider,parse,formula or cold reconstruction substantive failure;priority overrides scientific interpretations",
            "B_C_can_coexist":True,"no_retuning":"no further CLOW,SSPsteps,frequencies or grids in this stage"},
        "source_depth":"archived R1/R2 limits retained;no new source derivative accuracy certification;NO_ALIGNMENT_OR_C1e_INFORMATION_EVALUATION",
        "budget":{"new_KRAKEN_max":60,"per_call_s":90,"solver_total_s":5400,"analysis_s":900,"cold_review_s":900,
                  "overall_s":7200,"delivery_bytes":256000000,"new_FIELD":0,"new_MC":0,"new_audio":0,
                  "solver_retries":0,"run_once":True,"stop_after_provider_failure":True},
        "stop":"CommitB,push,remoteverify,STOP;A allows only recommendation for separate full-band research-lead design review;B/C/D stop current propagation derivative route",
        "original_P0":"H3_P0_NUMERICAL_OR_FORMULA_INCOMPLETE","source_depth_true_accuracy":"ACCURACY_UNCERTIFIED",
        "C1e_P0_depth_information":"NOT_EVALUATED","complete_RC2_support":"NOT_ESTABLISHED","R4_percent":0,
        "bindings":bindings})
    (O/"SOURCE_DEPTH_DERIVATIVE_LIMITS.md").write_text(
        "# Source-depth derivative limits\n\n"
        "Archived R1/R2 are retained as uncertainty context only. Their maximum projected-direction difference0.239% and effective-information difference0.355% do not certify a true error bound. They share samples and use complex-float32 mode shapes;source200m is a C-linear interpolation knot.\n\n"
        "The current stage evaluates SSP secants only. No source-depth/SSP angle, joint Fisher information or C1e-P0 depth survival is computed. Old1500m/s nominal source stencils cannot substitute for new-window SSP nominal fields. No true source-derivative error is set to zero.\n",encoding="utf-8")
    print("FROZEN",len(envrows),"envs",len(bindings),"bindings")
if __name__=="__main__":main()
