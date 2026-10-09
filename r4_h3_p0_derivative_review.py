"""Saved-data deterministic stencil review. No solver, FIELD, RNG or optimizer."""
from pathlib import Path
import time,json
import numpy as np
from scipy import linalg
import r4_h3_p0 as E
import r4_h3_g0 as G
O=E.ROOT/"results/R4_H3_P0_DERIVATIVE_REVIEW"
PARENT="ea021e88344167fc1fafcefe4f367631ff67efbe"
DEPTH=[198.,199.,199.5,200.,200.5,201.,202.]
METHODS=["D0.5","D1","D2","R1","R2"]
COEF=np.zeros((5,7))
for i,h in enumerate([.5,1.,2.]):
    COEF[i,DEPTH.index(200.+h)]=1/(2*h)
    COEF[i,DEPTH.index(200.-h)]=-1/(2*h)
COEF[3]=(4*COEF[0]-COEF[1])/3
COEF[4]=(4*COEF[1]-COEF[2])/3
PAIRS=[(i,j) for i in range(5) for j in range(i+1,5)]
def dump(p,x):E.dump(p,x)
def write(p,x):E.write(p,x)
def verify():
    for r in E.read(O/"DESIGN_FREEZE.json")["bindings"]:
        assert E.sha(E.ROOT/r["path"])==r["sha256"],r["path"]
def finite(x):return float(x) if np.isfinite(x) else "INF"
def rel(a,b):return float(np.linalg.norm(a-b)/max(np.linalg.norm(a),np.linalg.norm(b),1e-30))
def pressure_and_derivatives(mods,state,res):
    rr,_,zr=G.geometry(state,res);pr=[];dr=[];pc=[];dc=[];sc=[]
    for depths,phi,k in mods:
        ix={float(z):i for i,z in enumerate(depths)}
        assert all(z in ix for z in DEPTH+list(zr))
        recv=np.tile(phi[[ix[z] for z in zr]],(3,1)).T
        r=rr.ravel()
        prop=recv*np.sqrt(2*np.pi/(k[:,None]*r))*np.exp(-1j*k[:,None]*r-1j*np.pi/4)
        src=phi[[ix[z] for z in DEPTH]]
        p=src@prop
        dsrc=COEF@src;d=dsrc@prop
        pr.append(p.reshape(7,*rr.shape));dr.append(d.reshape(5,*rr.shape))
        pc.append((abs(src)@abs(prop)/np.maximum(abs(p),1e-300)).reshape(7,*rr.shape))
        dc.append((abs(dsrc)@abs(prop)/np.maximum(abs(d),1e-300)).reshape(5,*rr.shape))
        sc.append((abs(COEF)@abs(p)/np.maximum(abs(d),1e-300)).reshape(5,*rr.shape))
    return (np.stack(pr,-1),np.stack(dr,-1),np.stack(pc,-1),np.stack(dc,-1),np.stack(sc,-1))

def environment(path):
    lines=path.read_text(encoding="utf-8").splitlines()
    rows=[]
    for line in lines:
        a=line.split()
        if len(a)==6:
            try:
                if 0<=float(a[0])<=5000 and float(a[1])>1000:rows.append([float(a[0]),float(a[1])])
            except ValueError:pass
    low,high=map(float,lines[lines.index("'R' 0.0")+1].split())
    return np.array(rows),low,high,lines[3]

def modal_window():
    envs=[];margins=[];slopes=[]
    for g in E.GRIDS:
        for f in E.F:
            for delta in [-1,0,1]:
                path=E.modpath(f,g,delta).with_suffix(".env")
                ssp,cl,ch,opt=environment(path)
                zix=int(np.where(ssp[:,0]==200)[0][0])
                left=(ssp[zix,1]-ssp[zix-1,1])/(ssp[zix,0]-ssp[zix-1,0])
                right=(ssp[zix+1,1]-ssp[zix,1])/(ssp[zix+1,0]-ssp[zix,0])
                slopes.append({"frequency":f,"mesh":g,"offset":delta,"options":opt,
                    "source_at_SSP_knot":True,"left_c_slope":left,"right_c_slope":right,
                    "slope_jump":right-left,"smooth_high_order_assumption_certified":False})
                envs.append({"frequency":f,"mesh":g,"offset":delta,
                    "env_sha256":E.sha(path),"c_min":float(ssp[:,1].min()),"c_max":float(ssp[:,1].max()),
                    "CLOW":cl,"CHIGH":ch,"lower_window_above_minimum":cl>ssp[:,1].min(),
                    "actual_perturbed_mode_count":"NOT_EVALUATED","new_calls":0})
            z,phi,k=E.M.parse_mod(E.modpath(f,g,0))
            cp=2*np.pi*f/k.real
            for i,v in enumerate(cp):
                margins.append({"frequency":f,"mesh":g,"nominal_mode":i+1,"phase_speed":float(v),
                    "lower_margin":float(v-cl),"upper_margin":float(ch-v),
                    "within_1mps_lower_boundary":bool(v-cl<=1),
                    "predicted_perturbed_count":"NOT_PREDICTED"})
    write(O/"SSP_ENVIRONMENT_WINDOW_AUDIT.csv",envs)
    write(O/"NOMINAL_MODAL_WINDOW_MARGINS.csv",margins)
    write(O/"SOURCE_KNOT_SMOOTHNESS.csv",slopes)
    return envs,margins

def execute():
    start=time.monotonic();verify()
    assert not (O/"EXECUTION_STARTED.json").exists(),"ONE_FROZEN_REVIEW_ONLY"
    assert G.git("rev-parse","HEAD")==G.git("ls-remote","origin","refs/heads/main").split()[0]
    assert G.git("log","-1","--pretty=%s")=="R4 H3-P0: freeze depth derivative accuracy review"
    dump(O/"EXECUTION_STARTED.json",{"design_SHA":G.git("rev-parse","HEAD"),"new_KRAKEN":0,
                                   "new_FIELD":0,"new_MC":0,"new_received_audio":0})
    cmp=[];info=[];law=[];cancel=[];checks=[];saved={}
    old={(r["scene_id"],int(r["mesh"]),r["resource"],r["noise_model"],float(r["sigma"])):r
         for r in G.rows(G.O/"VERTICAL_INFORMATION_BY_SCENE.csv") if r["calibration"]=="C1b"}
    for grid in E.GRIDS:
        mods=[E.M.parse_mod(E.modpath(f,grid,0)) for f in E.F]
        for scene in E.scenes():
            for res in E.RES:
                ps,ds,pc,dc,sc=pressure_and_derivatives(mods,scene["state"],res)
                p=ps[3]
                with np.load(G.O/f"INPUT_{scene['scene_id']}_{res}_n{grid}.npz") as z:
                    original=z["pressure"];j=z["J"];jh=z["J_half"]
                for name,v in [("pressure",rel(p,original)),
                    ("D1",rel(ds[1],original*j[...,2]/G.SCALE[2])),
                    ("D0.5",rel(ds[0],original*jh[...,2]/G.SCALE[2]))]:
                    E.ck(checks,f"frozen_{name}_{scene['scene_id']}_{res}_{grid}",v,1e-9)
                np.savez(O/f"STENCILS_{scene['scene_id']}_{res}_n{grid}.npz",
                         source_depths=DEPTH,pressure=ps,derivatives=ds,
                         pressure_modal_cancellation=pc,derivative_modal_cancellation=dc,
                         pressure_stencil_cancellation=sc)
                for mi,method in enumerate(METHODS):
                    cancel.append({"scene":scene["scene_id"],"resource":res,"mesh":grid,"method":method,
                        "source_sample_depths":json.dumps([DEPTH[i] for i in range(7) if COEF[mi,i]!=0]),
                        "min_nominal_pressure":float(abs(p).min()),
                        "min_pressure_index_t_element_frequency":json.dumps(list(map(int,np.unravel_index(np.argmin(abs(p)),p.shape)))),
                        "derivative_cancellation_index_t_element_frequency":json.dumps(list(map(int,np.unravel_index(np.argmax(dc[mi]),p.shape)))),
                        "frequency_order":json.dumps(E.F),"snapshot_order_s":"[0,600,1200]",
                        "pressure_modal_cancellation_max":finite(pc.max()),
                        "derivative_modal_cancellation_max":finite(dc[mi].max()),
                        "derivative_stencil_cancellation_max":finite(sc[mi].max()),
                        "serialized_phi":"COMPLEX_FLOAT32;NOT_INDEPENDENT_PRECISION_CERTIFICATE",
                        "sample_dependence":"R1_AND_R2_SHARE_D1;NOT_INDEPENDENT_VALIDATION"})
                for noise in E.NOISE:
                    dz=G.projected(p,np.moveaxis(ds,0,-1)/p[...,None],noise,"C1b")
                    horiz=G.projected(p,j[...,[0,1,3,4]],noise,"C1b")
                    eff=G.remove(dz,horiz)
                    saved[scene["scene_id"],res,grid,noise]=(dz,eff)
                    for stage,a in [("GAIN_SOURCE_PROFILED",dz),("HORIZONTAL_PROFILED",eff)]:
                        for i,q in PAIRS:
                            cmp.append({"scene":scene["scene_id"],"resource":res,"mesh":grid,
                               "noise_model":noise,"projection_stage":stage,"method_a":METHODS[i],
                               "method_b":METHODS[q],"relative_direction_difference":rel(a[:,i],a[:,q]),
                               "diagnostic":"STENCIL_TRUNCATION_DIAGNOSTIC;NOT_TRUE_DERIVATIVE_ERROR"})
                        e10=a[:,1]-a[:,0];e21=a[:,2]-a[:,1]
                        ratio=float(e10@e21/max(e10@e10,1e-300))
                        law.append({"scene":scene["scene_id"],"resource":res,"mesh":grid,
                            "noise_model":noise,"projection_stage":stage,"fitted_h2_ratio":ratio,
                            "h2_expected_ratio_if_smooth_asymptotic":4,
                            "relative_h2_vector_residual":float(np.linalg.norm(e21-4*e10)/max(np.linalg.norm(e21),1e-30)),
                            "R1_R2_difference":rel(a[:,3],a[:,4]),
                            "true_error_bound":"NOT_ESTABLISHED"})
                    for mi,method in enumerate(METHODS):
                        for sigma in E.SIGMAS:
                            full=np.column_stack([horiz[:,0:2],dz[:,mi:mi+1]*200,horiz[:,2:4]])/sigma
                            met=G.metrics(full*sigma,sigma)
                            iz=met["depth_information_horizontal_profiled"]
                            prior=float(old[scene["scene_id"],grid,res,noise,sigma]["depth_information_horizontal_profiled"])
                            if method=="D1":E.ck(checks,f"frozen_Iz_{scene['scene_id']}_{res}_{grid}_{noise}_{sigma}",
                                               abs(iz-prior)/max(iz,prior,1e-30),1e-6)
                            # Rank of full 5-state matrix; depth column scaled as in original G0.
                            full=np.column_stack([horiz[:,0:2],dz[:,mi:mi+1]*200,horiz[:,2:4]])/sigma
                            sv=linalg.svdvals(full);rank=int(np.sum(sv>sv[0]*1e-10))
                            info.append({"scene":scene["scene_id"],"resource":res,"mesh":grid,"noise_model":noise,
                                "sigma":sigma,"method":method,"effective_depth_information":iz,
                                "relative_to_D1":iz/prior,"state_rank":rank,"depth_null_loading":met["depth_null_loading"],
                                "local_scale_m":1/np.sqrt(iz) if iz>0 else "INF",
                                "scope":"METHOD_SENSITIVITY_ONLY;NO_SSP_INFORMATION;NO_ACHIEVED_ACCURACY"})
                if time.monotonic()-start>900:raise TimeoutError("FROZEN_900S_CAP")
            print(f"REVIEW {grid} {scene['scene_id']}",flush=True)
    grids=[]
    for scene in E.scenes():
        for res in E.RES:
            for noise in E.NOISE:
                a,ea=saved[scene["scene_id"],res,80001,noise]
                b,eb=saved[scene["scene_id"],res,160001,noise]
                for mi,m in enumerate(METHODS):
                    ia=float(ea[:,mi]@ea[:,mi]);ib=float(eb[:,mi]@eb[:,mi])
                    grids.append({"scene":scene["scene_id"],"resource":res,"noise_model":noise,"method":m,
                        "gain_source_direction_grid_relative":rel(a[:,mi],b[:,mi]),
                        "horizontal_direction_grid_relative":rel(ea[:,mi],eb[:,mi]),
                        "effective_information_grid_relative":abs(ia-ib)/max(ia,ib,1e-30),
                        "diagnostic":"GRID_NUMERICAL_DIAGNOSTIC"})
    write(O/"STENCIL_COMPARISONS.csv",cmp);write(O/"TRUNCATION_LAW_DIAGNOSTIC.csv",law)
    write(O/"GRID_NUMERICAL_DIAGNOSTIC.csv",grids);write(O/"EFFECTIVE_INFORMATION_SENSITIVITY.csv",info)
    write(O/"CANCELLATION_AND_LOW_AMPLITUDE.csv",cancel);write(O/"NOMINAL_REPLAY_CHECKS.csv",checks)
    env,margin=modal_window()
    ri=[r for r in cmp if r["method_a"]=="R1" and r["method_b"]=="R2"]
    gi=[r for r in grids if r["method"] in ["R1","R2"]]
    paired=[]
    for s in E.scenes():
        for r in E.RES:
            for g in E.GRIDS:
                for n in E.NOISE:
                    xx=[v for v in info if v["scene"]==s["scene_id"] and v["resource"]==r and v["mesh"]==g and v["noise_model"]==n and v["sigma"]==.01 and v["method"] in ["R1","R2"]]
                    ia,ib=[x["effective_depth_information"] for x in xx]
                    paired.append(abs(ia-ib)/max(ia,ib,1e-30))
    consistency=all(r["PASS"] for r in checks) and max(r["relative_direction_difference"] for r in ri)<=.01 and max(paired)<=.05 and all(r["gain_source_direction_grid_relative"]<=.02 and r["effective_information_grid_relative"]<=.1 for r in gi)
    statuses=[]
    if consistency:statuses.append("H3_P0_DERIVATIVE_STENCIL_CONSISTENCY_SUPPORTED")
    statuses.append("H3_P0_DEPTH_DERIVATIVE_ACCURACY_UNCERTIFIED")
    window=any(r["offset"]==-1 and r["lower_window_above_minimum"] for r in env)
    if window:statuses.append("H3_P0_MODAL_WINDOW_REDESIGN_REQUIRED")
    dump(O/"REVIEW_DECISION.json",{"classifications":statuses,"stencil_internal_consistency":consistency,
        "true_derivative_accuracy":"NOT_CERTIFIED;NO_INDEPENDENT_REFERENCE_OR_REMAINDER_BOUND",
        "max_R1_R2_direction_relative":max(r["relative_direction_difference"] for r in ri),
        "max_R1_R2_effective_information_relative":max(paired),
        "max_high_order_grid_direction_relative":max(r["gain_source_direction_grid_relative"] for r in gi),
        "max_high_order_grid_information_relative":max(r["effective_information_grid_relative"] for r in gi),
        "current_window_ssp_coverage_certificate":"NOT_ESTABLISHED" if window else "REVIEW_ONLY",
        "perturbed_modes_observed":False,"original_P0":"NUMERICAL_OR_FORMULA_INCOMPLETE",
        "SSP_depth_information":"NOT_EVALUATED","R4_percent":0,"new_KRAKEN":0,"new_FIELD":0,
        "new_MC":0,"new_received_audio":0,"elapsed_s":time.monotonic()-start,
        "local_Fisher_environment_route":"STOP_ACCURACY_UNCERTIFIED",
        "next_stage":"NOT_AUTHORIZED;STOP"})
    assert sum(p.stat().st_size for p in O.rglob("*") if p.is_file())<=64000000
    print(json.dumps(E.read(O/"REVIEW_DECISION.json"),indent=2))
if __name__=="__main__":execute()
