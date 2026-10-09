"""H3-P1 once-only preregistered propagation provider and local SSP response screen."""
from pathlib import Path
import os,time,json,subprocess
os.environ["OPENBLAS_NUM_THREADS"]="1";os.environ["OMP_NUM_THREADS"]="1"
import numpy as np
from scipy import optimize
import r4_h3_p0 as E
import r4_h3_g0 as G
import r4_h3_p0_derivative_review as R
O=E.ROOT/"results/R4_H3_P1_UNIFIED_WINDOW_SSP"
PARENT="7aaab2b99150c11b6d9a18092d4a97076d1c64af"
F=[150.,200.,250.];GRIDS=[80001,160001];LOWS=[1490.,1480.];DELTA=[0.,-.01,.01,-.02,.02]
RES=E.RES;NOISE=E.NOISE
def path(c,g,f,d):return O/"modal"/f"p1_c{int(c)}_n{g}_f{int(f)}_d{int(round(d*1000)):+04d}.mod"
def configs():return [(c,g,f,d) for c in LOWS for g in GRIDS for f in F for d in DELTA]
def verify():
    for b in E.read(O/"DESIGN_FREEZE.json")["bindings"]:
        assert E.sha(E.ROOT/b["path"])==b["sha256"],b["path"]
def rel(a,b):return float(np.linalg.norm(a-b)/max(np.linalg.norm(a),np.linalg.norm(b),1e-30))
def ck(rows,name,v,lim):E.ck(rows,name,v,lim)
def field(mods,state,res):
    rr,_,zr=G.geometry(state,res);out=[]
    for z,phi,k in mods:
        ix={float(v):i for i,v in enumerate(z)}
        rec=np.tile(phi[[ix[v] for v in zr]],(3,1)).T
        prop=rec*np.sqrt(2*np.pi/(k[:,None]*rr.ravel()))*np.exp(-1j*k[:,None]*rr.ravel()-1j*np.pi/4)
        out.append((phi[ix[200.]]@prop).reshape(rr.shape))
    return np.stack(out,-1)
def subset_field(mods,state,res,sets):
    sub=[(z,phi[:,ix],k[ix]) for (z,phi,k),ix in zip(mods,sets)]
    return field(sub,state,res)
def prof(p,j,n,cal="C1b"):return G.projected(p,j,n,cal)
def compare(p,q):
    raw=rel(p,q);response=G.response_distance(p,q)
    alpha=np.sum(q.conj()*p,axis=1,keepdims=True)/np.sum(abs(q)**2,axis=1,keepdims=True)
    gauge=rel(p,alpha*q)
    values={"raw_field_relative":raw,"source_removed_field_relative":gauge,"relative_projector_distance":response}
    for n in NOISE:
        values["response_noise1_ratio_"+n]=response/G.response_noise(p,.01,n)
        zz=prof(p,(q/p-1)[...,None],n)
        # Whitened residual per original real observation and 1% noise assumption.
        values["C1b_noise1_RMS_"+n]=float(np.linalg.norm(zz)/(.01*np.sqrt(2*p.size)))
    return values
def match(a,b):
    za,pa,ka=a;zb,pb,kb=b
    corr=abs(pa.conj().T@pb)/np.maximum(np.linalg.norm(pa,axis=0)[:,None]*np.linalg.norm(pb,axis=0)[None,:],1e-300)
    gap=np.minimum(np.r_[np.inf,abs(np.diff(ka.real))],np.r_[abs(np.diff(ka.real)),np.inf])
    cost=1-corr+.05*np.minimum(abs(ka.real[:,None]-kb.real[None,:])/np.maximum(gap[:,None],1e-12),20)
    ia,ib=optimize.linear_sum_assignment(cost)
    # This is a sampled-shape assignment diagnostic, not an exact modal identity proof.
    good=corr[ia,ib]>=.95
    bad_a=ia[~good];bad_b=ib[~good]
    un_a=np.setdiff1d(np.arange(len(ka)),ia);un_b=np.setdiff1d(np.arange(len(kb)),ib)
    return ia,ib,corr[ia,ib],un_a,un_b,bad_a,bad_b
def generate():
    started=time.monotonic();rows=[];available={};stop=False
    for number,(c,g,f,d) in enumerate(configs(),1):
        p=path(c,g,f,d);begin=time.monotonic()
        row={"number":number,"CLOW":c,"mesh":g,"frequency":f,"delta_c":d,"attempted":False,
             "timeout_s":90,"elapsed_s":0.,"status":"NOT_ATTEMPTED_AFTER_PROVIDER_FAILURE","modes":0,"modal_sha256":""}
        if not stop:
            row["attempted"]=True
            try:
                remain=5400-(begin-started);assert remain>0,"SOLVER_TOTAL_CAP"
                result=subprocess.run([str(E.M.AT/"kraken.exe"),p.stem],cwd=p.parent,capture_output=True,timeout=min(90,remain))
                p.with_suffix(".stdout.txt").write_bytes(result.stdout+result.stderr)
                assert result.returncode==0,"KRAKEN_EXIT_"+str(result.returncode)
                m=E.M.parse_mod(p)
                assert np.all(np.isfinite(m[1])) and np.all(np.isfinite(m[2]))
                available[c,g,f,d]=m
                row.update(status="GENERATED",modes=len(m[2]),modal_sha256=E.sha(p))
            except Exception as exc:
                if isinstance(exc,subprocess.TimeoutExpired):p.with_suffix(".stdout.txt").write_bytes((exc.stdout or b"")+(exc.stderr or b""))
                row["status"]=type(exc).__name__+":"+str(exc);stop=True
        row["elapsed_s"]=time.monotonic()-begin;rows.append(row)
        E.write(O/"SOLVER_CALLS.csv",rows)
        print("SOLVER",number,c,g,f,d,row["status"],flush=True)
    return rows,available
def evaluate(mods):
    begin=time.monotonic();mode=[];assignment=[];sets={}
    for c,g,f,d in configs():
        z,phi,k=mods[c,g,f,d];speed=2*np.pi*f/k.real
        mode.append({"CLOW":c,"mesh":g,"frequency":f,"delta_c":d,"modes":len(k),
            "phase_speed_min":float(speed.min()),"phase_speed_max":float(speed.max()),
            "CLOW_margin_min":float(speed.min()-c),"CHIGH_margin_min":float(1800-speed.max()),
            "exact_depth_count":len(z),"all_requested_depths":bool(np.array_equal(z,E.M.DEPTHS)),
            "mode_count_equality_required":False,"CHIGH_certification":"NOT_INDEPENDENTLY_CERTIFIED"})
    # Compare nominal vs every offset, and compare windows at each offset.
    pairs=[("OFFSET",c,g,f,0.,c,g,f,d) for c in LOWS for g in GRIDS for f in F for d in DELTA if d!=0]
    pairs += [("WINDOW",1490.,g,f,d,1480.,g,f,d) for g in GRIDS for f in F for d in DELTA]
    for kind,c,g,f,d,cc,gg,ff,dd in pairs:
        ia,ib,corr,ua,ub,ba,bb=match(mods[c,g,f,d],mods[cc,gg,ff,dd])
        sets[kind,c,g,f,d,cc,gg,ff,dd]=(ua,ub,ba,bb)
        assignment.append({"kind":kind,"CLOW_a":c,"CLOW_b":cc,"mesh":g,"frequency":f,
            "delta_a":d,"delta_b":dd,"matched":len(ia),"min_sampled_shape_correlation":float(corr.min()),
            "unmatched_a":len(ua),"unmatched_b":len(ub),"low_correlation_a":len(ba),"low_correlation_b":len(bb),
            "unmatched_mode_ids_a":json.dumps(ua.tolist()),"unmatched_mode_ids_b":json.dumps(ub.tolist()),
            "matching_scope":"SAMPLED_13_DEPTH_ASSIGNMENT;NOT_EXACT_IDENTITY_CERTIFICATE"})
    E.write(O/"MODE_BOUNDARY_AND_CONTINUITY.csv",mode);E.write(O/"MODE_MATCHING.csv",assignment)
    fields={};checks=[];local=[];stable=[];contrib=[];cache={}
    for c in LOWS:
        for g in GRIDS:
            for sc in E.scenes():
                for res in RES:
                    for d in DELTA:
                        p=field([mods[c,g,f,d] for f in F],sc["state"],res)
                        assert np.isfinite(p).all() and np.all(abs(p)>0)
                        fields[c,g,sc["scene_id"],res,d]=p
                    p=fields[c,g,sc["scene_id"],res,0.]
                    for n in NOISE:
                        cols=[]
                        for h in [.01,.02]:
                            hi=fields[c,g,sc["scene_id"],res,h];lo=fields[c,g,sc["scene_id"],res,-h]
                            cols.extend([(hi-lo)/(2*h)/p,(hi-p)/h/p,(p-lo)/h/p])
                        zz=prof(p,np.stack(cols,-1),n);cache[c,g,sc["scene_id"],res,n]=zz
                        zz0=prof(p,np.stack(cols,-1),n,"C0")
                        for h,i in [(.01,0),(.02,3)]:
                            centered=zz[:,i];fb=rel(zz[:,i+1],zz[:,i+2])
                            common=rel(zz0[:,i+1],zz0[:,i+2])
                            hi=fields[c,g,sc["scene_id"],res,h];lo=fields[c,g,sc["scene_id"],res,-h]
                            # Source-invariant pressure/reference quotient phase (wrapped principal value).
                            spatial_hi=(hi/hi[:,0:1,:])/(p/p[:,0:1,:])
                            spatial_lo=(lo/lo[:,0:1,:])/(p/p[:,0:1,:])
                            maxphase=float(max(abs(np.angle(spatial_hi)).max(),abs(np.angle(spatial_lo)).max()))
                            row={"CLOW":c,"mesh":g,"scene":sc["scene_id"],"resource":res,"noise_model":n,
                                "step":h,"two_step_C1b_direction_relative":rel(zz[:,0],zz[:,3]),
                                "forward_backward_C1b_relative":fb,"forward_backward_source_only_relative":common,
                                "source_only_two_step_relative":rel(zz0[:,0],zz0[:,3]),
                                "C1b_to_source_only_norm":float(np.linalg.norm(centered)/max(np.linalg.norm(zz0[:,i]),1e-30)),
                                "max_source_removed_principal_phase_rad":maxphase,
                                "phase_wrap_control":"PRINCIPAL_PHASE_DIAGNOSTIC;NO_UNWRAPPED_DERIVATIVE_CLAIM",
                                "min_nominal_pressure":float(abs(p).min()),
                                "min_pressure_index":json.dumps(list(map(int,np.unravel_index(np.argmin(abs(p)),p.shape)))),
                                "SSP_step_scope":"NUMERICAL_LOCALITY_ONLY;NOT_OCEAN_UNCERTAINTY_RANGE"}
                            local.append(row)
                    np.savez_compressed(O/f"FIELDS_c{int(c)}_n{g}_{sc['scene_id']}_{res}.npz",
                                         delta_c=DELTA,pressure=np.stack([fields[c,g,sc["scene_id"],res,d] for d in DELTA]))
                print("FIELD",c,g,sc["scene_id"],flush=True)
    for sc in E.scenes():
        for res in RES:
            for c in LOWS:
                for d in DELTA:
                    a=fields[c,80001,sc["scene_id"],res,d];b=fields[c,160001,sc["scene_id"],res,d]
                    stable.append({"kind":"GRID","scene":sc["scene_id"],"resource":res,"CLOW":c,"mesh":"80001_VS_160001","delta_c":d,**compare(a,b)})
            for g in GRIDS:
                for d in DELTA:
                    a=fields[1490.,g,sc["scene_id"],res,d];b=fields[1480.,g,sc["scene_id"],res,d]
                    stable.append({"kind":"WINDOW","scene":sc["scene_id"],"resource":res,"CLOW":"1490_VS_1480","mesh":g,"delta_c":d,**compare(a,b)})
            for n in NOISE:
                for c in LOWS:
                    a=cache[c,80001,sc["scene_id"],res,n];b=cache[c,160001,sc["scene_id"],res,n]
                    for h,i in [(.01,0),(.02,3)]:
                        ck(checks,"GRID_SECANT:"+str((sc["scene_id"],res,n,c,h)),rel(a[:,i],b[:,i]),.02)
                for g in GRIDS:
                    a=cache[1490.,g,sc["scene_id"],res,n];b=cache[1480.,g,sc["scene_id"],res,n]
                    for h,i in [(.01,0),(.02,3)]:
                        ck(checks,"WINDOW_SECANT:"+str((sc["scene_id"],res,n,g,h)),rel(a[:,i],b[:,i]),.02)
            for kind,c,g,f,d,cc,gg,ff,dd in pairs:
                ua,ub,ba,bb=sets[kind,c,g,f,d,cc,gg,ff,dd]
                ma=mods[c,g,f,d];mb=mods[cc,gg,ff,dd]
                pa=field([ma],sc["state"],res);pb=field([mb],sc["state"],res)
                def fraction(m,ix,p):return float(np.linalg.norm(subset_field([m],sc["state"],res,[ix]))/max(np.linalg.norm(p),1e-30))
                contrib.append({"kind":kind,"scene":sc["scene_id"],"resource":res,"CLOW_a":c,"CLOW_b":cc,
                    "mesh":g,"frequency":f,"delta_a":d,"delta_b":dd,
                    "unmatched_field_fraction_a":fraction(ma,ua,pa),"unmatched_field_fraction_b":fraction(mb,ub,pb),
                    "low_correlation_field_fraction_a":fraction(ma,ba,pa),"low_correlation_field_fraction_b":fraction(mb,bb,pb),
                    "full_field_difference":rel(pa,pb),"all_modes_retained":True,
                    "interpretation":"MODE_SET_DIAGNOSTIC;COUNT_DIFFERENCE_ALONE_IS_NOT_FAILURE"})
    for row in stable:
        name=str((row["kind"],row["scene"],row["resource"],row["CLOW"],row["mesh"],row["delta_c"]))
        ck(checks,row["kind"]+"_RAW:"+name,row["raw_field_relative"],.01)
        for n in NOISE:
            ck(checks,row["kind"]+"_RESPONSE:"+name+n,row["response_noise1_ratio_"+n],.5)
            ck(checks,row["kind"]+"_GAIN_RESPONSE:"+name+n,row["C1b_noise1_RMS_"+n],.5)
    for row in local:
        name=str((row["CLOW"],row["mesh"],row["scene"],row["resource"],row["noise_model"],row["step"]))
        ck(checks,"LOCAL_STEP:"+name,row["two_step_C1b_direction_relative"],.05)
        ck(checks,"LOCAL_FB:"+name,row["forward_backward_C1b_relative"],.2)
        ck(checks,"LOCAL_PHASE:"+name,row["max_source_removed_principal_phase_rad"],np.pi/2)
    for row in contrib:
        if row["kind"]=="WINDOW":
            ck(checks,"WINDOW_UNMATCHED:"+str(len(checks)),max(row["unmatched_field_fraction_a"],row["unmatched_field_fraction_b"]),.005)
    for row in mode:
        ck(checks,"WINDOW_MARGIN:"+str(len(checks)),1 if row["CLOW_margin_min"]<=1 or row["CHIGH_margin_min"]<=1 else 0,0)
    E.write(O/"RELATIVE_RESPONSE_NUMERICAL_STABILITY.csv",stable)
    E.write(O/"SSP_SECANT_LOCALITY.csv",local);E.write(O/"WINDOW_FIELD_CONTRIBUTION.csv",contrib)
    E.write(O/"PREREGISTERED_CHECKS.csv",checks)
    assert time.monotonic()-begin<=900,"ANALYSIS_900S_CAP"
    window=all(r["PASS"] for r in checks if not r["check"].startswith("LOCAL_"))
    locality=all(r["PASS"] for r in checks if r["check"].startswith("LOCAL_"))
    classes=[]
    if not window:classes.append("H3_P1_MODAL_WINDOW_NOT_CERTIFIED")
    if not locality:classes.append("H3_P1_SSP_SECANT_LOCALITY_NOT_ESTABLISHED")
    if window and locality:classes.append("H3_P1_LOCAL_SSP_RESPONSE_NUMERICALLY_SUPPORTED")
    return classes,window,locality
def execute():
    verify();assert not (O/"EXECUTION_STARTED.json").exists()
    assert G.git("log","-1","--pretty=%s")=="R4 H3-P1: freeze unified modal window and SSP local response"
    assert G.git("rev-parse","HEAD")==G.git("ls-remote","origin","refs/heads/main").split()[0]
    E.dump(O/"EXECUTION_STARTED.json",{"design_SHA":G.git("rev-parse","HEAD"),"once_only":True})
    calls,mods=generate()
    classes=["H3_P1_PROVIDER_OR_FORMULA_INVALID"];window=locality=False
    if len(mods)==60:
        try:classes,window,locality=evaluate(mods)
        except Exception as exc:
            (O/"ANALYSIS_FAILURE.txt").write_text(type(exc).__name__+":"+str(exc),encoding="utf-8")
    E.dump(O/"H3_P1_DECISION.json",{"classifications":classes,"window_numerical_support":window,
        "SSP_secant_locality":locality,"new_KRAKEN":sum(r["attempted"] for r in calls),
        "generated_modal_files":len(mods),"new_FIELD":0,"new_MC":0,"new_audio":0,
        "original_P0":"H3_P0_NUMERICAL_OR_FORMULA_INCOMPLETE","source_depth_true_accuracy":"ACCURACY_UNCERTIFIED",
        "C1e_P0_depth_information":"NOT_EVALUATED","complete_RC2_support":"NOT_ESTABLISHED",
        "R4_percent":0,"CHIGH_independent_truncation_certification":"NOT_ESTABLISHED",
        "scope":"THREE_FREQUENCY_FROZEN_WINDOW_AND_STEP_NUMERICAL_LOCALITY_ONLY","next":"STOP"})
    assert sum(p.stat().st_size for p in O.rglob("*") if p.is_file())<=256000000
if __name__=="__main__":execute()
