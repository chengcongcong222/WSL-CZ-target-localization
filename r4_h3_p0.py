"""Once-only H3-P0 deterministic SSP sensitivity; no stochastic source/detector."""
from pathlib import Path
import os,json,csv,time,subprocess,hashlib
os.environ["OPENBLAS_NUM_THREADS"]="1";os.environ["OMP_NUM_THREADS"]="1"
import numpy as np
from scipy import linalg
import r4_h3_g0 as G
import r4_e2_g0 as M
import r4_e2_paired_recovery as P
import r4_h3_g0_audit as COLD
import r4_e2_nine_pair_pilot_audit as BIN
ROOT=Path(__file__).resolve().parent
O=ROOT/"results/R4_H3_P0_PROPAGATION_SENSITIVITY"
PARENT="1ab3c291dbf62bc4195aee79d665bac47a604a71"
F=list(P.RAW); GRIDS=[80001,160001]; RES=["B0","B1","B2"]
NOISE=["RELATIVE","ABSOLUTE_FLOOR"]; SIGMAS=[.01,.05]
START=0.; CALLS=0
def dump(path,data): G.dump(path,data)
def write(path,rows): G.write(path,rows)
def read(path): return G.read(path)
def sha(path):
    p=Path(path); b=p.read_bytes()
    if p.suffix in [".py",".json",".md",".env",".csv",".txt"]:
        b=b.replace(b"\r\n",b"\n")
    return hashlib.sha256(b).hexdigest()
def scenes():
    return [s for s in read(G.O/"DESIGN_FREEZE.json")["scenes"] if s["scene_id"].endswith("_M")]
def modpath(f,g,d):
    return P.path_for(f,g) if d==0 else O/"modal"/f"p0_f{int(f*2):04d}_s{d:+d}_n{g}.mod"
def guard():
    if time.monotonic()-START>5400: raise TimeoutError("FROZEN_TOTAL_5400S_CAP")
    if sum(p.stat().st_size for p in O.rglob("*") if p.is_file())>256000000:
        raise RuntimeError("FROZEN_DELIVERY_256MB_CAP")
def verify():
    for r in read(O/"DESIGN_FREEZE.json")["bindings"]:
        assert sha(ROOT/r["path"])==r["sha256"],r["path"]
def ck(rows,name,value,limit=1e-9):
    rows.append({"check":name,"value":float(value),"limit":limit,
                 "PASS":bool(np.isfinite(value) and value<=limit)})
def regression():
    controls=[]; fields={}
    old={ (r["scene_id"],int(r["mesh"]),r["resource"],r["noise_model"],float(r["sigma"])):r
          for r in G.rows(G.O/"VERTICAL_INFORMATION_BY_SCENE.csv") if r["calibration"]=="C1b"}
    for g in GRIDS:
        mods=[M.parse_mod(modpath(f,g,0)) for f in F]
        for sc in scenes():
            for res in RES:
                file=G.O/f"INPUT_{sc['scene_id']}_{res}_n{g}.npz"
                with np.load(file) as z: p,j,h=z["pressure"],z["J"],z["J_half"]
                pp,jj,hh,bound,free=G.field(mods,sc["state"],res)
                ck(controls,f"nominal_pressure_{sc['scene_id']}_{res}_{g}",
                   np.linalg.norm(pp-p)/np.linalg.norm(p))
                ck(controls,f"nominal_J_{sc['scene_id']}_{res}_{g}",
                   np.linalg.norm(jj-j)/np.linalg.norm(j))
                ck(controls,f"nominal_half_depth_J_{sc['scene_id']}_{res}_{g}",
                   np.linalg.norm(hh-h)/np.linalg.norm(h))
                ck(controls,f"actual_free_modal_gain_zero_{sc['scene_id']}_{res}_{g}",free,1e-10)
                fields[g,sc["scene_id"],res]=(p,j,h)
                # Cartesian path for perturbed polar-state parameters, no bearing/state truth prior.
                st=np.array(sc["state"]); r=np.linalg.norm(st[:2]);v=np.linalg.norm(st[2:])
                theta=np.arctan2(st[1],st[0]);psi=np.arctan2(st[3],st[2])
                steps=[.05,1e-6,1.,.00005,1e-5]
                for col in [0,1,3,4]:
                    vals=[r,theta,200.,v,psi]; ds=steps[col]; vals0=vals.copy(); vals1=vals.copy()
                    vals0[col]-=ds; vals1[col]+=ds
                    def cart(a): return [a[0]*np.cos(a[1]),a[0]*np.sin(a[1]),
                                        a[3]*np.cos(a[4]),a[3]*np.sin(a[4])]
                    lo,_=G.field(mods,cart(vals0),res,200.,False)
                    hi,_=G.field(mods,cart(vals1),res,200.,False)
                    fd=(hi-lo)/(2*ds)*G.SCALE[col]
                    ck(controls,f"analytic_FD_chain_{sc['scene_id']}_{res}_{g}_{col}",
                       np.linalg.norm(fd-p*j[...,col])/np.linalg.norm(p*j[...,col]),.002)
                for noise in NOISE:
                    z=G.projected(p,j,noise,"C1b")
                    zh=G.projected(p,h,noise,"C1b")
                    ck(controls,f"source_depth_steps_{sc['scene_id']}_{res}_{g}_{noise}",
                       np.linalg.norm(z[:,2]-zh[:,2])/np.linalg.norm(zh[:,2]),.02)
                    for sig in SIGMAS:
                        iz=G.metrics(z,sig)["depth_information_horizontal_profiled"]
                        prior=float(old[sc["scene_id"],g,res,noise,sig]["depth_information_horizontal_profiled"])
                        ck(controls,f"prior_information_{sc['scene_id']}_{res}_{g}_{noise}_{sig}",
                           abs(iz-prior)/max(iz,prior,1e-30),1e-7)
    write(O/"NOMINAL_REGRESSION_AND_CHAIN_CONTROLS.csv",controls)
    return controls,fields

def generate():
    global CALLS
    logs=[];mods={}
    for d in [-1,1]:
        for g in GRIDS:
            for f in F:
                guard()
                path=modpath(f,g,d);beg=time.monotonic();status="NOT_STARTED"
                assert CALLS<52
                CALLS+=1
                try:
                    remaining=5400-(time.monotonic()-START)
                    p=subprocess.run([str(M.AT/"kraken.exe"),path.stem],cwd=path.parent,
                                     capture_output=True,timeout=min(90,remaining))
                    path.with_suffix(".stdout.txt").write_bytes(p.stdout+p.stderr)
                    if p.returncode!=0: raise RuntimeError(f"KRAKEN_EXIT_{p.returncode}")
                    mods[d,g,f]=M.parse_mod(path);status="GENERATED"
                except Exception as e:
                    if isinstance(e,subprocess.TimeoutExpired):
                        path.with_suffix(".stdout.txt").write_bytes((e.stdout or b"")+(e.stderr or b""))
                    status=type(e).__name__+":"+str(e)
                    logs.append({"offset_mps":d,"mesh":g,"frequency_hz":f,
                     "call_number":CALLS,"timeout_s":90,"elapsed_s":time.monotonic()-beg,
                     "status":status,"mode_count":"","modal_sha256":sha(path) if path.exists() else ""})
                    write(O/"SOLVER_CALLS.csv",logs)
                    return logs,mods,False
                z,phi,k=mods[d,g,f]
                logs.append({"offset_mps":d,"mesh":g,"frequency_hz":f,
                     "call_number":CALLS,"timeout_s":90,"elapsed_s":time.monotonic()-beg,
                     "status":status,"mode_count":len(k),"modal_sha256":sha(path)})
                write(O/"SOLVER_CALLS.csv",logs)
                print(f"KRAKEN {CALLS}/52 offset={d:+} n={g} f={f} modes={len(k)} elapsed={time.monotonic()-START:.1f}",flush=True)
    return logs,mods,True

def source_chart(p,j,noise):
    n=p.shape[1];f=p.shape[2]
    D=np.zeros((n-1,n));D[:,0]=-1
    for i in range(n-1): D[i,i+1]=1
    L=linalg.block_diag(np.kron(D,np.eye(f)),np.kron(D,np.eye(f)))
    yy=[];gg=[]
    for t in range(3):
        var=np.ones(2*n*f)/2 if noise=="RELATIVE" else np.tile((G.REF/abs(p[t]).ravel())**2/2,2)
        W=linalg.solve_triangular(linalg.cholesky((L*var)@L.T,lower=True),L,lower=True)
        yy.append(W@G.rvec(j[t]));gg.append(W@np.eye(2*n*f))
    return G.remove(np.concatenate(yy),np.concatenate(gg))

def angular(a,b):
    na=np.linalg.norm(a);nb=np.linalg.norm(b)
    return abs(float(a@b)/(na*nb)) if na>0 and nb>0 else None

def forward_and_secants(mods,nominal):
    fidelity=[];alignment=[];numeric=[];zero=[];cache={}
    counts={}
    for d in [-1,0,1]:
        for f in F:
            z0,p0,k0=(M.parse_mod(modpath(f,80001,0)) if d==0 else mods[d,80001,f])
            z1,p1,k1=(M.parse_mod(modpath(f,160001,0)) if d==0 else mods[d,160001,f])
            same=len(k0)==len(k1);counts[d,f]=(len(k0),len(k1))
            corr=phase=None
            if same:
                corr=float(np.min(abs(np.sum(p0.conj()*p1,axis=0))/(np.linalg.norm(p0,axis=0)*np.linalg.norm(p1,axis=0))))
                phase=float(np.max(abs(k0.real-k1.real))*60000)
            good=same and corr is not None and corr>=.99 and phase<=.02
            fidelity.append({"kind":"MODAL","scene":"ALL","resource":"ALL","offset_mps":d,
              "frequency_hz":f,"source_receiver_depths_complete":bool(np.array_equal(z0,M.DEPTHS) and np.array_equal(z1,M.DEPTHS)),
              "modes_80001":len(k0),"modes_160001":len(k1),"min_phi_correlation":corr,
              "max_k_phase_difference_60km":phase,"raw_field_relative":"","response_projector_difference":"",
              "response_to_noise1_ratio":"","status":"PASS" if good else "INCOMPLETE"})
    # Fixed spectral window must not change included branch set over the registered difference stencil.
    for f in F:
        for gidx,g in enumerate(GRIDS):
            equal=len({counts[d,f][gidx] for d in [-1,0,1]})==1
            ck(numeric,f"fixed_phase_window_branch_set_{g}_{f}",0 if equal else 1,0)
    for sc in scenes():
        for res in RES:
            pp={}
            for g in GRIDS:
                p,j,h=nominal[g,sc["scene_id"],res]
                pp[0,g]=p
                for d in [-1,1]:
                    model=[mods[d,g,f] for f in F]
                    q,jq,hq,bounds,free=G.field(model,sc["state"],res)
                    pp[d,g]=q
                    # Independent binary parse and Cartesian field chain.
                    im=[BIN.mod(modpath(f,g,d)) for f in F]
                    qc,jc,hc=COLD.cold(im,sc["state"],res)
                    ck(numeric,f"cold_perturbed_pressure_{d}_{g}_{sc['scene_id']}_{res}",
                       np.linalg.norm(q-qc)/np.linalg.norm(q),1e-10)
                    ck(numeric,f"cold_perturbed_source_J_{d}_{g}_{sc['scene_id']}_{res}",
                       np.linalg.norm(q*jq[...,2]-jc[...,2])/np.linalg.norm(jc[...,2]),1e-9)
                    zero.append({"scene":sc["scene_id"],"resource":res,"mesh":g,"offset_mps":d,
                      "control":"C1g_FREE_SOURCE_DEPTH_ABSORPTION","residual":free,"limit":1e-10,
                      "PASS":free<=1e-10,"depth_information":0,"uncertainty":"INF"})
                    ck(numeric,f"free_modal_zero_{d}_{g}_{sc['scene_id']}_{res}",free,1e-10)
                eta=(pp[1,g]-pp[-1,g])/2
                forward=pp[1,g]-p;backward=p-pp[-1,g]
                aug=np.concatenate([j,(eta/p)[...,None]],axis=-1)
                for noise in NOISE:
                    zz=G.projected(p,aug,noise,"C1b")
                    # Positive and negative secants must agree AFTER source/fixed gain removal.
                    fb=np.stack([forward/p,backward/p],axis=-1)
                    fb=G.projected(p,fb,noise,"C1b")
                    locality=float(np.linalg.norm(fb[:,0]-fb[:,1])/max(np.linalg.norm(zz[:,5]),1e-30))
                    ck(numeric,f"SSP_secant_locality_{g}_{sc['scene_id']}_{res}_{noise}",locality,.2)
                    depth=zz[:,2]/G.SCALE[2];env=zz[:,5]
                    depthh=G.remove(depth[:,None],zz[:,[0,1,3,4]])[:,0]
                    envh=G.remove(env[:,None],zz[:,[0,1,3,4]])[:,0]
                    alignment.append({"scene":sc["scene_id"],"resource":res,"mesh":g,
                        "noise_model":noise,"absolute_cosine_after_gain_source":angular(depth,env),
                        "absolute_cosine_after_horizontal":angular(depthh,envh),
                        "positive_negative_secant_disagreement":locality,
                        "interpretation":"FINITE_SSP_SECANT_DIAGNOSTIC_NOT_CERTIFIED_LOCAL_DERIVATIVE"})
                    chart=source_chart(p,aug,noise)
                    ck(numeric,f"relative_chart_covariance_{g}_{sc['scene_id']}_{res}_{noise}",
                       np.linalg.norm(zz.T@zz-chart.T@chart)/max(np.linalg.norm(zz.T@zz),1e-30),1e-7)
                    raw=G.projected(p,aug,noise,"C0",True)
                    df=raw.T@raw-zz.T@zz
                    ck(numeric,f"DPI_same_raw_noise_{g}_{sc['scene_id']}_{res}_{noise}",
                       max(0,-np.linalg.eigvalsh((df+df.T)/2)[0]/np.linalg.norm(raw.T@raw)),1e-7)
                    C2=np.zeros_like(zz)
                    zero.append({"scene":sc["scene_id"],"resource":res,"mesh":g,"offset_mps":0,
                        "control":"C2_FREE_ALL_RESPONSES","residual":float(np.linalg.norm(C2)),
                        "limit":0,"PASS":True,"depth_information":0,"uncertainty":"INF"})
                    cache[g,sc["scene_id"],res,noise]=(zz,aug)
            for d in [-1,0,1]:
                a,b=pp[d,80001],pp[d,160001]
                raw=float(np.linalg.norm(a-b)/np.linalg.norm(b));resp=G.response_distance(a,b)
                ratio=max(resp/G.response_noise(b,.01,n) for n in NOISE)
                good=ratio<=.5
                fidelity.append({"kind":"FIELD","scene":sc["scene_id"],"resource":res,"offset_mps":d,
                  "frequency_hz":"ALL13","source_receiver_depths_complete":True,"modes_80001":"",
                  "modes_160001":"","min_phi_correlation":"","max_k_phase_difference_60km":"",
                  "raw_field_relative":raw,"response_projector_difference":resp,
                  "response_to_noise1_ratio":ratio,"status":"PASS" if good else "INCOMPLETE"})
            for noise in NOISE:
                za,_=cache[80001,sc["scene_id"],res,noise]
                zb,_=cache[160001,sc["scene_id"],res,noise]
                e=float(np.linalg.norm(za[:,5]-zb[:,5])/max(np.linalg.norm(zb[:,5]),1e-30))
                ck(numeric,f"SSP_secant_two_grid_{sc['scene_id']}_{res}_{noise}",e,.02)
                for sig in SIGMAS:
                    aa=G.remove(za[:,:5],za[:,5:6])
                    bb=G.remove(zb[:,:5],zb[:,5:6])
                    ma=G.metrics(aa,sig);mb=G.metrics(bb,sig)
                    ia=ma["depth_information_horizontal_profiled"];ib=mb["depth_information_horizontal_profiled"]
                    ck(numeric,f"effective_info_two_grid_{sc['scene_id']}_{res}_{noise}_{sig}",
                       abs(ia-ib)/max(ia,ib,1e-30),.1)
                    ck(numeric,f"profiled_rank_two_grid_{sc['scene_id']}_{res}_{noise}_{sig}",
                       0 if ma["rank"]==mb["rank"] else 1,0)
            guard()
    for row in fidelity: ck(numeric,"fidelity_"+str(len(numeric)),0 if row["status"]=="PASS" else 1,0)
    write(O/"FORWARD_NUMERICAL_FIDELITY.csv",fidelity)
    write(O/"DEPTH_SSP_DERIVATIVE_ALIGNMENT.csv",alignment)
    write(O/"NUMERICAL_AND_FORMULA_GATES.csv",numeric)
    write(O/"UNKNOWN_GAIN_ZERO_CONTROLS.csv",zero)
    return cache,numeric,fidelity

def information(nominal,cache,admitted):
    out=[]
    for sc in scenes():
        for res in RES:
            for noise in NOISE:
                for g in GRIDS:
                    p,j,h=nominal[g,sc["scene_id"],res]
                    z=G.projected(p,j,noise,"C1b")
                    for sig in SIGMAS:
                        b=G.metrics(z,sig)
                        row={"scene":sc["scene_id"],"resource":res,"mesh":g,"noise_model":noise,"sigma":sig,
                          "C1b_depth_information":b["depth_information_horizontal_profiled"],
                          "C1e_P0_depth_information":"","retention_ratio":"",
                          "depth_uncertainty":"NOT_EVALUATED_GATE_FAILED","rank":"",
                          "weak_direction":"","condition_identifiable":"",
                          "science_admitted":admitted,"horizontal_truth_prior":False,
                          "scope":"LOCAL_DETERMINISTIC_MEAN_MODEL_ONLY"}
                        if admitted:
                            aug,_=cache[g,sc["scene_id"],res,noise]
                            eff=G.remove(aug[:,:5],aug[:,5:6]);m=G.metrics(eff,sig)
                            iz=m["depth_information_horizontal_profiled"]
                            row.update(C1e_P0_depth_information=iz,retention_ratio=iz/b["depth_information_horizontal_profiled"],
                                depth_uncertainty=m["depth_local_scale_m"],rank=m["rank"],
                                weak_direction=m["weak_direction"],condition_identifiable=m["condition_identifiable"])
                        out.append(row)
    write(O/"PROFILED_DEPTH_INFORMATION.csv",out)
    match=[]
    ix={(r["scene"],r["resource"],r["mesh"],r["noise_model"],r["sigma"]):r for r in out}
    for r in out:
        if r["resource"]!="B1":continue
        q=ix[r["scene"],"B2",r["mesh"],r["noise_model"],r["sigma"]]
        match.append({"scene":r["scene"],"mesh":r["mesh"],"noise_model":r["noise_model"],"sigma":r["sigma"],
          "B1_C1b_information":r["C1b_depth_information"],"B2_C1b_information":q["C1b_depth_information"],
          "B1_to_B2_nominal":r["C1b_depth_information"]/q["C1b_depth_information"],
          "B1_C1e_P0_information":r["C1e_P0_depth_information"],"B2_C1e_P0_information":q["C1e_P0_depth_information"],
          "B1_to_B2_with_SSP":r["C1e_P0_depth_information"]/q["C1e_P0_depth_information"] if admitted and q["C1e_P0_depth_information"]>0 else "",
          "science_admitted":admitted})
    write(O/"RESOURCE_MATCHED_INFORMATION.csv",match)
    return out

def execute():
    global START
    verify()
    if (O/"EXECUTION_STARTED.json").exists():raise RuntimeError("ONE_FROZEN_EXECUTION_ONLY")
    assert G.git("rev-parse","HEAD")==G.git("ls-remote","origin","refs/heads/main").split()[0]
    assert G.git("log","-1","--pretty=%s")=="R4 H3-P0: freeze structured propagation sensitivity"
    START=time.monotonic()
    dump(O/"EXECUTION_STARTED.json",{"design_SHA":G.git("rev-parse","HEAD"),"epoch":time.time(),
                                    "once_only":True,"new_random_draws":0,"total_seconds_cap":5400})
    reg,nominal=regression()
    complete=False;numeric=[];fields=[];cache={};logs=[]
    if all(r["PASS"] for r in reg):
        logs,mods,complete=generate()
        if complete:cache,numeric,fields=forward_and_secants(mods,nominal)
    admitted=complete and all(r["PASS"] for r in reg+numeric)
    info=information(nominal,cache,admitted)
    decision="H3_P0_NUMERICAL_OR_FORMULA_INCOMPLETE"
    if admitted:
        main=[r for r in info if r["mesh"]==160001]
        ratios=[r["retention_ratio"] for r in main]
        if min(ratios)>=.5:decision="H3_P0_DEPTH_INFORMATION_RETAINS_UNDER_TESTED_SSP_PARAMETER"
        elif max(ratios)<=.1:decision="H3_P0_DEPTH_SSP_COUPLING_DOMINANT"
        else:decision="H3_P0_CONDITIONAL_OR_HETEROGENEOUS"
    dump(O/"H3_P0_DECISION.json",{"decision":decision,"design_SHA":G.git("rev-parse","HEAD"),
        "parent_SHA":PARENT,"new_KRAKEN_calls":CALLS,"provider_complete":complete,
        "nominal_regression_PASS":all(r["PASS"] for r in reg),
        "numerical_and_locality_gate_PASS":admitted,"failed_numerical_checks":[r["check"] for r in numeric if not r["PASS"]],
        "science_information_evaluated":admitted,"delta_c_mps":[-1,0,1],
        "delta_c_interpretation":"REGISTERED_MATHEMATICAL_DIAGNOSTIC_NOT_MEASURED_OCEAN_RANGE",
        "model":"G0_DETERMINISTIC_UNKNOWN_COMPLEX_SOURCE_AND_FIXED_ELEMENT_FREQUENCY_GAIN",
        "C1g_FREE":"DEPTH_ABSORBED_ZERO","C2":"ALL_STATE_ABSORBED_ZERO",
        "full_C1e":"NOT_ESTABLISHED","full_RC2_support":"NOT_ESTABLISHED",
        "actual_UUV_source":"NOT_ESTABLISHED","G2E_detector":"FROZEN_NOT_MODIFIED",
        "new_MC":0,"new_time_domain":0,"R4_percent":0,
        "elapsed_s":time.monotonic()-START,"next_stage":"NOT_AUTHORIZED;STOP"})
    guard()
    print(json.dumps(read(O/"H3_P0_DECISION.json"),indent=2))

if __name__=="__main__":execute()
