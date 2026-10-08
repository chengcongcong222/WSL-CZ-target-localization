"""Frozen conditional frequency-domain experiment; no propagation solver calls."""
import csv, hashlib, json, sys, time
from pathlib import Path
import numpy as np
from scipy.linalg import cho_factor, cho_solve
from threadpoolctl import threadpool_limits
import r4_h3_g2e_interface as I

ROOT=Path(__file__).resolve().parent
OUT=ROOT/"results/R4_H3_G2E_CONDITIONAL_EXTRACTION"
G0=ROOT/"results/R4_H3_G0_FINITE_VERTICAL_RESPONSE"
F=np.array([150,152,155,160,195,198,200,205,240,244,246,249,250])
PACK={"F1":[6],"F3":[0,6,12],"F13":list(range(13))}
CH={"B0":list(range(8)),"B1":list(range(12)),"B2":list(range(8))+list(range(12,16))}
KS=[1,8,32]
GRIDS=[160001,80001]
REF=.00010824612978073539
START=0

def write_json(p,x):
    p.write_text(json.dumps(x,ensure_ascii=False,indent=2,allow_nan=False)+"\n",encoding="utf-8")

def write_csv(p,rows):
    with p.open("w",encoding="utf-8",newline="") as f:
        w=csv.DictWriter(f,fieldnames=list(rows[0]))
        w.writeheader(); w.writerows(rows)

def sha(p):
    b=p.read_bytes()
    if p.suffix.lower() in (".py",".md",".json",".csv",".txt"):
        b=b.replace(b"\r\n",b"\n")
    return hashlib.sha256(b).hexdigest()

def array_sha(a):
    return hashlib.sha256(np.ascontiguousarray(a).tobytes()).hexdigest()

def gram():
    n=np.arange(2000)
    w=.5-.5*np.cos(2*np.pi*n/2000)
    e=np.exp(-2j*np.pi*F[:,None]*n[None,:]/1000)
    return I.herm((e*w)@(e*w).conj().T/(w@w))

def fields(scene,grid):
    b1=np.load(G0/f"INPUT_{scene}_M_B1_n{grid}.npz")["pressure"]
    b2=np.load(G0/f"INPUT_{scene}_M_B2_n{grid}.npz")["pressure"]
    b0=np.load(G0/f"INPUT_{scene}_M_B0_n{grid}.npz")["pressure"]
    assert np.array_equal(b0,b1[:,:8]) and np.array_equal(b0,b2[:,:8])
    return np.concatenate([b1,b2[:,8:]],axis=1)

def complex_draw(rng,shape):
    return (rng.standard_normal(shape)+1j*rng.standard_normal(shape))/np.sqrt(2)

def make_pool(scene_i,rep):
    # Separate substreams for source, observed noise, independent background.
    ss=np.random.SeedSequence([20261008,32051,scene_i,rep])
    rs,rn,rb=[np.random.Generator(np.random.PCG64(x)) for x in ss.spawn(3)]
    return {"source":complex_draw(rs,(3,32,13)),
            "noise":complex_draw(rn,(3,32,16,13)),
            "background":complex_draw(rb,(3,256,16,13))}

def generate(p,pool,kind,sigma,cf):
    d=sigma*np.abs(p) if kind=="RELATIVE" else np.full(p.shape,sigma*REF)
    # Row samples: multiply L.T so E[n_i conj(n_j)] = CF_ij.
    l=np.linalg.cholesky(cf)
    nn=np.einsum("tkcf,jf->tkcj",pool["noise"],l)
    bb=np.einsum("tkcf,jf->tkcj",pool["background"],l)
    y=p[:,None,:,:]*pool["source"][:,:,None,:]+d[:,None,:,:]*nn
    bg=d[:,None,:,:]*bb
    return y,bg,d

def candidate_response(samples,n,nf,t):
    # Deliberately non-truth, registered likelihood test parameters, not field fit.
    r=I.covariance(samples)
    scale=np.sqrt(max(float(np.trace(r).real)/(n*nf),np.finfo(float).tiny))
    c=np.arange(n)[:,None]; f=np.arange(nf)[None,:]
    h=(1+.01*c+.02*t)*np.exp(1j*(.03*c+.02*t*c+.005*f*c))
    g=(1+.03*c)*np.exp(1j*.04*c*(f+1))  # fixed across three templates
    q=.7+np.arange(nf)/20+t/10
    return scale*h*g,q

def guard():
    elapsed=time.perf_counter()-START
    size=sum(p.stat().st_size for p in OUT.rglob("*") if p.is_file())
    if elapsed>3600 or size>1000000000:
        raise RuntimeError(f"FROZEN_BUDGET_STOP seconds={elapsed} bytes={size}")

def verify_inputs(d):
    for r in d["bindings"]:
        assert sha(ROOT/r["path"])==r["sha256"], r["path"]

def extract_base(y,bg,d,p,cf):
    # Interface receives only observed Y/background, except separately labeled oracle.
    # All frequencies are extracted once; packages then select registered components.
    shape=(3,3,3,13)
    arrays={k:np.full(shape,np.nan) for k in (
        "before_error","after_error","oracle_error","background_error","energy",
        "threshold","oracle_energy","oracle_threshold","grid_unit",
        "noise_to_min_amplitude","signal_eigenvalue","rank")}
    arrays["detected"]=np.zeros(shape,bool)
    arrays["oracle_detected"]=np.zeros(shape,bool)
    arrays["u"]=np.zeros(shape+(12,),complex)
    arrays["u_before"]=np.zeros(shape+(12,),complex)
    likelihood=[]
    for ri,(res,ch) in enumerate(CH.items()):
        n=len(ch)
        for t in range(3):
            for fi in range(13):
                b=bg[t,:,ch,fi].T
                # numpy advanced indexing gives channels x samples here.
                rn=I.estimate_noise(b)
                oracle=np.diag(d[t,ch,fi]**2)
                bgerr=np.linalg.norm(rn-oracle,"fro")/np.linalg.norm(oracle,"fro")
                nominal=p[t,ch,fi]
                for ki,k in enumerate(KS):
                    x=y[t,:k,ch,fi].T
                    a=I.response_block(x,b)
                    o=I.response_block(x,b,supplied_oracle=oracle)
                    idx=(ri,ki,t,fi)
                    arrays["before_error"][idx]=I.response_error(a["before"],nominal)
                    ae=I.response_error(a["after"],nominal)
                    oe=I.response_error(o["after"],nominal)
                    arrays["after_error"][idx]=np.nan if ae is None else ae
                    arrays["oracle_error"][idx]=np.nan if oe is None else oe
                    arrays["detected"][idx]=a["detected"]
                    arrays["oracle_detected"][idx]=o["detected"]
                    arrays["background_error"][idx]=bgerr
                    arrays["energy"][idx]=a["energy"]
                    arrays["threshold"][idx]=a["threshold"]
                    arrays["oracle_energy"][idx]=o["energy"]
                    arrays["oracle_threshold"][idx]=o["threshold"]
                    arrays["signal_eigenvalue"][idx]=a["largest_signal_eigenvalue"]
                    arrays["noise_to_min_amplitude"][idx]=float(max(d[t,ch,fi]/np.abs(nominal)))
                    arrays["rank"][idx]=a["sample_rank"]
                    arrays["u_before"][idx][:n]=a["before"]
                    if a["after"] is not None:
                        arrays["u"][idx][:n]=a["after"]
            for pack,fs in PACK.items():
                xall=y[t][:,ch][:,:,fs].reshape(32,n*len(fs))
                ball=bg[t][:,ch][:,:,fs].reshape(256,n*len(fs))
                rn=I.estimate_noise(ball)
                resp,q=candidate_response(xall,n,len(fs),t)
                model=I.joint_model(rn,resp,q)
                # Background likelihood is common to both representations.
                rawb,csdb=I.likelihoods(ball,rn)
                dd=d[t][ch][:,fs].reshape(n*len(fs))
                oracle_full=np.kron(np.eye(n),cf[np.ix_(fs,fs)])*dd[:,None]*dd[None,:]
                full_noise_error=float(np.linalg.norm(rn-oracle_full,"fro")/np.linalg.norm(oracle_full,"fro"))
                for k in KS:
                    raw,csd=I.likelihoods(xall[:k],model)
                    denom=max(1,abs(raw),abs(csd))
                    likelihood.append({"resource":res,"template_s":[0,600,1200][t],
                        "pack":pack,"K":k,"joint_dimension":n*len(fs),
                        "raw_nll":raw,"csd_nll":csd,
                        "difference":raw-csd,"relative_difference":abs(raw-csd)/denom,
                        "background_nll_difference":rawb-csdb,
                        "model_cholesky":"PASS","sample_rank_upper_bound":min(k,n*len(fs)),
                        "sample_rank":int(np.linalg.matrix_rank(xall[:k])),
                        "full_background_relative_error":full_noise_error,
                        "sample_inverse_used":False,"same_nuisance_parameters":True})
    return arrays,likelihood

def cell_summary(arrays,ri,ki,fs):
    def sub(key):
        return arrays[key][ri,ki][:,fs]
    def vmax(key):
        v=sub(key); v=v[np.isfinite(v)]
        return float(v.max()) if len(v) else None
    fails=int((~sub("detected")).sum())
    err=vmax("after_error")
    return {"blocks":3*len(fs),"failed_blocks":fails,
        "oracle_failed_blocks":int((~sub("oracle_detected")).sum()),
        "before_max_projector_error":vmax("before_error"),
        "after_max_projector_error":err,"oracle_max_projector_error":vmax("oracle_error"),
        "background_relative_error_max":vmax("background_error"),
        "max_noise_to_channel_amplitude":vmax("noise_to_min_amplitude"),
        "accepted_quality":fails==0 and err is not None and err<=.1}

def run():
    global START
    START=time.perf_counter()
    design=json.loads((OUT/"DESIGN_FREEZE.json").read_text(encoding="utf-8"))
    verify_inputs(design)
    pre=json.loads((OUT/"PRE_RUN_CONTROLS.json").read_text(encoding="utf-8"))
    if not pre["PASS"]:
        raise RuntimeError("H3_G2E_STATISTICAL_IMPLEMENTATION_INVALID")
    if (OUT/"RUN_ONCE.json").exists():
        raise RuntimeError("FROZEN_EXPERIMENT_ALREADY_STARTED; no rerun")
    write_json(OUT/"RUN_ONCE.json",{"started_utc":__import__("datetime").datetime.now(
        __import__("datetime").timezone.utc).isoformat(),"design_sha":sha(OUT/"DESIGN_FREEZE.json"),
        "repeat_execution_allowed":False})
    (OUT/"random_records").mkdir(exist_ok=True)
    (OUT/"features").mkdir(exist_ok=True)
    cf=gram()
    cells=[]; likes=[]; records=[]; stability=[]
    for si,scene in enumerate(["H01","H06","H12"]):
        ps={g:fields(scene,g) for g in GRIDS}
        for rep in range(16):
            pool=make_pool(si,rep)
            poolfile=OUT/"random_records"/f"{scene}_rep{rep:02d}.npz"
            np.savez(poolfile,**pool)
            for kind in ["RELATIVE","ABSOLUTE_FLOOR"]:
                for sigma in [.01,.05]:
                    base=f"{scene}_{kind}_s{sigma:.2f}_rep{rep:02d}"
                    generated={}
                    extracted={}
                    for grid in GRIDS:
                        y,bg,d=generate(ps[grid],pool,kind,sigma,cf)
                        generated[grid]=(y,bg,d)
                        records.append({"base":base,"grid":grid,
                            "pool_path":str(poolfile.relative_to(ROOT)).replace("\\","/"),
                            "raw_y_sha256":array_sha(y),"observed_background_sha256":array_sha(bg),
                            "nominal_gain":"UNITY_GENERATOR_ONLY","source_power":"ONE_GENERATOR_ONLY"})
                        a,ll=extract_base(y,bg,d,ps[grid],cf)
                        extracted[grid]=a
                        np.savez(OUT/"features"/f"{base}_n{grid}.npz",**a)
                        for r in ll:
                            likes.append({"base":base,"grid":grid,**r})
                        guard()
                    for ri,res in enumerate(CH):
                        for ki,k in enumerate(KS):
                            for pack,fs in PACK.items():
                                main=cell_summary(extracted[160001],ri,ki,fs)
                                sec=cell_summary(extracted[80001],ri,ki,fs)
                                key={"scene":scene,"resource":res,"pack":pack,"K":k,
                                    "noise_model":kind,"sigma":sigma,"realization":rep}
                                cells.append({**key,**main,
                                    "secondary_failed_blocks":sec["failed_blocks"],
                                    "secondary_after_max_projector_error":sec["after_max_projector_error"],
                                    "same_random_pool":True,
                                    "interpretation":"GAIN_DISTORTED_RESPONSE_NOT_CALIBRATED_FIELD"})
                                aa=extracted[160001]; ab=extracted[80001]
                                diffs=[]
                                for t in range(3):
                                    for f in fs:
                                        if aa["detected"][ri,ki,t,f] and ab["detected"][ri,ki,t,f]:
                                            n=len(CH[res])
                                            diffs.append(I.response_error(aa["u"][ri,ki,t,f,:n],
                                                                          ab["u"][ri,ki,t,f,:n]))
                                stability.append({**key,
                                    "max_accepted_projector_grid_difference":max(diffs) if diffs else None,
                                    "detection_disagreements":int(np.sum(
                                        aa["detected"][ri,ki][:,fs]!=ab["detected"][ri,ki][:,fs])),
                                    "nominal_response_grid_error_max":max(
                                        I.response_error(ps[160001][t,CH[res],f],
                                                         ps[80001][t,CH[res],f])
                                        for t in range(3) for f in fs),
                                    "comparison":"PAIRED_IDENTICAL_INNOVATIONS_NOT_NEW_SAMPLES"})
            print(f"{scene} realization {rep+1}/16 elapsed={time.perf_counter()-START:.1f}s",flush=True)
    write_csv(OUT/"EXTRACTION_BY_CELL.csv",cells)
    write_csv(OUT/"RAW_VS_CSD_LIKELIHOOD_AUDIT.csv",likes)
    write_csv(OUT/"STOCHASTIC_RECORD_MANIFEST.csv",records)
    write_csv(OUT/"NUMERICAL_EXTRACTION_STABILITY.csv",stability)
    low=[r for r in cells if r["failed_blocks"] or r["max_noise_to_channel_amplitude"]>=.1]
    write_csv(OUT/"LOW_ENERGY_FAILURES.csv",low)
    lookup={(r["scene"],r["pack"],r["K"],r["noise_model"],r["sigma"],r["realization"],r["resource"]):r for r in cells}
    matched=[]
    for r in cells:
        if r["resource"]!="B1": continue
        key=tuple(r[x] for x in ["scene","pack","K","noise_model","sigma","realization"])
        b=lookup[key+("B2",)]
        ea=r["after_max_projector_error"]; eb=b["after_max_projector_error"]
        matched.append({**{k:r[k] for k in ["scene","pack","K","noise_model","sigma","realization"]},
            "B1_failed_blocks":r["failed_blocks"],"B2_failed_blocks":b["failed_blocks"],
            "B1_max_error":ea,"B2_max_error":eb,
            "B1_minus_B2_error":None if ea is None or eb is None else ea-eb,
            "B1_quality":r["accepted_quality"],"B2_quality":b["accepted_quality"],
            "comparison":"PAIRED_SOURCE_AND_ORIGINAL_EIGHT_NOISE;NOT_INDEPENDENT"})
    write_csv(OUT/"RESOURCE_MATCHED_EXTRACTION.csv",matched)
    write_json(OUT/"EXECUTION_METADATA.json",{"execution_seconds":time.perf_counter()-START,
        "registered_cells":len(cells),"new_independent_random_pools":48,
        "realizations_per_cell":16,"three_template_draws_per_pool":True,
        "K_nested":True,"frequency_packages_paired":True,"noise_conditions_paired":True,
        "two_meshes_paired":True,"new_propagation_calls":0,
        "localization_error_Monte_Carlo":"NOT_EXECUTED"})
    guard()

if __name__=="__main__":
    with threadpool_limits(limits=1):
        run()
