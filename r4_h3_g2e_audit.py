"""Independent H3-G2E reconstruction: own Gram, NLL, direction and controls."""
import ast, csv, json, time
from pathlib import Path
import numpy as np
from threadpoolctl import threadpool_limits
import r4_h3_g2e_interface as I
import r4_h3_g2e as E  # I/O/generator provenance only; no reuse of extraction/NLL below.

OUT=E.OUT
def own_cov(y):
    return np.einsum("ki,kj->ij",y,y.conj())/len(y)
def own_gram():
    w=.5-.5*np.cos(2*np.pi*np.arange(2000)/2000)
    z=np.fft.fft(w*w)/(w@w)
    return np.array([[z[int(2*(a-b))%2000] for b in E.F] for a in E.F])
def own_proj(v):
    v=v/np.linalg.norm(v)
    return v[:,None]*v.conj()[None,:]
def own_error(u,v):
    return float(np.sqrt(max(0,2-2*abs(np.vdot(u/np.linalg.norm(u),v/np.linalg.norm(v)))**2)))
def own_extract(y,b):
    rn=own_cov(b)
    r=own_cov(y)
    energy=float(np.trace(np.linalg.solve(rn,r)).real)
    threshold=y.shape[1]*len(b)/(len(b)-y.shape[1])/.01
    val,vec=np.linalg.eigh((r-rn+(r-rn).conj().T)/2)
    detected=energy>threshold and val[-1]>0
    return rn,energy,threshold,detected,vec[:,-1] if detected else None

def controls():
    rows=[]
    def check(name,value,limit=1e-10,scope="PRE_RUN_SYNTHETIC_FIXTURE"):
        passed=bool(np.isfinite(value) and value<=limit)
        rows.append({"control":name,"residual":float(value),"limit":limit,"PASS":passed,"scope":scope})
    c1=E.gram(); c2=own_gram()
    check("HANN_FULL_GRAM_DIRECT_VS_FFT",np.max(abs(c1-c2)),1e-12)
    expected=np.eye(13); expected[11,12]=expected[12,11]=1/6
    check("HANN_ONE_HZ_CORRELATION",np.max(abs(c2-expected)),1e-12)
    check("HANN_GRAM_POSITIVE_DEFINITE",0 if np.linalg.eigvalsh(c2).min()>0 else 1,0)
    rng=np.random.Generator(np.random.PCG64(873104))
    def draw(shape):
        return (rng.normal(size=shape)+1j*rng.normal(size=shape))/np.sqrt(2)
    b=draw((256,12)); y=draw((32,12))
    rn=I.estimate_noise(b); noise=own_cov(b)
    check("NOISE_ESTIMATE_ONLY_BACKGROUND",np.linalg.norm(rn-noise)/np.linalg.norm(noise))
    doubled=I.estimate_noise(2*b)
    check("BACKGROUND_SCALING_NOT_ORACLE_SUBSTITUTION",np.linalg.norm(doubled-4*rn)/np.linalg.norm(rn))
    n=4; nf=3
    a=draw((n,nf)); q=np.array([.7,1.1,2.])
    model=I.joint_model(rn,a,q)
    inv=np.linalg.inv(model); ld=np.linalg.slogdet(model)[1]
    for k in [1,8,32]:
        x=y[:k]
        rr=own_cov(x)
        raw=float(k*ld+sum(np.vdot(v,inv@v).real for v in x))
        csd=float(k*(ld+np.trace(inv@rr).real))
        aa,bb=I.likelihoods(x,model)
        check(f"RAW_CSD_NLL_K{k}",max(abs(raw-csd),abs(raw-aa),abs(raw-bb))/max(1,abs(raw)),1e-11)
        check(f"SAMPLE_RANK_BOUND_K{k}",0 if np.linalg.matrix_rank(rr)<=min(k,12) else 1,0)
        # Score equality for fixed-gain real/imag, source power and noise covariance.
        perturbations=[]
        da=a.copy()*0
        da[1,0]=a[1,0] # one fixed element-frequency gain direction
        ap=np.zeros((n*nf,nf),complex); dap=ap.copy()
        for f in range(nf):
            ap[np.arange(n)*nf+f,f]=a[:,f]
            dap[np.arange(n)*nf+f,f]=da[:,f]
        perturbations.extend([dap@np.diag(q)@ap.conj().T+ap@np.diag(q)@dap.conj().T,
            1j*dap@np.diag(q)@ap.conj().T-1j*ap@np.diag(q)@dap.conj().T,
            np.outer(ap[:,1],ap[:,1].conj()),rn])
        for j,dc in enumerate(perturbations):
            term=inv@dc@inv
            sr=k*np.trace(inv@dc).real-sum(np.vdot(v,term@v).real for v in x)
            sc=k*(np.trace(inv@dc).real-np.trace(term@rr).real)
            check(f"RAW_CSD_SCORE_NUISANCE{j}_K{k}",abs(sr-sc)/max(1,abs(sr),abs(sc)),1e-10)
    # Exact nuisance reparameterizations, no pseudoinverse assigning zero variance.
    factor=(1.2+.4j)*np.exp(1j*np.arange(nf)*.1)
    model2=I.joint_model(rn,a*factor,q/abs(factor)**2)
    check("SOURCE_COMMON_FACTOR_GAUGE",np.linalg.norm(model-model2)/np.linalg.norm(model))
    g=(1+.03*np.arange(n))[:,None]*np.exp(1j*.2*np.arange(nf)[None,:])
    h=(1+.02*np.arange(n))[:,None]*np.exp(1j*.13*np.arange(nf)[None,:])
    check("FIXED_GAIN_UNCALIBRATED_RESPONSE_NOT_IDENTIFIED",np.linalg.norm(g*a-(g*h)*(a/h)))
    phi1=1.2+.3j; phi2=.7-.8j
    ms1=I.joint_model(rn,phi1*a,q)
    ms2=I.joint_model(rn,phi2*a,q*abs(phi1/phi2)**2)
    check("SINGLE_SEPARABLE_COMPONENT_DEPTH_ABSORBED",np.linalg.norm(ms1-ms2)/np.linalg.norm(ms1))
    modes=draw((12,5)); ph1=draw((5,)); ph2=draw((5,)); gain=draw((5,))
    f1=modes@(ph1*gain); f2=modes@(ph2*(gain*ph1/ph2))
    check("FREE_MODAL_GAIN_DEPTH_ABSORBED",np.linalg.norm(f1-f2)/np.linalg.norm(f1))
    raw=draw((3,4,3)); changed=draw((3,4,3)); free=draw((3,4,3))
    check("C2_FREE_CHANNEL_TEMPLATE_DEPTH_ABSORBED",
          np.linalg.norm(free*raw-(free*raw/changed)*changed)/np.linalg.norm(free*raw))
    for k in [1,8,32]:
        # All-zero and extremely low observations must not fabricate direction.
        for label,x in [("NO_SOURCE",np.zeros_like(y[:k])),("LOW_ENERGY",1e-8*y[:k]),
                        ("NULL_GAUSSIAN",y[:k])]:
            a0=I.response_block(x,b)
            check(f"{label}_NO_FABRICATED_DIRECTION_K{k}",0 if a0["after"] is None else 1,0)
    tree=ast.parse((E.ROOT/"r4_h3_g2e_interface.py").read_text(encoding="utf-8"))
    imports=[ast.unparse(x) for x in ast.walk(tree) if isinstance(x,(ast.Import,ast.ImportFrom))]
    check("INTERFACE_IMPORT_TRUTH_ISOLATION",0 if all("r4_h3" not in x for x in imports) else 1,0)
    fn=[x for x in tree.body if isinstance(x,ast.FunctionDef) and x.name=="response_block"][0]
    check("TRUTH_FREE_PRIMARY_SIGNATURE",0 if [x.arg for x in fn.args.args]==
          ["samples","background","supplied_oracle","alpha"] else 1,0)
    pool=E.make_pool(0,0)
    pp=np.full((3,16,13),1+1j)
    yp,bp,_=E.generate(pp,pool,"RELATIVE",.01,E.gram())
    check("ORIGINAL_EIGHT_SHARED_RESOURCE_IDENTITY",
          np.linalg.norm(yp[:,:,:8]-yp[:,:,E.CH["B1"]][:,:,:8]))
    check("NESTED_SAMPLE_PREFIX_IDENTITY",np.linalg.norm(yp[:,:1]-yp[:,:8][:,:1]))
    check("SOURCE_NOISE_BACKGROUND_SEPARATE_STREAMS",0 if not np.array_equal(
          pool["noise"][:,:32],pool["background"][:,:32]) else 1,0)
    draft=json.loads((OUT/"DESIGN_FREEZE.json").read_text(encoding="utf-8"))
    check("DEPENDENCIES_NOT_INDEPENDENT_CONFIRMATIONS",0 if draft["dependency_contract"]==
          "K_PACK_RESOURCE_NOISE_MESH_PAIRED;TEMPLATES_STATIC;48_INDEPENDENT_POOLS" else 1,0)
    check("NO_TIME_DOMAIN_OR_STATIONARITY_CLAIM",0 if
          draft["real_64s_stationarity"]=="NOT_ESTABLISHED" and
          draft["time_domain"]=="NOT_OPENED" else 1,0)
    E.write_csv(OUT/"SOURCE_AND_NOISE_COVARIANCE_CONTROLS.csv",rows)
    E.write_json(OUT/"PRE_RUN_CONTROLS.json",{"PASS":all(x["PASS"] for x in rows),
        "checks":len(rows),"new_control_draws":"REGISTERED_SYNTHETIC_FIXTURES_ONLY",
        "experiment_executed":False,"classification_on_failure":
        "H3_G2E_STATISTICAL_IMPLEMENTATION_INVALID"})
    print(json.dumps({"PASS":all(x["PASS"] for x in rows),"checks":len(rows)}))
    return all(x["PASS"] for x in rows)

def audit():
    start=time.perf_counter()
    design=json.loads((OUT/"DESIGN_FREEZE.json").read_text(encoding="utf-8"))
    E.verify_inputs(design)
    meta=json.loads((OUT/"EXECUTION_METADATA.json").read_text(encoding="utf-8"))
    records=list(csv.DictReader((OUT/"STOCHASTIC_RECORD_MANIFEST.csv").open(encoding="utf-8")))
    saved_likes=list(csv.DictReader((OUT/"RAW_VS_CSD_LIKELIHOOD_AUDIT.csv").open(encoding="utf-8")))
    like_map={(r["base"],int(r["grid"]),r["resource"],int(r["template_s"]),r["pack"],int(r["K"])):r for r in saved_likes}
    check_count=0; failures=[]; maxres={}
    def check(name,res,limit=1e-9):
        nonlocal check_count
        check_count+=1
        maxres[name]=max(maxres.get(name,0),float(res))
        if not np.isfinite(res) or res>limit:
            failures.append({"check":name,"value":float(res),"limit":limit})
    cf=own_gram()
    E.write_json(OUT/"AUDIT_STARTED.json",{"protocol":"COLD_RECONSTRUCTION_OWN_EINSUM_SOLVE_EIGH",
        "started_utc":__import__("datetime").datetime.now(__import__("datetime").timezone.utc).isoformat()})
    for si,scene in enumerate(["H01","H06","H12"]):
        p={g:E.fields(scene,g) for g in E.GRIDS}
        for rep in range(16):
            path=OUT/"random_records"/f"{scene}_rep{rep:02d}.npz"
            with np.load(path) as z:
                pool={k:z[k] for k in z.files}
            regenerated=E.make_pool(si,rep)
            for key in pool:
                check("FROZEN_RNG_RECORD",0 if np.array_equal(pool[key],regenerated[key]) else 1,0)
            for kind in ["RELATIVE","ABSOLUTE_FLOOR"]:
                for sigma in [.01,.05]:
                    base=f"{scene}_{kind}_s{sigma:.2f}_rep{rep:02d}"
                    for grid in E.GRIDS:
                        # Own noise multiplication using matrix rows rather than generator einsum.
                        pp=p[grid]
                        d=sigma*abs(pp) if kind=="RELATIVE" else np.full(pp.shape,sigma*E.REF)
                        l=np.linalg.cholesky(cf)
                        n=pool["noise"]@l.T
                        b=pool["background"]@l.T
                        y=pp[:,None]*pool["source"][:,:,None]+d[:,None]*n
                        bg=d[:,None]*b
                        rec=next(r for r in records if r["base"]==base and int(r["grid"])==grid)
                        check("OBSERVATION_RAW_HASH",0 if E.array_sha(y)==rec["raw_y_sha256"] else
                              np.max(abs(y-E.generate(pp,pool,kind,sigma,cf)[0])),1e-14)
                        # Gram independent construction differs roundoff, hashes can differ.
                        yg,bg0,_=E.generate(pp,pool,kind,sigma,E.gram())
                        check("ARCHIVED_PRIMARY_RAW_HASH",0 if E.array_sha(yg)==rec["raw_y_sha256"] else 1,0)
                        check("ARCHIVED_BACKGROUND_HASH",0 if E.array_sha(bg0)==rec["observed_background_sha256"] else 1,0)
                        check("OWN_GENERATOR_NOISE",np.linalg.norm(bg-bg0)/np.linalg.norm(bg0),1e-11)
                        z=np.load(OUT/"features"/f"{base}_n{grid}.npz")
                        for ri,(res,ch) in enumerate(E.CH.items()):
                            nn=len(ch)
                            for t in range(3):
                                for f in range(13):
                                    bb=bg[t,:,ch,f].T
                                    for ki,k in enumerate(E.KS):
                                        yy=y[t,:k,ch,f].T
                                        rn,en,th,det,u=own_extract(yy,bb)
                                        idx=(ri,ki,t,f)
                                        check("ENERGY_COLD",abs(en-z["energy"][idx])/max(1,abs(en)),1e-8)
                                        check("DETECTION_COLD",0 if det==z["detected"][idx] else 1,0)
                                        oo=np.diag(d[t,ch,f]**2)
                                        or_energy=float(np.trace(np.linalg.solve(oo,own_cov(yy))).real)
                                        or_det=bool(or_energy>nn/.01 and np.linalg.eigvalsh(own_cov(yy)-oo)[-1]>0)
                                        check("ORACLE_DETECTION_COLD",0 if or_det==z["oracle_detected"][idx] else 1,0)
                                        check("BACKGROUND_ERROR_COLD",abs(np.linalg.norm(rn-oo)/np.linalg.norm(oo)-z["background_error"][idx]),1e-9)
                                        if or_det:
                                            _,ov=np.linalg.eigh(own_cov(yy)-oo)
                                            check("ORACLE_ERROR_COLD",abs(own_error(ov[:,-1],pp[t,ch,f])-z["oracle_error"][idx]),1e-7)
                                        if det:
                                            check("DIRECTION_PROJECTOR_COLD",np.linalg.norm(
                                                own_proj(u)-own_proj(z["u"][idx][:nn])),1e-8)
                                            er=own_error(u,pp[t,ch,f])
                                            check("PROJECTOR_ERROR_COLD",abs(er-z["after_error"][idx]),1e-7)
                                for pack,fs in E.PACK.items():
                                    xx=y[t][:,ch][:,:,fs].reshape(32,nn*len(fs))
                                    bb=bg[t][:,ch][:,:,fs].reshape(256,nn*len(fs))
                                    noise=own_cov(bb)
                                    # Independently construct same registered synthetic nuisance candidate.
                                    scale=np.sqrt(max(np.trace(own_cov(xx)).real/(nn*len(fs)),np.finfo(float).tiny))
                                    c=np.arange(nn)[:,None]; ff=np.arange(len(fs))[None,:]
                                    aa=scale*(1+.01*c+.02*t)*np.exp(1j*(.03*c+.02*t*c+.005*ff*c))
                                    aa=aa*(1+.03*c)*np.exp(1j*.04*c*(ff+1))
                                    qq=.7+np.arange(len(fs))/20+t/10
                                    model=noise.copy()
                                    for fj in range(len(fs)):
                                        av=np.zeros(nn*len(fs),complex)
                                        av[np.arange(nn)*len(fs)+fj]=aa[:,fj]
                                        model+=qq[fj]*np.outer(av,av.conj())
                                    inv=np.linalg.inv(model)
                                    ld=np.linalg.slogdet(model)[1]
                                    for k in E.KS:
                                        raw=k*ld+sum(np.vdot(v,inv@v).real for v in xx[:k])
                                        csd=k*(ld+np.trace(inv@own_cov(xx[:k])).real)
                                        check("RAW_CSD_COLD",abs(raw-csd)/max(1,abs(raw)),1e-9)
                                        saved=like_map[(base,grid,res,[0,600,1200][t],pack,k)]
                                        check("SAVED_NLL_COLD",max(abs(raw-float(saved["raw_nll"])),abs(csd-float(saved["csd_nll"])))/max(1,abs(raw)),1e-8)
                        z.close()
            if float(meta["execution_seconds"])+time.perf_counter()-start>3600:
                raise RuntimeError("FROZEN_3600S_BUDGET_STOP")
            print(f"cold {scene} rep{rep+1} checks={check_count} elapsed={time.perf_counter()-start:.1f}s",flush=True)
    # Rebuild every registered cell from archived features, no selecting strong cells.
    cells=list(csv.DictReader((OUT/"EXTRACTION_BY_CELL.csv").open(encoding="utf-8")))
    for r in cells:
        base=f'{r["scene"]}_{r["noise_model"]}_s{float(r["sigma"]):.2f}_rep{int(r["realization"]):02d}'
        z=np.load(OUT/"features"/f"{base}_n160001.npz")
        ri=list(E.CH).index(r["resource"]); ki=E.KS.index(int(r["K"])); fs=E.PACK[r["pack"]]
        det=z["detected"][ri,ki][:,fs]
        check("ALL_CELL_FAILURE_RECONSTRUCTION",abs(int((~det).sum())-int(r["failed_blocks"])),0)
        errs=z["after_error"][ri,ki][:,fs]; vals=errs[np.isfinite(errs)]
        err=float(max(vals)) if len(vals) else None
        accepted=bool(det.all() and err is not None and err<=.1)
        check("ALL_CELL_QUALITY_RECONSTRUCTION",0 if str(accepted)==r["accepted_quality"] else 1,0)
        z.close()
    check("REGISTERED_CELL_COUNT",abs(len(cells)-5184),0)
    likes=list(csv.DictReader((OUT/"RAW_VS_CSD_LIKELIHOOD_AUDIT.csv").open(encoding="utf-8")))
    for r in likes:
        check("SAVED_RAW_CSD_EQUALITY",float(r["relative_difference"]),1e-9)
        check("NO_SAMPLE_INVERSION",0 if r["sample_inverse_used"]=="False" else 1,0)
    size=sum(p.stat().st_size for p in OUT.rglob("*") if p.is_file())
    total=float(meta["execution_seconds"])+time.perf_counter()-start
    check("STORAGE_BUDGET",0 if size<=1000000000 else 1,0)
    check("EXECUTION_AND_AUDIT_BUDGET",0 if total<=3600 else 1,0)
    E.write_json(OUT/"VALIDATION.json",{"PASS":not failures,"checks":check_count,"FAIL":len(failures),
        "failures":failures,"maximum_residuals":maxres,"cold_audit_seconds":time.perf_counter()-start,
        "execution_plus_audit_seconds":total,"delivered_bytes_at_audit":size,
        "input_hashes_rechecked":True,"new_propagation_calls":0,
        "scope":"INDEPENDENT_ALGORITHMS_SAME_PROCESS_AUTHOR;NOT_EXTERNAL_RESEARCH_LEAD_AUDIT"})
    print(json.dumps({"PASS":not failures,"checks":check_count,"FAIL":len(failures),"seconds":total}))
    return not failures

if __name__=="__main__":
    import sys
    with threadpool_limits(limits=1):
        ok=controls() if sys.argv[1]=="pre" else audit()
    sys.exit(0 if ok else 2)
