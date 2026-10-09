"""Cold binary parsing, Cartesian/modal summation, independent covariance/chart/QR review."""
import time,json
import numpy as np
from scipy import linalg
import r4_h3_p1 as P
import r4_h3_g0_audit as C
import r4_e2_nine_pair_pilot_audit as B
import r4_h3_p0_derivative_review_audit as A
E=P.E;G=P.G;O=P.O
def field(mods,state,res):
    t=np.array([0.,600.,1200.]);post=np.maximum(t-600.,0)
    main=np.column_stack([2*np.minimum(t,600)+2*post*np.cos(np.pi/12),2*post*np.sin(np.pi/12)])
    x=np.array(state[:2])+t[:,None]*np.array(state[2:])
    dx=np.arange(-7.,8.,2.);zr=np.full(8,200.)
    if res!="B0":dx=np.r_[dx,[-1.]*4];zr=np.r_[zr,[198.,199.,201.,202.] if res=="B1" else [200.]*4]
    delta=x[:,None,:]-main[:,None,:]-np.column_stack([dx,np.zeros(len(dx))])[None,:,:]
    r=np.sqrt(np.sum(delta**2,axis=-1));values=[]
    for z,phi,k in mods:
        ix={float(v):i for i,v in enumerate(z)}
        recv=np.tile(phi[[ix[v] for v in zr]],(3,1)).T
        terms=phi[ix[200.],:,None]*recv*np.sqrt(2*np.pi/(k[:,None]*r.ravel()))*np.exp(-1j*k[:,None]*r.ravel()-1j*np.pi/4)
        values.append(np.sum(terms,axis=0).reshape(r.shape))
    return np.stack(values,-1)
def project(p,j,n):
    y,nu=A.system(p,j,n);return C.qrremove(y,nu)
def compare(p,q):
    raw=P.rel(p,q);uu=p/np.linalg.norm(p,axis=1,keepdims=True);vv=q/np.linalg.norm(q,axis=1,keepdims=True)
    pp=np.einsum("tnf,tmf->tfnm",uu,uu.conj());qq=np.einsum("tnf,tmf->tfnm",vv,vv.conj())
    response=float(np.sqrt(np.mean(np.sum(abs(pp-qq)**2,axis=(-2,-1)))))
    alpha=np.sum(q.conj()*p,axis=1,keepdims=True)/np.sum(abs(q)**2,axis=1,keepdims=True)
    out={"raw_field_relative":raw,"source_removed_field_relative":P.rel(p,alpha*q),"relative_projector_distance":response}
    for n in P.NOISE:
        if n=="RELATIVE":
            weights=abs(p)**2/np.sum(abs(p)**2,axis=1,keepdims=True)
            noise=float(np.sqrt(np.mean(2*.01**2*np.sum(weights*(1-weights),axis=1))))
        else:noise=float(np.sqrt(np.mean(2*.01**2*G.REF**2*(p.shape[1]-1)/np.sum(abs(p)**2,axis=1))))
        out["response_noise1_ratio_"+n]=response/noise
        out["C1b_noise1_RMS_"+n]=float(np.linalg.norm(project(p,(q/p-1)[...,None],n))/(.01*np.sqrt(2*p.size)))
    return out
def main():
    start=time.monotonic();P.verify();checks=[];fields={};mods={}
    def ck(name,error,limit=1e-7):E.ck(checks,name,error,limit)
    calls=G.rows(O/"SOLVER_CALLS.csv");decision=E.read(O/"H3_P1_DECISION.json")
    ck("attempted_cap",max(0,sum(r["attempted"]=="True" for r in calls)-60),0)
    for r in calls:
        if r["status"]!="GENERATED":continue
        key=float(r["CLOW"]),int(r["mesh"]),float(r["frequency"]),float(r["delta_c"])
        m=B.mod(P.path(*key));mods[key]=m
        ck("exact_depths:"+str(key),0 if np.array_equal(m[0],E.M.DEPTHS) else 1,0)
    if len(mods)==60 and (O/"SSP_SECANT_LOCALITY.csv").exists():
        # independent full-field reconstruction for every saved field, including new delta zero.
        for c in P.LOWS:
            for g in P.GRIDS:
                for sc in E.scenes():
                    for res in P.RES:
                        with np.load(O/f"FIELDS_c{int(c)}_n{g}_{sc['scene_id']}_{res}.npz") as z:
                            assert np.array_equal(z["delta_c"],P.DELTA)
                            for i,d in enumerate(P.DELTA):
                                pp=field([mods[c,g,f,d] for f in P.F],sc["state"],res)
                                fields[c,g,sc["scene_id"],res,d]=pp
                                ck("cold_pressure:"+str((c,g,sc["scene_id"],res,d)),P.rel(pp,z["pressure"][i]),1e-9)
                    print("COLD",c,g,sc["scene_id"],flush=True)
                    assert time.monotonic()-start<=900
        for row in G.rows(O/"RELATIVE_RESPONSE_NUMERICAL_STABILITY.csv"):
            sc=row["scene"];res=row["resource"];d=float(row["delta_c"])
            if row["kind"]=="GRID":
                c=float(row["CLOW"]);p=fields[c,80001,sc,res,d];q=fields[c,160001,sc,res,d]
            else:
                g=int(row["mesh"]);p=fields[1490.,g,sc,res,d];q=fields[1480.,g,sc,res,d]
            for key,val in compare(p,q).items():ck("response:"+str((row["kind"],sc,res,row["CLOW"],row["mesh"],d,key)),abs(val-float(row[key])))
        for row in G.rows(O/"SSP_SECANT_LOCALITY.csv"):
            c=float(row["CLOW"]);g=int(row["mesh"]);sc=row["scene"];res=row["resource"];n=row["noise_model"];h=float(row["step"])
            p=fields[c,g,sc,res,0.];cols=[]
            for hh in [.01,.02]:
                hi=fields[c,g,sc,res,hh];lo=fields[c,g,sc,res,-hh]
                cols.extend([(hi-lo)/(2*hh)/p,(hi-p)/hh/p,(p-lo)/hh/p])
            z=project(p,np.stack(cols,-1),n);i=0 if h==.01 else 3
            for key,val in [("two_step_C1b_direction_relative",P.rel(z[:,0],z[:,3])),("forward_backward_C1b_relative",P.rel(z[:,i+1],z[:,i+2]))]:
                ck("locality:"+str((c,g,sc,res,n,h,key)),abs(val-float(row[key])))
        matchrows={(r["kind"],float(r["CLOW_a"]),float(r["CLOW_b"]),int(r["mesh"]),float(r["frequency"]),float(r["delta_a"]),float(r["delta_b"])):r for r in G.rows(O/"MODE_MATCHING.csv")}
        scene={r["scene_id"]:r["state"] for r in E.scenes()}
        for row in G.rows(O/"WINDOW_FIELD_CONTRIBUTION.csv"):
            c=float(row["CLOW_a"]);cc=float(row["CLOW_b"]);g=int(row["mesh"]);f=float(row["frequency"]);d=float(row["delta_a"]);dd=float(row["delta_b"])
            ma=mods[c,g,f,d];mb=mods[cc,g,f,dd];mr=matchrows[row["kind"],c,cc,g,f,d,dd]
            for side,m in [("a",ma),("b",mb)]:
                ix=np.array(json.loads(mr["unmatched_mode_ids_"+side]),dtype=int)
                denom=field([m],scene[row["scene"]],row["resource"])
                z,phi,k=m
                numerator=field([(z,phi[:,ix],k[ix])],scene[row["scene"]],row["resource"])
                value=float(np.linalg.norm(numerator)/max(np.linalg.norm(denom),1e-30))
                ck("unmatched:"+str((row["kind"],c,cc,g,f,d,dd,row["scene"],row["resource"],side)),abs(value-float(row["unmatched_field_fraction_"+side])))
        for row in G.rows(O/"MODE_BOUNDARY_AND_CONTINUITY.csv"):
            key=float(row["CLOW"]),int(row["mesh"]),float(row["frequency"]),float(row["delta_c"]);m=mods[key]
            speed=2*np.pi*key[2]/m[2].real
            ck("mode_count:"+str(key),abs(len(m[2])-int(row["modes"])),0)
            ck("mode_low:"+str(key),abs(float(speed.min())-float(row["phase_speed_min"])),1e-10)
            ck("mode_high:"+str(key),abs(float(speed.max())-float(row["phase_speed_max"])),1e-10)
        registered=G.rows(O/"PREREGISTERED_CHECKS.csv")
        window=all(r["PASS"]=="True" for r in registered if not r["check"].startswith("LOCAL_"))
        locality=all(r["PASS"]=="True" for r in registered if r["check"].startswith("LOCAL_"))
        ck("classification_window",int(window!=decision["window_numerical_support"]),0)
        ck("classification_locality",int(locality!=decision["SSP_secant_locality"]),0)
    else:
        ck("provider_failure_scope",int("H3_P1_PROVIDER_OR_FORMULA_INVALID" not in decision["classifications"]),0)
    # Repeat the complete frozen environment contract checks, before interpreting any propagation.
    for c,g,f,d in P.configs():
        s,cl,ch,opt=P.R.environment(P.path(c,g,f,d).with_suffix(".env"))
        base,_,_,_=P.R.environment(E.modpath(f,g,0).with_suffix(".env"))
        ck("SSP_offset:"+str((c,g,f,d)),float(np.max(abs(s[:,1]-base[:,1]-d))),1e-10)
        ck("depth_positions:"+str((c,g,f,d)),float(np.max(abs(s[:,0]-base[:,0]))),0)
        ck("common_window:"+str((c,g,f,d)),abs(cl-c)+abs(ch-1800),0)
    P.verify()
    E.write(O/"INDEPENDENT_VALIDATION_CHECKS.csv",checks)
    value={"checks":len(checks),"PASS":sum(r["PASS"] for r in checks),"FAIL":sum(not r["PASS"] for r in checks),
        "cold_pressure_relative_max":max((r["value"] for r in checks if r["check"].startswith("cold_pressure")),default=0.),
        "reconstruction":"SEPARATE_BINARY_PARSER_CARTESIAN_TERMWISE_SUM_NUISANCE_CHART_PIVOTED_QR",
        "new_KRAKEN":0,"new_FIELD":0,"new_MC":0,"elapsed_s":time.monotonic()-start}
    E.dump(O/"VALIDATION.json",value)
    decision["independent_review"]="PASS" if value["FAIL"]==0 else "FAIL"
    if value["FAIL"]:decision["classifications"]=["H3_P1_PROVIDER_OR_FORMULA_INVALID"]
    E.dump(O/"H3_P1_DECISION.json",decision)
    assert time.monotonic()-start<=900,"AUDIT_CAP"
if __name__=="__main__":main()
