"""Independent cold modal parse/QR reconstruction, no forward solver execution."""
import json,time
from pathlib import Path
import numpy as np
from scipy import linalg
import r4_h3_p0 as E
import r4_h3_g0 as G
import r4_h3_g0_audit as C
import r4_e2_nine_pair_pilot_audit as B
def audit():
    start=time.monotonic();E.verify(); checks=[]
    def ck(name,v,lim=1e-8): E.ck(checks,name,v,lim)
    decision=E.read(E.O/"H3_P0_DECISION.json")
    reg=G.rows(E.O/"NOMINAL_REGRESSION_AND_CHAIN_CONTROLS.csv")
    num=G.rows(E.O/"NUMERICAL_AND_FORMULA_GATES.csv") if (E.O/"NUMERICAL_AND_FORMULA_GATES.csv").exists() else []
    admitted=decision["provider_complete"] and all(r["PASS"]=="True" for r in reg+num)
    ck("admission_not_rewritten",0 if admitted==decision["science_information_evaluated"] else 1,0)
    logs=G.rows(E.O/"SOLVER_CALLS.csv") if (E.O/"SOLVER_CALLS.csv").exists() else []
    ck("hard_call_cap",0 if len(logs)<=52 and len(logs)==decision["new_KRAKEN_calls"] else 1,0)
    for row in logs:
        path=E.modpath(float(row["frequency_hz"]),int(row["mesh"]),int(row["offset_mps"]))
        if row["status"]=="GENERATED": ck("modal_SHA_"+path.stem,0 if E.sha(path)==row["modal_sha256"] else 1,0)
        ck("per_call90_"+path.stem,0 if float(row["elapsed_s"])<92 else 1,0)
    # Every environment has exactly the inherited inputs except SSP sound speed.
    for d in [-1,1]:
        for grid in E.GRIDS:
            for f in E.F:
                nominal=E.P.path_for(f,grid).with_suffix(".env").read_text(encoding="utf-8").splitlines()
                shifted=E.modpath(f,grid,d).with_suffix(".env").read_text(encoding="utf-8").splitlines()
                ck("environment_line_count",0 if len(nominal)==len(shifted) else 1,0)
                for a,b in zip(nominal,shifted):
                    x=a.split();y=b.split()
                    if len(x)==6:
                        try:
                            zz=float(x[0]);cc=float(x[1])
                            ssp=0<=zz<=5000 and cc>1000
                        except ValueError:ssp=False
                    else:ssp=False
                    if ssp:
                        ck("offset_exact",abs(float(y[1])-float(x[1])-d),1e-10)
                        ck("other_SSP_properties_same",0 if x[:1]+x[2:]==y[:1]+y[2:] else 1,0)
                    else:ck("non_SSP_inputs_same",0 if a==b else 1,0)
    align={(r["scene"],r["resource"],int(r["mesh"]),r["noise_model"]):r
           for r in G.rows(E.O/"DEPTH_SSP_DERIVATIVE_ALIGNMENT.csv")} if decision["provider_complete"] else {}
    infos=G.rows(E.O/"PROFILED_DEPTH_INFORMATION.csv")
    lookup={(r["scene"],r["resource"],int(r["mesh"]),r["noise_model"],float(r["sigma"])):r for r in infos}
    for grid in E.GRIDS:
        nominal=[B.mod(E.modpath(f,grid,0)) for f in E.F]
        if decision["provider_complete"]:
            plus=[B.mod(E.modpath(f,grid,1)) for f in E.F]
            minus=[B.mod(E.modpath(f,grid,-1)) for f in E.F]
        for sc in E.scenes():
            for res in E.RES:
                p,j,h=C.cold(nominal,sc["state"],res);j=j/p[...,None]
                if decision["provider_complete"]:
                    qp,_,_=C.cold(plus,sc["state"],res)
                    qm,_,_=C.cold(minus,sc["state"],res)
                    ep=(qp-qm)/(2*p)
                    aug=np.concatenate([j,ep[...,None]],axis=-1)
                for noise in E.NOISE:
                    y,n=G.system(p,j,noise,"C1b")
                    z=C.qrremove(y,n)
                    if decision["provider_complete"]:
                        ya,na=G.system(p,aug,noise,"C1b")
                        za=C.qrremove(ya,na)
                        ss=np.stack([(qp-p)/p,(p-qm)/p],axis=-1)
                        yf,nf=G.system(p,ss,noise,"C1b");zf=C.qrremove(yf,nf)
                        loc=np.linalg.norm(zf[:,0]-zf[:,1])/np.linalg.norm(za[:,5])
                        row=align[sc["scene_id"],res,grid,noise]
                        ck("secant_locality_cold",abs(loc-float(row["positive_negative_secant_disagreement"]))/max(1,abs(loc)),1e-7)
                        depth=za[:,2]/G.SCALE[2];env=za[:,5]
                        cosine=abs(depth@env)/(np.linalg.norm(depth)*np.linalg.norm(env))
                        ck("secant_alignment_cold",abs(cosine-float(row["absolute_cosine_after_gain_source"])),1e-7)
                    for sig in E.SIGMAS:
                        r=lookup[sc["scene_id"],res,grid,noise,sig]
                        w=z/sig
                        rem=C.qrremove(w[:,2:3],w[:,[0,1,3,4]])
                        iz=float(np.sum(rem**2))/G.SCALE[2]**2
                        ck("nominal_Iz_QR_cold",abs(iz-float(r["C1b_depth_information"]))/max(iz,1e-30),1e-6)
                        if admitted:
                            prof=C.qrremove(za[:,:5]/sig,za[:,5:6]/sig)
                            rem=C.qrremove(prof[:,2:3],prof[:,[0,1,3,4]])
                            epiz=float(np.sum(rem**2))/G.SCALE[2]**2
                            ck("SSP_profile_Iz_QR_cold",abs(epiz-float(r["C1e_P0_depth_information"]))/max(epiz,1e-20),1e-5)
                        else:ck("blocked_science_not_filled",0 if not r["C1e_P0_depth_information"] and not r["retention_ratio"] else 1,0)
            print(f"COLD n={grid} {sc['scene_id']}",flush=True)
    ck("information_row_count",abs(len(infos)-72),0)
    if not admitted:ck("blocked_classification",0 if decision["decision"]=="H3_P0_NUMERICAL_OR_FORMULA_INCOMPLETE" else 1,0)
    ck("no_project_credit",decision["R4_percent"],0)
    # Numerical Gate failures are scientific admission blockers, not reconstruction bugs.
    duration=decision["elapsed_s"]+time.monotonic()-start
    ck("total_5400s_cap",0 if duration<=5400 else 1,0)
    size=sum(p.stat().st_size for p in E.O.rglob("*") if p.is_file())
    ck("256MB_cap",0 if size<=256000000 else 1,0)
    G.write(E.O/"INDEPENDENT_RECONSTRUCTION_CHECKS.csv",checks)
    E.dump(E.O/"VALIDATION.json",{"PASS":all(x["PASS"] for x in checks),
        "checks":len(checks),"FAIL":sum(not x["PASS"] for x in checks),
        "max_normalized_residual":max(x["value"]/max(x["limit"],1e-300) for x in checks if x["limit"]>0),
        "execution_plus_cold_audit_s":duration,"delivered_bytes":size,
        "new_calls_during_audit":0,"new_random_draws":0,
        "numerical_admission":"PASS" if admitted else "FAIL_RETAINED",
        "scope":"INTERNAL_INDEPENDENT_PARSE_CARTESIAN_SUM_QR;NOT_EXTERNAL_CERTIFICATION"})
    if not all(x["PASS"] for x in checks):
        decision.update(decision="H3_P0_NUMERICAL_OR_FORMULA_INCOMPLETE",science_information_evaluated=False,
                        reconstruction="FAIL",next_stage="NOT_AUTHORIZED;STOP")
        E.dump(E.O/"H3_P0_DECISION.json",decision)
    print(json.dumps(E.read(E.O/"VALIDATION.json"),indent=2))
    return all(x["PASS"] for x in checks)
if __name__=="__main__":
    raise SystemExit(0 if audit() else 2)
