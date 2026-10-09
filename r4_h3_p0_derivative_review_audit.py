"""Independent termwise sums, separate nuisance chart, pivoted QR and saved-row reconstruction."""
import time,json
import numpy as np
from scipy import linalg
import r4_h3_p0_derivative_review as R
import r4_h3_g0_audit as C
import r4_e2_nine_pair_pilot_audit as B
E=R.E;G=R.G;O=R.O
def system(p,j,noise):
    T,N,F=p.shape;K=2*N*F
    yy=[];nn=[]
    for t in range(T):
        weight=np.full(K,np.sqrt(2.)) if noise=="RELATIVE" else np.tile(np.sqrt(2.)*abs(p[t]).reshape(-1)/G.REF,2)
        y=np.concatenate([j[t].real.reshape(N*F,-1),j[t].imag.reshape(N*F,-1)],axis=0)
        n=np.zeros((K,2*T*F+K))
        for c in range(N):
            for f in range(F):
                for part in [0,1]:
                    row=part*N*F+c*F+f
                    n[row,2*t*F+part*F+f]=1
                    n[row,2*T*F+row]=1
        yy.append(y*weight[:,None]);nn.append(n*weight[:,None])
    return np.concatenate(yy),np.concatenate(nn)
def pressure(mods,st,res):
    tt=np.array([0.,600.,1200.]);post=np.maximum(tt-600,0)
    target=np.array(st[:2])[None,:]+tt[:,None]*np.array(st[2:])
    main=np.column_stack([2*np.minimum(tt,600)+2*post*np.cos(np.pi/12),2*post*np.sin(np.pi/12)])
    offsets=np.arange(-7.,8.,2.);zr=np.full(8,200.)
    if res!="B0":
        offsets=np.r_[offsets,[-1.]*4];zr=np.r_[zr,[198.,199.,201.,202.] if res=="B1" else [200.]*4]
    delta=target[:,None,:]-main[:,None,:]-np.column_stack([offsets,np.zeros(len(offsets))])[None,:,:]
    rr=np.sqrt(np.sum(delta**2,axis=-1))
    fields=[];derivs=[]
    for z,phi,k in mods:
        ix={float(v):i for i,v in enumerate(z)}
        recv=np.tile(phi[[ix[v] for v in zr]],(3,1)).T
        r=rr.ravel()
        prop=recv*np.sqrt(2*np.pi/(k[:,None]*r))*np.exp(-1j*k[:,None]*r-1j*np.pi/4)
        src=phi[[ix[v] for v in R.DEPTH]]
        pp=np.array([np.sum(v[:,None]*prop,axis=0) for v in src])
        d05=src[4]-src[2];d1=(src[5]-src[1])/2;d2=(src[6]-src[0])/4
        dd=np.array([np.sum(v[:,None]*prop,axis=0) for v in [d05,d1,d2,(4*d05-d1)/3,(4*d1-d2)/3]])
        fields.append(pp.reshape(7,*rr.shape));derivs.append(dd.reshape(5,*rr.shape))
    return np.stack(fields,-1),np.stack(derivs,-1)
def main():
    start=time.monotonic();R.verify();checks=[];saved={};maxp=maxd=maxi=0.
    def ck(name,v,limit=1e-7):
        E.ck(checks,name,v,limit)
    info=G.rows(O/"EFFECTIVE_INFORMATION_SENSITIVITY.csv")
    idx={(r["scene"],r["resource"],int(r["mesh"]),r["noise_model"],float(r["sigma"]),r["method"]):r for r in info}
    cmps=G.rows(O/"STENCIL_COMPARISONS.csv");ci={(r["scene"],r["resource"],int(r["mesh"]),r["noise_model"],r["projection_stage"],r["method_a"],r["method_b"]):r for r in cmps}
    laws=G.rows(O/"TRUNCATION_LAW_DIAGNOSTIC.csv");li={(r["scene"],r["resource"],int(r["mesh"]),r["noise_model"],r["projection_stage"]):r for r in laws}
    oldstep={r["check"]:r for r in G.rows(E.O/"NOMINAL_REGRESSION_AND_CHAIN_CONTROLS.csv") if r["check"].startswith("source_depth_steps")}
    step=[]
    for g in E.GRIDS:
        mods=[B.mod(E.modpath(f,g,0)) for f in E.F]
        for sc in E.scenes():
            for res in E.RES:
                ps,ds=pressure(mods,sc["state"],res);p=ps[3]
                with np.load(O/f"STENCILS_{sc['scene_id']}_{res}_n{g}.npz") as z:
                    pe=R.rel(ps,z["pressure"]);de=R.rel(ds,z["derivatives"])
                    ck("cold_pressure:"+str((sc["scene_id"],res,g)),pe,1e-9)
                    ck("cold_derivative:"+str((sc["scene_id"],res,g)),de,1e-9)
                    maxp=max(maxp,pe);maxd=max(maxd,de)
                cp,cj,ch=C.cold(mods,sc["state"],res)
                ck("cold_center:"+str((sc["scene_id"],res,g)),R.rel(cp,p),1e-9)
                for noise in E.NOISE:
                    y,n=system(p,np.moveaxis(ds,0,-1)/p[...,None],noise);dz=C.qrremove(y,n)
                    y,n=system(p,(cj/p[...,None])[...,[0,1,3,4]],noise);h=C.qrremove(y,n)
                    eff=C.qrremove(dz,h);saved[sc["scene_id"],res,g,noise]=(dz,eff)
                    val=np.linalg.norm(dz[:,1]-dz[:,0])/np.linalg.norm(dz[:,0])
                    key=f"source_depth_steps_{sc['scene_id']}_{res}_{g}_{noise}"
                    ck("old_2pct_reconstruction:"+key,abs(val-float(oldstep[key]["value"])),1e-9)
                    step.append({"check":key,"relative_difference":float(val),"original_2pct_PASS":bool(val<=.02),"original_gate_unchanged":True})
                    for stage,a in [("GAIN_SOURCE_PROFILED",dz),("HORIZONTAL_PROFILED",eff)]:
                        for i,j in R.PAIRS:
                            prior=ci[sc["scene_id"],res,g,noise,stage,R.METHODS[i],R.METHODS[j]]
                            ck("pair:"+str((sc["scene_id"],res,g,noise,stage,i,j)),abs(R.rel(a[:,i],a[:,j])-float(prior["relative_direction_difference"])),1e-7)
                        u=a[:,1]-a[:,0];v=a[:,2]-a[:,1];prior=li[sc["scene_id"],res,g,noise,stage]
                        ck("h2ratio:"+str((sc["scene_id"],res,g,noise,stage)),abs(float(u@v/(u@u))-float(prior["fitted_h2_ratio"])),1e-5)
                    for i,m in enumerate(R.METHODS):
                        for sig in E.SIGMAS:
                            full=np.column_stack([h[:,:2],dz[:,i:i+1]*200,h[:,2:]])/sig
                            _,sv,vh=linalg.svd(full,full_matrices=False)
                            rank=int(np.sum(sv>sv[0]*1e-10));null=vh[rank:]
                            nz=np.linalg.norm(null[:,2]) if len(null) else 0
                            iv=float(eff[:,i]@eff[:,i])/sig**2
                            iz=iv if nz<=1e-8 and iv*200**2>(sv[0]*1e-10)**2 else 0.
                            row=idx[sc["scene_id"],res,g,noise,sig,m];prior=float(row["effective_depth_information"])
                            err=abs(iz-prior)/max(iz,prior,1e-30);maxi=max(maxi,err)
                            ck("information:"+str((sc["scene_id"],res,g,noise,sig,m)),err,1e-6)
                            ck("rank:"+str((sc["scene_id"],res,g,noise,sig,m)),abs(rank-int(row["state_rank"])),0)
            assert time.monotonic()-start<=900,"INDEPENDENT_900S_CAP"
            print("COLD_REVIEW",g,sc["scene_id"],flush=True)
    for row in G.rows(O/"GRID_NUMERICAL_DIAGNOSTIC.csv"):
        a,ea=saved[row["scene"],row["resource"],80001,row["noise_model"]]
        b,eb=saved[row["scene"],row["resource"],160001,row["noise_model"]];i=R.METHODS.index(row["method"])
        for key,val in [("gain_source_direction_grid_relative",R.rel(a[:,i],b[:,i])),("horizontal_direction_grid_relative",R.rel(ea[:,i],eb[:,i])),("effective_information_grid_relative",abs(ea[:,i]@ea[:,i]-eb[:,i]@eb[:,i])/max(ea[:,i]@ea[:,i],eb[:,i]@eb[:,i],1e-30))]:
            ck("grid:"+str((row["scene"],row["resource"],row["noise_model"],i,key)),abs(val-float(row[key])),1e-7)
    envs=G.rows(O/"SSP_ENVIRONMENT_WINDOW_AUDIT.csv")
    for row in envs:
        f=float(row["frequency"]);g=int(row["mesh"]);d=int(row["offset"])
        lines=E.modpath(f,g,d).with_suffix(".env").read_text(encoding="utf-8").splitlines()
        s=[]
        for line in lines:
            a=line.split()
            try:
                if len(a)==6 and 0<=float(a[0])<=5000 and float(a[1])>1000:s.append(list(map(float,a[:2])))
            except ValueError:pass
        s=np.array(s);windowline=next(i for i,v in enumerate(lines) if v.strip()=="'R' 0.0")
        cl,ch=map(float,lines[windowline+1].split());opt=lines[3]
        base,_,_,_=R.environment(E.modpath(f,g,0).with_suffix(".env"))
        ck("SSP_exact_offset:"+str((f,g,d)),float(np.max(abs(s[:,1]-base[:,1]-d))),1e-10)
        ck("SSP_locations:"+str((f,g,d)),float(np.max(abs(s[:,0]-base[:,0]))),0)
        for key,v in [("c_min",s[:,1].min()),("c_max",s[:,1].max()),("CLOW",cl),("CHIGH",ch)]:
            ck("environment:"+str((f,g,d,key)),abs(float(row[key])-v),1e-10)
        if d:ck("no_perturbed_modal:"+str((f,g,d)),int(E.modpath(f,g,d).exists()),0)
    for row in G.rows(O/"NOMINAL_MODAL_WINDOW_MARGINS.csv"):
        f=float(row["frequency"]);g=int(row["mesh"]);i=int(row["nominal_mode"])-1
        # Cache binary parsing so the complete mode-by-mode review stays bounded.
        if not hasattr(main,"mcache"):main.mcache={}
        if (f,g) not in main.mcache:main.mcache[f,g]=B.mod(E.modpath(f,g,0))[2]
        speed=2*np.pi*f/main.mcache[f,g][i].real
        ck("modal_phase:"+str((f,g,i)),abs(speed-float(row["phase_speed"])),1e-10)
    ck("complete_stencil_rows",abs(len(cmps)-720),0)
    ck("complete_information_rows",abs(len(info)-360),0)
    ck("old_failures_retained",abs(sum(not r["original_2pct_PASS"] for r in step)-20),0)
    decision=E.read(O/"REVIEW_DECISION.json")
    ck("no_credit",decision["R4_percent"],0)
    ck("SSP_not_evaluated",int(decision["SSP_depth_information"]!="NOT_EVALUATED"),0)
    R.verify()
    E.write(O/"ORIGINAL_DEPTH_GATE_RECONSTRUCTION.csv",step)
    E.write(O/"INDEPENDENT_VALIDATION_CHECKS.csv",checks)
    E.dump(O/"VALIDATION.json",{"checks":len(checks),"PASS":sum(r["PASS"] for r in checks),
        "FAIL":sum(not r["PASS"] for r in checks),"cold_pressure_relative_max":maxp,
        "cold_derivative_relative_max":maxd,"cold_information_relative_max":maxi,
        "new_KRAKEN":0,"new_FIELD":0,"new_MC":0,"scope":"INDEPENDENT_SUMMATION_BINARY_PARSER_NUISANCE_CHART_QR_RECONSTRUCTION",
        "elapsed_s":time.monotonic()-start})
    assert all(r["PASS"] for r in checks),"INDEPENDENT_VALIDATION_FAILED;STOP"
if __name__=="__main__":main()
