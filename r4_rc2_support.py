"""Once-only saved-observation horizontal support; no truth in support constructor."""
from pathlib import Path
import time,json,csv
import numpy as np
import r4_h3_p0 as E
import r4_rc2_interval as I
O=E.ROOT/"results/R4_RC2_OBSERVATION_SUPPORT"
PARENT="08da59d7cc1b1b01b1dc27e8c6551257ab7ff530"
SOURCE=E.ROOT/"results/R4_A1_FIX_CONTINUOUS_SEARCH"
LOW=np.array([45.,-5.,1.,-15.]);HIGH=np.array([60.,5.,3.,15.])
RESOLUTION=np.array([.25,.05,.05,1.])
CAP=2047
def write(p,rows,headers=None):
    with Path(p).open("w",encoding="utf-8",newline="") as f:
        w=csv.DictWriter(f,fieldnames=list(rows[0]) if rows else headers);w.writeheader();w.writerows(rows)
def verify():
    for b in E.read(O/"DESIGN_FREEZE.json")["bindings"]:assert E.sha(E.ROOT/b["path"])==b["sha256"],b["path"]
def solve(obs,t,platform):
    # Receives only registered observation and platform arrays, no evaluation truth.
    begin=time.monotonic();queue=[(1,0,np.column_stack([LOW,HIGH]))];head=0;records=[];kept=[];count=0
    sigma=float(I.UP(I.RAD[1]*.1))
    while head<len(queue) and count<CAP and time.monotonic()-begin<30:
        node,parent,box=queue[head];head+=1;count+=1
        s,z,lo,hi=I.bound(box,t,platform,obs,sigma)
        decision="RETAIN_LEAF";dim=-1;mid=0.
        # A fixed roundoff security gap prevents threshold-tie rejections.
        if s>195+1e-6 or z>5+1e-6:decision="REJECT"
        elif np.all(box[:,1]-box[:,0]<=RESOLUTION):
            kept.append((node,box,"RESOLUTION_UNRESOLVED"))
        elif count>=CAP or time.monotonic()-begin>=30:
            kept.append((node,box,"BUDGET_UNRESOLVED"))
        else:
            decision="SPLIT";dim=int(np.argmax((box[:,1]-box[:,0])/RESOLUTION))
            mid=float((box[dim,0]+box[dim,1])/2)
            a=box.copy();b=box.copy();a[dim,1]=mid;b[dim,0]=mid
            queue.append((node*2,node,a));queue.append((node*2+1,node,b))
        rec={"node":node,"parent":parent,"decision":decision,"split_dim":dim,"mid":mid,
            "SSE_lower":s,"max_abs_standardized_lower":z}
        for j in range(4):rec[f"lo{j}"]=box[j,0];rec[f"hi{j}"]=box[j,1]
        records.append(rec)
    for node,parent,box in queue[head:]:kept.append((node,box,"BUDGET_UNRESOLVED"))
    return kept,records,{"evaluated_boxes":count,"retained_boxes":len(kept),
        "unresolved_boxes":len(kept),"rejected_boxes":sum(r["decision"]=="REJECT" for r in records),
        "budget_exhausted":count>=CAP or time.monotonic()-begin>=30,"elapsed_s":time.monotonic()-begin}
def execute():
    verify();assert not (O/"EXECUTION_STARTED.json").exists()
    assert E.G.git("log","-1","--pretty=%s")=="R4 RC2: freeze observation-defined multibranch support pilot"
    assert E.G.git("rev-parse","HEAD")==E.G.git("ls-remote","origin","refs/heads/main").split()[0]
    E.dump(O/"EXECUTION_STARTED.json",{"design_SHA":E.G.git("rev-parse","HEAD"),"no_new_observations":True})
    with np.load(SOURCE/"CASE_OBSERVATIONS.npz") as z:
        ids=z["case_ids"].tolist();obs=z["bearing_rad"];t=z["times_s"]
    platform=E.read(O/"DESIGN_FREEZE.json")["platform_xy"]
    platform=np.array(platform);boxes=[];logs=[];diag=[]
    for i,cid in enumerate(ids):
        if "_s0.1_" not in cid:continue
        kept,records,stats=solve(obs[i],t,platform)
        for r in records:logs.append(dict(case_id=cid,**r))
        for node,box,status in kept:
            rec={"architecture":"U0","case_id":cid,"node":node,"status":status,
                 "sign_branch":"NEGATIVE_THETA" if box[1,1]<=0 else "POSITIVE_THETA" if box[1,0]>=0 else "CROSSES_ZERO",
                 "component_scope":"SIGN_SECTORS_NOT_CERTIFIED_CONNECTED_COMPONENTS"}
            for j in range(4):rec[f"lo{j}"]=box[j,0];rec[f"hi{j}"]=box[j,1]
            boxes.append(rec)
        if kept:
            aa=np.array([x[1] for x in kept]);bounds=np.column_stack([aa[:,:,0].min(axis=0),aa[:,:,1].max(axis=0)])
        else:bounds=np.full((4,2),np.nan)
        row={"architecture":"U0","case_id":cid,**stats}
        for j,name in enumerate(["range_km","theta_deg","speed_mps","heading_deg"]):
            row[name+"_low"]=float(bounds[j,0]) if kept else "EMPTY"
            row[name+"_high"]=float(bounds[j,1]) if kept else "EMPTY"
            row[name+"_width"]=float(bounds[j,1]-bounds[j,0]) if kept else "EMPTY"
        row["sign_sectors_retained"]=";".join(sorted(set(x["sign_branch"] for x in boxes if x["case_id"]==cid)))
        diag.append(row);print("SUPPORT",cid,stats,flush=True)
    write(O/"SUPPORT_BOXES.csv",boxes,["architecture","case_id","node","status"])
    write(O/"PARTITION_AND_REJECTION_CERTIFICATES.csv",logs)
    write(O/"BRANCH_AND_BUDGET_DIAGNOSTICS.csv",diag)
    # Truth files are read only AFTER all support construction has ended.
    truths={}
    for p in [E.ROOT/"results/R4_A1_OFFGRID_BEARING_BOUNDARY/OFFGRID_TRUTH_PANEL.csv",SOURCE/"HOLDOUT_TRUTH_PANEL.csv"]:
        for r in E.G.rows(p):truths[r["panel_id"]]=np.array([float(r[k]) for k in ["r_km","theta_deg","v_mps","psi_deg"]])
    evaluated=[]
    for d in diag:
        cid=d["case_id"];state=truths[cid.split("_")[0]]
        i=ids.index(cid);pred=I.predict(np.column_stack([state,state]),t,platform)
        beta=(pred[0]+pred[1])/2
        residual=np.minimum(abs(beta-obs[i]),abs(beta+obs[i]))
        scalar=np.maximum(0,residual-1e-10)/(I.RAD[1]*.1)
        compatible=float(scalar@scalar)<=195 and max(scalar)<=5
        actual=any(all(float(r[f"lo{j}"])<=state[j]<=float(r[f"hi{j}"]) for j in range(4)) for r in boxes if r["case_id"]==cid)
        evaluated.append({"case_id":cid,"truth_compatible_event":bool(compatible),"algorithm_kept_truth":actual,
            "compatible_but_deleted":bool(compatible and not actual),"truth_SSE":float(scalar@scalar),
            "scope":"POSTHOC_FIXED_FINITE_PANEL;NOT_EMPIRICAL_GLOBAL_COVERAGE"})
    write(O/"TRUTH_COVERAGE_EVALUATION.csv",evaluated)
    comparison=[{"architecture":"U0","input_status":"ARCHIVED_P_Q_FALLBACK","cases":len(diag),"H01_H06_H12_available":False,
        "coverage_scope":"FROZEN_GAUSSIAN_EXACT_NAV_CONDITIONAL;REFLECTION_OUTER_UNION","hardware":"ONE_MAIN_PLATFORM_HLA"},
        {"architecture":"U","input_status":"OBSERVATION_OR_NOISE_CONTRACT_INCOMPLETE","cases":0,"H01_H06_H12_available":False,
        "coverage_scope":"NOT_EVALUATED;NO_RAW_BEARING_OR_REPORTED_NODE_ARCHIVE","hardware":"MAIN_PLUS_5KM_AUXILIARY_SEPARATE_RESOURCE"}]
    write(O/"SINGLE_VS_AUXILIARY_SUPPORT.csv",comparison)
    informative=all(any(r[name+"_width"]!="EMPTY" and float(r[name+"_width"])<.9*width for name,width in zip(["range_km","theta_deg","speed_mps","heading_deg"],HIGH-LOW)) for r in diag)
    broad=any(r["range_km_width"]=="EMPTY" or float(r["range_km_width"])>1 or float(r["speed_mps_width"])>.2 for r in diag)
    E.dump(O/"RC2_SUPPORT_DECISION.json",{"architecture_classifications":{
        "U0":"RC2_SUPPORT_COMPUTATION_INCOMPLETE" if not informative else "RC2_OUTER_SUPPORT_VALID_BUT_BROAD" if broad else "RC2_CONDITIONAL_OUTER_SUPPORT_ESTABLISHED",
        "U":"RC2_OBSERVATION_OR_NOISE_CONTRACT_INCOMPLETE"},
        "implementation_status":"PENDING_INDEPENDENT_CERTIFICATE_REVIEW","case_count":len(diag),
        "source_panel":"ARCHIVED_P_Q;NOT_H01_H06_H12","candidate_builder_used_truth":False,
        "connected_component_completeness":"OUTER_UNION_COVERS_ALL_BRANCHES;NO_EXACT_COMPONENT_COUNT",
        "new_KRAKEN":0,"new_FIELD":0,"new_MC":0,"new_bearing_draws":0,"new_audio":0,
        "H3_C1e":"NOT_EVALUATED","complete_RC2_engineering_support":"NOT_ESTABLISHED",
        "H3_G0":"CONDITIONAL_LOCAL_DEPTH_INFORMATION_ACCEPTED","R4_percent":0,"next":"STOP"})
if __name__=="__main__":execute()
