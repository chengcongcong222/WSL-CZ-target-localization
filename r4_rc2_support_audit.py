"""Decimal independent rejection replay and complete dyadic-tree coverage validation."""
from decimal import Decimal as D,localcontext
from fractions import Fraction
from math import factorial
import time,json,itertools
import numpy as np
import r4_rc2_support as P
E=P.E;O=P.O
PAD=D("1e-35")
def frac(q):return D(q.numerator)/D(q.denominator)
def atanfrac(x,n=35):
    s=Fraction()
    for i in range(n):s+=(-1)**i*x**(2*i+1)/Fraction(2*i+1)
    t=(-1)**n*x**(2*n+1)/Fraction(2*n+1)
    return min(s,s+t),max(s,s+t)
al,ah=atanfrac(Fraction(1,5));bl,bh=atanfrac(Fraction(1,239))
PIL=16*al-4*bh;PIH=16*ah-4*bl
def dec(x):return D.from_float(float(x))
def poly(x,kind):
    if kind=="s":s=sum(((-1)**i*x**(2*i+1)/D(factorial(2*i+1)) for i in range(7)),D());rem=frac(Fraction(1,3)**15/Fraction(factorial(15)))
    elif kind=="c":s=sum(((-1)**i*x**(2*i)/D(factorial(2*i)) for i in range(7)),D());rem=frac(Fraction(1,3)**14/Fraction(factorial(14)))
    else:s=sum(((-1)**i*x**(2*i+1)/D(2*i+1) for i in range(9)),D());rem=frac(Fraction(1,3)**19/Fraction(19))
    return s-rem-PAD,s+rem+PAD
def mul(a,b):
    xs=[a[i]*b[j] for i in [0,1] for j in [0,1]]
    return min(xs)-PAD,max(xs)+PAD
def trig(lo,hi):
    rad=(frac(PIL)/D(180)-PAD,frac(PIH)/D(180)+PAD)
    xl=mul((lo,lo),rad);xh=mul((hi,hi),rad)
    sl=poly(xl[0],"s");sh=poly(xh[1],"s")
    cr=[poly(x,"c") for x in [xl[0],xl[1],xh[0],xh[1]]]
    return (sl[0],sh[1]),(min(x[0] for x in cr),D(1) if lo<=0<=hi else max(x[1] for x in cr))
def cold(box,t,platform,obs):
    r0,r1=map(dec,box[0]);sn,cs=trig(*map(dec,box[1]));sv,cv=trig(*map(dec,box[3]))
    rx=mul((r0*1000,r1*1000),cs);ry=mul((r0*1000,r1*1000),sn)
    vx=mul(tuple(map(dec,box[2])),cv);vy=mul(tuple(map(dec,box[2])),sv)
    rad=(frac(PIL)/D(180)-PAD,frac(PIH)/D(180)+PAD)
    sigma=dec(float(P.I.UP(P.I.RAD[1]*.1)))
    ss=D();maxz=D()
    for time_,node,bearing in zip(t,platform,obs):
        tt=dec(time_);px=dec(node[0]);py=dec(node[1]);b=dec(bearing)
        x=(rx[0]+vx[0]*tt-px-PAD,rx[1]+vx[1]*tt-px+PAD)
        y=(ry[0]+vy[0]*tt-py-PAD,ry[1]+vy[1]*tt-py+PAD)
        assert x[0]>0
        ratios=[y[i]/x[j] for i in [0,1] for j in [0,1]]
        lo=poly(min(ratios)-PAD,"a")[0];hi=poly(max(ratios)+PAD,"a")[1]
        delta=min(max(D(),lo-b,b-hi),max(D(),lo+b,-b-hi))
        z=max(D(),delta-D("1e-10"))/sigma
        ss+=z*z;maxz=max(maxz,z)
    return float(ss),float(maxz)
def main():
    start=time.monotonic();P.verify();checks=[]
    def ck(n,ok,detail=""):checks.append({"check":n,"PASS":bool(ok),"detail":str(detail)})
    records=E.G.rows(O/"PARTITION_AND_REJECTION_CERTIFICATES.csv");boxes=E.G.rows(O/"SUPPORT_BOXES.csv")
    f=E.read(O/"DESIGN_FREEZE.json");platform=np.array(f["platform_xy"])
    with np.load(P.SOURCE/"CASE_OBSERVATIONS.npz") as z:
        ids=z["case_ids"].tolist();allobs=z["bearing_rad"].copy();t=z["times_s"].copy()
    cold_count=0
    for cid in f["case_ids"]:
        rr={int(r["node"]):r for r in records if r["case_id"]==cid}
        kept={int(r["node"]):r for r in boxes if r["case_id"]==cid}
        ck("root:"+cid,1 in rr or 1 in kept)
        def getbox(r):return np.array([[float(r[f"lo{j}"]),float(r[f"hi{j}"])] for j in range(4)])
        ck("root_domain:"+cid,np.array_equal(getbox(rr.get(1,kept.get(1))),np.column_stack([P.LOW,P.HIGH])))
        for node,r in rr.items():
            box=getbox(r)
            if r["decision"]=="SPLIT":
                dim=int(r["split_dim"]);mid=float(r["mid"])
                for child in [node*2,node*2+1]:
                    row=rr.get(child,kept.get(child));ck("partition_child:"+cid+":"+str(child),row is not None)
                    expected=box.copy()
                    if child==node*2:expected[dim,1]=mid
                    else:expected[dim,0]=mid
                    ck("closed_partition:"+cid+":"+str(child),row is not None and np.array_equal(getbox(row),expected))
            elif r["decision"]=="REJECT":
                with localcontext() as ctx:
                    ctx.prec=50
                    ss,z=cold(box,t,platform,allobs[ids.index(cid)])
                cold_count+=1
                ck("Decimal_reject:"+cid+":"+str(node),ss>195 or z>5,(ss,z))
                ck("lower_bound_reconstruction:"+cid+":"+str(node),float(r["SSE_lower"])<=ss+1e-7 and float(r["max_abs_standardized_lower"])<=z+1e-7)
                ck("reject_not_exported:"+cid+":"+str(node),node not in kept)
            else:ck("retained_leaf:"+cid+":"+str(node),node in kept)
        for node,r in kept.items():
            ck("kept_parent_exists:"+cid+":"+str(node),node==1 or (node//2 in rr and rr[node//2]["decision"]=="SPLIT"))
            ck("kept_inside_domain:"+cid+":"+str(node),np.all(getbox(r)[:,0]>=P.LOW) and np.all(getbox(r)[:,1]<=P.HIGH))
        print("DECIMAL_TREE",cid,"cumulative rejects",cold_count,flush=True)
        assert time.monotonic()-start<=900,"INDEPENDENT_900S_CAP"
    cov=E.G.rows(O/"TRUTH_COVERAGE_EVALUATION.csv")
    ck("no_compatible_truth_deleted",all(r["compatible_but_deleted"]=="False" for r in cov))
    # Analytic geometry controls use prescribed mathematical states, not panel-truth observations.
    controls=[]
    state=np.array([50000.,0.,2.,0.]);times=np.arange(121)*10.
    straight=np.column_stack([2*times,np.zeros(121)])
    for name,track,dual in [("COLLINEAR_STRAIGHT",straight,False),("TURNING_SINGLE",platform,False),("TURNING_AUX5KM",platform,True)]:
        nodes=np.stack([track,track+[0,5000]],axis=1) if dual else track[:,None,:]
        target=state[:2]+times[:,None]*state[2:]
        delta=target[:,None,:]-nodes
        bearing=np.arctan2(delta[:,:,1],delta[:,:,0]).ravel();tt=np.repeat(times,nodes.shape[1])
        mat=np.column_stack([-np.sin(bearing)*50000,np.cos(bearing)*50000,-tt*np.sin(bearing)*2,tt*np.cos(bearing)*2])
        sv=np.linalg.svd(mat,compute_uv=False);rank=int(np.sum(sv>sv[0]*1e-10))
        controls.append({"control":name,"rank":rank,"singular_values":json.dumps(sv.tolist()),"scope":"ANALYTIC_DESIGN_GEOMETRY;NO_NEW_NOISY_OBSERVATIONS"})
    ck("collinear_rank_deficient",controls[0]["rank"]<4)
    ck("turn_and_aux_rank4",controls[1]["rank"]==controls[2]["rank"]==4)
    # Fixed root corners, epoch endpoints, and midpoint verify the bound at evaluation-only points.
    root=np.column_stack([P.LOW,P.HIGH]);lo,hi=P.I.predict(root,t,platform)
    for bits in itertools.product([0,1],repeat=4):
        point=root[np.arange(4),bits]
        theta=np.deg2rad(point[1]);psi=np.deg2rad(point[3])
        x=point[0]*1000*np.cos(theta)+point[2]*t*np.cos(psi)-platform[:,0]
        y=point[0]*1000*np.sin(theta)+point[2]*t*np.sin(psi)-platform[:,1]
        b=np.arctan2(y,x)
        ck("root_corner:"+str(bits),np.all((lo<=b)&(b<=hi)))
    P.verify();P.write(O/"GEOMETRY_CONTROLS.csv",controls)
    P.write(O/"INDEPENDENT_VALIDATION_CHECKS.csv",checks)
    v={"checks":len(checks),"PASS":sum(r["PASS"] for r in checks),"FAIL":sum(not r["PASS"] for r in checks),
       "Decimal_rejection_certificates":cold_count,"complete_dyadic_partition_rebuilt":True,
       "statistical_scope":"ANALYTIC_CONDITIONAL_BOUND;NOT_VALIDATED_BY_EMPIRICAL_RETENTION",
       "new_observations":0,"new_solvers":0,"elapsed_s":time.monotonic()-start}
    E.dump(O/"VALIDATION.json",v)
    d=E.read(O/"RC2_SUPPORT_DECISION.json")
    d["implementation_status"]="CERTIFIED_CONSERVATIVE_UNDER_REGISTERED_ARITHMETIC_CONTRACT" if v["FAIL"]==0 else "IMPLEMENTATION_INVALID"
    if v["FAIL"]:d["architecture_classifications"]["U0"]="RC2_SUPPORT_COMPUTATION_INCOMPLETE"
    d["independent_review"]="PASS" if v["FAIL"]==0 else "FAIL"
    E.dump(O/"RC2_SUPPORT_DECISION.json",d)
if __name__=="__main__":main()
