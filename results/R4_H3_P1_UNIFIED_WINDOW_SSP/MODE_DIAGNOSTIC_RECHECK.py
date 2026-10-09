"""Independent saved-mode assignment diagnostic and full subset contribution replay."""
import json,time
import numpy as np
from scipy.optimize import linear_sum_assignment
import r4_h3_p1 as P
import r4_h3_p1_audit as I
import r4_e2_nine_pair_pilot_audit as B
E=P.E;G=P.G;O=P.O
begin=time.monotonic();mods={};cache={};checks=[]
def ck(n,v,lim=1e-7):E.ck(checks,n,v,lim)
for row in G.rows(O/"MODE_MATCHING.csv"):
    c=float(row["CLOW_a"]);cc=float(row["CLOW_b"]);g=int(row["mesh"]);f=float(row["frequency"]);d=float(row["delta_a"]);dd=float(row["delta_b"])
    a=(c,g,f,d);b=(cc,g,f,dd)
    for key in [a,b]:
        if key not in mods:mods[key]=B.mod(P.path(*key))
    za,pa,ka=mods[a];zb,pb,kb=mods[b]
    normalized_a=pa/np.maximum(np.linalg.norm(pa,axis=0),1e-300)
    normalized_b=pb/np.maximum(np.linalg.norm(pb,axis=0),1e-300)
    corr=abs(normalized_a.conj().T@normalized_b)
    gap=np.empty(len(ka));gap[0]=abs(ka[0].real-ka[1].real);gap[-1]=abs(ka[-1].real-ka[-2].real)
    gap[1:-1]=np.minimum(abs(ka.real[1:-1]-ka.real[:-2]),abs(ka.real[2:]-ka.real[1:-1]))
    cost=1-corr+.05*np.minimum(abs(ka.real[:,None]-kb.real[None,:])/np.maximum(gap[:,None],1e-12),20)
    ia,ib=linear_sum_assignment(cost)
    bad=corr[ia,ib]<.95
    key=(row["kind"],c,cc,g,f,d,dd);cache[key]=(ia[bad],ib[bad])
    ck("sampled_correlation:"+str(key),abs(float(corr[ia,ib].min())-float(row["min_sampled_shape_correlation"])))
    ck("low_correlation_count_a:"+str(key),abs(int(bad.sum())-int(row["low_correlation_a"])),0)
    ck("low_correlation_count_b:"+str(key),abs(int(bad.sum())-int(row["low_correlation_b"])),0)
states={s["scene_id"]:s["state"] for s in E.scenes()}
for row in G.rows(O/"WINDOW_FIELD_CONTRIBUTION.csv"):
    c=float(row["CLOW_a"]);cc=float(row["CLOW_b"]);g=int(row["mesh"]);f=float(row["frequency"]);d=float(row["delta_a"]);dd=float(row["delta_b"])
    key=(row["kind"],c,cc,g,f,d,dd);aa,bb=cache[key]
    for side,m,ix in [("a",mods[c,g,f,d],aa),("b",mods[cc,g,f,dd],bb)]:
        z,phi,k=m;st=states[row["scene"]];res=row["resource"]
        p=I.field([m],st,res);q=I.field([(z,phi[:,ix],k[ix])],st,res)
        val=float(np.linalg.norm(q)/max(np.linalg.norm(p),1e-30))
        ck("ambiguous_subset_contribution:"+str((key,row["scene"],res,side)),abs(val-float(row["low_correlation_field_fraction_"+side])))
E.write(O/"MODE_DIAGNOSTIC_INDEPENDENT_CHECKS.csv",checks)
E.dump(O/"MODE_DIAGNOSTIC_VALIDATION.json",{"checks":len(checks),"PASS":sum(r["PASS"] for r in checks),"FAIL":sum(not r["PASS"] for r in checks),"elapsed_s":time.monotonic()-begin,"new_KRAKEN":0,"scope":"INDEPENDENT_NORMALIZED_MODAL_SHAPE_COST_AND_TERMWISE_LOW_CORRELATION_SUBSET_RECONSTRUCTION"})
assert all(r["PASS"] for r in checks)
print(json.dumps(E.read(O/"MODE_DIAGNOSTIC_VALIDATION.json"),indent=2))
