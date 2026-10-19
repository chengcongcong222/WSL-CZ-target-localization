from common import *
def main():
 rows=[];mat={};maxerr=0.;t=np.array(read('CAP06_FROZEN_SPEC.json')['times_s']);scale=np.array([1000.,1000.,1.,1.])
 for en in read('P0_FROZEN_PANEL.json')['states']:
  for side in [-1,1]:
   s=np.array(en['state']);s[1]*=side;s[4]*=side;x=np.r_[s[0]*np.array([np.cos(s[1]),np.sin(s[1])]),s[3]*np.array([np.cos(s[4]),np.sin(s[4])])]
   Fs={};free={}
   for node,off in [('MAIN',[0,0]),('B0',[0,side*5000]),('B3',[3000,side*15000])]:
    pos=platform(t)+off;d=x[:2]+t[:,None]*x[2:]-pos;r2=np.sum(d*d,axis=1);alpha=np.arctan2(d[:,1],d[:,0])-axis(t)
    Jang=np.column_stack((-d[:,1]/r2,d[:,0]/r2,-t*d[:,1]/r2,t*d[:,0]/r2));J=-np.sin(alpha)[:,None]*Jang*scale
    def f(q):
     dd=q[:2]+t[:,None]*q[2:]-pos
     return np.cos(np.arctan2(dd[:,1],dd[:,0])-axis(t))
    D=np.column_stack([(f(x+np.eye(4)[j]*scale[j]*1e-4)-f(x-np.eye(4)[j]*scale[j]*1e-4))/(2e-4) for j in range(4)])
    err=float(np.linalg.norm(D-J)/max(np.linalg.norm(J),1e-12));maxerr=max(maxerr,err)
    # Unknown complex source independently removed: spatial phase information is proportional to centered aperture variance.
    weight=2*128*100*(2*np.pi*200/1500)**2*np.sum(np.arange(-7,8,2,dtype=float)**2)
    F=weight*J.T@J;B=np.sin(alpha);B=B/np.linalg.norm(B) if np.linalg.norm(B)>1e-12 else B;R=J-B[:,None]*(B@J)[None,:]
    Fs[node]=F;free[node]=weight*R.T@R
   for res,nodes in [('S0',['MAIN']),('S2_B0',['MAIN','B0']),('S2_B3',['MAIN','B3'])]:
    for mode,dictF in [('CAL_ORIENTATION_IDEAL',Fs),('FREE_NODE_ORIENTATION',free)]:
     F=sum(dictF[n] for n in nodes);eig=np.linalg.eigvalsh(F);Jp=F[:2,:2];couple=F[:2,2:];V=F[2:,2:]-couple.T@np.linalg.pinv(Jp,rcond=1e-12)@couple;ev=np.linalg.eigvalsh(V)
     key=f"{en['id']}_{side}_{res}_{mode}";mat[key]=F
     rows.append(dict(state_id=en['id'],split=en['split'],side=side,resource=res,nuisance=mode,scaled_rank=int(np.sum(eig>max(eig[-1],1)*1e-10)),weak_eigenvalue=float(eig[0]),velocity_profile_eigen_min=float(ev[0]),velocity_profile_trace=float(np.trace(V)),z='ZERO_INFORMATION',scope='LOCAL_DESIGN_UPPER_BOUND_NOT_MEASURED_PRECISION; source coefficient profiled; known relative channel gains; no real CZ forward'))
 tab('P2_LOCAL_INFORMATION.csv',rows);np.savez_compressed(RAW/'P2_LOCAL_INFORMATION_MATRICES.npz',**mat)
 good=maxerr<=1e-6 and all(r['weak_eigenvalue']>=-1e-7 for r in rows)
 write('P2_INFORMATION_GATE.json',dict(PASS=good,jacobian_relative_max=maxerr,threshold=1e-6,depth_zero=True,scope='No pseudoinverse covariance or zero variance claim; ideal calibrated local design only'))
 gate=read('P2_GATE.json');gate['independent_information_control']=good
 if not good:gate['P3_eligible']=[]
 write('P2_GATE.json',gate)
if __name__=='__main__':main()
