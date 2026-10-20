"""Independent mathematics. No production modules, model, generator or statistics imports."""
import hashlib,struct,math
import numpy as np

def rebuild_panel(cfg,registry):
 p=cfg['panel'];prior=[[float(x[k]) if x[k] else None for k in ['r_km','theta_deg','v_mps','psi_deg']] for x in registry];states=[];log=[]
 for i in range(p['n']):
  for a in range(p['max_attempts_per_stratum']):
   digest=hashlib.sha256(f"{p['seed_text']}/{i}/{a}".encode('ascii')).digest();n=struct.unpack('>QQQQ',digest);u=[((v>>12)+.5)/(2**52) for v in n]
   h=[p['lower'][0]+(i+u[0])*(p['upper'][0]-p['lower'][0])/12]
   for j in range(1,4):h.append(p['lower'][j]+u[j]*(p['upper'][j]-p['lower'][j]))
   why='ACCEPT'
   for j in range(4):
    loc=(h[j]-p['grid_origin'][j])/p['grid_step'][j]
    if abs(loc-round(loc))<=p['grid_fraction_guard']:why=f'GRID_{j}';break
   if why=='ACCEPT':
    for kind,old in [('HISTORY',prior),('NEW',states)]:
     for ii,z in enumerate(old):
      axes=range(2) if z[2] is None else range(4);w=p['static_position_guard'] if z[2] is None else p['old_and_new_four_state_neighborhood_half_width']
      for mirror in (1,-1):
       passed=True
       for j in axes:
        d=h[j]-(mirror*z[j] if j in (1,3) else z[j]);d=(d+180)%360-180 if j in (1,3) else d
        if abs(d)>w[j]:passed=False;break
       if passed:why=f'{kind}_{ii}_{mirror}';break
      if why!='ACCEPT':break
     if why!='ACCEPT':break
   log.append(dict(case_index=i,attempt=a,digest=digest.hex(),h=h,decision=why))
   if why=='ACCEPT':states.append(h);break
  else:raise ValueError('independent panel exhausted')
 return dict(states=states,attempts=log)

def residual_and_jac(h,t,P,y,quality):
 r,th,v,ps,c,d=map(float,h);r*=1000;th*=math.pi/180;ps*=math.pi/180
 loc=np.array([r*math.cos(th),r*math.sin(th)])+P[0,0];vel=v*np.array([math.cos(ps),math.sin(ps)]);D=(loc+t[:,None]*vel)[:,None,:]-P;sq=np.einsum('tni,tni->tn',D,D)
 mu=np.arctan2(D[...,1],D[...,0])+(c+d*np.array([-1.,1.]))*math.pi/180;e=(mu-y+math.pi)%(2*math.pi)-math.pi
 s0=(quality['bearing_sigma_deg']*math.pi/180)**2;var=s0+quality['nav_sigma_m']**2/sq
 derivatives=np.zeros((len(t),2,6));derivatives[:,:,0]=[1000*math.cos(th),1000*math.sin(th)];derivatives[:,:,1]=[-r*math.sin(th)*math.pi/180,r*math.cos(th)*math.pi/180];derivatives[:,:,2]=t[:,None]*[math.cos(ps),math.sin(ps)];derivatives[:,:,3]=t[:,None]*[-v*math.sin(ps)*math.pi/180,v*math.cos(ps)*math.pi/180]
 H=(D[...,0,None]*derivatives[:,None,1,:]-D[...,1,None]*derivatives[:,None,0,:])/sq[...,None];H[...,4]=math.pi/180;H[...,5]=np.array([-1.,1.])*math.pi/180
 dvar=-2*quality['nav_sigma_m']**2*np.einsum('tni,tij->tnj',D,derivatives)/sq[...,None]**2
 rr=e/np.sqrt(var);ll=np.sqrt(np.log1p((var-s0)/s0));J1=H/np.sqrt(var[...,None])-rr[...,None]*dvar/(2*var[...,None]);J2=dvar/(2*var[...,None]*ll[...,None])
 return np.concatenate([rr.ravel(),ll.ravel()]),np.concatenate([J1.reshape(-1,6),J2.reshape(-1,6)])

def initial(t,P,y):
 # Independent 2D linear intersection solved in determinant form, then frozen IRLS.
 u0=np.c_[np.cos(y[:,0]),np.sin(y[:,0])];u1=np.c_[np.cos(y[:,1]),np.sin(y[:,1])];d=P[:,1]-P[:,0]
 det=u0[:,0]*u1[:,1]-u0[:,1]*u1[:,0]
 with np.errstate(divide='ignore',invalid='ignore'):
  first=(d[:,0]*u1[:,1]-d[:,1]*u1[:,0])/det;second=(d[:,0]*u0[:,1]-d[:,1]*u0[:,0])/det
 points=P[:,0]+first[:,None]*u0;keep=(abs(det)>1e-12)&(first>0)&(second>0)&np.isfinite(points).all(axis=1)
 if keep.sum()<4:return np.array([52.5,0,2,0,0,0]),True,int(keep.sum())
 A=np.c_[np.ones(keep.sum()),t[keep]/1200];Y=points[keep];fit=np.linalg.lstsq(A,Y,rcond=1e-12)[0]
 for _ in range(5):
  err=np.sqrt(np.sum((Y-A@fit)**2,axis=1));scale=max(10,float(np.median(err)));w=np.sqrt(np.minimum(1,1.5*scale/np.maximum(err,1e-15)));fit=np.linalg.lstsq(A*w[:,None],Y*w[:,None],rcond=1e-12)[0]
 a=fit[0]-P[0,0];b=fit[1]/1200
 return np.array([np.linalg.norm(a)/1000,np.degrees(np.arctan2(a[1],a[0])),np.linalg.norm(b),np.degrees(np.arctan2(b[1],b[0])),0,0]),False,int(keep.sum())

def fisher_design(state,t,P,q,scales):
 D=(np.asarray(state[:2])+t[:,None]*np.asarray(state[2:4]))[:,None,:]-P;r2=np.einsum('tni,tni->tn',D,D);sigma=(q['bearing_sigma_deg']*np.pi/180)**2;var=sigma+q['nav_sigma_m']**2/r2
 gx=-D[:,:,1]/r2;gy=D[:,:,0]/r2
 H=np.zeros((len(t),2,6));H[:,:,0]=gx;H[:,:,1]=gy;H[:,:,2]=gx*t[:,None];H[:,:,3]=gy*t[:,None];H[:,:,4]=1;H[:,:,5]=[-1,1]
 V=np.zeros_like(H);base=-2*q['nav_sigma_m']**2/r2**2;V[:,:,0]=base*D[:,:,0];V[:,:,1]=base*D[:,:,1];V[:,:,2]=V[:,:,0]*t[:,None];V[:,:,3]=V[:,:,1]*t[:,None]
 A=np.concatenate([(H*scales/np.sqrt(var[:,:,None])).reshape(-1,6),(V*scales/(np.sqrt(2)*var[:,:,None])).reshape(-1,6)])
 return A,H.reshape(-1,6),V.reshape(-1,6),var
