"""Frozen paired generator. Truth is never an estimator argument."""
import numpy as np

def platform():
    t=np.arange(121,dtype=float)*10;post=np.maximum(t-600,0)
    p=np.column_stack([2*np.minimum(t,600)+2*post*np.cos(np.pi/12),2*post*np.sin(np.pi/12)])
    return t,p

def cartesian(truth):
    r,th,v,ps=[float(truth[k]) for k in ['r0_km','theta0_deg','v_mps','psi_deg']];th,ps=np.deg2rad([th,ps])
    return np.array([r*1000*np.cos(th),r*1000*np.sin(th),v*np.cos(ps),v*np.sin(ps)])

def node_geometry(scheme,side):
    t,main=platform();off=np.array([scheme['forward_km']*1000,side*scheme['lateral_km']*1000])
    return t,np.stack([main,main+off],axis=1)

def observations(truth,scheme,quality,normals,uniforms,sides):
    t,main=platform();n=len(sides);true=np.broadcast_to(main[None,:,None,:],(n,121,2,2)).copy()
    true[:,:,1,0]+=scheme['forward_km']*1000;true[:,:,1,1]+=sides[:,None]*scheme['lateral_km']*1000
    state=cartesian(truth);target=state[:2]+t[:,None]*state[2:];d=target[None,:,None,:]-true
    c=(uniforms[:,0]*2-1)*quality['common_bound_deg'];half=(uniforms[:,1]*2-1)*quality['half_diff_bound_deg'];bias=np.deg2rad(c[:,None]+half[:,None]*np.array([-1.,1.]))
    y=np.arctan2(d[:,:,:,1],d[:,:,:,0])+bias[:,None,:]+np.deg2rad(quality['bearing_sigma_deg'])*normals[:,:,:2]
    y=(y+np.pi)%(2*np.pi)-np.pi
    reported=true+quality['nav_sigma_m']*normals[:,:,2:].reshape(n,121,2,2)
    return t,reported,y,np.column_stack([c,half])
