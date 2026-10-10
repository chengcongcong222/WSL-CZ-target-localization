"""Independent Cartesian geometry and all-sign one-dimensional bias audit."""
import numpy as np
T=np.arange(121,dtype=float)*10.
IDX=np.arange(7)*20
def wrap(a):return np.arctan2(np.sin(a),np.cos(a))
def geometry(states):
    h=np.atleast_2d(np.asarray(states,dtype=float));theta=h[:,1]*np.pi/180;psi=h[:,3]*np.pi/180
    initial=np.stack((1000*h[:,0]*np.cos(theta),1000*h[:,0]*np.sin(theta)),axis=1)
    velocity=np.stack((h[:,2]*np.cos(psi),h[:,2]*np.sin(psi)),axis=1)
    obs=np.empty((121,2));mask=T<=600
    obs[mask]=np.stack((2*T[mask],np.zeros(mask.sum())),axis=1)
    obs[~mask]=np.array([1200.,0.])+2*(T[~mask]-600)[:,None]*np.array([np.cos(np.pi/12),np.sin(np.pi/12)])
    relative=initial[:,None,:]+T[None,:,None]*velocity[:,None,:]-obs[None,:,:]
    distance=np.sqrt(np.einsum('nti,nti->nt',relative,relative))
    beta=np.arctan2(relative[:,:,1],relative[:,:,0])
    q=(distance[:,IDX[1:]]-distance[:,IDX[:-1]])*.005
    return beta,q,distance
def profile(q,u):
    q=np.atleast_2d(q);u=np.atleast_2d(u)
    signs=np.where((np.arange(64)[:,None]>>np.arange(6))&1,1.,-1.)
    left=np.maximum(-.2,np.max(np.where(signs[None,:,:]>0,-q[:,None,:],-np.inf),axis=-1))
    right=np.minimum(.2,np.min(np.where(signs[None,:,:]<0,-q[:,None,:],np.inf),axis=-1))
    b=np.clip(np.mean(signs[None,:,:]*u[:,None,:]-q[:,None,:],axis=-1),left,right)
    cost=np.square(np.abs(q[:,None,:]+b[:,:,None])-u[:,None,:]).sum(axis=-1)
    cost=np.where(left<=right,cost,np.inf)
    ix=cost.argmin(axis=-1)
    return b[np.arange(len(q)),ix],cost[np.arange(len(q)),ix]
def residual(state,beta_obs,q_obs,sigma):
    beta,q,_=geometry(state)
    bearing=wrap(beta[0]-beta_obs)/np.radians(.1)
    return np.r_[bearing,(q[0]-q_obs)/sigma] if sigma is not None else bearing
