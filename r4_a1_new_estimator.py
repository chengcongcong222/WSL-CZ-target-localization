"""Observation-only horizontal estimator. No generator or truth dependencies."""
import numpy as np
from scipy.optimize import least_squares
SCALE=np.array([50000.,50000.,2.,2.])
def wrap(x):return (x+np.pi)%(2*np.pi)-np.pi

def initialize(times,node_positions,bearings):
 p1=node_positions[:,0];p2=node_positions[:,1];u1=np.column_stack([np.cos(bearings[:,0]),np.sin(bearings[:,0])]);u2=np.column_stack([np.cos(bearings[:,1]),np.sin(bearings[:,1])]);d=p2-p1;det=u1[:,0]*u2[:,1]-u1[:,1]*u2[:,0]
 with np.errstate(divide='ignore',invalid='ignore'):
  l1=(d[:,0]*u2[:,1]-d[:,1]*u2[:,0])/det;l2=(d[:,0]*u1[:,1]-d[:,1]*u1[:,0])/det;xy=p1+l1[:,None]*u1
 valid=(np.abs(det)>1e-12)&(l1>0)&(l2>0)&np.isfinite(xy).all(axis=1)
 if valid.sum()<4 or np.ptp(times[valid])<=0:raise ValueError('insufficient mathematically valid intersections')
 t=times[valid];points=xy[valid];T=np.column_stack([np.ones(len(t)),t/1200]);coef=np.linalg.lstsq(T,points,rcond=None)[0];weights=np.ones(len(t))
 for _ in range(5):
  norm=np.linalg.norm(points-T@coef,axis=1);scale=max(10.,float(np.median(norm)));weights=np.minimum(1.,1.5*scale/np.maximum(norm,1e-15));sw=np.sqrt(weights);coef=np.linalg.lstsq(T*sw[:,None],points*sw[:,None],rcond=None)[0]
 state=np.r_[coef[0],coef[1]/1200]
 return state,dict(valid_intersections=int(valid.sum()),invalid_intersections=int((~valid).sum()),min_abs_det=float(np.min(np.abs(det))),initializer_weight_min=float(weights.min()))

def residual(q,times,node_positions,bearings,sigma_deg):
 state=q*SCALE;target=state[:2]+times[:,None]*state[2:];d=target[:,None,:]-node_positions;angles=np.arctan2(d[:,:,1],d[:,:,0]);return (wrap(angles-bearings)/np.deg2rad(sigma_deg)).ravel()
def jacobian(q,times,node_positions,bearings,sigma_deg):
 state=q*SCALE;d=(state[:2]+times[:,None]*state[2:])[:,None,:]-node_positions;rho=np.maximum((d*d).sum(axis=2),1e-24);gx=-d[:,:,1]/rho;gy=d[:,:,0]/rho
 return np.stack([gx*SCALE[0],gy*SCALE[1],gx*times[:,None]*SCALE[2],gy*times[:,None]*SCALE[3]],axis=-1).reshape(-1,4)/np.deg2rad(sigma_deg)

def estimate(times,node_positions,bearings,sigma_deg):
 times=np.asarray(times);node_positions=np.asarray(node_positions);bearings=np.asarray(bearings)
 if node_positions.shape!=(len(times),2,2) or bearings.shape!=(len(times),2):raise ValueError('shape')
 x,diag=initialize(times,node_positions,bearings);fit=least_squares(residual,x/SCALE,jac=jacobian,args=(times,node_positions,bearings,sigma_deg),loss='soft_l1',f_scale=1.5,method='trf',max_nfev=100,ftol=1e-10,xtol=1e-10,gtol=1e-10)
 state=fit.x*SCALE;rr=residual(fit.x,times,node_positions,bearings,sigma_deg);J=jacobian(fit.x,times,node_positions,bearings,sigma_deg);sv=np.linalg.svd(J,compute_uv=False);variance=float(rr@rr/max(1,len(rr)-4));covq=np.linalg.pinv(J.T@J,rcond=1e-12)*variance;se=np.sqrt(np.maximum(np.diag(covq),0))*SCALE
 diag.update(initial_x_m=float(x[0]),initial_y_m=float(x[1]),initial_vx_mps=float(x[2]),initial_vy_mps=float(x[3]),solver_success=bool(fit.success),solver_status=int(fit.status),nfev=int(fit.nfev),cost=float(fit.cost),residual_RMS_deg=float(np.sqrt(np.mean(rr**2))*sigma_deg),jacobian_condition_scaled=float(sv[0]/sv[-1]),approx_se_x_m=float(se[0]),approx_se_y_m=float(se[1]),approx_se_vx_mps=float(se[2]),approx_se_vy_mps=float(se[3]),covariance_scope='NAIVE_LOCAL_DIAGNOSTIC; SYSTEMATIC/NAV_ERROR_NOT_ACCOUNTED')
 if not fit.success or not np.isfinite(state).all():raise ValueError('solver nonconvergence/nonfinite')
 return state,diag
