"""Bounded, unregularized six-state observation-only candidate."""
import numpy as np
from scipy.optimize import least_squares
from r4_a1_new_estimator import initialize,wrap,jacobian as four_jacobian
SCALE4=np.array([50000.,50000.,2.,2.])
def residual(q,times,node_positions,bearings,sigma_deg,common_bound_deg,diff_bound_deg):
 state=q[:4]*SCALE4;d=(state[:2]+times[:,None]*state[2:])[:,None,:]-node_positions;pred=np.arctan2(d[:,:,1],d[:,:,0])+np.deg2rad(q[4]*common_bound_deg+q[5]*diff_bound_deg*np.array([-1.,1.]));return (wrap(pred-bearings)/np.deg2rad(sigma_deg)).ravel()
def jacobian(q,times,node_positions,bearings,sigma_deg,common_bound_deg,diff_bound_deg):
 J4=four_jacobian(q[:4],times,node_positions,bearings,sigma_deg);return np.column_stack([J4,np.full(2*len(times),common_bound_deg/sigma_deg),np.tile([-diff_bound_deg/sigma_deg,diff_bound_deg/sigma_deg],len(times))])
def structure(q,times,node_positions,bearings,sigma_deg,common_bound_deg,diff_bound_deg):
 out=[]
 for segment,mask in [('STRAIGHT',times<=600),('POST_TURN',times>600),('FULL',np.ones(len(times),dtype=bool))]:
  J=jacobian(q,times[mask],node_positions[mask],bearings[mask],sigma_deg,common_bound_deg,diff_bound_deg);sv=np.linalg.svd(J,compute_uv=False);ratio=float(sv[-1]/sv[0]);rank=int(np.sum(sv/sv[0]>1e-10));row=dict(segment=segment,n_epochs=int(mask.sum()),numerical_rank=rank,full_rank=rank==6,smallest_largest_ratio=ratio,condition_number=float(sv[0]/sv[-1]) if sv[-1]>0 else 'INF')
  row.update({f'singular_value_{i+1}':float(v) for i,v in enumerate(sv)});out.append(row)
 return out

def estimate_bias_aware(times,node_positions,bearings,sigma_deg,common_bound_deg,diff_bound_deg):
 x,init=initialize(times,node_positions,bearings);q0=np.r_[x/SCALE4,0.,0.];args=(times,node_positions,bearings,sigma_deg,common_bound_deg,diff_bound_deg)
 fit=least_squares(residual,q0,jac=jacobian,args=args,bounds=(np.r_[np.full(4,-np.inf),-1.,-1.],np.r_[np.full(4,np.inf),1.,1.]),method='trf',loss='soft_l1',f_scale=1.5,max_nfev=300,ftol=1e-10,xtol=1e-10,gtol=1e-10)
 state=fit.x[:4]*SCALE4;diag=dict(**init,initial_x_m=float(x[0]),initial_y_m=float(x[1]),initial_vx_mps=float(x[2]),initial_vy_mps=float(x[3]),initial_common_bias_deg=0.,initial_diff_bias_deg=0.,solver_success=bool(fit.success and np.isfinite(fit.x).all()),solver_status=int(fit.status),nfev=int(fit.nfev),cost=float(fit.cost),optimality=float(fit.optimality),residual_RMS_deg=float(np.sqrt(np.mean(residual(fit.x,*args)**2))*sigma_deg),common_hat_deg=float(fit.x[4]*common_bound_deg),diff_HALF_hat_deg=float(fit.x[5]*diff_bound_deg),common_at_bound=bool(abs(fit.x[4])>=1-1e-6),diff_at_bound=bool(abs(fit.x[5])>=1-1e-6))
 return state,diag,structure(fit.x,*args)
