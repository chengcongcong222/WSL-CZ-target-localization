"""Registered E0 geometry/format controls. No target inference."""
from env_provider import *
def main():
 D=W/'e0_controls';D.mkdir(exist_ok=True)
 marker=D/'SCIENCE_STARTED.json'
 if marker.exists():raise RuntimeError('No rerun permitted')
 js(D/'CONTROL_FREEZE.json',{'straight_ray_z_tol_m':0.1,'constant_medium_TL_slope_abs_tol_dB':0.5,'two_beam_counts_TL_max_dB':0.5,'gradient_integral_relative_tol':1e-9,'calls':3,'no_boundary_contact':True,'control_range_km':[10,20],'H_m':10000,'source_z_m':5000})
 js(marker,{'code_sha256':sha(__file__),'provider_sha256':sha(Path(__file__).with_name('env_provider.py')),'purpose':'pre-freeze deterministic analytic implementation controls; counted in campaign solver budget'})
 checks=[]
 curve=return_curve(np.array([0.,1000.]),np.array([1500.,1520.]),0.,[5.])
 exact=2*1500*np.tan(np.deg2rad(5))/.02/1000
 checks.append({'check':'C_linear_exact_turn_integral','PASS':abs(curve[0]['first_return_km']/exact-1)<1e-9,'relative_error':abs(curve[0]['first_return_km']/exact-1)})
 assert return_curve(np.array([0.,5000.]),np.array([1500.,1500.]),200.,[1,5,10])==[]
 checks.append({'check':'constant_speed_no_deep_refractive_return','PASS':True})
 p=write_env(D,'straight',np.array([0.,10000.]),np.array([1500.,1500.]),10000,5000,'R',3,25,(-1,1),20,[5000]);execute(p,'E0_CONTROLS')
 rays=readrays(p.with_suffix('.ray'))
 for a,ns,nb,x in rays:
  err=np.max(np.abs(x[:,-1]-(5000+x[:,0]*np.tan(np.deg2rad(a)))))
  checks.append({'check':f'straight_ray_{a}','PASS':err<=.1 and ns==0 and nb==0,'max_z_error_m':err})
 curves=[]
 for b in [801,1601]:
  p=write_env(D,f'open_{b}',np.array([0.,10000.]),np.array([1500.,1500.]),10000,5000,'I',b,25,(-5,5),20,[5000]);execute(p,'E0_CONTROLS')
  rr,zz,pp=readshd(p.with_suffix('.shd'));tl=-20*np.log10(abs(pp[0]));curves.append(tl)
  delta=(tl[-1]-tl[0])-20*np.log10(rr[-1]/rr[0]);checks.append({'check':f'open_TL_spherical_spreading_{b}','PASS':abs(delta)<=.5,'error_dB':delta})
 delta=np.max(abs(curves[0]-curves[1]));checks.append({'check':'beam_control_convergence','PASS':delta<=.5,'max_dB':delta})
 js(O/'E0_PROVIDER_CONTROL_VALIDATION.json',{'PASS':all(c['PASS'] for c in checks),'checks':checks,'external_solver_calls':3})
 print(json.dumps(checks,indent=2))
 if not all(c['PASS'] for c in checks):sys.exit(2)
if __name__=='__main__':main()
