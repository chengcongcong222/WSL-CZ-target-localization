"""E2 archived paper controls + analytic source gauges, not CZ signal generation."""
from env_provider import *
from scipy.integrate import trapezoid
import struct

def main():
 D=W/'E2';D.mkdir(exist_ok=True);marker=D/'SCIENCE_STARTED.json'
 if marker.exists():raise RuntimeError('No repeat of E2 unit')
 freeze=json.loads((O/'E2_UNIT_FREEZE.json').read_text());assert all(sha(R/a['path'])==a['sha256'] for a in freeze['inputs'])
 js(marker,{'config_sha256':sha(O/'E2_UNIT_FREEZE.json'),'code_sha256':sha(__file__),'source_inputs':freeze['inputs'],'seed':'DETERMINISTIC','purpose':'stored shallow oracle control / analytic necessary conditions, no CZ pressure'})
 p=R/'results/R3_C2_Yang_SA_depth/R3_C2_2C';raw=(p/'_kraken_paper/yang2014_f350.mod').read_bytes();recl=4*struct.unpack_from('<i',raw)[0];nz=struct.unpack_from('<i',raw,92)[0];M=struct.unpack_from('<i',raw,5*recl)[0]
 zz=np.frombuffer(raw,dtype='<f4',count=nz,offset=4*recl).astype(float);phi=np.stack([np.frombuffer(raw,dtype='<c8',count=nz,offset=(7+m)*recl) for m in range(M)],axis=1).astype(complex)
 km=np.frombuffer(raw,dtype='<c8',count=M,offset=(7+M)*recl).astype(complex);assert M==23 and nz==89 and np.array_equal(zz,np.arange(89))
 peaks=list(csv.DictReader((p/'paper_mode_identification.csv').open()));ref=list(csv.DictReader((p/'paper_depth_ambiguity.csv').open()));rows=[];checks=[];scorecurves=[]
 for zs in [4.,50.]:
  for zr in [18.,70.]:
   saved=[v for v in peaks if v['method']=='ORACLE_MODE_ID' and float(v['zs'])==zs and float(v['zr'])==zr];saved=sorted(saved,key=lambda a:int(a['mode_id']));assert len(saved)==M
   gm=np.array([complex(float(a['g_real']),float(a['g_imag'])) for a in saved]);pr=phi[int(zr)];
   for reg in [.05,.1,.2]:
    delta=reg*abs(phi).max();weights=pr/(pr**2+delta**2);s=abs(phi@(gm*weights))**2;s/=s.sum();zhat=float(zz[np.argmax(s)])
    r0=next(a for a in ref if a['route']=='THEORY_CONTROL' and float(a['zs_true'])==zs and float(a['zr'])==zr and float(a['Delta_ratio'])==reg)
    checks.append({'check':'archived_shallow_depth_score','zs':zs,'zr':zr,'reg':reg,'PASS':abs(zhat-float(r0['z_hat']))<1e-9 and abs(s.max()-float(r0['global_peak_value']))<1e-9,'max_score_difference':abs(s.max()-float(r0['global_peak_value']))})
    for z,val in zip(zz,s):scorecurves.append({'zs':zs,'zr':zr,'regularization_ratio':reg,'z_candidate_m':z,'score':val})
   rows.append({'env_id':'YANG_SHALLOW_H88_NON_CZ','method':'YANG_KNOWN_MODE_DELTA0_STORED_CONTROL','zs_m':zs,'zr_m':zr,'Lr_m':4990,'total_modes':M,'evidence':'L4_ORACLE_INFORMATION_HISTORICAL_RECONSTRUCTION','status':'RECONFIRMED_ON_STORED_MODAL_OUTPUT','nonoracle_extracted_modes':'NOT_EVALUATED','z_hat_m':zhat,'why_not_CZ_demo':'Shallow model; known true mode IDs, CW frequency, range increments and delta=0'})
 # Eq(1)->single-mode Eq(5), independently stable exponential interval integral including alpha=0.
 rr=np.arange(5010.,10000.+.0001,.625)
 for m in range(8):
  k=float(km[m].real);a=float(abs(km[m].imag));pressure=np.sqrt(2*np.pi/(k*rr))*np.exp(-1j*k*rr-a*rr-1j*np.pi/4)
  numerical=trapezoid(pressure*np.exp(1j*k*rr)*np.sqrt(rr),rr)*np.exp(1j*np.pi/4)/np.sqrt(2*np.pi*k)
  analytic=np.exp(-a*5010)*(-np.expm1(-a*4990))/(a*k) if a else 4990/k
  err=abs(numerical-analytic)/abs(analytic);checks.append({'check':'single_mode_Eq1_to_Eq5','mode':m+1,'PASS':err<=1e-6,'relative_error':err})
 # Unknown common complex source: normalized plane-wave shape has no r/z dependence.
 x=np.arange(8)*2-7;k=2*np.pi*200/1500;u=np.exp(1j*k*x*.1);u/=np.linalg.norm(u)
 pp=np.outer(u,u.conj());other=np.exp(.31j)*u;err=np.linalg.norm(pp-np.outer(other,other.conj()));checks.append({'check':'plane_wave_unknown_source_range_depth_zero','PASS':err<1e-12,'error':err})
 # Reflection ambiguity for any coplanar x-z array in horizontally isotropic field.
 for ay in [0.]:
  dx=50000-x;d1=np.sqrt(dx*dx+5000**2+15**2);d2=np.sqrt(dx*dx+(-5000)**2+15**2);err=np.max(abs(d1-d2));checks.append({'check':'HLA_left_right_mirror_zero','PASS':err==0,'error_m':err})
 # Exact frequency gauge at constant radial rate, source frequency unknown.
 f=200.;q=1.;q2=1.4;f2=f*(1-q/1500)/(1-q2/1500);err=abs(f*(1-q/1500)-f2*(1-q2/1500));checks.append({'check':'unknown_source_frequency_constant_q_gauge','PASS':err<1e-12,'error_Hz':err})
 # Independent arbitrary modal complex gains absorb source-depth changes (nonzero modes).
 m=np.arange(1,21);phi1=np.sin((m-.5)*np.pi*185/5000);phi2=np.sin((m-.5)*np.pi*215/5000);g=np.exp(.2j*m)*phi1/phi2;err=np.max(abs(g*phi2-np.exp(.2j*m)*phi1));checks.append({'check':'free_modal_gain_depth_absorption','PASS':err<1e-12,'error':err})
 # Erratum phase: density unchanged, relative modal phases generally change.
 g0=np.array([1.,.7]);ks=np.array([.8,.81]);gd=g0*np.exp(1j*ks*100);err=np.max(abs(abs(gd)-abs(g0)));checks.append({'check':'Erratum_density_invariant_phase_not','PASS':err<1e-12 and abs(abs(gd.sum())**2-abs(g0.sum())**2)>.1,'density_error':err,'coherent_score_delta':abs(abs(gd.sum())**2-abs(g0.sum())**2)})
 # Five overlapping differences share six window errors; preserve covariance and DPI.
 A=np.diff(np.eye(6),axis=0);cov=A@A.T;j=np.arange(6.)[:,None]**2;rawinfo=float(j.T@j);info=float((A@j).T@np.linalg.solve(cov,A@j));checks.append({'check':'shared_window_noise_and_DPI','PASS':np.allclose(np.diag(cov,1),-1) and info<=rawinfo+1e-12,'raw_information':rawinfo,'difference_information':info,'neighbor_correlation':-.5})
 for env in sorted(set(a['env_id'] for a in json.loads((O/'E0_FREEZE.json').read_text())['E1_cases'])):
  for method in ['YANG_CZ_A_ORACLE','YANG_CZ_B_FINITE_APERTURE','YANG_CZ_C_NONORACLE','LIANG_MMAR_CZ_AR_HANKEL']:
   rows.append({'env_id':env,'method':method,'evidence':'NOT_EVALUATED','status':'NOT_EVALUATED','why_not_CZ_demo':'No L3 complete coherent modal provider for new SSP/bottom identity; E1 intensity/rays cannot synthesize admitted modal-pressure truth','nonoracle_extracted_modes':'NOT_EVALUATED'})
 table(O/'E2_YANG_MMAR_COMPARISON.csv',rows);table(O/'E2_PAPER_DEPTH_CONTROL_CURVES.csv',scorecurves);table(O/'E2_ANALYTIC_AND_PAPER_CONTROLS.csv',checks)
 js(O/'E2_METHOD_GATE.json',{'analytic_and_stored_control_PASS':all(a['PASS'] for a in checks),'checks':len(checks),'Yang_shallow_control':'RECONFIRMED_STORED_ORACLE_4_CASES','MMAR_identity':'PAPER_IMPLEMENTATION_PARTIAL_HISTORICAL','CZ_L3_coherent_provider':'NOT_ADMITTED','CZ_L5_nonoracle_extraction':'NOT_EVALUATED','CZ_L6_estimator_performance':'NOT_EVALUATED','allow_E3_CZ_performance':False,'reasons':['New E1 provider computes rays and incoherent band energy only','No complete verified complex eigenfunctions/wavenumbers plus finite boundary/mode-completeness certificate at new SSP identity','Prior R3/E-STD oracle outputs and T2 failures cannot certify new environmental modal truth'],'new_solver_calls':0,'new_MC':0})
 print('E2 controls:',len(checks),'PASS=',all(a['PASS'] for a in checks),'CZ extraction NOT_EVALUATED',flush=True)
if __name__=='__main__':main()
