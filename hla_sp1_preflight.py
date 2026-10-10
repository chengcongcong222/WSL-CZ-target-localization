"""Once bounded static provider preparation and pre-A algebra controls."""
from hla_sp1_core import *
import subprocess,shutil

def main():
 if (OUT/'PREFLIGHT.json').exists():raise RuntimeError('preflight already finished')
 assert subprocess.check_output(['git','rev-parse','HEAD'],text=True).strip()==PARENT
 assert subprocess.check_output(['git','ls-remote','origin','refs/heads/main'],text=True).split()[0]==PARENT
 assert sha(OUT/'DELEGATED_TASK.md')=='ac32c7f9c62c6cc19513b670bf1398e0f73736ca95c19b296fc1e182cd490efb'
 md=LOCAL/'modal';md.mkdir(exist_ok=True);template=ROOT/'results/R3_C2_Yang_SA_depth/R3_C2_1/_kraken_zgrid/zgrid_f201.env';lines=template.read_text(encoding='utf-8-sig').splitlines();boundary=next(i for i,l in enumerate(lines) if l.strip().startswith("'R'"));at=ROOT/'tools/acoustics_toolbox/atWin10/at/bin';plans=[]
 for mesh in [80001,160001]:
  for f in FREQ:
   stem=f'sp1_f{int(f)}_n{mesh}';env=["'HLA_SP1_ESTD'",str(f),'1',lines[3],f'{mesh} 0.0 5000.0']+lines[5:boundary]+["'R' 0.0",'1500.0 1800.0','60.0','1','200.0',str(len(ZS))]+[str(z) for z in ZS]+['R'];(md/(stem+'.env')).write_text('\n'.join(env)+'\n',encoding='utf-8');plans.append((mesh,f,stem))
 dump(OUT/'PRECOMPUTE_CONTRACT.json',{'time_utc':stamp(),'six_static_KRAKEN_max':6,'FIELD_max':3,'timeout_KRAKEN_s':90,'timeout_FIELD_s':60,'env_sha256':{s+'.env':sha(md/(s+'.env')) for _,_,s in plans},'kraken_sha256':sha(at/'kraken.exe'),'field_sha256':sha(at/'field.exe'),'physics_gate':'exact source depths, nonfinite/attenuation checks; independent FIELD shape <=0.002; algebra <=1e-10; waveform relative <=1e-10','no_dynamic_provider':True})
 calls=[]
 for mesh,f,stem in plans:
  guard();beg=time.monotonic();result=subprocess.run([str(at/'kraken.exe'),stem],cwd=md,capture_output=True,timeout=90);(md/(stem+'.stdout.txt')).write_bytes(result.stdout+result.stderr);m=parse_mod(md/(stem+'.mod'));assert result.returncode==0;calls.append(dict(mesh=mesh,frequency_hz=f,new_KRAKEN=True,returncode=result.returncode,seconds=time.monotonic()-beg,modes=len(m[2]),mod_sha256=sha(md/(stem+'.mod')),env_sha256=sha(md/(stem+'.env'))));table(OUT/'SOLVER_CALLS.csv',calls);print('STATIC KRAKEN',len(calls),'/6',flush=True)
 field=[]
 for f in FREQ:
  stem=f'sp1_check_f{int(f)}';base=md/f'sp1_f{int(f)}_n160001';dest=md/stem;dest.with_suffix('.mod').write_bytes(base.with_suffix('.mod').read_bytes());dest.with_suffix('.env').write_bytes(base.with_suffix('.env').read_bytes());dest.with_suffix('.flp').write_text('\n'.join(["'SP1_FIELD'","'RA'",'9999','1','0.0','8','54.993 55.007 /','1','200.0','1','200.0','1','0.0 /'])+'\n',encoding='utf-8');res=subprocess.run([str(at/'field.exe'),stem],cwd=md,capture_output=True,timeout=60);dest.with_suffix('.stdout.txt').write_bytes(res.stdout+res.stderr);assert res.returncode==0;b=dest.with_suffix('.shd').read_bytes();rec=4*int(np.frombuffer(b[:4],'<i4')[0]);hdr=np.frombuffer(b[2*rec:3*rec],'<i4');ns,nrd,nr=map(int,hdr[4:7]);assert (ns,nrd,nr)==(1,1,8);rr=np.frombuffer(b[9*rec:10*rec],'<f4')[:nr].astype(float);obs=np.frombuffer(b[10*rec:10*rec+nr*8],'<c8').astype(complex);pred=pressure(parse_mod(base.with_suffix('.mod')),rr,200);common=np.vdot(pred,obs)/np.vdot(pred,pred);e=float(np.linalg.norm(obs-common*pred)/np.linalg.norm(obs));field.append(dict(frequency_hz=f,new_FIELD=True,shape_relative_error=e,limit=.002,PASS=e<=.002,actual_ranges_m=json.dumps(rr.tolist())));table(OUT/'FIELD_CONTROLS.csv',field)
 rng=np.random.default_rng(2026101017);panel=[dict(geometry=i,r_km=rng.uniform(45.2,59.8),theta_deg=rng.uniform(-4.8,4.8),z_label=float(rng.choice(ZS))) for i in range(16)];table(OUT/'CALIBRATION_PANEL.csv',panel)
 high=Provider(160001);low=Provider(80001);states=np.array([[p['r_km'],p['theta_deg'],p['z_label']] for p in panel]);hp=high.field(states,exact=True);ref=float(np.median(np.mean(abs(hp)**2,axis=(1,2))));checks=[]
 def add(n,e,t=1e-10):checks.append(dict(control=n,error=float(e),limit=t,PASS=bool(e<=t)))
 rng=np.random.default_rng(np.random.SeedSequence(2026101020));h=rng.normal(size=(3,8))+1j*rng.normal(size=(3,8));source=(.7+rng.random((8,3,1)))*np.exp(1j*rng.uniform(-np.pi,np.pi,(8,3,1)));y=h[None]*source;gain=10**(rng.uniform(-1,1,8)/20)*np.exp(1j*np.radians(rng.uniform(-30,30,8)))
 for method in range(3):
  add('unknown_source_P'+str(method),abs(float(score(h,y,method))-1));c=moment(y,method);ft=realfeatures(transform(h,method));sc=float(np.sum(ft*realmoment(c))/3);add('moment_formula_P'+str(method),abs(sc-float(score(h,y,method))))
 add('P2_crossfreq_gain_invariance',abs(float(score(h,y*gain,2))-1));add('P1_amplitude_NONinvariance_expected',0 if float(score(h,y*gain,1))<1-1e-8 else 1)
 pw=np.exp(2j*np.pi*FREQ[:,None]*OFF[None,:]/1500);a=pw*np.array([.2+1j,3-1j,.3-.5j])[:,None]
 for m in range(3):add('plane_wave_range_depth_null_P'+str(m),abs(float(score(pw,a[None],m))-1))
 # single-element direct normalized identity
 add('single_element_no_resolution',abs(abs(np.conj(2+3j)*(4-1j))**2/(abs(2+3j)**2*abs(4-1j)**2)-1))
 n=np.arange(2048);tone=np.exp(2j*np.pi*FREQ[:,None]*n[None,:]/1024);x=2*np.real(np.einsum('kfm,fn->knm',y,tone));dft=np.einsum('knm,fn->kfm',x,tone.conj())/2048;add('real_waveform_DFT_scale',np.linalg.norm(dft-y)/np.linalg.norm(y));basis=np.stack([np.cos(2*np.pi*f*n/1024) for f in FREQ]+[np.sin(2*np.pi*f*n/1024) for f in FREQ],axis=1);fit=np.linalg.lstsq(basis,x[0],rcond=None)[0];coeff=(fit[:3]-1j*fit[3:])/2;add('independent_sine_fit',np.linalg.norm(coeff-y[0])/np.linalg.norm(y[0]));noise=rng.normal(size=(4096,2048))*np.sqrt(2048*ref/100);nb=noise@tone[0].conj()/2048;add('noise_bin_variance_empirical',abs(float(np.mean(abs(nb)**2))/(ref/100)-1),.08)
 d=[y[:,p]*y[:,q].conj() for p,q in PAIRS];add('shared_pair_cycle',np.linalg.norm(d[0]*d[1]-d[2]*abs(y[:,1])**2)/np.linalg.norm(d[2]*abs(y[:,1])**2))
 # fixed off-grid distances, all21 depths and8elements, plus same score checks
 tests=np.array([[r,t,z] for r,t in [(45.271,.17),(50.123,2.31),(55.321,4.73),(59.817,-1.29)] for z in ZS]);interp=[]
 for mesh,p in [(80001,low),(160001,high)]:
  direct=p.field(tests,exact=True);approx=p.field(tests);e=np.linalg.norm(approx-direct,axis=-1)/np.linalg.norm(direct,axis=-1);add('interpolation_full_vector_n'+str(mesh),e.max(),1e-5);add('theta_mirror_n'+str(mesh),np.max(abs(p.field(tests*np.array([1,-1,1]))-approx)),1e-12)
  for m in range(3):
   for i in range(len(tests)):
    diff=abs(float(score(approx[i],direct[i][None],m))-1);interp.append(dict(mesh=mesh,test=i,method=m,max_raw_vector_relative=float(e[i].max()),score_error=diff))
 table(OUT/'INTERPOLATION_PREFLIGHT.csv',interp);table(OUT/'ALGEBRA_AND_WAVEFORM_CONTROLS.csv',checks)
 success=all(r['PASS'] for r in checks+field);dump(OUT/'PREFLIGHT.json',{'PASS':success,'checks':checks,'FIELD':field,'new_KRAKEN':len(calls),'new_FIELD':len(field),'reference_power':ref,'nu_by_reference_SNR_dB':{'20':ref/100,'5':ref/10**.5},'budgets':guard(),'failure_action':'no test execution if false'})
 dump(OUT/'PROVIDER_PROVENANCE.json',{'template':str(template.relative_to(ROOT)),'template_sha256':sha(template),'KRAKEN_sha256':sha(at/'kraken.exe'),'FIELD_sha256':sha(at/'field.exe'),'calls':calls,'exact_depths_m':ZS.tolist(),'receiver_z_m':200,'complex_pressure':'sum phi(zs)*phi(zr)*sqrt(2pi/(complex k*r))*exp(-i k*r-i pi/4), all modes','time_convention':'exp(+i omega t), attenuation k.imag<=0','window':[1500,1800],'CHIGH_truncation_UNCERTIFIED':True,'old_E2_FAIL_UNCHANGED':True,'grid_spline_step_m':1,'carrier':'midpoint returned min/max real k; restore exp(-i carrier*r)','physical_model':'STATIC_HELMHOLTZ_E0_MATCHED_NOMINAL_ENVIRONMENT','archive_qualification':'existing201/235/283 modes not at both registered meshes or exact21 depths; no qualified full six cache'})
 print('PREFLIGHT',success,'reference power',ref,flush=True)
 if not success:raise RuntimeError('SP1 preflight failed; do not execute test')
if __name__=='__main__':main()
