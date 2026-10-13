"""Cold N3 verifier: independent struct reader, brent roots, J0-iY0 radial kernel.
Does NOT import t2_n3_common or main field/shape/reference functions.
"""
from pathlib import Path
import csv,json,hashlib,struct,math,subprocess
import numpy as np
from scipy.optimize import brentq
from scipy.special import jv,yv
ROOT=Path(__file__).resolve().parent;O=ROOT/'results/R4_T2_N3_BOUNDARY_PHYSICS_ANALYTIC_REFERENCE';L=Path('D:/ProjectStorage/WSL-CZ/T2_N3')
checks=[];max_field_delta=0.;cache={};rootcache={}
def rows(p):
 with Path(p).open(encoding='utf-8-sig',newline='') as f:return list(csv.DictReader(f))
def js(p):return json.loads(Path(p).read_text(encoding='utf-8-sig'))
def hashfile(p):
 h=hashlib.sha256()
 with Path(p).open('rb') as f:
  for x in iter(lambda:f.read(1048576),b''):h.update(x)
 return h.hexdigest()
def check(kind,name,value,limit=0.):checks.append({'kind':kind,'identity':name,'difference':float(value),'tolerance':limit,'PASS':bool(math.isfinite(value) and value<=limit)})
def num(a,b):return abs(float(a)-float(b))
def modal(p):
 p=Path(p)
 if str(p) in cache:return cache[str(p)]
 bb=p.read_bytes();r=4*struct.unpack_from('<i',bb)[0];nf,nmedia,nz,np_=struct.unpack_from('<4i',bb,84);frequency=struct.unpack_from('<d',bb,3*r)[0];count=struct.unpack_from('<i',bb,5*r)[0]
 depth=np.array(struct.unpack_from('<'+str(nz)+'f',bb,4*r));pm=[]
 for n in range(count):
  v=np.array(struct.unpack_from('<'+str(2*np_)+'f',bb,(7+n)*r));pm.append(v[::2]+1j*v[1::2])
 vv=np.array(struct.unpack_from('<'+str(2*count)+'f',bb,(7+count)*r));k=vv[::2]+1j*vv[1::2]
 result=(depth,np.array(pm).T,k,frequency);cache[str(p)]=result;return result

def roots(f,bc,high):
 key=(f,bc,high)
 if key in rootcache:return rootcache[key]
 offset=.5 if bc=='PR_R' else 1.;q=2*f*5000/1500
 if high!='ALL_PROPAGATING':q*=math.sqrt(1-(1500/float(high))**2)
 n=max(0,math.ceil(q-offset));fn=math.cos if bc=='PR_R' else math.sin
 x=np.array([brentq(fn,(i+offset-.24)*math.pi,(i+offset+.24)*math.pi,xtol=1e-13) for i in range(n)])
 mu=x/5000;w=2*math.pi*f;k=np.sqrt((w/1500)**2-mu**2);rootcache[key]=(mu,k);return mu,k

a={}
for row in rows(O/'package/T2_TILT_G0_ARRAY_GEOMETRY_20261013_v1.csv'):
 if row['resource'] in ('A0_0','A_P4M','B1_V4','B2_DUP'):a.setdefault(row['resource'],[]).append((int(row['channel']),*[float(row[k]) for k in ('x_m','y_m','z_m')]))
for key in a:a[key]=np.array([r[1:] for r in sorted(a[key])])
probes={'P01':(45.5,0.,180.),'P02':(50.,0.,200.),'P03':(55.,2.5,200.),'P04':(58.,4.5,210.),'P05':(59.5,5.,220.)}
def near(d,z):
 n=int(np.argmin(abs(d-z)))
 if abs(d[n]-z)>1e-5:raise ValueError('cold depth mismatch')
 return n

def pressure(f,bc,high,pos,coords,far=False):
 mu,k=roots(f,bc,high);r,t,z=pos;xy=(1000*r*math.cos(math.radians(t)),1000*r*math.sin(math.radians(t)));out=[]
 for x,y,d in coords:
  rho=math.hypot(xy[0]-x,xy[1]-y);argument=k*rho
  radial=np.sqrt(2*math.pi/argument)*np.exp(-1j*argument-1j*math.pi/4) if far else -1j*math.pi*(jv(0,argument)-1j*yv(0,argument))
  out.append(np.dot(.0004*np.sin(mu*z)*np.sin(mu*d),radial))
 return np.array(out)
def mpressure(mod,pos,coords,mask=None):
 d,pm,k,f=mod
 if mask is None:mask=np.ones(len(k),bool)
 r,t,z=pos;xy=(1000*r*math.cos(math.radians(t)),1000*r*math.sin(math.radians(t)));kk=k[mask];out=[]
 for x,y,zr in coords:
  rho=math.hypot(xy[0]-x,xy[1]-y);w=pm[near(d,z),mask]*pm[near(d,zr),mask];rad=np.sqrt(2*math.pi/(kk*rho))*np.exp(-1j*kk*rho-1j*math.pi/4);out.append(np.dot(w,rad))
 return np.array(out)
def error(v,w):
 if np.linalg.norm(v)==0 or np.linalg.norm(w)==0:raise ValueError('undefined norm')
 Q=np.linalg.qr(w.reshape(-1,1),mode='reduced')[0];return float(np.linalg.norm(v-Q@(Q.conj().T@v))/np.linalg.norm(v))

def main():
 global max_field_delta
 for r in js(O/'DESIGN_FREEZE.json')['frozen_files']:check('frozen_input',r['path'],0 if hashfile(r['path'])==r['sha256'] else 1)
 for r in js(O/'OLD_N2_SOURCE_IDENTITIES.json')['records']:
  check('old_source',r['path'],0 if hashfile(r['path'])==r['sha256'] else 1)
  if 'parent_blob_sha256' in r:
   path=Path(r['path']).relative_to(ROOT).as_posix();blob=subprocess.check_output(['git','show',js(O/'DESIGN_FREEZE.json')['parent']+':'+path]);check('parent_blob',path,0 if hashlib.sha256(blob).hexdigest()==r['parent_blob_sha256'] else 1)
 seal=js(O/'SCIENCE_OUTPUT_FREEZE.json')
 for section in ('raw_files','scientific_tables'):
  for p,h in seal[section].items():check('sealed_science',p,0 if hashfile(p)==h else 1)
 for r in rows(O/'package/ANALYTIC_WAVEGUIDE_MODE_COUNTS.csv'):
  count=len(roots(float(r['frequency_hz']),r['boundary'],r['chigh_mps'])[0]);check('preview_modes',str(r),abs(count-int(r['analytic_modes'])))
 for r in rows(O/'package/ANALYTIC_WAVEGUIDE_WINDOW_SHAPES.csv'):
  p=probes[r['probe']];rx=a[r['resource']];f=float(r['frequency_hz']);v=pressure(f,r['boundary'],r['left_chigh_mps'],p,rx);w=pressure(f,r['boundary'],r['right_chigh_mps'],p,rx);delta=num(error(v,w),r['shape_difference']);check('preview_window',str(r),delta,1e-8)
 for r in rows(O/'package/ANALYTIC_EXACT_RADIAL_VS_FARFIELD.csv'):
  f=float(r['frequency_hz']);rx=a[r['resource']];p=probes[r['probe']];v=pressure(f,r['boundary'],4000,p,rx);w=pressure(f,r['boundary'],4000,p,rx,True);check('preview_radial',str(r),num(error(v,w),r['exact_Hankel_vs_farfield_shape']),1e-8)
 old=lambda v,f:Path('D:/ProjectStorage/WSL-CZ/T2_TILT_G0/modal')/f't2_f{f}_n160001.mod' if v=='K0' else Path('D:/ProjectStorage/WSL-CZ/T2_N2/K')/v/f'n2_{v}_f{f}.mod'
 env=old('K0',201).with_suffix('.env').read_text().splitlines();bottom=next(i for i in range(5,len(env)) if env[i].startswith("'R'"));ssp=np.array([[float(v) for v in r.split()] for r in env[5:bottom]]);cmax=ssp[:,1].max()
 for r in rows(O/'SSP_AND_TURNING_CLASSIFICATION.csv'):
  f=float(r['frequency_hz']);d,pm,k,mf=modal(old(r['variant'],int(f)));j=int(r['mode_index'])-1;cp=2*math.pi*f/k[j].real;kz2=(2*math.pi*f/ssp[:,1])**2-(2*math.pi*f/cp)**2;check('cp_classification',r['variant']+'_'+str(f)+'_'+str(j),max(num(cp,r['cp_mps']),num(kz2.min(),r['min_kz2_at_SSP_nodes']),num(kz2.max(),r['max_kz2_at_SSP_nodes'])),1e-10);check('no_turning_boolean',str(r['mode_index'])+r['variant']+str(f),0 if str(bool(cp>cmax and np.all(kz2>0)))==r['no_internal_water_turning'] else 1)
 bandfields={};sents={r['sentinel_id']:tuple(float(r[k]) for k in ('r_km','abs_theta_deg','source_z_m')) for r in rows(O/'package/T2_N2_SENTINELS_20261013_v1.csv')}
 stored=np.load(L/'OLD_MODAL_BAND_FIELDS.npz')
 for r in rows(O/'OLD_N2_MODAL_BAND_CONTRIBUTIONS.csv'):
  key=f"{r['variant']}_{r['frequency_hz']}_{r['sentinel_id']}_{r['resource']}"
  if key not in bandfields:
   mod=modal(old(r['variant'],int(r['frequency_hz'])));cp=2*math.pi*mod[3]/mod[2].real;groups=np.where(cp<=cmax,0,np.where(cp<=1800,1,np.where(cp<=2000,2,3)));full=mpressure(mod,sents[r['sentinel_id']],a[r['resource']]);parts=[mpressure(mod,sents[r['sentinel_id']],a[r['resource']],groups==i) for i in range(4)];bandfields[key]=(full,parts)
   allp=np.stack([full]+parts);check('old_complex_band_fields',key,float(np.linalg.norm(allp-stored[key])/np.linalg.norm(allp)),1e-10)
  full,parts=bandfields[key];part=parts[int(r['band_index'])];fn=np.linalg.norm(full);pn=np.linalg.norm(part);cancel=fn/sum(np.linalg.norm(x) for x in parts);check('old_band_metrics',key+'_'+r['band_index'],max(num(fn,r['full_norm']),num(pn,r['band_norm']),num(pn/fn,r['band_over_full_norm']),num(cancel,r['complex_cancellation_full_over_sum_band_norms'])),1e-8)
 calls=rows(O/'SOLVER_CALLS.csv');by={r['stem']:r for r in calls};nconfig=0;crossrows=rows(O/'ANALYTIC_KRAKEN_MODE_AND_FIELD_CROSSCHECK.csv');cv={};mainfields={}
 for r in crossrows:
  if 'stem' not in r:continue
  c=by[r['stem']];mod=modal(Path(c['cwd'])/(r['stem']+'.mod'));f=float(r['frequency_hz']);bc=r['boundary'];high=float(r['CHIGH']);mu,k=roots(f,bc,high);n=min(len(k),len(mod[2]));krerr=float(np.max(abs(mod[2].real[:n]-k[:n])/k[:n]));pos=probes[r['probe']];rx=a[r['resource']];aa=pressure(f,bc,high,pos,rx);af=pressure(f,bc,high,pos,rx,True);bb=mpressure(mod,pos,rx);e=error(aa,bb);far=error(aa,af);delta=num(e,r['spatial_shape_error']);max_field_delta=max(max_field_delta,delta);check('new_analytic_field',r['stem']+'_'+r['probe']+'_'+r['resource'],delta,1e-8);check('new_kr_max',r['stem']+'_'+r['probe']+'_'+r['resource'],num(krerr,r['max_kr_relative_error']),1e-10)
  actualpass=(len(k)==len(mod[2]) and krerr<=2e-5 and e<=.002 and far<=1e-4);check('analytic_gate_reporting',r['stem']+'_'+r['probe']+'_'+r['resource'],0 if str(actualpass)==r['PASS'] else 1)
  cv.setdefault(r['stem'],[]).append(actualpass);key=(bc,int(f),int(high),int(r['mesh']),r['probe']+'_'+r['resource']);mainfields[key]=bb;npz=np.load(Path(c['cwd'])/'FIELDS.npz');check('new_saved_field',str(key),float(np.linalg.norm(bb-npz[r['probe']+'_'+r['resource']])/np.linalg.norm(bb)),1e-10)
 for stem,flags in cv.items():
  nconfig+=1;check('config_gate',stem,0 if str(all(flags))==by[stem]['analytic_gate_PASS'] else 1)
 for r in rows(O/'ANALYTIC_MODAL_ERRORS.csv'):
  if 'stem' not in r or r.get('matched')!='True':continue
  c=by[r['stem']];mod=modal(Path(c['cwd'])/(r['stem']+'.mod'));mu,k=roots(int(c['frequency_hz']),c['boundary'],int(c['CHIGH']));j=int(r['mode'])-1;v=.02*np.sin(mu[j]*mod[0]);w=mod[1][:,j];fac=np.vdot(v,w)/np.vdot(v,v);err=abs(w-fac*v);check('modal_depth_function',r['stem']+'_'+r['mode'],max(num(err.max(),r['modal_function_max_abs_error']),num(error(v,w),r['modal_function_shape_diagnostic'])),1e-8)
 for r in rows(O/'ANALYTIC_MESH_COMPARISON.csv'):
  if 'boundary' not in r:continue
  key=(r['boundary'],int(r['frequency_hz']),int(r['CHIGH']));v=mainfields[(*key,160001,r['probe_resource'])];w=mainfields[(*key,80001,r['probe_resource'])];check('analytic_mesh',str(r),num(error(v,w),r['shape_difference']),1e-8)
 for r in rows(O/'E_STD_BOTTOM_EXTREMES.csv'):
  if 'frequency_hz' not in r:continue
  f=int(r['frequency_hz']);c=next(c for c in calls if c['kind']=='E_STD_V' and int(c['frequency_hz'])==f);mod=modal(Path(c['cwd'])/(c['stem']+'.mod'));v=mpressure(mod,sents[r['sentinel_id']],a[r['resource']]);w=mpressure(modal(old('K2',f)),sents[r['sentinel_id']],a[r['resource']]);check('E_STD_bottom',str(r),num(error(w,v),r['spatial_shape_difference']),1e-8)
 check('attempted_calls',str(len(calls)),abs(len(calls)-seal['calls_attempted']))
 planned=rows(O/'ALL_PLANNED_CALL_STATUS.csv');check('full_call_denominator','19',abs(len(planned)-19))
 failed=[r for r in calls if r.get('analytic_gate_PASS')=='False' or r['status']!='SUCCESS']
 if failed:check('stop_immediately_after_failure',failed[0]['stem'],0 if int(failed[0]['call_order'])==len(calls) else 1)
 for r in calls:
  check('solver_raw_MOD_hash',r['stem'],0 if r.get('MOD_sha256')==hashfile(Path(r['cwd'])/(r['stem']+'.mod')) else 1)
 fail=sum(not r['PASS'] for r in checks)
 with (O/'COLD_AUDIT_CHECKS.csv').open('w',encoding='utf-8',newline='') as f:
  w=csv.DictWriter(f,fieldnames=list(checks[0]),lineterminator='\n');w.writeheader();w.writerows(checks)
 result={'PASS':fail==0,'checks':len(checks),'FAIL':fail,'independent_root_method':'brentq cos/sin; independent struct MOD reader; J0-iY0; QR projection','preview_records':336,'old_modal_classification_records':len(rows(O/'SSP_AND_TURNING_CLASSIFICATION.csv')),'old_band_field_vectors':len(bandfields),'analytic_configs_cold_rebuilt':nconfig,'analytic_field_vectors_cold_rebuilt':len(crossrows),'E_STD_bottom_vectors_cold_rebuilt':len(rows(O/'E_STD_BOTTOM_EXTREMES.csv')),'max_independent_analytic_shape_discrepancy':max_field_delta,'scope':'REBUILD_AND_REPORTED_GATE_CORRECTNESS; cold PASS does not turn a scientific FAIL into PASS'}
 (O/'COLD_AUDIT.json').write_text(json.dumps(result,indent=2)+'\n',encoding='utf-8');print(json.dumps(result),flush=True)
 if fail:raise RuntimeError('COLD_AUDIT_FAIL')
if __name__=='__main__':main()