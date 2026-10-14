#!/usr/bin/env python3
"""T2-N3-PREC1: frozen ideal-waveguide MOD precision attribution.

No solvers; never modifies historical files. Requires only numpy/scipy and two
committed PR_R 201 Hz MOD files, geometry CSV and N3 recorded field comparisons.
C0: stored MOD k, stored MOD phi (historical result)
C1: analytic FP64 k, stored MOD phi
C2: stored MOD k, analytic FP64 phi
C3: analytic FP64 k, analytic FP64 phi (far-field reference)

The analytic ideal-waveguide k cannot be used as a repair for E-STD variable-SSP fields.
"""
from __future__ import annotations
import argparse,csv,datetime,hashlib,json,math,struct,time
from pathlib import Path
import numpy as np
from scipy.special import hankel2

SCOPE='IDEAL_ISOVELOCITY_PRESSURE_RELEASE_SURFACE_RIGID_BOTTOM_ONLY'
REPO_BASE='results/R4_T2_N3_BOUNDARY_PHYSICS_ANALYTIC_REFERENCE'
NAMES=('A0_0','A_P4M','B1_V4','B2_DUP')
PROBES=(('P01',45.5,0.0,180.0),('P02',50.0,0.0,200.0),('P03',55.0,2.5,200.0),('P04',58.0,4.5,210.0),('P05',59.5,5.0,220.0))
H=5000.0; C=1500.0; F=201.0
WINDOWS=(1800,2400)
ORIGINAL_GATE=.002
PARENT='8ee2bd53fe1b7a50725833e6c7b5c5dea35fbe23'
EXPECTED_SOURCE_SHA256={
 f'{REPO_BASE}/modal_archive/n3_PR_R_f201_h1800_n160001/n3_PR_R_f201_h1800_n160001.mod':'cf5745862a3b367283c9ea9deef723ff194ca354c7f389fa5b185b992f028d0f',
 f'{REPO_BASE}/modal_archive/n3_PR_R_f201_h2400_n160001/n3_PR_R_f201_h2400_n160001.mod':'d9e476d0ed84c08c30f7551ad090e3823067d36918ed41375eeb3ec720afae42',
 f'{REPO_BASE}/ANALYTIC_KRAKEN_MODE_AND_FIELD_CROSSCHECK.csv':'2de1541cdfcc93c5f3dc76310fb3c6b8b47bf2e382518d8534d5ede2f56860bd',
 f'{REPO_BASE}/package/T2_TILT_G0_ARRAY_GEOMETRY_20261013_v1.csv':'bc2ccaebfa52893ba7bc7f217791ec980c47cb37d743abcd20cae055cf732ffd',
}

def sha(p):
 h=hashlib.sha256()
 with Path(p).open('rb') as fp:
  for b in iter(lambda:fp.read(2**20),b''): h.update(b)
 return h.hexdigest()

def rows(p):
 with Path(p).open(newline='',encoding='utf-8-sig') as f:return list(csv.DictReader(f))

def save_csv(p,data,columns=None):
 data=list(data);cols=columns or list(data[0])
 with Path(p).open('w',newline='',encoding='utf-8') as f:
  w=csv.DictWriter(f,fieldnames=cols,lineterminator='\n');w.writeheader();w.writerows(data)

def save_json(p,data):
 Path(p).write_text(json.dumps(data,indent=2,ensure_ascii=False,allow_nan=False)+'\n',encoding='utf-8')

def make_geometry(repopath):
 tab=rows(repopath/REPO_BASE/'package/T2_TILT_G0_ARRAY_GEOMETRY_20261013_v1.csv')
 data={n:np.array([[float(r[k]) for k in ('x_m','y_m','z_m')] for r in sorted([r for r in tab if r['resource']==n],key=lambda q:int(q['channel']))]) for n in NAMES}
 assert tuple(len(data[n]) for n in NAMES)==(8,8,12,12)
 assert all(np.max(abs(g[:,1]))==0 for g in data.values())
 return data

def mod_read(path):
 b=Path(path).read_bytes()
 if len(b)<1024:raise RuntimeError('INVALID_MOD_LENGTH')
 rv=4*struct.unpack_from('<i',b,0)[0]
 nf,nmedia,nt,nz=struct.unpack_from('<4i',b,84)
 if rv<8*nt or nf!=1 or nt!=nz:raise RuntimeError('INVALID_MOD_HEADER')
 freq=struct.unpack_from('<d',b,3*rv)[0]
 depths=np.frombuffer(b,dtype='<f4',count=nt,offset=4*rv).astype(float)
 m=struct.unpack_from('<i',b,5*rv)[0]
 if m<=0 or 7+m>=len(b)//rv:raise RuntimeError('INVALID_MOD_COUNT')
 phi=np.stack([np.frombuffer(b,dtype='<c8',count=nz,offset=(7+j)*rv).astype(np.complex128) for j in range(m)],axis=1)
 k=np.frombuffer(b,dtype='<c8',count=m,offset=(7+m)*rv).astype(np.complex128)
 if len(k)!=m or not np.isfinite(phi).all() or not np.isfinite(k).all():raise RuntimeError('INVALID_MODAL_ARRAY')
 if abs(freq-F)>1e-8 or np.any(k.real<=0) or np.any(abs(k.imag)>1e-12):raise RuntimeError('FREQUENCY_OR_WAVENUMBER_CONTRACT_FAIL')
 return {'path':str(path),'depths':depths,'phi':phi,'k':k,'mode_count':m,'sha256':sha(path),'freq':freq,'record_bytes':rv}

def depth_indices(depths,req):
 out=[]
 for z in req:
  i=int(np.argmin(abs(depths-z)))
  if abs(depths[i]-z)>1e-5:raise RuntimeError(f'MOD_EXACT_DEPTH_REQUIRED: {z}')
  out.append(i)
 return out

def analytic_spectrum(m):
 mu=(np.arange(m,dtype=np.float64)+.5)*np.pi/H
 return mu,np.sqrt((2*np.pi*F/C)**2-mu**2)

def field(mod,position,coord,use_analytic_k=False,use_analytic_phi=False):
 """Complex 64-bit sums; modal radial factor is historical large-range kernel."""
 pid,r_km,theta,zs=position
 z=mod['depths'];m=mod['mode_count'];phi=mod['phi']
 mu,ana_k=analytic_spectrum(m)
 k=ana_k if use_analytic_k else mod['k'].real
 omega=2*np.pi*F
 if np.any(~np.isfinite(k)) or np.any(k<=0):raise RuntimeError('K_INVALID')
 if use_analytic_phi:
  ps=np.sqrt(2/H)*np.sin(mu*zs)
  pr=np.sqrt(2/H)*np.sin(mu[:,None]*coord[None,:,2])
 else:
  ix=depth_indices(z,[zs])[0];rec=depth_indices(z,coord[:,2]);ps=phi[ix];pr=phi[rec].T
 weights=ps[:,None]*pr
 theta=math.radians(theta);xy=r_km*1000*np.array([math.cos(theta),math.sin(theta)])
 rho=np.linalg.norm(xy[None,:]-coord[:,:2],axis=1)
 phase=-1j*k[:,None]*rho[None,:]-1j*np.pi/4
 ker=np.sqrt(2*np.pi/(k[:,None]*rho[None,:]))*np.exp(phase)
 return np.sum(weights*ker,axis=0,dtype=np.complex128)

def stable_shape(u,v):
 if not np.isfinite(u).all() or not np.isfinite(v).all() or np.linalg.norm(u)==0 or np.linalg.norm(v)==0: raise RuntimeError('FIELD_UNDEFINED')
 a=u/np.linalg.norm(u); b=v/np.linalg.norm(v)
 return float(np.linalg.norm(a-b*np.vdot(b,a)))

def run(repo,output,verify_frozen=True):
 repo=Path(repo).resolve();output=Path(output).resolve();output.mkdir(parents=True,exist_ok=True)
 start=time.perf_counter()
 if (output/'PREC1_EXECUTION_STARTED.json').exists():raise RuntimeError('PREC1_ONCE_ONLY_EXECUTION_ALREADY_STARTED')
 geo=make_geometry(repo)
 history=rows(repo/REPO_BASE/'ANALYTIC_KRAKEN_MODE_AND_FIELD_CROSSCHECK.csv')
 modpaths={high:repo/REPO_BASE/f'modal_archive/n3_PR_R_f201_h{high}_n160001/n3_PR_R_f201_h{high}_n160001.mod' for high in WINDOWS}
 input_paths=list(modpaths.values())+[repo/REPO_BASE/'ANALYTIC_KRAKEN_MODE_AND_FIELD_CROSSCHECK.csv',repo/REPO_BASE/'package/T2_TILT_G0_ARRAY_GEOMETRY_20261013_v1.csv']
 for p in input_paths:
  if not p.is_file():raise FileNotFoundError(p)
 if verify_frozen:
  for source,expected in EXPECTED_SOURCE_SHA256.items():
   candidate=repo/source
   if sha(candidate)!=expected:raise RuntimeError('FROZEN_SOURCE_SHA256_MISMATCH '+source)
 save_json(output/'PREC1_EXECUTION_STARTED.json',{'stage':'T2_N3_PREC1','utc':datetime.datetime.now(datetime.timezone.utc).isoformat(),'parent':PARENT,'inputs_checked':verify_frozen,'solver_calls':0})
 mds={high:mod_read(modpaths[high]) for high in WINDOWS}
 for high,expected_modes in ((1800,741),(2400,1046)):
  if mds[high]['mode_count']!=expected_modes:raise RuntimeError('FROZEN_MODE_COUNT_CHANGED '+str(high))
 samples=[];kmodes=[];exceed=[]
 for high in WINDOWS:
  mod=mds[high];count=mod['mode_count'];mu,k64=analytic_spectrum(count);k32=np.float32(k64).astype(np.float64)
  dk=mod['k'].real-k64;dr=mod['k'].real-k32
  if np.max(abs(dr))!=0.:raise RuntimeError('MOD_K_NOT_EXACT_FP32_ROUNDED_ANALYTIC: '+str(high))
  if np.max(abs(mod['k'].imag))!=0.:raise RuntimeError('MODE_IMAGINARY_NONZERO')
  kr_stats=dict(CHIGH=high,mode_count=count,all_stored_k_equal_float32_analytic=True,max_abs_k_difference=float(np.max(abs(dk))),max_relative_k_difference=float(np.max(abs(dk)/k64)),max_phase_error_58km_rad=float(np.max(abs(dk))*58000),max_analytic_cp_mps=float(np.max(2*np.pi*F/k64)),nonzero_im_k_count=int(np.sum(abs(mod['k'].imag)>0)))
  kmodes.append(kr_stats)
  for p in PROBES:
   for res in NAMES:
    xyz=geo[res]
    fields={(ka,pa):field(mod,p,xyz,ka,pa) for ka,pa in [(False,False),(True,False),(False,True),(True,True)]}
    C0=fields[(False,False)];C1=fields[(True,False)];C2=fields[(False,True)];C3=fields[(True,True)]
    original=next((q for q in history if q['stem']==f'n3_PR_R_f201_h{high}_n160001' and q['probe']==p[0] and q['resource']==res),None)
    if not original:raise RuntimeError('MISSING_HISTORY_ROW')
    e0=stable_shape(C0,C3);e1=stable_shape(C1,C3);e2=stable_shape(C2,C3)
    diff=abs(e0-float(original['spatial_shape_error']))
    record=dict(CHIGH=high,probe=p[0],resource=res,r_km=p[1],abs_theta_deg=p[2],source_z_m=p[3],channels=len(xyz),C0_stored_k_stored_phi_vs_C3=e0,C1_analytic_k_stored_phi_vs_C3=e1,C2_stored_k_analytic_phi_vs_C3=e2,C3_analytic_k_analytic_phi_vs_C3=0.0,C0_to_C1_direct=stable_shape(C0,C1),C0_to_C2_direct=stable_shape(C0,C2),C1_to_C2_direct=stable_shape(C1,C2),history_analytic_shape=original['spatial_shape_error'],history_reconstruction_error=diff,history_gate_pass=original['PASS'],C0_below_0p002=e0<=ORIGINAL_GATE,C1_below_0p002=e1<=ORIGINAL_GATE,C2_below_0p002=e2<=ORIGINAL_GATE,real_environment_generalization='NOT_EVALUATED')
    if diff>2e-6:raise RuntimeError('N3_RECORDED_FIELD_NOT_REPRODUCED '+str((high,p[0],res,diff)))
    samples.append(record)
  print('READ_ONLY_ATTRIBUTION',high,'modes',count,'rows_so_far',len(samples),flush=True)
 assert len(samples)==40
 save_csv(output/'PREC1_ALL_40_CASES.csv',samples)
 save_csv(output/'PREC1_MODE_WAVENUMBER_STORAGE.csv',kmodes)
 primary=[r for r in samples if r['CHIGH']==2400 and r['probe']=='P04']
 save_csv(output/'PREC1_P04_FROZEN_FOUR_RESOURCE.csv',primary)
 for s in primary:exceed.append(dict(resource=s['resource'],original_exceeds=float(s['C0_stored_k_stored_phi_vs_C3'])>.002,k_only_replacement_clears=float(s['C1_analytic_k_stored_phi_vs_C3'])<=.002,phi_only_replacement_clears=float(s['C2_stored_k_analytic_phi_vs_C3'])<=.002))
 save_csv(output/'PREC1_GATE_ATTRIBUTION.csv',exceed)
 decision={
 'stage':'T2_N3_PREC1','scope':SCOPE,'status':'PREC1_K_STORAGE_ROUNDING_MECHANISM_SUPPORTED_CONDITIONALLY' if all(x['original_exceeds'] and x['k_only_replacement_clears'] for x in exceed) else 'PREC1_ATTRIBUTION_INCOMPLETE',
 'original_N3_scientific_status':'T2_N3_PARTIAL_NUMERIC_NOT_ADMITTED_UNCHANGED',
 'all_40_control_rows':len(samples),'K_float32_exact_modes':sum(q['mode_count'] for q in kmodes),'P04_k_replacement_cleared_count':sum(x['k_only_replacement_clears'] for x in exceed),'P04_phi_replacement_cleared_count':sum(x['phi_only_replacement_clears'] for x in exceed),
 'all_40_k_replacement_below_original_threshold':all(s['C1_below_0p002'] for s in samples),
 'max_history_delta':max(s['history_reconstruction_error'] for s in samples),
 'source_hashes':{str(p.relative_to(repo)):sha(p) for p in input_paths},
 'binary_solver_calls':0,'new_KRAKEN':0,'new_FIELD':0,'new_MC_waveforms':0,'R4_percent':0,
 'E_STD_CHIGH_real_seabed_and_A_B_performance':'NOT_EVALUATED',
 'wall_s':round(time.perf_counter()-start,6),
 'stop':'STOP_ONCE_ONLY',
 'research_lead_audit':'PENDING',
 }
 save_json(output/'PREC1_DECISION.json',decision)
 print(json.dumps({k:v for k,v in decision.items() if k not in ('source_hashes',)},ensure_ascii=False,indent=2))
 return decision


def selftest():
 a=np.array([1+1j,2-1j,3+2j]);assert stable_shape(a,5j*a)<1e-15
 z=np.array([.2,.3,.4]);mu,k=analytic_spectrum(6); assert np.all(k>0)
 from tempfile import TemporaryDirectory
 with TemporaryDirectory() as td:
  p=Path(td)/'synthetic.mod'; rec=128; M=5; N=4
  b=bytearray((7+M+1)*rec)
  struct.pack_into('<i',b,0,rec//4);struct.pack_into('<4i',b,84,1,1,N,N)
  struct.pack_into('<d',b,3*rec,F);struct.pack_into('<4f',b,4*rec,180,190,200,210);struct.pack_into('<i',b,5*rec,M)
  for j in range(M):
   for t,d in enumerate((180,190,200,210)):
    val=np.sqrt(2/H)*np.sin((j+.5)*np.pi*d/H)
    struct.pack_into('<2f',b,(7+j)*rec+8*t,val,0.)
   struct.pack_into('<2f',b,(7+M)*rec+8*j,np.float32(k[j]),0.)
  p.write_bytes(b); mod=mod_read(p); assert mod['mode_count']==5
  assert np.all(mod['k'].real==np.float32(k[:M]).astype(np.float64))
  x=np.array([[-1,0,200],[1,0,210]])
  v=field(mod,('XX',55,2.5,180),x,True,True)
  assert v.shape==(2,) and np.isfinite(v).all()
 print('SELFTEST_PASS')

if __name__=='__main__':
 ap=argparse.ArgumentParser();ap.add_argument('--repo',type=Path);ap.add_argument('--out',type=Path);ap.add_argument('--selftest',action='store_true');args=ap.parse_args()
 if args.selftest:selftest()
 else:
  if not args.repo or not args.out:ap.error('--repo and --out are required')
  run(args.repo,args.out)
