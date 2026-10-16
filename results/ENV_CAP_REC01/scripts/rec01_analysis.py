"""REC01 deterministic archived-data analysis; never imports old scientific code."""
from pathlib import Path
import csv,json,hashlib,struct,math,time,sys
import numpy as np
from scipy.ndimage import uniform_filter1d
from scipy.signal import find_peaks
O=Path(__file__).resolve().parents[1];R=O.parents[1];OLD=R/'results/ENV_CAP_FULL_CAMPAIGN';W=Path('D:/ProjectStorage/WSL-CZ/ENV_CAP_REC01');RAW=W.parent/'ENV_CAP/E1'
EPS=np.finfo(float).eps

def sha(p):
 h=hashlib.sha256()
 with Path(p).open('rb') as f:
  for b in iter(lambda:f.read(1048576),b''):h.update(b)
 return h.hexdigest()
def js(p,x):
 Path(p).write_text(json.dumps(x,ensure_ascii=False,indent=2,allow_nan=False,default=lambda v:v.item() if hasattr(v,'item') else str(v))+'\n',encoding='utf-8',newline='\n')
def table(p,rows):
 keys=list(dict.fromkeys(k for a in rows for k in a))
 with Path(p).open('w',encoding='utf-8',newline='') as f:
  d=csv.DictWriter(f,keys);d.writeheader();d.writerows(rows)
def readshd(p):
 b=p.read_bytes();rec=4*struct.unpack_from('<i',b)[0];dims=struct.unpack_from('<7i',b,2*rec);nf,nt,nx,ny,ns,nz,nr=dims
 assert dims[:5]==(1,1,1,1,1) and nz==81 and nr==881 and len(b)>=(10+nz)*rec
 z=np.frombuffer(b,'<f4',nz,8*rec).astype(float);r=np.frombuffer(b,'<f4',nr,9*rec).astype(float)/1000
 p=np.array([np.frombuffer(b,'<c8',nr,(10+j)*rec).astype(np.complex128) for j in range(nz)])
 assert np.isfinite(p).all() and np.array_equal(r,np.arange(10,120.001,.125)) and np.array_equal(z,np.arange(100,300.01,2.5))
 return r,z,p

def positive_energy(p,z,zs):
 q=p[np.abs(z-zs)<=25];e=np.mean(q.real*q.real+q.imag*q.imag,axis=0,dtype=np.float64)
 ep=np.pad(e,4,mode='edge');out=np.array([np.sum(ep[i:i+9],dtype=np.float64)/9 for i in range(len(e))]);assert (out>=0).all() and np.isfinite(out).all()
 return e,out

def tl(e,floor=1e-30):return -10*np.log10(np.maximum(e,floor))

def segment(z,c,p):
 u=p*c;assert np.max(u)<=1+64*EPS
 s=np.sqrt(np.maximum(0,(1-u)*(1+u)));den=s[:-1]+s[1:];assert (den>0).all()
 return float(np.sum(p*(c[:-1]+c[1:])*np.diff(z)/den))

def path(z,c,zs,H,a):
 cs=float(np.interp(zs,z,c));p=math.cos(math.radians(a))/cs;cp=1/p
 deep=np.flatnonzero((z[1:]>zs)&(c[:-1]<cp)&(c[1:]>=cp)&(np.diff(c)>0))
 if len(deep):
  j=deep[0];bottom=False;turn=z[j]+(cp-c[j])/(c[j+1]-c[j])*(z[j+1]-z[j])
  if turn>=H-1:bottom=True;turn=H
 else:bottom=True;turn=H
 zz=np.r_[zs,z[(z>zs)&(z<turn)],turn];cc=np.interp(zz,z,c)
 if not bottom:cc[-1]=cp
 if np.max(p*cc)>1+64*EPS:return None
 down=2*segment(zz,cc,p)/1000
 upper=np.flatnonzero((z[:-1]<zs)&(c[:-1]>=cp)&(c[1:]<cp)&(np.diff(c)<0))
 if len(upper):
  j=upper[-1];zt=z[j]+(cp-c[j])/(c[j+1]-c[j])*(z[j+1]-z[j]);surface=False
 else:zt=0.;surface=True
 zz=np.r_[zt,z[(z>zt)&(z<zs)],zs];cc=np.interp(zz,z,c)
 if not surface:cc[0]=cp
 up=2*segment(zz,cc,p)/1000
 return {'angle_deg':float(a),'deep_return_km':down,'turn_z_m':float(turn),'phase_speed_mps':cp,'deep_bottom_bounce':bottom,'upper_return_km':up,'upper_surface_bounce':surface,'upper_then_deep_km':up+down,'turn_nearest_SSP_node_m':float(np.min(np.abs(z-turn)))}

def curve(z,c,zs,H,n):return [b for a in np.linspace(.1,20,n) if (b:=path(z,c,zs,H,a)) is not None]

def choose_fold(rows):
 zero=[a for a in rows if not a['deep_bottom_bounce']];y=np.array([a['deep_return_km'] for a in zero]);ix=find_peaks(-y)[0]
 ii=[i for i in ix if 10<=y[i]<=120]
 if not ii:return None
 return zero[min(ii,key=lambda i:y[i])]

def rayreturns(p,zs,H,rows):
 lines=p.read_text().splitlines();assert 'rz' in lines[6].lower();i=7;errs=[];out=[]
 good=[a for a in rows if not a['deep_bottom_bounce']];aa=np.array([a['angle_deg'] for a in good]);rr=np.array([a['deep_return_km'] for a in good])
 while i<len(lines):
  if not lines[i].strip():i+=1;continue
  angle=float(lines[i]);n,ns,nb=map(int,lines[i+1].split());i+=2;x=np.array([[float(v) for v in s.split()] for s in lines[i:i+n]]);i+=n
  hit=np.flatnonzero((x[:-1,1]>zs)&(x[1:,1]<=zs))
  if not len(hit):continue
  j=int(hit[0]);seg=x[:j+2,1]
  if seg.min()<1 or seg.max()>=H-1:continue
  f=(zs-x[j,1])/(x[j+1,1]-x[j,1]);r=(x[j,0]+f*(x[j+1,0]-x[j,0]))/1000
  if len(aa) and aa[0]<=angle<=aa[-1] and 10<=r<=120:
   err=abs(r-float(np.interp(angle,aa,rr)));errs.append(err);out.append({'angle_deg':angle,'range_km':r,'Snell_delta_km':err,'trace_total_surface_bounces_not_prefix':ns,'trace_total_bottom_bounces_not_prefix':nb})
 return {'n':len(errs),'max_error_km':max(errs) if errs else None,'PASS':bool(len(errs)>=10 and max(errs)<=.1)},out

def focus(r,t,expected):
 ids=np.flatnonzero(abs(r-expected)<=5);ii=int(ids[np.argmin(t[ids])]);rp=float(r[ii]);flank=(abs(r-rp)>=5)&(abs(r-rp)<=10);prom=float(np.median(t[flank])-t[ii]);half=t[ii]+max(prom,0)/2;lo=hi=ii
 while lo>0 and t[lo]<=half:lo-=1
 while hi<len(r)-1 and t[hi]<=half:hi+=1
 return {'peak_km':rp,'width_km':float(r[hi]-r[lo]),'prominence_dB':prom,'boundary':bool(ii in (ids[0],ids[-1]))}

def mirror_geometry(src,array):
 src=np.asarray(src,dtype=np.float64);array=np.asarray(array,dtype=np.float64)
 if src.shape!=(3,) or array.ndim!=2 or array.shape[1]!=3 or not np.isfinite(src).all() or not np.isfinite(array).all():raise ValueError('finite xyz metres required')
 if np.max(np.abs(src))>1e7 or np.max(np.abs(array))>1e7:raise ValueError('registered coordinate bound exceeded')
 d=src-array;return np.hypot(np.hypot(d[:,0],d[:,1]),d[:,2])

def arrays():
 s=np.arange(8,dtype=float)*2-7;zero=np.zeros(8)
 h=np.c_[s,zero,np.full(8,200.)];a=math.asin(4/14)
 return {'HLA8':h,'TILT_POS8':np.c_[s*math.cos(a),zero,200+s*math.sin(a)],'TILT_NEG8':np.c_[s*math.cos(a),zero,200-s*math.sin(a)],'HLA8_V4':np.vstack((h,[[-1,0,198],[-1,0,199],[-1,0,201],[-1,0,202]]))}

def main():
 cfg=json.loads((O/'DESIGN_FREEZE.json').read_text(encoding='utf-8'));old=json.loads((OLD/'E0_FREEZE.json').read_text(encoding='utf-8'))
 if (W/'FORMAL_STARTED.json').exists():raise RuntimeError('formal science cannot be repeated')
 for a in cfg['execution_inputs']:assert sha(a['path'])==a['sha256'],a['path']
 js(W/'FORMAL_STARTED.json',{'time_unix':time.time(),'freeze_sha256':sha(O/'DESIGN_FREEZE.json'),'inputs_verified':len(cfg['execution_inputs']),'new_solver_calls':0,'new_MC':0})
 full=[];summ=[];dec=[];rayrows=[];factors=[]
 for case in old['E1_cases']:
  name=case['case_id'];print('START',name,flush=True);zs=case['source_z_m'];H=case['H_m'];D=RAW/name
  prof=list(csv.DictReader((OLD/case['profile']).open(encoding='utf-8-sig')));z=np.array([float(a['depth_m']) for a in prof]);c=np.array([float(a['sound_speed_mps']) for a in prof]);z=np.r_[z[z<H],H];c=np.interp(z,[float(a['depth_m']) for a in prof],[float(a['sound_speed_mps']) for a in prof])
  paths=curve(z,c,zs,H,2001);finepaths=curve(z,c,zs,H,4001);fold=choose_fold(paths);foldfine=choose_fold(finepaths)
  table(O/'inputs'/(name+'_PATHS.csv'),paths)
  d={'case_id':name,'source_family':case['source_family'],'nominal_C_linear_only':True,'old_stage_unchanged':True,'fold':fold,'refined_angle_fold':foldfine,'classification':None,'geometry_supported':False,'ray_TL_conditional_supported':False,'boundary_contributions_isolated':False}
  if not (D/'I0.shd').exists():
   d.update(classification='CZ_NOT_PRESENT_IN_TESTED_DOMAIN',reason='Constant-speed negative: no refracted turning return; old fields were never generated. Other boundary returns not claimed absent.');dec.append(d);summ.append({'case_id':name,'level':'ALL','status':'NO_ARCHIVED_SHD_NEGATIVE_CONTROL','full_range_bins':0});continue
  allcurves=[];levsumm=[]
  for lev in ['I0','I1','I2']:
   r,rz,p=readshd(D/(lev+'.shd'));raw,e=positive_energy(p,rz,zs);newtl=tl(e)
   legacy_raw=np.mean(np.abs(p[np.abs(rz-zs)<=25])**2,axis=0);legacy=uniform_filter1d(legacy_raw,9,mode='nearest');saved=np.load(D/(lev+'_ENERGY_CURVE.npz'))['TL_dB'];assert np.max(abs(tl(legacy)-saved))<=1e-9
   maxe=float(e.max());floorcurves=[tl(e,x) for x in cfg['display_energy_floors']];gamma=(2*len(r)+9)*EPS/(1-(2*len(r)+9)*EPS);roll_bound=4*gamma*float(max(legacy_raw.max(),raw.max()))
   row={'case_id':name,'level':lev,'status':'RECONSTRUCTED_ALL_BINS','full_range_bins':len(r),'zero_positive_windows':int(np.sum(e==0)),'negative_legacy_windows':int(np.sum(legacy<0)),'max_legacy_linear_error':float(np.max(abs(legacy-e))),'rolling_global_absolute_roundoff_diagnostic_bound':roll_bound,'max_legacy_TL_delta_dB':float(np.max(abs(saved-newtl))),'max_energy':maxe,'median_energy':float(np.median(e))}
   row['rolling_error_within_global_diagnostic_bound']=bool(np.max(abs(legacy-e))<=roll_bound)
   for dyn in cfg['relative_power_dynamic_range_dB']:row[f'windows_below_{dyn}dB_assumed_dynamic_range']=int(np.sum(e<maxe*10**(-dyn/10)))
   for k,rr in enumerate(r):
    rec={'case_id':name,'level':lev,'range_km':float(rr),'raw_band_energy':float(raw[k]),'positive_window_energy':float(e[k]),'legacy_rolling_energy_signed':float(legacy[k]),'legacy_minus_positive_energy':float(legacy[k]-e[k]),'positive_TL_floor1e30_dB':float(newtl[k]),'old_saved_TL_dB':float(saved[k]),'exact_zero_window':bool(e[k]==0),'max_display_floor_TL_spread_dB':float(max(a[k] for a in floorcurves)-min(a[k] for a in floorcurves))}
    for dyn in cfg['relative_power_dynamic_range_dB']:rec[f'within_assumed_{dyn}dB_range']=bool(e[k]>=maxe*10**(-dyn/10))
    full.append(rec)
   levsumm.append(row);allcurves.append(newtl)
  summ.extend(levsumm)
  base=[a for a in paths if not a['deep_bottom_bounce']]
  rays=[]
  for lev in ['R0','R1']:
   met,rrs=rayreturns(D/(lev+'.ray'),zs,H,base);met['level']=lev;rays.append(met)
   rayrows.extend(dict(case_id=name,level=lev,**a) for a in rrs)
  d['raw_ray_checks']=rays
  if fold:
   expected=fold['deep_return_km'];mask=abs(r-expected)<=5;ff=[focus(r,t,expected) for t in allcurves];metrics=[]
   for i,j in [(0,1),(1,2)]:
    delta=abs(allcurves[i][mask]-allcurves[j][mask]);m={'comparison':f'I{i}_I{j}','median_dB':float(np.median(delta)),'P95_dB':float(np.quantile(delta,.95)),'peak_delta_km':abs(ff[i]['peak_km']-ff[j]['peak_km']),'width_delta_km':abs(ff[i]['width_km']-ff[j]['width_km'])};m['PASS']=bool(m['median_dB']<=1 and m['P95_dB']<=2 and m['peak_delta_km']<=1 and m['width_delta_km']<=2);metrics.append(m)
   stablefold=foldfine is not None and abs(expected-foldfine['deep_return_km'])<=.1
   geometry=stablefold and all(a['PASS'] for a in rays);enhance=not ff[-1]['boundary'] and ff[-1]['prominence_dB']>max(.5,2*max(a['P95_dB'] for a in metrics))
   support=geometry and enhance and all(a['PASS'] for a in metrics)
   d.update(geometry_supported=bool(geometry),ray_TL_conditional_supported=bool(support),focus_metrics=ff,field_refinement=metrics,fold_angle_refinement_PASS=bool(stablefold),enhancement_PASS=bool(enhance),classification='CZ_RAY_TL_CONDITIONAL_SUPPORTED' if support else ('CZ_RAY_GEOMETRY_SUPPORTED' if geometry else 'CZ_NUMERIC_UNCERTIFIED'))
   nearby=[]
   for a in paths:
    for kind,key,bounce in [('FIRST_BOTTOM','deep_return_km','deep_bottom_bounce'),('FIRST_SURFACE','upper_return_km','upper_surface_bounce'),('UPPER_THEN_DEEP_SURFACE','upper_then_deep_km','upper_surface_bounce')]:
     if a[bounce] and abs(a[key]-expected)<=5:nearby.append({'path_type':kind,'angle_deg':a['angle_deg'],'return_km':a[key]})
   d['possible_boundary_path_overlap_in_focus']=bool(nearby);d['boundary_path_overlap_samples']=len(nearby)
   factors.extend(dict(case_id=name,**a) for a in nearby)
   d['physical_scope']='Zero-boundary first-return fold plus total incoherent energy enhancement; boundary energy not separated. Conditional nominal C-linear evidence, not exclusive refractive intensity or detection.'
  else:
   d.update(classification='CZ_PHYSICAL_IDENTITY_UNRESOLVED',reason='No interior zero-boundary first-return fold in frozen downward-angle/range domain. Archived field retained; no adaptive focus selection.')
  dec.append(d);js(W/(name+'_EXIT.json'),d);print('EXIT',name,d['classification'],flush=True)
 table(O/'REC01_FULL_RANGE_ENERGY.csv',full);table(O/'REC01_E1_NUMERICAL_RECONSTRUCTION.csv',summ);table(O/'REC01_RAY_RECONSTRUCTION.csv',rayrows);table(O/'REC01_POSSIBLE_BOUNDARY_RETURNS.csv',factors)
 mirror=[]
 for name,ar in arrays().items():
  for s in [[45000,5000,185],[60000,5000,215],[1000000,200000,250],[0,1,200]]:
   x=np.array(s,dtype=np.int64);xm=x.copy();xm[1]*=-1;a=mirror_geometry(x,ar);b=mirror_geometry(xm,ar)
   mirror.append({'array':name,'source_m':s,'mirror_source_m':xm.tolist(),'max_distance_delta_m':float(np.max(abs(a-b))),'PASS':bool(np.array_equal(a,b)),'distances_m':a.tolist(),'dtype':str(a.dtype),'coordinates_m':ar.tolist()})
 guard=[]
 for s in [[float('nan'),0,0],[float('inf'),0,0],[1e10,0,0]]:
  try:mirror_geometry(s,arrays()['HLA8']);ok=False
  except ValueError:ok=True
  guard.append(ok)
 js(O/'REC01_MIRROR_CONTROL.json',{'PASS':all(a['PASS'] for a in mirror) and all(guard),'cases':mirror,'invalid_input_guards':guard,'units':'m; xyz; radians not used in distances','old_25_of_26_unchanged':True})
 conditional=[a['case_id'] for a in dec if a['ray_TL_conditional_supported']]
 js(O/'REC01_CZ_PHYSICS_ADJUDICATION.json',{'stage':'REC01_ONLY_PRE_COLD_PROVISIONAL','cases':dec,'conditional_cases':conditional,'distinct_CZ_source_families':sorted(set(a['source_family'] for a in dec if a['ray_TL_conditional_supported'])),'global_energy_reconstruction_pending_cold':True,'coherent_provider':'NOT_CERTIFIED','new_external_solver_calls':0,'new_MC':0,'old_E1_final_unchanged':True})
 js(O/'REC01_ANALYSIS_EXIT.json',{'execution_completed':True,'R1_cases':len(dec),'positive_energy_bins':len(full),'R2_PASS':all(a['PASS'] for a in mirror) and all(guard),'R3_new_numeric_execution':'NOT_REQUESTED_NO_CERTIFIED_COHERENT_PROVIDER','new_solver_calls':0,'new_MC':0})
if __name__=='__main__':main()
