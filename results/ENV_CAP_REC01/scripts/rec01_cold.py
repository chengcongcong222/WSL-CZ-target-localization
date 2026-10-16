"""Independent REC01 cold audit: struct, fsum, Decimal and quadrature; no main import."""
from pathlib import Path
import csv,json,math,struct,hashlib,time
from decimal import Decimal,localcontext
import numpy as np
from scipy.integrate import quad
O=Path(__file__).resolve().parents[1];OLD=O.parent/'ENV_CAP_FULL_CAMPAIGN';W=Path('D:/ProjectStorage/WSL-CZ/ENV_CAP_REC01');RAW=W.parent/'ENV_CAP/E1'
def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def main():
 if (W/'COLD_STARTED.json').exists():raise RuntimeError('cold audit cannot be repeated')
 cfg=json.loads((O/'DESIGN_FREEZE.json').read_text(encoding='utf-8'))
 for a in cfg['execution_inputs']:assert sha(a['path'])==a['sha256'],a['path']
 (W/'COLD_STARTED.json').write_text(json.dumps({'time':time.time(),'code_sha':sha(__file__),'freeze_sha':sha(O/'DESIGN_FREEZE.json')})+'\n',encoding='utf-8')
 checks=[];summ=[]
 def check(name,ok,detail=None):checks.append({'check':name,'PASS':bool(ok),'detail':detail})
 old=json.loads((OLD/'E0_FREEZE.json').read_text(encoding='utf-8'));main=json.loads((O/'REC01_CZ_PHYSICS_ADJUDICATION.json').read_text(encoding='utf-8'))
 saved=list(csv.DictReader((O/'REC01_FULL_RANGE_ENERGY.csv').open(encoding='utf-8-sig')));aggregate=list(csv.DictReader((O/'REC01_E1_NUMERICAL_RECONSTRUCTION.csv').open(encoding='utf-8-sig')));savedrays=list(csv.DictReader((O/'REC01_RAY_RECONSTRUCTION.csv').open(encoding='utf-8-sig')))
 for case in old['E1_cases']:
  name=case['case_id'];zs=case['source_z_m'];H=case['H_m'];D=RAW/name;decision=next(a for a in main['cases'] if a['case_id']==name);series=[]
  prof=list(csv.DictReader((OLD/case['profile']).open(encoding='utf-8-sig')));zz=[float(a['depth_m']) for a in prof];cc=[float(a['sound_speed_mps']) for a in prof]
  paths=list(csv.DictReader((O/'inputs'/(name+'_PATHS.csv')).open(encoding='utf-8-sig')))
  if not (D/'I0.shd').exists():
   check(name+'_negative_no_fabricated_field',case['env_id']=='NO_CZ_CONSTANT' and decision['classification']=='CZ_NOT_PRESENT_IN_TESTED_DOMAIN' and all(a['deep_bottom_bounce']=='True' for a in paths));continue
  for lev in ['I0','I1','I2']:
   b=(D/(lev+'.shd')).read_bytes();stride=4*struct.unpack_from('<i',b)[0];nf,nt,nx,ny,ns,nz,nr=struct.unpack_from('<7i',b,2*stride)
   check(name+lev+'_layout', (nf,nt,nx,ny,ns,nz,nr)==(1,1,1,1,1,81,881) and len(b)>=(10+nz)*stride)
   depths=struct.unpack_from('<'+str(nz)+'f',b,8*stride);r=np.array(struct.unpack_from('<'+str(nr)+'f',b,9*stride))/1000;ix=[j for j,v in enumerate(depths) if abs(v-zs)<=25];rows=[]
   for j in ix:
    row=struct.unpack_from('<'+str(nr*2)+'f',b,(10+j)*stride);rows.append(list(zip(row[::2],row[1::2])))
   new=[];bounds=[];maxratio=0.;exact=0
   for i in range(nr):
    indices=[max(0,min(nr-1,j)) for j in range(i-4,i+5)];terms=[re*re+im*im for row in rows for j in indices for re,im in [row[j]]];en=math.fsum(terms)/(9*len(rows));new.append(en);bound=16*len(terms)*np.finfo(float).eps*en+np.nextafter(0.,1.);bounds.append(bound)
   output=[a for a in saved if a['case_id']==name and a['level']==lev];check(name+lev+'_all_range_and_rows',len(output)==nr and all(float(a['range_km'])==r[i] for i,a in enumerate(output)))
   diffs=np.array([abs(float(a['positive_window_energy'])-new[i]) for i,a in enumerate(output)]);check(name+lev+'_linear_energy_fsum',np.all(diffs<=bounds),{'max_absolute':float(diffs.max()),'max_bound_ratio':float(np.max(diffs/np.array(bounds)))})
   check(name+lev+'_exact_zero_preserved',all((float(a['positive_window_energy'])==0)==(new[i]==0) for i,a in enumerate(output)))
   # Fixed endpoints/midpoint, plus every exact-zero bin: Decimal products and sums independent of grouping.
   ids=sorted(set([0,nr//2,nr-1]+[i for i,v in enumerate(new) if v==0]))
   with localcontext() as ctx:
    ctx.prec=80
    for i in ids:
     js=[max(0,min(nr-1,j)) for j in range(i-4,i+5)];ds=sum((Decimal.from_float(re)**2+Decimal.from_float(im)**2 for row in rows for j in js for re,im in [row[j]]),Decimal(0))/Decimal(9*len(rows));err=abs(float(ds)-new[i]);check(name+lev+f'_decimal_bin{i}',err<=bounds[i],err)
   arr=np.array(new);t=-10*np.log10(np.maximum(arr,1e-30));series.append(t)
   check(name+lev+'_TL_and_dynamic_labels',all(abs(float(a['positive_TL_floor1e30_dB'])-t[i])<1e-10 and all((a[f'within_assumed_{d}dB_range']=='True')==(arr[i]>=arr.max()*10**(-d/10)) for d in cfg['relative_power_dynamic_range_dB']) for i,a in enumerate(output)))
   agg=next(a for a in aggregate if a['case_id']==name and a['level']==lev)
   check(name+lev+'_summary_zero_and_legacy',int(agg['zero_positive_windows'])==int(np.sum(arr==0)) and abs(float(agg['max_legacy_TL_delta_dB'])-max(abs(float(a['old_saved_TL_dB'])-t[i]) for i,a in enumerate(output)))<1e-9)
   summ.append({'case_id':name,'level':lev,'bins':nr,'max_fsum_absolute_difference':float(diffs.max()),'exact_zero_bins':int(np.sum(arr==0))})
  # Independent endpoint-regularized numerical quadrature vs rational C-linear segment integrals.
  for target in [1.,3.,5.]:
   possible=[a for a in paths if a['deep_bottom_bounce']=='False']
   if not possible:continue
   hit=min(possible,key=lambda a:abs(float(a['angle_deg'])-target));ang=float(hit['angle_deg']);zs_turn=float(hit['turn_z_m']);p=math.cos(math.radians(ang))/float(np.interp(zs,zz,cc));cp=1/p
   j=next(j for j in range(len(zz)-1) if zz[j]<=zs_turn<=zz[j+1] and cc[j+1]>cc[j]);g=(cc[j+1]-cc[j])/(zz[j+1]-zz[j]);top=math.sqrt(zs_turn-zs)
   def fun(u):
    if u<1e-4:return math.sqrt(2*cp/g)
    c=float(np.interp(zs_turn-u*u,zz,cc));return 2*u*p*c/math.sqrt(max(1-(p*c)**2,1e-30))
   cuts=sorted(math.sqrt(zs_turn-z) for z in zz if zs<z<zs_turn);q=2*quad(fun,0,top,points=cuts,epsabs=1e-5,epsrel=1e-8,limit=500)[0]/1000
   check(name+f'_Snell_independent_quadrature_{target}',abs(q-float(hit['deep_return_km']))<=.001,abs(q-float(hit['deep_return_km'])))
  valid=[a for a in paths if a['deep_bottom_bounce']=='False'];aa=[float(a['angle_deg']) for a in valid];rs=[float(a['deep_return_km']) for a in valid];rayg=[]
  for lev in ['R0','R1']:
   lines=(D/(lev+'.ray')).read_text().splitlines();i=7;errs=[];hits=[]
   while i<len(lines):
    if not lines[i].strip():i+=1;continue
    a=float(lines[i]);n,ns,nb=map(int,lines[i+1].split());i+=2;x=[[float(v) for v in line.split()] for line in lines[i:i+n]];i+=n
    for j in range(n-1):
     if x[j][1]>zs and x[j+1][1]<=zs:
      zprefix=[p[1] for p in x[:j+2]]
      if min(zprefix)>=1 and max(zprefix)<H-1:
       f=(zs-x[j][1])/(x[j+1][1]-x[j][1]);rr=(x[j][0]+f*(x[j+1][0]-x[j][0]))/1000
       if aa and aa[0]<=a<=aa[-1] and 10<=rr<=120:errs.append(abs(rr-float(np.interp(a,aa,rs))));hits.append(rr)
      break
   recorded=[a for a in savedrays if a['case_id']==name and a['level']==lev];check(name+lev+'_ray_prefix_replay',len(recorded)==len(errs) and all(abs(float(a['Snell_delta_km'])-v)<1e-9 for a,v in zip(recorded,errs)))
   rayg.append(len(errs)>=10 and max(errs)<=.1)
  fold=decision['fold']
  # Cold discrete local minima from raw analytic path table, no call to scipy peak finder.
  minima=[i for i in range(1,len(valid)-1) if rs[i]<rs[i-1] and rs[i]<rs[i+1] and 10<=rs[i]<=120]
  check(name+'_fold_identity',bool(minima)==bool(fold) and (not minima or abs(min(rs[i] for i in minima)-fold['deep_return_km'])<1e-8))
  if fold:
   expected=fold['deep_return_km'];mask=abs(r-expected)<=5;focus=[]
   for t in series:
    choices=[i for i in range(len(r)) if mask[i]];ix=min(choices,key=lambda i:t[i]);prom=float(np.median([t[i] for i in range(len(r)) if 5<=abs(r[i]-r[ix])<=10])-t[ix]);half=t[ix]+max(prom,0)/2;l=u=ix
    while l>0 and t[l]<=half:l-=1
    while u<len(r)-1 and t[u]<=half:u+=1
    focus.append((float(r[ix]),float(r[u]-r[l]),prom,ix in (choices[0],choices[-1])))
   met=[]
   for i,j in [(0,1),(1,2)]:
    delta=abs(series[i][mask]-series[j][mask]);med=float(np.median(delta));p95=float(np.quantile(delta,.95));okay=med<=1 and p95<=2 and abs(focus[i][0]-focus[j][0])<=1 and abs(focus[i][1]-focus[j][1])<=2;met.append((okay,p95))
   rf=decision['refined_angle_fold'];geometry=bool(rf and abs(rf['deep_return_km']-expected)<=.1 and all(rayg));supports=geometry and all(x[0] for x in met) and not focus[-1][3] and focus[-1][2]>max(.5,2*max(x[1] for x in met))
   check(name+'_conditional_CZ_classification',geometry==decision['geometry_supported'] and supports==decision['ray_TL_conditional_supported'])
   check(name+'_focus_and_refinement_metrics',all(abs(f[0]-g['peak_km'])<1e-9 and abs(f[1]-g['width_km'])<1e-9 and abs(f[2]-g['prominence_dB'])<1e-9 for f,g in zip(focus,decision['focus_metrics'])) and all(abs(m[1]-v['P95_dB'])<1e-9 for m,v in zip(met,decision['field_refinement'])))
 mirror=json.loads((O/'REC01_MIRROR_CONTROL.json').read_text(encoding='utf-8'))
 with localcontext() as ctx:
  ctx.prec=80
  for k,case in enumerate(mirror['cases']):
   for j,(xyz,dd) in enumerate(zip(case['coordinates_m'],case['distances_m'])):
    s=case['source_m'];dist=sum(((Decimal.from_float(float(s[i]))-Decimal.from_float(float(xyz[i])))**2 for i in range(3)),Decimal(0)).sqrt();ref=float(dist);check(f'mirror_{k}_channel{j}_decimal',abs(ref-dd)<=8*math.ulp(ref),abs(ref-dd))
   check(f'mirror_{k}_symmetry',case['max_distance_delta_m']==0 and case['PASS'])
 check('all_units_and_full_bins',len(main['cases'])==15 and len(saved)==12*3*881)
 check('old_audit_12_failures_retained',json.loads((OLD/'E1_INDEPENDENT_COLD_AUDIT.json').read_text(encoding='utf-8'))['FAIL_count']==12)
 check('new_MC_and_solver_zero',json.loads((O/'REC01_ANALYSIS_EXIT.json').read_text(encoding='utf-8'))['new_solver_calls']==0)
 check('mirror_result_and_guard',mirror['PASS'] and all(mirror['invalid_input_guards']))
 result={'PASS':all(a['PASS'] for a in checks),'checks':len(checks),'FAIL_count':sum(not a['PASS'] for a in checks),'scope':'All 31716 stored-field bins reconstructed independently; Decimal sentinels/zeros; ray prefix and path quadrature; conditional classifications; mirror Decimal distances. No absolute real-ocean or coherent modal accuracy certificate.','checks_detail':checks,'per_curve':summ}
 (O/'REC01_INDEPENDENT_COLD_AUDIT.json').write_text(json.dumps(result,ensure_ascii=False,indent=2,allow_nan=False)+'\n',encoding='utf-8',newline='\n');print(json.dumps({k:v for k,v in result.items() if k not in ('checks_detail','per_curve')},ensure_ascii=False));return 0 if result['PASS'] else 2
if __name__=='__main__':raise SystemExit(main())
