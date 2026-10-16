"""Independent CZROB02 audit: own Hermite slopes/root bisection, adaptive quadrature, raw-ray parser. No physics/main import."""
from pathlib import Path
import csv,json,math,hashlib,time,warnings
import numpy as np
from scipy.integrate import quad,IntegrationWarning
O=Path(__file__).resolve().parents[1];R=O.parents[1];OLD=R/'results/ENV_CAP_FULL_CAMPAIGN';REC=R/'results/ENV_CAP_REC01';W=Path('D:/ProjectStorage/WSL-CZ/ENV_CAP_CZROB02')
def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
class IndependentProfile:
 def __init__(self,z,c,kind):
  self.z=np.array(z,float);self.c=np.array(c,float);self.kind=kind;h=np.diff(self.z);d=np.diff(self.c)/h
  if kind=='PCHIP':
   m=np.zeros(len(z))
   for i in range(1,len(z)-1):
    if d[i-1]*d[i]>0:
     w1=2*h[i]+h[i-1];w2=h[i]+2*h[i-1];m[i]=(w1+w2)/(w1/d[i-1]+w2/d[i])
   def end(h0,h1,d0,d1):
    v=((2*h0+h1)*d0-h0*d1)/(h0+h1)
    if v*d0<=0:return 0.
    if d0*d1<0 and abs(v)>3*abs(d0):return 3*d0
    return v
   m[0]=end(h[0],h[1],d[0],d[1]);m[-1]=end(h[-1],h[-2],d[-1],d[-2]);self.coeff=np.c_[(m[:-1]+m[1:]-2*d)/h**2,(3*d-2*m[:-1]-m[1:])/h,m[:-1],self.c[:-1]]
  else:self.coeff=np.c_[np.zeros(len(h)),np.zeros(len(h)),d,self.c[:-1]]
 def ix(self,x):return min(len(self.z)-2,max(0,int(np.searchsorted(self.z,x,side='right')-1)))
 def val(self,x):
  j=self.ix(x);a,b,c,d=self.coeff[j];u=x-self.z[j];return ((a*u+b)*u+c)*u+d
 def derivatives(self,x,j):
  a,b,c,d=self.coeff[j];u=x-self.z[j];return (3*a*u*u+2*b*u+c,6*a*u+2*b,6*a)
 def turn(self,zs,cp,down):
  candidates=[j for j in range(len(self.z)-1) if (self.z[j+1]>zs and self.c[j]<cp<=self.c[j+1]) or (not down and self.z[j]<zs and self.c[j]>=cp>self.c[j+1])]
  if down:candidates=[j for j in candidates if self.c[j+1]>self.c[j] and self.z[j+1]>zs]
  else:candidates=[j for j in candidates if self.c[j]>self.c[j+1] and self.z[j]<zs]
  if not candidates:return (float(self.z[-1]) if down else 0.),True
  j=candidates[0] if down else candidates[-1];a=max(zs,self.z[j]) if down else self.z[j];b=self.z[j+1] if down else min(zs,self.z[j+1])
  for k in range(64):
   mid=(a+b)/2;f=self.val(mid)-cp
   if (down and f<0) or (not down and f>0):a=mid
   else:b=mid
  return float((a+b)/2),False
 def leg(self,a,b,cp,where):
  nodes=[a]+[float(v) for v in self.z if a<v<b]+[b];total=0.;maxerr=0.
  for k,(aa,bb) in enumerate(zip(nodes[:-1],nodes[1:])):
   endpoint=(where=='HIGH' and k==len(nodes)-2) or (where=='LOW' and k==0)
   if endpoint:
    zt=bb if where=='HIGH' else aa;sgn=-1 if where=='HIGH' else 1;j=self.ix((aa+bb)/2);d1,d2,d3=self.derivatives(zt,j);g=-sgn*d1
    if g<=1e-12:raise ArithmeticError('critical turning endpoint')
    def fun(u):
     cv=self.val(zt+sgn*u*u);factor=g-.5*d2*u*u-sgn*d3/6*u**4;return 2*cv/math.sqrt(factor*(cp+cv))
    val,err=quad(fun,0,math.sqrt(bb-aa),epsabs=1e-5,epsrel=1e-9,limit=150)
   else:
    def fun(x):
     cv=self.val(x);return cv/math.sqrt((cp-cv)*(cp+cv))
    val,err=quad(fun,aa,bb,epsabs=1e-5,epsrel=1e-9,limit=150)
   total+=val;maxerr=max(maxerr,err)
  return total,maxerr
 def calc(self,zs,a):
  cp=self.val(zs)/math.cos(math.radians(a));down,bd=self.turn(zs,cp,True);up,sf=self.turn(zs,cp,False);d,ed=self.leg(zs,down,cp,None if bd else 'HIGH');u,eu=self.leg(up,zs,cp,None if sf else 'LOW');return 2*d/1000,2*u/1000,down,bd,sf,max(ed,eu)/500

def main():
 if (W/'COLD_STARTED.json').exists():raise RuntimeError('No audit retry')
 cfg=json.loads((O/'DESIGN_FREEZE.json').read_text(encoding='utf-8'))
 for a in cfg['execution_inputs']:assert sha(a['path'])==a['sha256'],a['path']
 (W/'COLD_STARTED.json').write_text(json.dumps({'time_unix':time.time(),'freeze_sha':sha(O/'DESIGN_FREEZE.json'),'code_sha':sha(__file__)})+'\n',encoding='utf-8')
 allchecks=[];comparisons=[]
 coeffs=list(csv.DictReader((O/'CZROB02_INTERPOLATION_COEFFICIENTS.csv').open(encoding='utf-8-sig')))
 def ck(n,ok,detail=None):allchecks.append({'check':n,'PASS':bool(ok),'detail':detail})
 rows=list(csv.DictReader((O/'CZROB02_RETURN_CURVES.csv').open(encoding='utf-8-sig')));res=json.loads((O/'CZROB02_PRIMARY_RESULTS.json').read_text(encoding='utf-8'));old=json.loads((OLD/'E0_FREEZE.json').read_text(encoding='utf-8'));ode=list(csv.DictReader((O/'CZROB02_ODE_COMPARISON.csv').open(encoding='utf-8-sig')));sheets=list(csv.DictReader((O/'CZROB02_SSP_NODE_DIAGNOSTICS.csv').open(encoding='utf-8-sig')))
 for case in res['cases']:
  name=case['case_id'];print('AUDIT',name,flush=True);conf=next(a for a in old['E1_cases'] if a['case_id']==name);data=list(csv.DictReader((OLD/conf['profile']).open(encoding='utf-8-sig')));z=[float(a['depth_m']) for a in data];c=[float(a['sound_speed_mps']) for a in data];zs=conf['source_z_m'];H=conf['H_m']
  for method in ['C_LINEAR','PCHIP']:
   p=IndependentProfile(z,c,method)
   for jj in range(len(z)-1):
    co=next(a for a in coeffs if a['case_id']==name and a['SSP_method']==method and int(a['interval_index'])==jj);ck(name+method+f'_independent_Hermite_coeff{jj}',max(abs(float(co[k])-p.coeff[jj,j]) for j,k in enumerate(['cubic','quadratic','linear','constant']))<=1e-10)
   curve=[a for a in rows if a['case_id']==name and a['SSP_method']==method];ck(name+method+'_all_angle_rows',len(curve)==cfg['angle_samples'] and all(int(a['angle_index'])==i for i,a in enumerate(curve)))
   for j in range(len(z)):
    sheet=next(a for a in sheets if a['case_id']==name and a['SSP_method']==method and int(a['node_index'])==j);ck(name+method+f'_node{j}_original_samples',abs(p.val(z[j])-c[j])<=1e-10 and abs(float(sheet['interpolated_c_mps'])-p.val(z[j]))<=1e-10)
   # All C-linear angles and registered PCHIP sentinel stride; every discovered fold added by fixed rule.
   ids=set(range(len(curve))) if method=='C_LINEAR' else set(range(0,len(curve),cfg['cold_PCHIP_stride']))|{len(curve)-1}
   for i in range(1,len(curve)-1):
    if all(curve[j]['status']=='FINITE' and curve[j]['deep_bottom_bounce']=='False' and curve[j]['deep_refracted_clearance_PASS']=='True' for j in [i-1,i,i+1]) and float(curve[i]['deep_return_km'])<min(float(curve[i-1]['deep_return_km']),float(curve[i+1]['deep_return_km'])):ids.update([i-1,i,i+1])
   # Exact C-linear primitive gives independent, fast full-curve check.
   if method=='C_LINEAR':
    for i in sorted(ids):
     a=curve[i]
     if a['status']!='FINITE':continue
     angle=float(a['angle_deg']);cp=p.val(zs)/math.cos(math.radians(angle));t,bb=p.turn(zs,cp,True);zz=np.r_[zs,np.array(z)[(np.array(z)>zs)&(np.array(z)<t)],t];cc=np.interp(zz,z,c);u=cc/cp;ss=np.sqrt(np.maximum(0,1-u*u));slope=np.diff(cc)/np.diff(zz)
     terms=[(ss[k]-ss[k+1])*cp/slope[k] if abs(slope[k])>1e-12 else (zz[k+1]-zz[k])*u[k]/ss[k] for k in range(len(slope))];rr=2*math.fsum(terms)/1000
     ck(name+f'_C_LINEAR_angle{i}',abs(rr-float(a['deep_return_km']))<=cfg['cold_quadrature_tolerance_km'] and bb==(a['deep_bottom_bounce']=='True'),abs(rr-float(a['deep_return_km'])))
   else:
    for i in sorted(ids):
     a=curve[i]
     if a['status']!='FINITE':continue
     try:
      with warnings.catch_warnings(record=True) as warns:
       warnings.simplefilter('always',IntegrationWarning);v,u,t,bb,sf,err=p.calc(zs,float(a['angle_deg']))
      delta=max(abs(v-float(a['deep_return_km'])),abs(u-float(a['upper_return_km'])));good=delta<=cfg['cold_quadrature_tolerance_km'] and err<=cfg['cold_quadrature_tolerance_km'] and bb==(a['deep_bottom_bounce']=='True') and sf==(a['upper_surface_bounce']=='True');comparisons.append({'case_id':name,'method':method,'angle_index':i,'angle_deg':a['angle_deg'],'cold_deep_km':v,'cold_upper_km':u,'max_delta_km':delta,'quad_error_estimate_km':err,'quad_warnings':len(warns)});ck(name+f'_PCHIP_adaptive_angle{i}',good,delta)
     except Exception as e:ck(name+f'_PCHIP_adaptive_angle{i}',False,repr(e))
   def folds(stride):
    a=curve[::stride];f=[]
    for i in range(1,len(a)-1):
     if all(t['status']=='FINITE' and t['deep_bottom_bounce']=='False' and t['deep_refracted_clearance_PASS']=='True' for t in a[i-1:i+2]):
      vv=float(a[i]['deep_return_km'])
      if vv<float(a[i-1]['deep_return_km']) and vv<float(a[i+1]['deep_return_km']) and 10<=vv<=120:f.append(a[i])
    return f
   ff=folds(1);cf=folds(2);first=min(ff,key=lambda a:float(a['deep_return_km'])) if ff else None;coarse=min(cf,key=lambda a:float(a['deep_return_km'])) if cf else None;mr=case['methods'][method]
   ck(name+method+'_first_fold_replay',len(ff)==mr['folds_count'] and bool(first)==bool(mr['first_fold']) and (first is None or abs(float(first['deep_return_km'])-mr['first_fold']['deep_return_km'])<1e-10))
   finite=[a for a in curve if a['status']=='FINITE'];qok=len(finite)==len(curve) and max(float(a['max_order_difference_km']) for a in finite)<=cfg['quadrature_tolerance_km'];aok=bool(first and coarse and abs(float(first['deep_return_km'])-float(coarse['deep_return_km']))<=cfg['angle_fold_tolerance_km']);ck(name+method+'_order_and_angle_Gates',qok==mr['quadrature_PASS'] and aok==mr['angle_refinement_PASS'])
   orows=[a for a in ode if a['case_id']==name and a['SSP_method']==method]
   for k,a in enumerate(orows):
    buf=np.load(a['path_file']);rx=buf['range_m'];zz=buf['depth_m'];theta=buf['angle_rad'];iv=np.cos(theta)/np.array([p.val(float(v)) for v in zz]);p0=math.cos(math.radians(float(a['angle_deg'])))/p.val(zs);inv=float(np.max(abs(iv/p0-1)));ck(name+method+f'_ODE_path{k}_invariant',abs(inv-float(a['max_relative_Snell_invariant_error']))<=1e-10 and len(rx)==int(a['n_points']) and abs(zz[0]-zs)<1e-10)
    v,u,t,bb,sf,qe=p.calc(zs,float(a['angle_deg']));rd=float(a['range_km']) if a['range_km'] else None;dif=abs(rd-v) if rd is not None else None;passed=a['success']=='True' and a['returned']=='True' and dif is not None and dif<=cfg['ODE_Snell_tolerance_km'] and inv<=cfg['ODE_invariant_relative_tolerance'];ck(name+method+f'_ODE{k}_independent_Snell_and_gate',passed==(a['PASS']=='True') and (dif is None or abs(dif-float(a['Snell_difference_km']))<=cfg['cold_quadrature_tolerance_km']))
   ck(name+method+'_ODE_count',len(orows)==mr['ODE_cases'])
   odepass=None if not orows else all(a['PASS']=='True' and a['step_refinement_PASS']=='True' for a in orows)
   for a in orows:
    other=next(t for t in orows if t['angle_deg']==a['angle_deg'] and t['step_max_m']!=a['step_max_m']);delta=abs(float(a['range_km'])-float(other['range_km'])) if a['range_km'] and other['range_km'] else None;ck(name+method+'_ODE_step_'+a['angle_deg']+'_'+a['step_max_m'],(delta is not None and delta<=cfg['ODE_step_tolerance_km'])==(a['step_refinement_PASS']=='True'))
   ck(name+method+'_combined_numeric_gate',mr['ODE_PASS']==odepass and mr['numeric_PASS']==bool(qok and aok and odepass is not False))
  # Independent archive parser: all source-depth crossings, prefix contacts, focus membership.
  pathsummary=list(csv.DictReader((O/'CZROB02_PATH_ATTRIBUTION.csv').open(encoding='utf-8-sig')))
  for lev in ['R0','R1']:
   sheet=next(a for a in pathsummary if a['case_id']==name and a['data']=='ARCHIVED_C_LINEAR_'+lev);lo=float(sheet['focus_low_km']);hi=float(sheet['focus_high_km']);lines=(W.parent/'ENV_CAP/E1'/name/(lev+'.ray')).read_text().splitlines();i=7;total=zero=boundary=surface=bottom=0
   while i<len(lines):
    if not lines[i].strip():i+=1;continue
    ang=float(lines[i]);n,ns,nb=map(int,lines[i+1].split());i+=2;xy=[[float(v) for v in s.split()] for s in lines[i:i+n]];i+=n
    for j in range(n-1):
     a,b=xy[j],xy[j+1]
     if not ((a[1]>zs and b[1]<=zs) or (a[1]<zs and b[1]>=zs)):continue
     r=(a[0]+(zs-a[1])/(b[1]-a[1])*(b[0]-a[0]))/1000
     if 10<=r<=120:
      total+=1
      if lo<=r<=hi:
       hasS=any(v[1]<=1 for v in xy[:j+2]);hasB=any(v[1]>=H-1 for v in xy[:j+2]);zero+=not(hasS or hasB);boundary+=(hasS or hasB);surface+=hasS;bottom+=hasB
   ck(name+lev+'_boundary_prefix_summary',total==int(sheet['all_10_120km_crossings']) and zero==int(sheet['zero_boundary_focus_crossings']) and boundary==int(sheet['boundary_focus_crossings']) and surface==int(sheet['surface_focus_crossings']) and bottom==int(sheet['bottom_focus_crossings']))
  cl=case['methods']['C_LINEAR'];pc=case['methods']['PCHIP'];numeric=cl['numeric_PASS'] and pc['numeric_PASS'];both=cl['first_fold'] and pc['first_fold'];stable=bool(numeric and both and abs(cl['first_fold']['deep_return_km']-pc['first_fold']['deep_return_km'])<=cfg['representation_range_tolerance_km'] and abs(cl['first_fold']['deep_turn_z_m']-pc['first_fold']['deep_turn_z_m'])<=cfg['representation_turn_tolerance_m']);classification='GEOMETRY_NUMERIC_UNCERTIFIED' if not numeric else ('SSP_GEOMETRY_ROBUST_FOR_TESTED_REPRESENTATIONS' if stable else 'SSP_INTERPOLATION_SENSITIVE');ck(name+'_classification',stable==case['geometry_stable'] and classification==case['classification']);ck(name+'_no_energy_claim',case['energy_attribution']=='EXCLUSIVE_ENERGY_ATTRIBUTION_UNRESOLVED' and case['coherent_provider_investment']=='DEFER')
 ck('all9_and_one_family',len(res['cases'])==9 and {a['source_family'] for a in res['cases']}=={'WOA23_17N118E'})
 ck('no_external_calls_MC_or_new_fields',res['new_external_solver_calls']==0 and res['new_MC']==0 and res['PCHIP_energy']=='NOT_EVALUATED')
 with (O/'CZROB02_COLD_QUADRATURE.csv').open('w',encoding='utf-8',newline='') as f:
  d=csv.DictWriter(f,list(dict.fromkeys(k for a in comparisons for k in a)));d.writeheader();d.writerows(comparisons)
 out={'PASS':all(a['PASS'] for a in allchecks),'checks':len(allchecks),'FAIL_count':sum(not a['PASS'] for a in allchecks),'scope':'Own Hermite slopes/root bisection; all C-linear return ranges; registered PCHIP adaptive-quadrature sentinels and all folds; ODE saved-path invariants and reference returns; all raw-ray prefix summaries. No independent energy decomposition.','checks_detail':allchecks};(O/'CZROB02_INDEPENDENT_COLD_AUDIT.json').write_text(json.dumps(out,ensure_ascii=False,indent=2)+'\n',encoding='utf-8',newline='\n');print(json.dumps({k:v for k,v in out.items() if k!='checks_detail'}));return 0 if out['PASS'] else 2
if __name__=='__main__':raise SystemExit(main())
