"""Independent read-only E1 auditor: no import of env_provider/e1_campaign."""
from pathlib import Path
import numpy as np,csv,json,struct,hashlib,math
from scipy.integrate import quad
R=Path(__file__).resolve().parents[3];O=R/'results/ENV_CAP_FULL_CAMPAIGN';W=Path('D:/ProjectStorage/WSL-CZ/ENV_CAP')
def main():
 out=O/'E1_INDEPENDENT_COLD_AUDIT.json'
 if out.exists():raise RuntimeError('Cold audit already executed')
 rows=list(csv.DictReader((O/'E1_CZ_ENVIRONMENTS.csv').open()));checks=[];diffs=[]
 def test(name,ok,detail=''):
  checks.append({'check':name,'PASS':bool(ok),'detail':str(detail)})
 frozen=json.loads((O/'E0_FREEZE.json').read_text());cfg={c['case_id']:c for c in frozen['E1_cases']};summary=list(csv.DictReader((O/'E1_NUMERICAL_STABILITY.csv').open()));rayrows=list(csv.DictReader((O/'E1_RAY_INDEPENDENT_COMPARISON.csv').open()))
 for row in rows:
  c=cfg[row['case_id']];D=W/'E1'/row['case_id'];prof=list(csv.DictReader((O/c['profile']).open()));z=np.array([float(a['depth_m']) for a in prof]);speeds=np.array([float(a['sound_speed_mps']) for a in prof]);zs=c['source_z_m'];H=c['H_m']
  for x in c['input_sha']:test(row['case_id']+'_input_'+Path(x['path']).stem,hashlib.sha256(Path(x['path']).read_bytes()).hexdigest()==x['sha256'])
  ar=list(csv.DictReader((D/'SNELL_RETURN_CURVE.csv').open()))
  # Independent numerical quadrature after z=turn-u^2 removes endpoint singularity.
  for target_a in [1.,3.,5.]:
   if not ar:continue
   hit=min(ar,key=lambda a:abs(float(a['angle_deg'])-target_a));angle=float(hit['angle_deg']);p=math.cos(math.radians(angle))/float(np.interp(zs,z,speeds));cp=1/p
   candidates=[i for i in range(len(z)-1) if z[i+1]>zs and speeds[i]<cp<=speeds[i+1] and speeds[i+1]>speeds[i]]
   if not candidates:continue
   i=candidates[0];g=(speeds[i+1]-speeds[i])/(z[i+1]-z[i]);turn=z[i]+(cp-speeds[i])/g;top=math.sqrt(turn-zs)
   cuts=sorted(math.sqrt(turn-zz) for zz in z if zs<zz<turn)
   def fun(u):
    if u<1e-4:return math.sqrt(2*cp/g)
    cc=float(np.interp(turn-u*u,z,speeds));den=max(1-(p*cc)**2,1e-30);return 2*u*p*cc/math.sqrt(den)
   got=2*quad(fun,0,top,points=cuts,epsabs=1e-5,epsrel=1e-8,limit=500)[0]/1000;error=abs(got-float(hit['first_return_km']));test(row['case_id']+f'_Snell_quadrature_a{target_a}',error<=.001,error)
  if row['execution_status']!='COMPLETED':test(row['case_id']+'_failed_preserved',row['CZ_class']=='CZ_NUMERIC_UNCERTIFIED');continue
  if row['CZ_class']=='CZ_NOT_PRESENT_IN_TESTED_DOMAIN':
   test(row['case_id']+'_no_return_domain',not any(10<=float(a['first_return_km'])<=120 for a in ar));continue
  arrays=[];ff=[];expected=float(row['ray_first_caustic_km'])
  for level in ['I0','I1','I2']:
   buf=(D/(level+'.shd')).read_bytes();record=4*struct.unpack_from('<i',buf)[0];nf,nt,nx,ny,ns,nrd,nrr=struct.unpack_from('<7i',buf,2*record)
   depths=np.array(struct.unpack_from('<'+str(nrd)+'f',buf,8*record));ranges=np.array(struct.unpack_from('<'+str(nrr)+'f',buf,9*record))/1000
   pressure=[]
   for d in range(nrd):
    x=np.array(struct.unpack_from('<'+str(nrr*2)+'f',buf,(10+d)*record));pressure.append(x[::2]+1j*x[1::2])
   pressure=np.array(pressure);energy=np.mean(abs(pressure[abs(depths-zs)<=25])**2,axis=0);pad=np.pad(energy,(4,4),mode='edge');energy=np.array([sum(pad[i:i+9])/9 for i in range(nrr)]);tl=-10*np.log10(np.maximum(energy,1e-30));saved=np.load(D/(level+'_ENERGY_CURVE.npz'));err=float(np.max(abs(tl-saved['TL_dB'])));test(row['case_id']+'_'+level+'_TL_from_raw_SHD',err<=1e-9,err);test(row['case_id']+'_'+level+'_range_depth_units',np.allclose(ranges,saved['r_km']) and ranges[0]==10 and ranges[-1]==120 and nrd==81)
   ids=np.where(abs(ranges-expected)<=5)[0];ii=min(ids,key=lambda j:tl[j]);prom=float(np.median(tl[(abs(ranges-ranges[ii])>=5)&(abs(ranges-ranges[ii])<=10)])-tl[ii]);half=tl[ii]+max(prom,0)/2;lo=hi=ii
   while lo>0 and tl[lo]<=half:lo-=1
   while hi<nrr-1 and tl[hi]<=half:hi+=1
   ff.append({'peak':ranges[ii],'width':ranges[hi]-ranges[lo],'prom':prom,'boundary':ii==ids[0] or ii==ids[-1]});arrays.append(tl)
  comparisons=[]
  for i,j in [(0,1),(1,2)]:
   mask=abs(ranges-expected)<=5;d=sorted(abs(arrays[i][mask]-arrays[j][mask]));med=np.quantile(d,.5);p95=np.quantile(d,.95);pk=abs(ff[i]['peak']-ff[j]['peak']);wi=abs(ff[i]['width']-ff[j]['width']);okay=med<=1 and p95<=2 and pk<=1 and wi<=2;comparisons.append((okay,p95))
   saved=next(a for a in summary if a['case_id']==row['case_id'] and a['comparison']==f'I{i}_I{j}');test(row['case_id']+f'_comparison_I{i}_I{j}',okay==(saved['PASS']=='True') and abs(p95-float(saved['P95_TL_delta_dB']))<1e-9)
  rays=[]
  for lev in ['R0','R1']:
   lines=(D/(lev+'.ray')).read_text().splitlines();ii=7;errors=[];aa=[float(v['angle_deg']) for v in ar];rs=[float(v['first_return_km']) for v in ar]
   while ii<len(lines):
    if not lines[ii].strip():ii+=1;continue
    angle=float(lines[ii]);n,ns,nb=map(int,lines[ii+1].split());ii+=2;pts=np.array([[float(v) for v in l.split()] for l in lines[ii:ii+n]]);ii+=n
    zz=pts[:,-1];hit=np.where((zz[:-1]>zs)&(zz[1:]<=zs))[0]
    if len(hit):
     jj=hit[0]
     if zz[:jj+2].min()>=1 and zz[:jj+2].max()<H-1:
      frac=(zs-zz[jj])/(zz[jj+1]-zz[jj]);r=(pts[jj,0]+frac*(pts[jj+1,0]-pts[jj,0]))/1000
      if aa[0]<=angle<=aa[-1] and 10<=r<=120:errors.append(abs(r-np.interp(angle,aa,rs)))
   saved=next(a for a in rayrows if a['case_id']==row['case_id'] and a['level']==lev)
   test(row['case_id']+'_'+lev+'_raw_ray_replay',len(errors)==int(saved['n_compared']) and abs(max(errors)-float(saved['max_return_range_difference_km']))<1e-8)
   rays.append({'n_compared':len(errors),'max_return_range_difference_km':max(errors) if errors else math.inf})
  raygood=all(a['n_compared']>=10 and a['max_return_range_difference_km']<=.1 for a in rays)
  passes=row['fold_interior']=='True' and raygood and all(a[0] for a in comparisons) and ff[-1]['prom']>max(.5,2*max(a[1] for a in comparisons)) and not ff[-1]['boundary']
  test(row['case_id']+'_classification',passes==(row['CZ_class']=='CZ2_RAY_TL_STABLE'))
  test(row['case_id']+'_fine_peak',abs(ff[-1]['peak']-float(row['peak_km']))<1e-9)
 # Outputs and deterministic no-retry record per attempted unit.
 test('all_15_units_present',len(rows)==15 and len({a['case_id'] for a in rows})==15)
 for c in frozen['E1_cases']:test(c['case_id']+'_unit_marker_and_exit',(W/'E1'/c['case_id']/'SCIENCE_STARTED.json').exists() and (W/'E1'/c['case_id']/'UNIT_EXIT.json').exists())
 result={'PASS':all(a['PASS'] for a in checks),'checks':len(checks),'FAIL_count':sum(not a['PASS'] for a in checks),'scope':'Every attempted E1 raw SHD reconstructed with struct parser and explicit energy smoothing; analytic integrals independently quadrature-checked; input identities/classifications recomputed. No new solver calls or science rerun. WOA reconstruction identity covered by source audit, not external CTD or coherent modal certification.','rows':checks}
 out.write_text(json.dumps(result,indent=2)+'\n',encoding='utf-8');print(json.dumps({k:v for k,v in result.items() if k!='rows'}));return 0 if result['PASS'] else 2
if __name__=='__main__':raise SystemExit(main())
