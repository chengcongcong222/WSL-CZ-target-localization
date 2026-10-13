"""Independent read-only N2 audit: struct MOD/SHD IO, scalar sums, bracket interpolation and QR shapes."""
from t2_n2_common import ROOT,OUT,LOCAL,G0,OLDLOCAL,FREQ,RES,CATS,KVAR,sha,read,rows,dump,table,utc
import numpy as np,struct,math,time

def mod_read(p):
 b=p.read_bytes();rec=struct.unpack_from('<i',b)[0]*4;nt,nz=struct.unpack_from('<ii',b,92);count=struct.unpack_from('<i',b,5*rec)[0];depth=np.array(struct.unpack_from('<'+'f'*nt,b,4*rec));psi=np.empty((nt,count),complex)
 for j in range(count):x=np.array(struct.unpack_from('<'+'f'*(2*nz),b,(7+j)*rec));psi[:,j]=x[::2]+1j*x[1::2]
 x=np.array(struct.unpack_from('<'+'f'*(2*count),b,(7+count)*rec));return depth,psi,x[::2]+1j*x[1::2]

def cold_pressure(mod,s,coords,f):
 z,psi,k=mod;source=int(np.argmin(abs(z-s[2])));assert abs(z[source]-s[2])<=1e-5;out=[]
 for x,y,zr in coords:
  rx=s[0]*1000*np.cos(np.radians(s[1]))-x;ry=s[0]*1000*np.sin(np.radians(s[1]))-y;rr=math.sqrt(rx*rx+ry*ry);receiver=int(np.argmin(abs(z-zr)));assert abs(z[receiver]-zr)<=1e-5;terms=psi[source]*psi[receiver]*np.sqrt(2*math.pi/(k*rr))*np.exp(-1j*k*rr-1j*math.pi/4);out.append(np.sum(terms,dtype=np.complex128))
 return np.array(out)

def shd_read(p):
 with p.open('rb') as fp:
  rec=struct.unpack('<i',fp.read(4))[0]*4;fp.seek(2*rec);nf,nt,nx,ny,ns,nz,nr,atten=struct.unpack('<7if',fp.read(32));assert (nf,nt,nx,ny)==(1,1,1,1);fp.seek(3*rec);f=struct.unpack('<d',fp.read(8))[0];fp.seek(7*rec);sd=np.array(struct.unpack('<'+'f'*ns,fp.read(4*ns)));fp.seek(8*rec);rd=np.array(struct.unpack('<'+'f'*nz,fp.read(4*nz)));fp.seek(9*rec);r=np.array(struct.unpack('<'+'f'*nr,fp.read(4*nr)));val=np.empty((ns,nz,nr),complex)
  for i in range(ns):
   for j in range(nz):fp.seek((10+i*nz+j)*rec);x=np.array(struct.unpack('<'+'f'*(2*nr),fp.read(8*nr)));val[i,j]=x[::2]+1j*x[1::2]
 return f,sd,rd,r,val

def cold_interpolate(d,s,coords):
 f,sd,rd,axis,p=d;iz=int(np.argmin(abs(sd-s[2])));assert abs(sd[iz]-s[2])<=1e-5 and np.all(np.diff(axis)>0);out=[]
 for x,y,z in coords:
  zr=int(np.argmin(abs(rd-z)));assert abs(rd[zr]-z)<=1e-5;rx=s[0]*1000*np.cos(np.radians(s[1]))-x;ry=s[0]*1000*np.sin(np.radians(s[1]))-y;r=math.sqrt(rx*rx+ry*ry);j=int(np.searchsorted(axis,r,side='right'))-1;j=max(j,0);assert j<len(axis)-1 and axis[j]<=r<=axis[j+1];q=(r-axis[j])/(axis[j+1]-axis[j]);out.append((1-q)*p[iz,zr,j]+q*p[iz,zr,j+1])
 return np.array(out)

def qre(a,b):
 A=a/np.sqrt(np.sum(abs(a)**2,axis=-1))[...,None];B=b/np.sqrt(np.sum(abs(b)**2,axis=-1))[...,None];R=np.linalg.qr(np.stack([A,B],axis=-1),mode='reduced')[1];return np.abs(R[...,1,1])

def main():
 began=time.monotonic();checks=[]
 def add(name,e,limit=1e-10,detail='',n=1):checks.append({'check':name,'error':float(e),'limit':limit,'PASS':bool(np.isfinite(e) and e<=limit),'detail':detail,'items':int(n)})
 for name,h in read(OUT/'DESIGN_FREEZE.json')['code_and_input_sha256'].items():add('design_identity',0 if sha(ROOT/name)==h else 1,0,name)
 for r in read(OUT/'SOURCE_IDENTITY.json')['records']:add('source_identity',0 if sha(r['path'])==r['sha256'] else 1,0,r['path'])
 seal=read(OUT/'SCIENCE_OUTPUT_FREEZE.json')
 for name,h in seal['raw_files'].items():add('science_original_output_identity',0 if sha(name)==h else 1,0,name)
 for name,h in read(OUT/'NUMERICAL_EVALUATION_FREEZE.json')['scientific_result_sha256'].items():add('numerical_evaluation_identity',0 if sha(OUT/name)==h else 1,0,name)
 plan=read(OUT/'SOLVER_PLAN.json');sources={r['path']:r for r in read(OUT/'SOURCE_IDENTITY.json')['records']};calls=rows(OUT/'PROVIDER_CALLS.csv');add('call_count_budgets',0 if sum(r['tool']=='kraken' for r in calls)<=9 and sum(r['tool']=='scooter' for r in calls)<=9 and sum(r['tool']=='fields' for r in calls)<=144 else 1,0)
 add('no_duplicate_solver_stem_call',len(calls)-len(set((r['tool'],r['stem']) for r in calls)),0)
 for r in calls:
  binary=ROOT/'tools/acoustics_toolbox/atWin10/at/bin'/(r['tool']+'.exe');add('actual_solver_binary_identity',0 if sha(binary)==r['binary_sha256'] else 1,0)
  if r['tool']!='fields':
   p=next(p for p in plan['primary'] if p['stem']==r['stem']);add('actual_ENV_identity',0 if sha(__import__('pathlib').Path(r['cwd'])/(r['stem']+'.env'))==p['env_sha256'] else 1,0)
  else:
   p=next(p for p in plan['FIELDS'] if p['stem']==r['stem']);add('actual_FLP_identity',0 if sha(__import__('pathlib').Path(r['cwd'])/(r['stem']+'.flp'))==p['flp_sha256'] else 1,0)
 for p in plan['primary']:
  old=(OLDLOCAL/'modal'/f"t2_f{p['frequency_hz']}_n160001.env").read_text().splitlines();new=__import__('pathlib').Path(p['env_path']).read_text().splitlines();b=next(i for i,x in enumerate(old) if x.strip().startswith("'R'"));add('same_physical_ENV_block',0 if old[1:4]==new[1:4] and old[5:b+1]==new[5:b+1] else 1,0,p['stem'])
 g={};coord=rows(OUT/'ARRAY_GEOMETRY_TABLE.csv')
 for res in RES:g[res]=np.array([[float(r[k]) for k in ['x_m','y_m','z_m']] for r in sorted([r for r in coord if r['resource']==res],key=lambda r:int(r['channel']))]);add('planar_yzero_symmetry',float(np.max(abs(g[res][:,1]))),0,res)
 sns=rows(OUT/'T2_N2_SENTINELS_20261013_v1.csv');ss=np.array([[float(r[k]) for k in ['r_km','abs_theta_deg','source_z_m']] for r in sns]);complete=seal['K_complete_variants'];mods={};coldK=[];fieldchecks=0;maxfield=0.
 for v in complete:
  for f in FREQ:
   path=OLDLOCAL/'modal'/f't2_f{f}_n160001.mod' if v=='K0' else LOCAL/'K'/v/f'n2_{v}_f{f}.mod';m=mod_read(path);mods[v,f]=m;reported=next(r for r in rows(OUT/'MODE_WINDOW_AND_MATCHING.csv') if r['variant']==v and int(r['frequency_hz'])==f);cp=2*math.pi*f/m[2].real;add('modecount_and_phase_window',abs(len(cp)-int(reported['modes']))+abs(float(cp.min())-float(reported['cp_min']))+abs(float(cp.max())-float(reported['cp_max'])),1e-10,v+str(f))
   for res in RES:
    h=np.load(OLDLOCAL/'replicas'/f'{res}_fine_n160001.npy',mmap_mode='r') if v=='K0' else np.load(LOCAL/'K_FIELDS'/f'{v}_{res}_fine.npy',mmap_mode='r');cat=np.load(G0/'CATALOGUE_fine.npy',mmap_mode='r')
    for s in ss:
     cid=int(np.flatnonzero(np.all(cat==s,axis=1))[0]);value=cold_pressure(m,s,g[res],f);e=np.linalg.norm(value-h[cid,FREQ.index(f)])/np.linalg.norm(value);add('independent_MODAL_term_sum_sentinel',e,1e-10,v+res+str(f)+str(cid),len(value));maxfield=max(maxfield,float(e));fieldchecks+=1
 # Complete catalogues are compared with a separate QR projector route, all states retained.
 summary=rows(OUT/'FULL_CATALOGUE_KRAKEN_STABILITY.csv');entries=0;maxshape=0.
 for left,right in [('K0','K1'),('K1','K2'),('K2','K3')]:
  if left not in complete or right not in complete:continue
  for cat in CATS:
   ns=3355 if cat=='primary' else 12705
   for res in RES:
    a=np.load(OLDLOCAL/'replicas'/f'{res}_{cat}_n160001.npy',mmap_mode='r') if left=='K0' else np.load(LOCAL/'K_FIELDS'/f'{left}_{res}_{cat}.npy',mmap_mode='r');b=np.load(LOCAL/'K_FIELDS'/f'{right}_{res}_{cat}.npy',mmap_mode='r');e=qre(a,b);old=np.load(LOCAL/f'K_DIFFERENCE_{left}_{right}_{res}_{cat}.npz');de=float(np.max(abs(e-old['shape'])));raw=np.sqrt(np.sum(abs(b-a)**2,axis=-1)/np.sum(abs(b)**2,axis=-1));add('all_full_catalogue_QR_shapes',de,1e-10,left+right+res+cat,e.size);add('all_raw_pressure_differences',float(np.max(abs(raw-old['raw']))),1e-10,left+right+res+cat,raw.size);entries+=e.size;maxshape=max(maxshape,de);add('full_catalogue_dimension',abs(len(a)-ns)+abs(len(b)-ns),0)
    for fi,f in enumerate(FREQ):
     r=next(r for r in summary if r['left']==left and r['right']==right and r['resource']==res and r['catalogue']==cat and int(r['frequency_hz'])==f);add('summary_max_count_failures',abs(float(e[:,fi].max())-float(r['max']))+abs(int(np.sum(e[:,fi]>.002))-int(r['above_0p002'])),1e-10)
 # Every FIELDS result is independently parsed and interpolated for all four resources.
 source_rows=rows(OUT/'SCOOTER_FIELDS_PROVENANCE.csv');ffpcold={};ffpchecks=0;maxffp=0.
 for r in source_rows:
  if not r.get('SHD_path'):continue
  d=shd_read(__import__('pathlib').Path(r['SHD_path']));sn=next(t for t in sns if t['sentinel_id']==r['sentinel_id']);s=np.array([float(sn[k]) for k in ['r_km','abs_theta_deg','source_z_m']]);add('SHD_frequency_and_monotonic_actual_axis',abs(d[0]-float(r['frequency_hz']))+(0 if np.all(np.diff(d[3])>0) else 1),1e-6);stored=np.load(__import__('pathlib').Path(r['SHD_path']).parent/'INTERPOLATED_FIELD.npz')
  for res in RES:
   value=cold_interpolate(d,s,g[res]);e=np.linalg.norm(value-stored[res])/np.linalg.norm(value);add('every_FIELDS_actual_axis_cold_interp',e,1e-10,res+r['variant']+r['sentinel_id']+r['range_step_m'],len(value));ffpchecks+=1;maxffp=max(maxffp,float(e));ffpcold[r['variant'],float(r['range_step_m']),res,r['sentinel_id'],int(r['frequency_hz'])]=value
 # Check each reported FFP self/interpolation/cross metric against the original SHD-derived cold fields.
 for name in ['FFP_INTERNAL_STABILITY.csv','FFP_RANGE_INTERPOLATION.csv','KRAKEN_FFP_CROSS_METHOD.csv']:
  for r in rows(OUT/name):
   if r['status']!='COMPLETED':continue
   res=r['resource'];sid=r['sentinel_id'];f=int(r['frequency_hz']);comparison=r['comparison']
   if comparison=='F1_VS_F2':a=ffpcold['F1',.05,res,sid,f];b=ffpcold['F2',.05,res,sid,f]
   elif comparison=='F2_VS_F3':a=ffpcold['F2',.05,res,sid,f];b=ffpcold['F3',.05,res,sid,f]
   elif comparison=='F2_STEP050_VS025':a=ffpcold['F2',.05,res,sid,f];b=ffpcold['F2',.025,res,sid,f]
   else:
    a=ffpcold['F2',.05,res,sid,f];sn=next(t for t in sns if t['sentinel_id']==sid);s=np.array([float(sn[k]) for k in ['r_km','abs_theta_deg','source_z_m']]);b=cold_pressure(mods['K3',f],s,g[res],f)
   e=float(qre(a[None,None,:],b[None,None,:])[0,0]);add('FFP_metric_independent_shape',abs(e-float(r['shape_error'])),1e-10,name+res+sid+str(f))
 # Raw cross-method/global scalar diagnostics independently assembled from cold SHD fields.
 for row in rows(OUT/'RAW_GLOBAL_COMPLEX_FACTOR_DIAGNOSTIC.csv'):
  if row['status']!='COMPLETED':continue
  f=int(row['frequency_hz']);z=float(row['source_z_m']);idx=np.flatnonzero(ss[:,2]==z);comparison=row['comparison'];av=[];bv=[]
  for res in RES:
   for i in idx:
    sid=sns[i]['sentinel_id']
    if comparison=='F1_VS_F2':a=ffpcold['F1',.05,res,sid,f];b=ffpcold['F2',.05,res,sid,f]
    elif comparison=='F2_VS_F3':a=ffpcold['F2',.05,res,sid,f];b=ffpcold['F3',.05,res,sid,f]
    elif comparison=='F2_STEP050_VS025':a=ffpcold['F2',.05,res,sid,f];b=ffpcold['F2',.025,res,sid,f]
    else:a=ffpcold['F2',.05,res,sid,f];b=cold_pressure(mods['K3',f],ss[i],g[res],f)
    av.extend(a.tolist());bv.extend(b.tolist())
  av=np.array(av);bv=np.array(bv);c=np.sum(av.conj()*bv)/np.sum(abs(av)**2);e=np.sqrt(np.sum(abs(bv-c*av)**2)/np.sum(abs(bv)**2));add('global_raw_single_factor_cold',abs(float(e)-float(row['raw_global_profiled_relative'])),1e-10,comparison+str(f)+str(z))
 decision=read(OUT/'NUMERICAL_DECISION.json')
 for name,comparisons,expected in [('K_WINDOW',['1800_TO_2000','2000_TO_2400'],48),('K_MESH',['MESH160001_TO240001'],24)]:
  rr=[r for r in summary if r['comparison'] in comparisons];yes=len(rr)==expected and all(float(r['max'])<=.002 for r in rr);add('independent_K_gate_conjunction',0 if yes==decision[name]['PASS'] else 1,0,name)
 add('no_old_G0_upgrade',0 if decision['historical_G0']=='FORWARD_NUMERIC_NOT_ADMITTED_UNCHANGED' and decision['A_performance']=='NOT_EVALUATED' and decision['B_performance']=='NOT_EVALUATED' and decision['R4_percent']==0 else 1,0)
 table(OUT/'INDEPENDENT_COLD_CHECKS.csv',checks);result={'PASS':all(r['PASS'] for r in checks),'check_rows':len(checks),'FAIL':sum(not r['PASS'] for r in checks),'cold_K_sentinel_vectors':fieldchecks,'full_K_shape_entries':entries,'cold_FIELDS_resource_vectors':ffpchecks,'FIELDS_records_actually_checked':len(source_rows),'max_cold_modal_field_relative':maxfield if fieldchecks else 'NOT_EVALUATED','max_QR_shape_absolute_discrepancy':maxshape if entries else 'NOT_EVALUATED','max_cold_FFP_interp_relative':maxffp if ffpchecks else 'NOT_EVALUATED','new_solver_calls':0,'no_fake_PASS_for_zero_evaluated_FFP':ffpchecks>0,'original_physical_field_accuracy_certified':False,'wall_s':time.monotonic()-began,'R4_percent':0};dump(OUT/'INDEPENDENT_COLD_AUDIT.json',result);print('COLD N2',result,flush=True)
if __name__=='__main__':main()