"""N3 read-only physical audit and solver-free design preparation."""
from t2_n3_common import *
import shutil,subprocess,sys,zipfile,ast,platform

def main():
 if OUT.exists() or LOCAL.exists():raise RuntimeError('N3_ALREADY_PREPARED')
 assert subprocess.check_output(['git','rev-parse','HEAD'],text=True).strip()==PARENT
 assert subprocess.check_output(['git','ls-remote','origin','refs/heads/main'],text=True).split()[0]==PARENT
 manifest=read(PKG/'T2_N3_FILE_MANIFEST_20261013_v1.json')
 assert sha(PKG/'T2_N3_CODEX_TASK_20261013_v1.md')==TASKSHA
 zip_path=PKG.with_suffix('.zip');assert sha(zip_path)=='c3a54c3eef949f75b760191aa22bd718ad547c957e21cd3520f540192141a423'
 with zipfile.ZipFile(zip_path) as z:
  for r in manifest['files']:
   p=PKG/r['name'];assert sha(p)==r['sha256'] and p.stat().st_size==r['size_bytes']
   matches=[n for n in z.namelist() if Path(n).name==r['name']];assert len(matches)==1 and hashlib.sha256(z.read(matches[0])).hexdigest()==r['sha256']
 OUT.mkdir();LOCAL.mkdir();(OUT/'package').mkdir();(OUT/'environments').mkdir()
 for p in PKG.iterdir():
  if p.is_file():shutil.copyfile(p,OUT/'package'/p.name)
 (OUT/'.gitattributes').write_text('* -text\n',encoding='utf-8')
 identities=[]
 def identity(p,role,expected=None):
  p=Path(p);h=sha(p)
  if expected is not None:assert h==expected,(str(p),h,expected)
  r={'path':str(p),'sha256':h,'bytes':p.stat().st_size,'role':role}
  if p.is_relative_to(ROOT):
   rel=p.relative_to(ROOT).as_posix();blob=subprocess.check_output(['git','show',PARENT+':'+rel]);r['parent_blob_sha256']=hashlib.sha256(blob).hexdigest()
  identities.append(r)
 oldident=read(N2/'SOURCE_IDENTITY.json')['records'];kraken=AT/'kraken.exe';identity(kraken,'frozen_solver_binary',next(r['sha256'] for r in oldident if r.get('name')=='kraken'))
 for stage in (G0,ROOT/'results/R4_T2_N1_ARCHIVED_FIELD_SCORE_STABILITY',N2):
  for p in sorted(stage.rglob('*')):
   if p.is_file() and p.suffix not in ('.npy','.npz','.mod','.grn','.shd') and p.stat().st_size<8*1024**2:identity(p,'historical_read_only_report_or_input')
 frozen=read(N2/'SCIENCE_OUTPUT_FREEZE.json')['raw_files'];gindex={r['path']:r for r in rows(G0/'LOCAL_ARCHIVE_INDEX.csv')};envs={}
 for v in ('K0','K1','K2'):
  for f in FREQ:
   p=mod_path(v,f)
   for ext in ('.mod','.env','.prt'):
    q=p.with_suffix(ext)
    expected=gindex.get(str(q),{}).get('SHA256') if v=='K0' else frozen[str(q)]
    if expected is None:raise RuntimeError('ARCHIVE_IDENTITY_MISSING '+str(q))
    identity(q,'old_'+v+ext,expected)
   d=env_data(p.with_suffix('.env'));envs[v,f]=d
   assert d['f']==f and d['option']=='CVWT' and d['bottom']=='R' and d['ssp'][-1,0]==5000 and np.all(d['ssp'][:,3]==1) and np.all(d['ssp'][:,2]==0) and np.all(d['ssp'][:,4:]==0)
   prt=p.with_suffix('.prt').read_text(encoding='utf-8-sig');assert 'THORP attenuation added' in prt and 'Perfectly RIGID' in prt and 'VACUUM' in prt
   assert np.array_equal(d['ssp'],envs['K0',201]['ssp'])
 assert sha(OUT/'package/T2_N2_SENTINELS_20261013_v1.csv')==sha(N2/'T2_N2_SENTINELS_20261013_v1.csv')
 a=arrays();old=rows(G0/'ARRAY_GEOMETRY_TABLE.csv')
 for res in RES:
  coords=np.array([[float(r[k]) for k in ('x_m','y_m','z_m')] for r in sorted([r for r in old if r['resource']==res],key=lambda r:int(r['channel']))]);assert np.array_equal(coords,a[res])
 assert np.array_equal(a['A0_0'],a['B1_V4'][:8]) and np.array_equal(a['A0_0'],a['B2_DUP'][:8])
 ssp=envs['K0',201]['ssp'];cmin=float(ssp[:,1].min());cmax=float(ssp[:,1].max());assert cmin==1500 and cmax==1635.1811
 classification=[];bands=[];modal_summary=[];raw={}
 for v in ('K0','K1','K2'):
  for f in FREQ:
   mod=parse_mod(mod_path(v,f));z,phi,k,mf=mod;assert mf==f;cp=2*np.pi*f/k.real
   group=np.where(cp<=cmax,0,np.where(cp<=1800,1,np.where(cp<=2000,2,3)));assert np.all(cp>1500) and np.all(cp<=envs[v,f]['chigh'])
   for i,c in enumerate(cp):
    kz2=(2*np.pi*f/ssp[:,1])**2-(2*np.pi*f/c)**2
    classification.append({'variant':v,'frequency_hz':f,'mode_index':i+1,'Re_kr':float(k[i].real),'Im_kr':float(k[i].imag),'cp_mps':float(c),'cmax_mps':cmax,'band_index':int(group[i]),'band':('(1500,cmax]','(cmax,1800]','(1800,2000]','(2000,2400]')[group[i]],'min_kz2_at_SSP_nodes':float(kz2.min()),'max_kz2_at_SSP_nodes':float(kz2.max()),'cp_above_cmax':bool(c>cmax),'no_internal_water_turning':bool(c>cmax and np.all(kz2>0)),'source_MOD_SHA256':sha(mod_path(v,f))})
   modal_summary.append({'variant':v,'frequency_hz':f,'modes':len(k),'cp_min':float(cp.min()),'cp_max':float(cp.max()),'above_cmax':int(np.sum(cp>cmax)),'bands':[int(np.sum(group==j)) for j in range(4)]})
   for sn in sentinels():
    pos=tuple(float(sn[k]) for k in ('r_km','abs_theta_deg','source_z_m'))
    for res,coords in a.items():
     full=modal_field(mod,pos,coords);parts=[modal_field(mod,pos,coords,group==j) for j in range(4)];error=float(np.linalg.norm(sum(parts)-full)/np.linalg.norm(full));assert error<=1e-10
     totalnorm=float(np.linalg.norm(full));sum_norm=sum(float(np.linalg.norm(x)) for x in parts);cancel=totalnorm/sum_norm
     raw[f'{v}_{f}_{sn["sentinel_id"]}_{res}']=np.stack([full]+parts)
     for j,p in enumerate(parts):bands.append({'variant':v,'frequency_hz':f,'sentinel_id':sn['sentinel_id'],'resource':res,'r_km':pos[0],'abs_theta_deg':pos[1],'source_z_m':pos[2],'band_index':j,'band_modes':int(np.sum(group==j)),'band_norm':float(np.linalg.norm(p)),'full_norm':totalnorm,'band_over_full_norm':float(np.linalg.norm(p)/totalnorm),'complex_cancellation_full_over_sum_band_norms':cancel,'relative_reconstruction_error':error,'weak_relative_to_same_variant_frequency_median':'PENDING','scope':'COHERENT_COMPLEX_CONTRIBUTION_NOT_INDEPENDENT_RAY_ENERGY'})
 for v in ('K0','K1','K2'):
  for f in FREQ:
   norms=[r['full_norm'] for r in bands if r['variant']==v and r['frequency_hz']==f];ref=float(np.median(norms))
   for r in bands:
    if r['variant']==v and r['frequency_hz']==f:r['weak_relative_to_same_variant_frequency_median']=bool(r['full_norm']<.01*ref);r['same_variant_frequency_reference_norm']=ref
 np.savez(LOCAL/'OLD_MODAL_BAND_FIELDS.npz',**raw)
 table(OUT/'SSP_AND_TURNING_CLASSIFICATION.csv',classification);table(OUT/'OLD_N2_MODAL_BAND_CONTRIBUTIONS.csv',bands)
 dump(OUT/'PHYSICAL_SOURCE_AUDIT.json',{'status':'PHYSICAL_SOURCE_AND_TURNING_AUDIT_COMPLETE','c_min':cmin,'c_max':cmax,'axis_depths':ssp[ssp[:,1]==cmin,0].tolist(),'SSP_nodes':len(ssp),'water_depth':5000,'option':'CVWT','top':'V','bottom':'R','water_rho':1,'modal_records':len(classification),'contribution_records':len(bands),'mode_summary':modal_summary,'no_real_seabed_parameters':True,'all_complex_band_reconstruction_PASS':True})
 dump(OUT/'OLD_N2_SOURCE_IDENTITIES.json',{'parent_SHA':PARENT,'records':identities})
 (OUT/'PHYSICAL_ENV_CONTRACT.md').write_text(f'# Original nominal E-STD physical contract\n\nRaw ENV/PRT/MOD identities verified against sealed history. H=5000m; surface V; bottom R; C-linear SSP; CVWT with actual PRT Thorp; water density1g/cm3. SSP has {len(ssp)} nodes, cmin={cmin}, cmax={cmax}, axis1300m. All high-cp new modes above cmax have positive kz^2 at every SSP node; C-linear interpolation cannot create a larger internal c. This establishes no internal water-column turning point in this declared model.\n\nIt does not establish a received bottom-hit count or real high-mode energy. Coherent modal-band norms are not independent ray energies. Actual seabed sound speed, density, attenuation and terrain are UNKNOWN. No measured client ocean parameters were substituted.\n\nBoundary meanings: [official KRAKEN manual](https://oalib-acoustics.org/website_resources/AcousticsToolbox/manual/node47.html), [boundary conditions](https://oalib-acoustics.org/website_resources/AcousticsToolbox/manual/node10.html). Current archived PRT defines runtime absorption, avoiding treating legacy manual option layouts as this executable build.\n',encoding='utf-8')
 ref=LOCAL/'reference_precheck';shutil.copytree(OUT/'package',ref)
 subprocess.run([sys.executable,'-B',str(ref/'t2_n3_analytic_reference.py'),'--outdir',str(ref)],check=True)
 subprocess.run([sys.executable,'-B',str(ref/'t2_n3_independent_reference_check.py')],check=True)
 checks=[]
 for name in ('ANALYTIC_WAVEGUIDE_MODE_COUNTS.csv','ANALYTIC_WAVEGUIDE_WINDOW_SHAPES.csv','ANALYTIC_EXACT_RADIAL_VS_FARFIELD.csv'):
  supplied=rows(OUT/'package'/name);computed=rows(ref/name);assert len(supplied)==len(computed)
  for i,(x,y) in enumerate(zip(supplied,computed)):
   delta=0.
   for key in x:
    try:d=abs(float(x[key])-float(y[key]));delta=max(delta,d);assert d<=1e-8
    except ValueError:assert x[key]==y[key]
   checks.append({'file':name,'row':i,'max_abs_difference':delta,'PASS':True})
 assert max(float(r['exact_Hankel_vs_farfield_shape']) for r in rows(ref/'ANALYTIC_EXACT_RADIAL_VS_FARFIELD.csv'))<=1e-4
 table(OUT/'ANALYTIC_REFERENCE_VALIDATION.csv',checks)
 shutil.copyfile(ref/'ANALYTIC_REFERENCE_INDEPENDENT_CHECK.json',OUT/'PRECHECK_INDEPENDENT_REFERENCE.json')
 plans=[];receiver=sorted(set(np.concatenate([x[:,2] for x in a.values()]).tolist()+[10.,500.,1250.,2500.,3750.,4990.]))
 for mesh,windows in ((160001,(1800,2400,4000)),(80001,(2400,))):
  for bc in ('PR_R','PR_PR'):
   for f in (201,283):
    for high in windows:
     stem=f'n3_{bc}_f{f}_h{high}_n{mesh}';p=OUT/'environments'/(stem+'.env')
     lines=[f"'{stem}'",f'{f}.0','1',"'CVW'",f'{mesh} 0.0 5000.0','0.0 1500.0 0.0 1.0 0.0 0.0','5000.0 1500.0 0.0 1.0 0.0 0.0',"'R' 0.0" if bc=='PR_R' else "'V' 0.0",f'1500.0 {high}.0','60.0','4','180.0','200.0','210.0','220.0',str(len(receiver))]+[f'{x:.12f}' for x in receiver]+['R']
     p.write_text('\n'.join(lines)+'\n',encoding='utf-8',newline='\n');d=env_data(p);assert d['option']=='CVW' and np.all(d['ssp'][:,4:]==0)
     plans.append({'call_order':len(plans)+1,'kind':'ANALYTIC','boundary':bc,'f':f,'high':high,'mesh':mesh,'stem':stem,'ENV':str(p),'sha256':sha(p),'timeout_s':180})
 for f in FREQ:
  d=envs['K2',f];lines=d['lines'].copy();lines[0]=f"'n3_ESTD_V_f{f}'";lines[d['b']]=lines[d['b']].replace("'R'","'V'");assert [i for i in range(len(lines)) if lines[i]!=d['lines'][i]]==[0,d['b']]
  p=OUT/'environments'/f'n3_ESTD_V_f{f}.env';p.write_text('\n'.join(lines)+'\n',encoding='utf-8',newline='\n')
  plans.append({'call_order':len(plans)+1,'kind':'E_STD_V','boundary':'V','f':f,'high':2400,'mesh':160001,'stem':p.stem,'ENV':str(p),'sha256':sha(p),'timeout_s':180,'allowed_changed_lines':[0,d['b']]})
 dump(OUT/'SOLVER_PLAN.json',plans)
 ledger=ROOT/'results/R4_MASTER/R4_EVIDENCE_LEDGER.csv';shutil.copyfile(ledger,LOCAL/'LEDGER_ORIGINAL.csv');dump(OUT/'LEDGER_BEFORE.json',{'working_sha256':sha(ledger),'parent_blob_sha256':hashlib.sha256(subprocess.check_output(['git','show',PARENT+':results/R4_MASTER/R4_EVIDENCE_LEDGER.csv'])).hexdigest()})
 finalize()

def finalize():
 if (OUT/'DESIGN_FREEZE.json').exists():raise RuntimeError('DESIGN_ALREADY_FROZEN')
 verify_sources()
 for p in ROOT.glob('t2_n3_*.py'):ast.parse(p.read_text(encoding='utf-8-sig'))
 files=[p for p in OUT.rglob('*') if p.is_file()]+list(ROOT.glob('t2_n3_*.py'))+[LOCAL/'OLD_MODAL_BAND_FIELDS.npz',LOCAL/'LEDGER_ORIGINAL.csv']
 dump(OUT/'DESIGN_FREEZE.json',{'parent':PARENT,'task_SHA256':TASKSHA,'spec':read(OUT/'package/T2_N3_EXPERIMENT_SPEC_20261013_v1.json'),'frozen_files':[{'path':str(p),'sha256':sha(p),'bytes':p.stat().st_size} for p in files],'noise_MC_receiver':0,'weak_field_flag':'norm<0.01 same_variant_frequency_median; descriptive only; never removed','analytic_mesh_gate':'same analytic count/kr/field hard gates applied also to4coarse runs; intermeshdifference diagnostic only','mode_shape_diagnostic':'single mode common complex phase only; near nodes|analytic phi|<1e-6maxabs flagged with absolute error; no gate substitution','stop_on_first_hard_gate':True,'budget':guard(),'python':sys.version,'numpy':np.__version__})
 print('FROZEN',len(files),'code/input files;19planned calls;0solver calls',flush=True)
if __name__=='__main__':main()