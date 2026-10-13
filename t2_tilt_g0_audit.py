"""Independent read-only cold reconstruction. No solver, MC generation or reselection."""
from t2_tilt_g0_core import ROOT,OUT,LOCAL,RES,FREQ,read,rows,dump,table,sha,utc
import numpy as np,math,struct,json,time

def mod_independent(p):
 b=p.read_bytes();record=struct.unpack_from('<i',b,0)[0]*4;total,points=struct.unpack_from('<ii',b,92);nz=total;nm=struct.unpack_from('<i',b,5*record)[0];z=np.array(struct.unpack_from('<'+'f'*nz,b,4*record));phi=np.empty((nz,nm),complex)
 for j in range(nm):
  vals=np.array(struct.unpack_from('<'+'f'*(2*points),b,(7+j)*record));phi[:,j]=vals[::2]+1j*vals[1::2]
 vals=np.array(struct.unpack_from('<'+'f'*(2*nm),b,(7+nm)*record));kk=vals[::2]+1j*vals[1::2];return z,phi,kk

def cold_field(mod,r,theta,zs,xyz):
 depths,phi,kk=mod;source=int(np.argmin(abs(depths-zs)));assert abs(depths[source]-zs)<=1e-5;out=[]
 for x,y,zr in xyz:
  rec=int(np.argmin(abs(depths-zr)));assert abs(depths[rec]-zr)<=1e-5
  rho=math.hypot(r*1000*math.cos(math.radians(theta))-x,r*1000*math.sin(math.radians(theta))-y)
  out.append(np.sum(phi[source]*phi[rec]*np.sqrt(2*math.pi/(kk*rho))*np.exp(-1j*kk*rho-1j*math.pi/4),dtype=np.complex128))
 return np.array(out)

def main():
 began=time.monotonic();check=[]
 def add(name,e,t=1e-10,detail=''):
  check.append({'check':name,'error':float(e),'limit':t,'PASS':bool(np.isfinite(e) and e<=t),'detail':detail})
 freeze=read(OUT/'DESIGN_FREEZE.json')
 for name,expected in freeze['code_and_input_SHA256'].items():add('frozen_input_and_code_hash',0 if sha(ROOT/name)==expected else 1,0,name)
 coords=rows(OUT/'ARRAY_GEOMETRY_TABLE.csv');g={}
 for res in RES:
  rr=sorted([r for r in coords if r['resource']==res],key=lambda r:int(r['channel']));a=[]
  for row in rr:
   x=float(row['x_m']);y=float(row['y_m']);z=float(row['z_m']);alpha=float(row['tilt_from_horizontal_deg'])
   if row['arc_s_m']:
    u=float(row['arc_s_m']);add('independent_axis_coordinates',max(abs(x-u*math.cos(math.radians(alpha))),abs(z-200-u*math.sin(math.radians(alpha)))),1e-10,res+row['channel'])
   a.append([x,y,z])
  g[res]=np.array(a)
  add('planar_y0_mirror_geometry',np.max(abs(g[res][:,1])),0,res)
 add('positive_primary_span4m',abs(np.ptp(g['A_P4M'][:,2])-4),1e-10)
 mods={(m,int(f)):mod_independent(LOCAL/'modal'/f't2_f{int(f)}_n{m}.mod') for m in [80001,160001] for f in FREQ} if (OUT/'FORWARD_DECISION.json').exists() else {}
 if mods:
  for row in rows(OUT/'PROVIDER_MODE_WINDOW_AND_MESH.csv'):
   mod=mods[int(row['mesh']),int(float(row['frequency_hz']))];add('mode_count_all_returned',abs(len(mod[2])-int(row['modes'])),0,row['stem']);add('mode_MOD_hash',0 if sha(LOCAL/'modal'/(row['stem']+'.mod'))==row['MOD_SHA256'] else 1,0,row['stem'])
  summary=rows(OUT/'FORWARD_CATALOGUE_NUMERICAL_STABILITY.csv');fd=read(OUT/'FORWARD_DECISION.json');actualmax=0.
  for cat in ['primary','fine']:
   states=np.load(OUT/f'CATALOGUE_{cat}.npy')
   for res in RES:
    hi=np.load(LOCAL/'replicas'/f'{res}_{cat}_n160001.npy');lo=np.load(LOCAL/'replicas'/f'{res}_{cat}_n80001.npy');err=np.empty((len(states),3))
    # Scalar elimination independently via normalized projection, full data, no frozen score helper.
    for fi in range(3):
     a=lo[:,fi];b=hi[:,fi];c=np.sum(np.conj(a)*b,axis=1)/np.sum(np.abs(a)**2,axis=1);err[:,fi]=np.sqrt(np.sum(np.abs(b-c[:,None]*a)**2,axis=1)/np.sum(np.abs(b)**2,axis=1))
     report=next(r for r in summary if r['resource']==res and r['catalogue']==cat and int(float(r['frequency_hz']))==int(FREQ[fi]));add('crossmesh_summary_rebuilt',abs(float(err[:,fi].max())-float(report['max_profiled_shape_error'])),1e-12,res+cat);actualmax=max(actualmax,float(err[:,fi].max()))
    ids=np.unique([0,len(states)//2,len(states)-1,int(np.argmax(err.max(axis=1)))])
    for mesh,cached in [(80001,lo),(160001,hi)]:
     for idx in ids:
      s=states[idx]
      for fi,f in enumerate(FREQ):
       h=cold_field(mods[mesh,int(f)],*s,g[res]);den=np.linalg.norm(h);e=np.linalg.norm(h-cached[idx,fi])/den;add('direct_full_modal_cold_field',e,1e-5,f'{res}/{cat}/{mesh}/{idx}/{int(f)}');mirror=cold_field(mods[mesh,int(f)],s[0],-s[1],s[2],g[res]);add('cold_theta_mirror',np.max(abs(h-mirror)),1e-12)
  add('global_crossmesh_max_rebuilt',abs(actualmax-fd['max_cross_mesh_profiled_error']),1e-12)
  for f in FREQ:
   v=np.load(OUT/f'FIELD_COMPARISON_f{int(f)}.npz');pp=rows(OUT/f'FIELD_POINTS_f{int(f)}.csv');fresh=[];plan=read(OUT/'PROVIDER_PROVENANCE.json')['FIELD_plan']
   for row in pp:
    s=plan['states'][int(row['state'])];fresh.append(cold_field(mods[160001,int(f)],*s,g[row['resource']])[int(row['channel'])])
   fresh=np.array(fresh);add('FIELD_direct_predictor_cold',np.linalg.norm(fresh-v['pred'])/np.linalg.norm(fresh),1e-5,str(int(f)))
   c=np.sum(np.conj(fresh)*v['FIELD'])/np.sum(abs(fresh)**2);e=np.linalg.norm(v['FIELD']-c*fresh)/np.linalg.norm(v['FIELD']);r=next(r for r in rows(OUT/'FIELD_CONTROLS.csv') if int(float(r['frequency_hz']))==int(f));add('FIELD_one_common_scalar_error_rebuilt',abs(e-float(r['error'])),1e-10,str(int(f)))
  expected=all(r['PASS']=='True' for r in rows(OUT/'FORWARD_CONTROLS.csv')) and all(r['PASS']=='True' for r in summary);add('forward_admission_conjunction',0 if expected==fd['admitted'] else 1,0)
 # Independent DFT of every saved waveform, plus full deterministic source/noise reconstruction.
 recs=0;paired=0;max_dft=0.;max_wave=0.
 for panel in ['calibration','test']:
  ip=OUT/f'{panel.upper()}_PUBLIC_RECORD_INDEX.csv'
  if not ip.exists():continue
  public=rows(ip);yp=OUT/f'{panel.upper()}_DFT_INDEX.csv';df={int(r['record_id']):r for r in rows(yp)} if yp.exists() else {};private=read(LOCAL/panel/'private/MANIFEST.json');tr={r['base_id']:r for r in private['private_records']};last=-1;innov=None;truefields={};wave_pair={}
  for r in public:
   i=int(r['record_id']);base=int(r['base_id']);res=r['resource'];fp=__import__('pathlib').Path(r['waveform_path']);x=np.load(fp);add('raw_wave_hash',0 if sha(fp)==r['waveform_SHA256'] else 1,0,f'{panel}/{i}')
   fft=np.fft.rfft(x,axis=-1)[:,:,np.array([402,470,566])].transpose(0,2,1)/2048
   if i in df:
    y=np.load(df[i]['DFT_path']);e=float(np.linalg.norm(fft-y)/np.linalg.norm(y));max_dft=max(max_dft,e);add('independent_FFT_every_record',e,1e-10,f'{panel}/{i}');add('DFT_hash',0 if sha(df[i]['DFT_path'])==df[i]['DFT_SHA256'] else 1,0)
   if base!=last:
    innov=np.load(tr[base]['innovation_file']);last=base;truefields={};wave_pair={};add('private_innovation_hash',0 if sha(tr[base]['innovation_file'])==tr[base]['innovation_SHA256'] else 1,0)
   t=tr[base]['truth'];m=x.shape[1]
   if res not in truefields:truefields[res]=np.array([cold_field(mods[160001,int(f)],t['r_km'],t['theta_deg'],t['z_m'],g[res]) for f in FREQ])
   gain=np.ones(m) if r['calibration']=='C0' else (innov['gain'][:8] if m==8 else np.concatenate([innov['gain'][:8],innov['gain'][8:12] if res=='B1_V4' else innov['gain'][12:16]]));spectra=innov['source'][:,:,None]*gain[None,None,:]*truefields[res][None];timepoints=np.arange(2048)/1024;rebuild=np.zeros_like(x)
   for fi,f in enumerate(FREQ):
    rebuild+=2*(spectra[:,fi,:,None].real*np.cos(2*math.pi*f*timepoints)-spectra[:,fi,:,None].imag*np.sin(2*math.pi*f*timepoints))
   nu=read(OUT/'NOISE_REFERENCE_FREEZE.json')['nu'][r['reference_SNR_dB']];rebuild+=innov['noise'][:,:m]*math.sqrt(2048*nu);e=np.linalg.norm(rebuild-x)/np.linalg.norm(x);max_wave=max(max_wave,e);add('independent_every_waveform_reconstruction',e,1e-10,f'{panel}/{i}');recs+=1
   k=(r['calibration'],r['reference_SNR_dB'])
   if res=='A0_0':wave_pair[k]=x[:,:8].copy()
   if res in ['B1_V4','B2_DUP']:
    add('common_geometry_first8_raw_bytes',0 if wave_pair[k].tobytes()==x[:,:8].tobytes() else 1,0);paired+=1
  add('record_denominator',abs(len(public)-(9216 if panel=='calibration' else 5184)),0,panel)
 if (OUT/'THRESHOLD_FREEZE.json').exists():
  ds=rows(OUT/'CALIBRATION_NEAREST_SCORES.csv')
  for r in rows(OUT/'CALIBRATION_THRESHOLDS.csv'):
   keys=['resource','catalogue','mesh','calibration','reference_SNR_dB'];a=[float(t['D']) for t in ds if all(t[k]==r[k] for k in keys)];add('calibration256',abs(len(a)-256),0);add('rank244_threshold',abs(np.sort(a)[243]-float(r['tau_D'])),0)
 # Reconstruct score and acceptance boundary strata from sealed raw DFT and cached replicas.
 flips=0;max_score=0.;boundary_checks=0
 if (OUT/'INFERENCE_FREEZE.json').exists():
  pp={int(r['record_id']):r for r in rows(OUT/'TEST_DFT_INDEX.csv')};results=rows(OUT/'INFERENCE_RESULTS_NO_TRUTH.csv')
  for res in RES:
   for cat in ['primary','fine']:
    for mesh in [80001,160001]:
     subset=[r for r in results if r['resource']==res and r['catalogue']==cat and int(r['mesh'])==mesh];h=np.load(LOCAL/'replicas'/f'{res}_{cat}_n{mesh}.npy');d=np.load(LOCAL/'inference'/f'SCORES_{res}_{cat}_n{mesh}.npy',mmap_mode='r');add('inference576per_resource',abs(len(subset)-576),0)
     for i in [0,1,128,257,575]:
      r=subset[i];y=np.load(pp[int(r['record_id'])]['DFT_path']);tau=float(r['tau_D']);ids=np.unique([int(np.argmin(d[i])),int(np.argmin(abs(d[i]-tau))),0,len(h)-1]);obs=y/np.sqrt(np.sum(abs(y)**2,axis=-1))[...,None]
      for idx in ids:
       vec=h[idx]/np.sqrt(np.sum(abs(h[idx])**2,axis=-1))[:,None];qv=0.
       for k in range(8):
        for fi in range(3):qv+=abs(np.vdot(vec[fi],obs[k,fi]))**2/24
       fresh=1-qv;e=abs(fresh-d[i,idx]);max_score=max(max_score,e);add('independent_best_boundary_score',e,1e-7,f'{res}/{cat}/{mesh}/{i}/{idx}');flips+=int((fresh<=tau)!=(d[i,idx]<=tau));boundary_checks+=1
  add('cold_acceptance_flips',flips,0)
 # Independent summary/decision reconstruction from full sealed terminal labels.
 if (OUT/'ACCEPTANCE_AND_POINT_BY_CASE.csv').exists():
  evidence=rows(OUT/'ACCEPTANCE_AND_POINT_BY_CASE.csv');summaries=rows(OUT/'COMPARISON_BY_RESOURCE.csv');signals=rows(OUT/'PRIMARY_SIGNAL_BY_GEOMETRY.csv');cold_counts={}
  for r in summaries:
   keys=['case_id','generated_z_m','resource','calibration','reference_SNR_dB','catalogue','mesh'];a=[t for t in evidence if all(str(t[k])==str(r[k]) for k in keys)];add('summary_full16_denominator',abs(len(a)-16),0)
   accepted=sum(t['accepted_output']=='True' for t in a);retained=sum(t['horizontal_label_retained']=='True' for t in a);add('summary_accepted_and_truth_retained',abs(accepted-int(r['accepted']))+abs(retained-int(r['horizontal_retained'])),0)
   width=np.sort([float(t['envelope_km']) for t in a]);median=float((width[7]+width[8])/2);reported=float(r['median_unconditional_envelope_km']);add('unconditional_median_width',0 if median==reported else abs(median-reported),1e-12)
   error=np.sort([float(t['unconditional_range_error_km']) for t in a]);reported=float(r['P95_unconditional_point_error_km']);add('nearest_rank_P95_with_INF',0 if error[15]==reported else abs(error[15]-reported),1e-12)
  for r in signals:
   def find(res):return next(t for t in summaries if t['case_id']==r['case_id'] and float(t['generated_z_m'])==200 and t['resource']==res and t['calibration']=='C0' and t['reference_SNR_dB']=='20' and t['catalogue']==r['catalogue'] and t['mesh']==r['mesh'])
   t=find(r['treated']);c=find(r['control']);width=float(t['median_unconditional_envelope_km']);cw=float(c['median_unconditional_envelope_km']);yes=int(t['accepted'])>=15 and int(t['horizontal_retained'])>=15 and width<=10 and cw>0 and np.isfinite(cw) and width/cw<=.7
   add('primary_rule_independent',0 if yes==(r['signal']=='True') else 1,0)
   key=f"{r['comparison']}_{r['catalogue']}_{r['mesh']}";cold_counts[key]=cold_counts.get(key,0)+int(yes)
  decision=read(OUT/'DECISION.json');add('all_numeric_geometry_counts',0 if cold_counts==decision['all_numeric_signal_counts'] else 1,0)
  for l in ['A','B']:add('decision_main_count',abs(cold_counts[f'{l}_primary_160001']-decision['primary_geometry_counts'][l]),0)
  reversal=any((cold_counts[f'{l}_primary_160001']>=2)!=(cold_counts[f'{l}_{cat}_{mesh}']>=2) for l in ['A','B'] for cat in ['primary','fine'] for mesh in [80001,160001]);add('numerical_reversal_classification',0 if reversal==(decision['numerical_status']=='NUMERICAL_OR_GRID_LIMITED') else 1,0)
 if (OUT/'CALIBRATION_NEAREST_SCORES.csv').exists():
  public={int(r['record_id']):r for r in rows(OUT/'CALIBRATION_DFT_INDEX.csv')};cache={};ys={}
  for r in rows(OUT/'CALIBRATION_NEAREST_SCORES.csv'):
   key=(r['resource'],r['catalogue'],r['mesh']);rid=int(r['record_id'])
   if key not in cache:cache[key]=np.load(LOCAL/'replicas'/f"{key[0]}_{key[1]}_n{key[2]}.npy")
   if rid not in ys:ys[rid]=np.load(public[rid]['DFT_path'])
   h=cache[key][int(r['nearest_candidate_id'])];y=ys[rid];numer=np.abs(np.sum(h.conj()[None]*y,axis=-1))**2;den=np.sum(abs(h)**2,axis=-1)[None]*np.sum(abs(y)**2,axis=-1);coldD=1-float(np.sum(numer/den)/24)
   add('calibration_score_independent',abs(coldD-float(r['D'])),1e-7)
 table(OUT/'INDEPENDENT_COLD_CHECKS.csv',check)
 result={'PASS':all(r['PASS'] for r in check),'checks':len(check),'FAIL':sum(not r['PASS'] for r in check),'cold_records_rebuilt':recs,'common8_raw_pairs':paired,'all_resources_source_noise_pairing':'same privatebaseinnovation; tilted raw signal differs by geometry, commongeometry raw8 byte-equal','max_DFT_relative':max_dft,'max_waveform_relative':max_wave,'max_cold_score_absolute':max_score,'cold_boundary_checks':boundary_checks,'acceptance_flips':flips,'new_solver_MC_or_inference':0,'audit_wall_s':time.monotonic()-began,'frozen_hashes_checked':len(freeze['code_and_input_SHA256']),'external_mode_or_derivative_accuracy_certified':False}
 dump(OUT/'VALIDATION.json',result);print('COLD AUDIT',result,flush=True)
if __name__=='__main__':main()