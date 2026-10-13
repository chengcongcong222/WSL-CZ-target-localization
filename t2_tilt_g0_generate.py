"""Private once-only waveform generator, with source/noise pairing and full raw D archives."""
from t2_tilt_g0_core import *
import sys

def main(panel):
 if (OUT/f'{panel.upper()}_GENERATION_STARTED.json').exists():raise RuntimeError('generator already started; no repeat')
 assert read(OUT/'FORWARD_DECISION.json')['admitted']
 if panel=='test':assert (OUT/'THRESHOLD_FREEZE.json').exists()
 dump(OUT/f'{panel.upper()}_GENERATION_STARTED.json',{'utc':utc(),'once_only':True})
 start=time.monotonic();g=geometries();p=Provider(160001);contract=read(OUT/'EXECUTION_CONTRACT.json');ref=read(OUT/'NOISE_REFERENCE_FREEZE.json');nr={'20':ref['nu']['20'],'5':ref['nu']['5']}
 folder=LOCAL/panel;folder.mkdir();(folder/'waveforms').mkdir();(folder/'private').mkdir()
 if panel=='calibration':
  cp=rows(OUT/'CALIBRATION_GEOMETRIES.csv');base=[{'id':i*16+j,'geometry':i,'replicate':j,'r_km':float(t['r_km']),'theta_deg':float(t['theta_deg']),'z_m':float(t['z_m'])} for i,t in enumerate(cp) for j in range(16)]
 else:
  truth=[r for r in rows(TRUTH) if r['case_id'] in ['H01','H06','H12']]
  base=[{'id':i*48+iz*16+j,'geometry':i,'case_id':t['case_id'],'replicate':j,'r_km':float(t['r0_km']),'theta_deg':float(t['theta0_deg']),'z_m':z} for i,t in enumerate(truth) for iz,z in enumerate([185,200,215]) for j in range(16)]
 source_seed=contract['source'][f'{panel}_seed'];noise_seed=contract['noise'][f'{panel}_seed'];gain_seed=contract['gain'][f'{panel}_seed'];public=[];private=[];idx=0;fieldcache={}
 for b in base:
  guard(start);i=b['id'];s_rng=np.random.default_rng(np.random.SeedSequence([source_seed,i]));n_rng=np.random.default_rng(np.random.SeedSequence([noise_seed,i]));g_rng=np.random.default_rng(np.random.SeedSequence([gain_seed,i]));source=10**(s_rng.uniform(-3,3,(K,3))/20)*np.exp(1j*s_rng.uniform(-np.pi,np.pi,(K,3)));noise=n_rng.normal(size=(K,12,N));gain=10**(g_rng.uniform(-1,1,16)/20)*np.exp(1j*np.radians(g_rng.uniform(-30,30,16)))
  ip=folder/'private'/f'innovation_{i:04d}.npz';np.savez(ip,source=source,noise=noise,gain=gain);private.append({'base_id':i,'truth':b,'innovation_file':str(ip),'innovation_SHA256':sha(ip),'shared_noise_first8_SHA256':hashlib.sha256(noise[:,:8].tobytes()).hexdigest(),'extra4_SHA256':hashlib.sha256(noise[:,8:].tobytes()).hexdigest()})
  shared_wave={}
  for res in RES:
   key=(res,b['r_km'],b['theta_deg'],b['z_m'])
   if key not in fieldcache:fieldcache[key]=p.field([[b['r_km'],b['theta_deg'],b['z_m']]],g[res])[0]
   h=fieldcache[key];m=len(g[res])
   for c in ['C0','C1']:
    for snr in ['20','5']:
     usedgain=np.ones(m) if c=='C0' else (gain[:8] if m==8 else np.concatenate([gain[:8],gain[8:12] if res=='B1_V4' else gain[12:16]]));x=waveform(h,source,usedgain,noise[:,:m],nr[snr]);
     if res=='A0_0':shared_wave[(c,snr)]=x.copy()
     if res in ['B1_V4','B2_DUP']:x[:,:8]=shared_wave[(c,snr)]
     fp=folder/'waveforms'/f'{panel}_{idx:05d}.npy';np.save(fp,x)
     public.append({'record_id':idx,'base_id':i,'resource':res,'calibration':c,'reference_SNR_dB':snr,'waveform_path':str(fp),'waveform_SHA256':sha(fp),'fs_Hz':FS,'N':N,'K':K,'channels':m});idx+=1
  if i%16==0:table(OUT/f'{panel.upper()}_PUBLIC_RECORD_INDEX.csv',public);print('WAVEFORM',panel,i+1,'/',len(base),'configs',idx,flush=True)
 table(OUT/f'{panel.upper()}_PUBLIC_RECORD_INDEX.csv',public)
 dump(folder/'private/MANIFEST.json',{'panel':panel,'private_records':private,'generated_configs':idx,'base_innovations':len(base),'source_seed':source_seed,'noise_seed':noise_seed,'gain_seed':gain_seed})
 dump(OUT/f'{panel.upper()}_PRIVATE_MANIFEST.json',{'path':str(folder/'private/MANIFEST.json'),'SHA256':sha(folder/'private/MANIFEST.json'),'private_generator_not_receiver_or_estimator':True,'configs':idx,'base_innovations':len(base),'utc':utc()})
 dump(OUT/f'{panel.upper()}_GENERATION_FREEZE.json',{'utc':utc(),'record_index_SHA256':sha(OUT/f'{panel.upper()}_PUBLIC_RECORD_INDEX.csv'),'configs':idx,'science_wall_s':time.monotonic()-start,'raw_records_complete_on_D':True})
 print('WAVEFORMS SEALED',panel,idx,flush=True)
if __name__=='__main__':main(sys.argv[1])