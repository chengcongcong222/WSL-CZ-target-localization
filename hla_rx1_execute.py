"""Once-only post-Commit-A controls, provider and receiver panel."""
import csv,time,json,os,subprocess
from pathlib import Path
import numpy as np
from scipy import signal
from hla_rx1_common import *
from hla_rx1_generate import solve,group_provider,clean_field
from hla_rx1_extract import extract
from hla_rx1_controls import run

def noise(rng,shape,fs):
 z=(rng.normal(size=shape)+1j*rng.normal(size=shape))/np.sqrt(2)
 ft=np.fft.fft(z,axis=0);ft[abs(np.fft.fftfreq(len(z),1/fs))>.75]=0
 return np.fft.ifft(ft,axis=0)/np.sqrt(1.5/fs)

def main():
 cfg=read(OUT/'DESIGN_FREEZE.json');start=time.monotonic()
 if (OUT/'EXECUTION_STARTED.json').exists():raise RuntimeError('one frozen execution only; record already exists')
 for p,s in cfg['input_and_code_sha256'].items():
  if sha(ROOT/p)!=s:raise RuntimeError('frozen input changed '+p)
 dump(OUT/'EXECUTION_STARTED.json',{'time_utc':__import__('datetime').datetime.now(__import__('datetime').timezone.utc).isoformat(),'design_sha':subprocess.check_output(['git','rev-parse','HEAD'],text=True).strip(),'single_execution':True})
 allrows=[];manifest=[];diag=[];records=[];partial='';resident=0
 try:
  controls=run(cfg)
  if not all(x['PASS'] for x in controls):raise RuntimeError('receiver control failure; main provider and panel not executed')
  mods=solve(cfg,start);groups={};drift_ready=True
  for line in range(3):
   for mesh in [80001,160001]:
    center=mods[(line,mesh,0.)]
    try:
     fine,q=group_provider(center,mods[(line,mesh,-.006)],mods[(line,mesh,.006)],.006)
     coarse,qq=group_provider(center,mods[(line,mesh,-.012)],mods[(line,mesh,.012)],.012)
     error=float(np.max(abs(fine-coarse)/np.maximum(abs(fine),1e-30)))
     banderror=0.
     # Frozen central modal amplitude / linear dispersion approximation check.
     z,ph,k=center
     for delta in [-.006,.006]:
      za,pha,ka=mods[(line,mesh,delta)]
      for zs in [180.,200.,220.]:
       w=ph[np.flatnonzero(z==zs)[0]]*ph[np.flatnonzero(z==200.)[0]]
       wa=pha[np.flatnonzero(za==zs)[0]]*pha[np.flatnonzero(za==200.)[0]]
       r=np.array([45000.,52500.,60000.])
       pp=w@ (np.sqrt(2*np.pi/(k[:,None]*r))*np.exp(-1j*(k+delta*fine)[:,None]*r-1j*np.pi/4))
       aa=wa@ (np.sqrt(2*np.pi/(ka[:,None]*r))*np.exp(-1j*ka[:,None]*r-1j*np.pi/4))
       banderror=max(banderror,float(np.linalg.norm(pp-aa)/np.linalg.norm(aa)))
     diag.append(dict(line=line,mesh=mesh,modes=len(fine),minimum_mode_match=min(q+qq),group_secant_difference=error,narrowband_field_error=banderror,limit=.02,drift_ready=error<.02 and banderror<.01))
     if error>=.02 or banderror>=.01:drift_ready=False
     groups[(line,mesh)]=fine
    except Exception as e:
     drift_ready=False;diag.append(dict(line=line,mesh=mesh,modes=len(center[2]),minimum_mode_match=None,group_secant_difference=None,narrowband_field_error=None,limit=.02,drift_ready=False))
     groups[(line,mesh)]=np.ones(len(center[2]))*2*np.pi/1500
  table(OUT/'PROVIDER_GROUP_AND_BANDWIDTH_CHECK.csv',diag)
  truth=list(csv.DictReader((ROOT/cfg['truth_file']).open(encoding='utf-8-sig')))
  cache={};stability=[];meta=cfg['receiver_metadata'];fs=meta['sample_rate_hz']
  for gi,row in enumerate(truth):
   for family in ['STABLE','DRIFT']:
    if family=='DRIFT' and not drift_ready:continue
    vals={}
    for mesh in [80001,160001]:
     guard(start);vals[mesh]=np.stack([clean_field(row,line,family,cfg,mods[(line,mesh,0.)],groups[(line,mesh)]) for line in range(3)],axis=1)
    cache[(gi,family)]=vals[160001]
    for line in range(3):
     e=np.linalg.norm(vals[80001][:,line]-vals[160001][:,line])/np.linalg.norm(vals[160001][:,line])
     stability.append(dict(geometry=gi,source_family=family,line=line,relative_complex_field_grid_error=float(e),descriptive_limit=.01,PASS=bool(e<.01)))
    print('CLEAN_RECEIVER',gi,family,flush=True)
  table(OUT/'RECEIVER_PROVIDER_STABILITY.csv',stability)
  # Input source family/geometry labels below remain evaluation-side.
  waves=LOCAL/'receiver_records';waves.mkdir(exist_ok=True)
  rng=np.random.Generator(np.random.PCG64(cfg['seed']));record_id=0
  for gi in range(8):
   for rep in range(8):
    phase=rng.uniform(0,2*np.pi,3);amps=np.exp(rng.uniform(-.4,.4,3));zz=noise(rng,(4800,3,8),fs)
    # Same source amplitudes/phases and basic noise across source family and SNR.
    for family in ['STABLE','DRIFT']:
     for snr in [20,5]:
      rid=f'R{record_id:04d}';record_id+=1
      if (gi,family) not in cache:
       records.append(dict(record_id=rid,geometry=gi,source_family=family,SNR_dB=snr,replicate=rep,status='PROVIDER_NOT_READY'));continue
      guard(start);clean=cache[(gi,family)]*amps[None,:,None]*np.exp(1j*phase[None,:,None]);power=np.mean(abs(clean[:,:,0])**2,axis=0)
      sigma=np.sqrt(power/10**(snr/10));x=clean+zz*sigma[None,:,None]
      dest=waves/f'{rid}.npy';np.save(dest,x)
      ref={'source_frequency_hz':cfg['private_source_frequency_hz'],'drift_hz_per_s':cfg['private_drift_hz_per_s'] if family=='DRIFT' else [0.,0.,0.]}
      rr,spectra=extract(x,meta,ref,keep_spectra=True)
      specdest=waves/f'{rid}_correlation_spectra.npz';np.savez(specdest,**{key+'_'+name:value for key,vals in spectra.items() for name,value in zip(['times','residual','correlation','frequencies','spectrum'],vals)})
      tags=dict(record_id=rid,geometry=gi,source_family=family,SNR_dB=snr,replicate=rep)
      allrows.extend([dict(**tags,**v) for v in rr]);records.append(dict(**tags,status='COMPLETED'))
      manifest.append(dict(**tags,path=str(dest),sha256=sha(dest),correlation_spectra_path=str(specdest),correlation_spectra_sha256=sha(specdest),shape=list(x.shape),dtype=str(x.dtype),sample_rate_hz=fs,reference_channel=0,source_phase_private_rad=phase.tolist(),source_amplitude_private=amps.tolist(),reference_signal_power=power.tolist(),constant_noise_scale=sigma.tolist(),realized_noise_power=np.mean(abs(zz*sigma[None,:,None])**2,axis=0).tolist()))
      if rid in ['R0000','R0002']:np.savez(OUT/f'RECOMPUTABLE_{family}_SAMPLE.npz',samples=x[:800],clean=clean[:800],base_innovation=zz[:800],sigma=sigma,receiver_meta=json.dumps(meta))
      if record_id%32==0:print('RECEIVER_RECORDS',record_id,'/256',flush=True)
  if not drift_ready:partial='DRIFT_PROVIDER_NOT_READY'
 except Exception as e:
  partial=repr(e)
 finally:
  table(OUT/'EXTRACTED_FEATURES.csv',allrows);table(OUT/'RECORD_EXECUTION.csv',records)
  dump(OUT/'RECEIVER_DATA_MANIFEST.json',{'records':manifest,'receiver_public_metadata':cfg['receiver_metadata'],'private_source_frequencies_hz':cfg['private_source_frequency_hz'],'generator_source_time_drift_hz_per_s':cfg['private_drift_hz_per_s'],'PRNG':'numpy PCG64','seed':cfg['seed'],'paired_noise_across_families_and_SNR':True,'sample_interval':'[0,1200); no extra coherence','phase_preservation':'complex128 npy lossless','partial':partial})
  dump(OUT/'EXECUTION_STATUS.json',{'elapsed_s':time.monotonic()-start,'receiver_records_completed':len(manifest),'registered_records':256,'feature_rows':len(allrows),'partial':partial,'new_KRAKEN_calls':len(list((LOCAL/'modal').glob('*.stdout.txt'))) if (LOCAL/'modal').exists() else 0,'new_FIELD_calls':0,'H2_feedback':False,'R4_percent':0})
  print(read(OUT/'EXECUTION_STATUS.json'),flush=True)
if __name__=='__main__':main()
