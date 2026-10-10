"""Supplementary independent review of frozen raw receiver records, not a new trial.
Written after A; does not replace or modify frozen generator/extractor/evaluator.
"""
import csv,json,time
import numpy as np
from scipy.optimize import linear_sum_assignment
from hla_rx1_common import ROOT,OUT,LOCAL,read,dump,table,sha

def raw_modes(p):
 b=p.read_bytes();stride=int.from_bytes(b[:4],'little',signed=True)*4
 ndepth=int.from_bytes(b[92:96],'little',signed=True);n=int.from_bytes(b[5*stride:5*stride+4],'little',signed=True)
 z=np.frombuffer(b, dtype='<f4',count=ndepth,offset=4*stride).astype(float)
 phi=np.column_stack([np.frombuffer(b,dtype='<c8',count=ndepth,offset=(7+i)*stride) for i in range(n)]).astype(complex)
 k=np.frombuffer(b,dtype='<c8',count=n,offset=(7+n)*stride).astype(complex)
 return z,phi,k

def match(ref,other):
 ph,k=ref[1:];po,ko=other[1:]
 cor=abs(np.conj(ph).T@po)/(np.linalg.norm(ph,axis=0)[:,None]*np.linalg.norm(po,axis=0)[None,:]+1e-300)
 cost=1-cor+.01*abs(k[:,None]-ko[None,:])/max(np.median(np.diff(np.sort(k.real))),1e-8)
 a,b=linear_sum_assignment(cost);ix=np.empty(len(k),int);ix[a]=b
 return ko[ix]

def main():
 start=time.monotonic();cfg=read(OUT/'DESIGN_FREEZE.json');manifest=read(OUT/'RECEIVER_DATA_MANIFEST.json');meta=cfg['receiver_metadata'];checks=[]
 def check(name,e,lim):checks.append(dict(check=name,error=float(e),limit=lim,PASS=bool(e<=lim)))
 truth=list(csv.DictReader((ROOT/cfg['truth_file']).open(encoding='utf-8-sig')))
 index={(c['line'],c['mesh'],c['offset_hz']):LOCAL/'modal'/f"{c['stem']}.mod" for c in cfg['propagation_plan']}
 models={}
 if all(p.exists() for p in index.values()):
  for l in range(3):
   a=raw_modes(index[(l,160001,0.)]);minus=raw_modes(index[(l,160001,-.006)]);plus=raw_modes(index[(l,160001,.006)])
   diag=list(csv.DictReader((OUT/'PROVIDER_GROUP_AND_BANDWIDTH_CHECK.csv').open(encoding='utf-8-sig')))
   registered=next(x for x in diag if int(x['line'])==l and int(x['mesh'])==160001)
   if registered['minimum_mode_match'] and minus[1].shape==a[1].shape==plus[1].shape:
    dk=(match(a,plus).real-match(a,minus).real)/.012
   else:dk=np.full(len(a[2]),2*np.pi/1500)
   models[l]=(a,dk)
 samples=[0,137,799,1801,2400,3199,4001,4799];forwardworst=0.;phaseworst=0.;stftworst=0.;base={}
 rng=np.random.Generator(np.random.PCG64(cfg['seed']))
 record_by={(r['geometry'],r['replicate'],r['source_family'],r['SNR_dB']):r for r in manifest['records']}
 for gi in range(8):
  row=truth[gi];r0=float(row['r_km'])*1000;theta=np.deg2rad(float(row['theta_deg']));p0=r0*np.array([np.cos(theta),np.sin(theta)]);v=float(row['v_mps'])*np.array([np.cos(np.deg2rad(float(row['psi_deg']))),np.sin(np.deg2rad(float(row['psi_deg'])))])
  for rep in range(8):
   phases=rng.uniform(0,2*np.pi,3);amps=np.exp(rng.uniform(-.4,.4,3));zz=(rng.normal(size=(4800,3,8))+1j*rng.normal(size=(4800,3,8)))/np.sqrt(2)
   zzf=np.fft.fft(zz,axis=0);zzf[abs(np.fft.fftfreq(4800,.25))>.75]=0;zz=np.fft.ifft(zzf,axis=0)/np.sqrt(1.5/4)
   for family in ['STABLE','DRIFT']:
    for snr in [20,5]:
     rec=record_by.get((gi,rep,family,snr))
     if rec is None:continue
     x=np.load(rec['path']);check('private phase replay '+rec['record_id'],max(abs(phases-np.array(rec['source_phase_private_rad']))),1e-15)
     # Independent Newton retarded emission versus saved receiver-noise subtraction.
     for l in range(3):
      if l not in models:continue
      (depths,phi,k),dk=models[l];vg=2*np.pi/dk;weights=phi[np.flatnonzero(depths==float(row['z_label']))[0]]*phi[np.flatnonzero(depths==200)[0]]
      f0=cfg['private_source_frequency_hz'][l];alpha=cfg['private_drift_hz_per_s'][l] if family=='DRIFT' else 0.
      for idx in samples:
       t=idx/4;post=max(t-600,0);center=np.array([2*min(t,600)+2*post*np.cos(np.pi/12),2*post*np.sin(np.pi/12)])
       rx=center+np.column_stack([meta['offsets_m'],np.zeros(8)]);u=np.full((8,len(k)),t)-np.linalg.norm(p0+t*v-rx,axis=1)[:,None]/vg
       for _ in range(6):
        d=p0+u[:,:,None]*v-rx[:,None,:];r=np.linalg.norm(d,axis=-1);dr=np.sum(d*v,axis=-1)/r
        u-=(u+r/vg-t)/(1+dr/vg)
       d=p0+u[:,:,None]*v-rx[:,None,:];r=np.linalg.norm(d,axis=-1);delay=t-u
       physical_phase=2*np.pi*f0*u+np.pi*alpha*u*u+2*np.pi*f0*delay-k[None,:]*r-2*np.pi*meta['digital_reference_hz'][l]*t
       pred=np.sum(weights*np.sqrt(2*np.pi/(k[None,:]*r))*np.exp(1j*physical_phase-1j*np.pi/4),axis=1)*amps[l]*np.exp(1j*phases[l])
       clean=x[idx,l]-zz[idx,l]*rec['constant_noise_scale'][l]
       e=np.linalg.norm(pred-clean)/max(np.linalg.norm(pred),1e-300);forwardworst=max(forwardworst,float(e));phaseworst=max(phaseworst,float(np.max(abs(u+r/vg-t))))
     if rep==0 and snr==20:
      with np.load(rec['correlation_spectra_path']) as spectra:
       # Manual fixed steering and explicit Hann samples, independent of extractor.
       w=np.exp(-2j*np.pi*np.array(meta['digital_reference_hz'])[:,None]*np.array(meta['offsets_m'])[None,:]/meta['c_ref_mps'])/8
       beam=np.sum(x*w[None,:,:],axis=2);hann=(1-np.cos(2*np.pi*np.arange(33)/32))/2;hann/=hann.sum()
       ys=np.stack([np.sum(beam[i:i+33]*hann[:,None],axis=0) for i in range(0,len(x)-32,4)]);ts=(np.arange(len(ys))*4+16)/4
       rows=list(csv.DictReader((OUT/'EXTRACTED_FEATURES.csv').open(encoding='utf-8-sig')))
       rr=[a for a in rows if a['record_id']==rec['record_id'] and a['lag_s']=='200']
       for a in rr:
        j=int(a['window']);l=int(a['line']);use=(ts>=200*j)&(ts<200*(j+1));tt=ts[use]
        if a['method']=='RX':phase=2*np.pi*float(a['rx_center_offset_hz'])*tt
        else:phase=2*np.pi*(cfg['private_source_frequency_hz'][l]-meta['digital_reference_hz'][l])*tt+np.pi*(cfg['private_drift_hz_per_s'][l] if family=='DRIFT' else 0)*tt**2
        y=ys[use,l]*np.exp(-1j*phase);saved=spectra[f"w{j}_f{l}_{a['method']}_L200_residual"]
        stftworst=max(stftworst,float(np.linalg.norm(y-saved)/max(np.linalg.norm(y),1e-300)))
 check('all states/SNR/source sampled Newton forward reconstruction',forwardworst,1e-7)
 check('all states independent retarded-time residual seconds',phaseworst,1e-8)
 check('all states/source raw waveform fixed beam/STFT/reference provenance',stftworst,1e-10)
 table(OUT/'COLD_SIGNAL_PROVENANCE.csv',checks);dump(OUT/'SIGNAL_PROVENANCE_VALIDATION.json',{'PASS':all(x['PASS'] for x in checks),'checks':len(checks),'FAIL':sum(not x['PASS'] for x in checks),'forward_relative_worst':forwardworst,'retarded_time_residual_s':phaseworst,'raw_record_to_residual_relative_worst':stftworst,'elapsed_s':time.monotonic()-start,'source_truth_not_returned_to_RX':True,'scope':'cold deterministic saved-record review only; no changed provider or experiment'})
 print('COLD_SIGNAL',len(checks),'checks',sum(not x['PASS'] for x in checks),'FAIL',flush=True)
if __name__=='__main__':main()
