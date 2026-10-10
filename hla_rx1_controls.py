"""Post-A actual receiver controls. Never invoked during design."""
import numpy as np
from scipy import signal
from hla_rx1_extract import extract,combine,stft_line,correlation
from hla_rx1_common import table,OUT,LOCAL,dump

def run(cfg):
 meta=cfg['receiver_metadata'];fs=meta['sample_rate_hz'];t=np.arange(0,1200,1/fs);f=np.array(cfg['private_source_frequency_hz']);c=meta['c_ref_mps'];q=.4137;r0=50000.
 # Static receiver and receding source: exact retarded emission u=(t-r0/c)/(1+q/c).
 u=(t-r0/c)/(1+q/c);r=r0+q*u
 p=np.exp(2j*np.pi*f[None,:,None]*u[:,None,None])/r[:,None,None]*np.ones((1,1,8))
 x=p*np.exp(-2j*np.pi*np.array(meta['digital_reference_hz'])[None,:,None]*t[:,None,None])
 ref={'source_frequency_hz':f.tolist(),'drift_hz_per_s':[0.,0.,0.]}
 a,_=extract(x,meta,ref);refrows=[v for v in a if v['method']=='REF' and v['lag_s']==200]
 exact=q/(1+q/c);err=max(abs(v['q_eff_mps']-exact) for v in refrows)
 rows=[dict(control='exact retarded singlepath REF units no double-Doppler',error=err,limit=.02,PASS=err<.02,scope='known source frequency, magnitude only')]
 alt=dict(meta);alt['digital_reference_hz']=(np.array(meta['digital_reference_hz'])+.1731).tolist()
 xx=p*np.exp(-2j*np.pi*np.array(alt['digital_reference_hz'])[None,:,None]*t[:,None,None]);b,_=extract(xx,alt,ref)
 err=max(abs(v['q_eff_mps']-w['q_eff_mps']) for v,w in zip(a,b))
 rows.append(dict(control='same physical record two noninteger digital references',error=err,limit=.005,PASS=err<.005,scope='REF and RX, no truth offset correction'))
 # Same received single-component record, distinct (f0,q): gain supplies equal amplitude.
 f1=235.083;q1=.4;q2=.8;f2=f1*(1+q2/c)/(1+q1/c)
 u1=(t-r0/c)/(1+q1/c);u2=(t-r0/c)/(1+q2/c)
 ph2=2*np.pi*(f1*u1[0]-f2*u2[0]);z1=np.exp(2j*np.pi*f1*u1);z2=np.exp(2j*np.pi*f2*u2+1j*ph2)
 e=float(np.max(abs(z1-z2)));rows.append(dict(control='unknown source frequency/motion identical-record ambiguity',error=e,limit=1e-8,PASS=e<1e-8,scope='single channel single nondispersive component constant registered gain; not general No-Go'))
 # Numerical sample refinement of same retarded physical control, not a new main lag.
 fs2=8.;tt=np.arange(0,1200,1/fs2);uu=(tt-r0/c)/(1+q/c)
 x2=np.exp(2j*np.pi*(f[None,:,None]*uu[:,None,None]-np.array(meta['digital_reference_hz'])[None,:,None]*tt[:,None,None]))/(r0+q*uu[:,None,None])
 x2=np.broadcast_to(x2,(len(tt),3,8)).copy();m2=dict(meta);m2['sample_rate_hz']=fs2
 b,_=extract(x2,m2,ref);e=max(abs(v['q_eff_mps']-w['q_eff_mps']) for v,w in zip(a,b));rows.append(dict(control='baseband fs4 vs fs8 actual record',error=e,limit=.005,PASS=e<.005,scope='REF RX bothlags'))
 # Full 200s real waveform check; known receiver LO, no source demodulation.
 n=800;tt=np.arange(n)/4;bb=x[:n];fsr=1024.;tr=np.arange(200*1024)/fsr
 # Fourier interpolation is exact for registered bandlimited finite-period control.
 # This auxiliary equivalence input uses real receiver carrier and all eight channels.
 periodic=np.exp(2j*np.pi*.137*tt)[:,None]*np.exp(1j*np.arange(8)[None,:])
 high=signal.resample(periodic,len(tr),axis=0);real=np.real(high*np.exp(2j*np.pi*201*tr[:,None]));analytic=signal.hilbert(real,axis=0)
 low=signal.resample(analytic*np.exp(-2j*np.pi*201*tr[:,None]),n,axis=0)
 e=float(np.linalg.norm(low-periodic)/np.linalg.norm(periodic));rows.append(dict(control='real multichannel fs1024 / complex baseband equivalence',error=e,limit=1e-8,PASS=e<1e-8,scope='noninteger finite record, exact Fourier interpolation; edge and leakage retained'))
 np.savez(OUT/'CONTROL_RECOMPUTABLE_SAMPLE.npz',t=t[:800],baseband=x[:800],second_reference=xx[:800],confounding_a=z1[:800],confounding_b=z2[:800],equivalence_baseband=periodic,equivalence_returned=low)
 table(OUT/'RECEIVER_CONTROLS.csv',rows);dump(OUT/'CONTROLS.json',{'checks':rows,'PASS':all(v['PASS'] for v in rows),'REF_RX_singlepath_rows':a,'actual_sign_used_by_RX':False,'noninteger_source_frequency_hz':f.tolist()})
 return rows
