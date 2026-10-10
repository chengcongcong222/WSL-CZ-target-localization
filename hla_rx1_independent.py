"""Cold verification after frozen once-only execution; no new solver or optimizer.
Independent direct-sum ACF, cosine spectrum and real receiver conversion.
"""
import csv,json,time
import numpy as np
from scipy import signal
from hla_rx1_common import *

def main():
 coldstart=time.monotonic();cfg=read(OUT/'DESIGN_FREEZE.json');status=read(OUT/'EXECUTION_STATUS.json');man=read(OUT/'RECEIVER_DATA_MANIFEST.json');checks=[]
 def check(name,e,lim):checks.append(dict(check=name,error=float(e),limit=lim,PASS=bool(e<=lim)))
 for p,h in cfg['input_and_code_sha256'].items():check('frozen hash '+p,0 if sha(ROOT/p)==h else 1,0)
 allrows=list(csv.DictReader((OUT/'EXTRACTED_FEATURES.csv').open(encoding='utf-8-sig'))) if (OUT/'EXTRACTED_FEATURES.csv').exists() else []
 by={}
 for x in allrows:by.setdefault(x['record_id'],[]).append(x)
 worstacf=0.;worstspectrum=0.;count=0
 for rec in man['records']:
  guard(coldstart,1800)
  check('receiver hash '+rec['record_id'],0 if sha(rec['path'])==rec['sha256'] else 1,0)
  check('spectra hash '+rec['record_id'],0 if sha(rec['correlation_spectra_path'])==rec['correlation_spectra_sha256'] else 1,0)
  with np.load(rec['correlation_spectra_path']) as data:
   for row in by[rec['record_id']]:
    key=f"w{row['window']}_f{row['line']}_{row['method']}_L{row['lag_s']}"
    y=data[key+'_residual'];c=data[key+'_correlation'];f=data[key+'_frequencies'];s=data[key+'_spectrum']
    direct=np.array([np.vdot(y[:len(y)-j],y[j:]).real for j in range(len(c))])/np.vdot(y,y).real
    worstacf=max(worstacf,float(np.max(abs(c-direct))))
    # Independently use cosine series at all saved detected peaks and fixed bins.
    ix=np.unique(np.r_[np.linspace(0,len(f)-1,33,dtype=int),[int(round(p['nu_hz']*16384)) for p in json.loads(row['all_peaks'])]])
    dd=direct.copy();dd[1:]*=2
    recomputed=np.maximum(np.cos(2*np.pi*f[ix,None]*np.arange(len(c))[None,:])@dd,0)
    worstspectrum=max(worstspectrum,float(np.max(abs(recomputed-s[ix]))/max(float(np.max(s)),1e-300)))
    check('mapping arithmetic '+rec['record_id']+' '+key,abs(float(row['q_eff_mps'])-cfg['receiver_metadata']['c_ref_mps']*float(row['nu_hz'])/float(row['f_rx_hz'])),1e-12)
    count+=1
 check('all direct-sum complex-record ACF',worstacf,1e-10);check('all independent cosine lag-spectrum peak bins',worstspectrum,1e-9)
 # Actual first stable/first drift multichannel records vs synthesized real receiver.
 equivalence=[]
 for family in ['STABLE','DRIFT']:
  rec=next((r for r in man['records'] if r['source_family']==family),None)
  if rec is None:continue
  x=np.load(rec['path'])[:800];n=200*1024;t=np.arange(n)/1024.;lo=np.array(cfg['receiver_metadata']['digital_reference_hz'])
  high=signal.resample(x,n,axis=0);real=np.real(np.sum(high*np.exp(2j*np.pi*t[:,None,None]*lo[None,:,None]),axis=1));analytic=signal.hilbert(real,axis=0)
  recovered=np.stack([signal.resample(analytic*np.exp(-2j*np.pi*f*t[:,None]),800,axis=0) for f in lo],axis=1)
  e=np.linalg.norm(recovered-x)/np.linalg.norm(x);check('main real/baseband equivalence '+family,e,1e-8);equivalence.append(dict(source_family=family,relative_error=float(e),record_id=rec['record_id'],scope='200s all8channels all3lines, exact Fourier receiver representation, same samples and edge leakage'))
 # Independently rebuild saved magnitude error summaries using csv literals.
 metrics=list(csv.DictReader((OUT/'FEATURE_ERROR_AND_COVERAGE.csv').open(encoding='utf-8-sig'))) if (OUT/'FEATURE_ERROR_AND_COVERAGE.csv').exists() else []
 maps=list(csv.DictReader((OUT/'WINDOW_FEATURE_MAPPING.csv').open(encoding='utf-8-sig'))) if (OUT/'WINDOW_FEATURE_MAPPING.csv').exists() else []
 for row in metrics:
  g=[a for a in maps if all(a[k]==row[k] for k in ['source_family','SNR_dB','method','lag_s'])];accepted=[a for a in g if a['accepted']=='True'];e=sorted(float(a['qeff_vs_abs_qbar_error_mps']) for a in accepted)
  for quant,name in [(.5,'conditional_median_abs_error_mps'),(.9,'conditional_P90_abs_error_mps'),(.95,'conditional_P95_abs_error_mps')]:
   v=e[int(np.ceil(len(e)*quant))-1] if e else float('inf');saved=float(row[name]);check('cold metric '+name+' '+str(tuple(row[k] for k in ['source_family','SNR_dB','method','lag_s'])),0 if v==saved else abs(v-saved),1e-12)
  check('cold denominator '+str(tuple(row[k] for k in ['source_family','SNR_dB','method','lag_s'])),abs(len(g)-int(row['evaluated_line_windows'])),0)
 table(OUT/'COLD_VALIDATION_CHECKS.csv',checks)
 dump(OUT/'VALIDATION.json',{'PASS':all(v['PASS'] for v in checks),'checks':len(checks),'FAIL':sum(not v['PASS'] for v in checks),'cold_line_window_lag_rows':count,'maximum_ACF_difference':worstacf,'maximum_relative_spectral_difference':worstspectrum,'real_baseband_equivalence':equivalence,'receiver_controls':read(OUT/'CONTROLS.json'),'scope':'algorithm arithmetic and receiver representation; no proof of radial qbar physical mapping','new_solver_calls':0,'R4_percent':0})
 print('COLD',len(checks),'checks',sum(not v['PASS'] for v in checks),'FAIL',flush=True)
if __name__=='__main__':main()
