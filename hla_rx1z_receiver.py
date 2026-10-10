"""Receive-only saved-data representation, separate process before truth evaluation.
Inputs: sanitized receiver columns, public metadata and raw record index only.
"""
import numpy as np
from scipy import signal
from hla_rx1z_common import *

def frequency_rows(x,meta):
 fs=meta['sample_rate_hz'];off=np.asarray(meta['offsets_m']);lo=np.asarray(meta['digital_reference_hz']);c=meta['c_ref_mps']
 weights=np.exp(-2j*np.pi*lo[:,None]*off[None,:]/c)/len(off)
 beam=np.einsum('tfm,fm->tf',x,weights)
 win=signal.windows.hann(int(8*fs)+1,sym=True);win/=win.sum()
 y=signal.convolve(beam,win[:,None],mode='valid',method='direct');ix=np.arange(0,len(y),int(fs));y=y[ix];t=(ix+(len(win)-1)/2)/fs
 out={}
 for j in range(6):
  use=(t>=j*200)&(t<(j+1)*200);tt=t[use]
  for line in range(3):
   a=abs(np.fft.fft(y[use,line]*signal.windows.hann(sum(use)),16384))**2
   f=np.fft.fftfreq(16384,tt[1]-tt[0]);usef=np.flatnonzero(abs(f)<=.6);k=usef[np.argmax(a[usef])];delta=0.
   if 0<k<16383:
    v=np.log(np.maximum(a[k-1:k+2],1e-300));den=v[0]-2*v[1]+v[2]
    if den!=0:delta=.5*(v[0]-v[2])/den
   offset=float(f[k]+delta/(16384*(tt[1]-tt[0])))
   out[(line,j)]=(float(lo[line]+offset),offset,float(a[k]/max(float(np.median(a[usef])),1e-300)))
 return out

def main():
 if (OUT/'EXECUTION_STARTED.json').exists():raise RuntimeError('same frozen RX1Z execution already started; no rerun')
 import subprocess
 dump(OUT/'EXECUTION_STARTED.json',{'time_utc':stamp(),'design_sha':subprocess.check_output(['git','rev-parse','HEAD'],text=True).strip(),'new_task_sha256':read(OUT/'NEW_TASK_AND_PARENT.json')['new_task_sha256'],'once_only':True})
 meta=read(OUT/'PUBLIC_RECEIVER_METADATA.json');index=rows(OUT/'PUBLIC_RECORD_INDEX.csv');saved=rows(OUT/'PUBLIC_SAVED_CARRIER_ROWS.csv');ledger=[];comp=[];memory=0
 for rec in index:
  _,rss=guard();memory=max(memory,rss);rr=[r for r in saved if r['record_id']==rec['record_id']]
  recomputed=None
  if rec['raw_available']=='True':
   if sha(rec['path'])!=rec['sha256']:raise RuntimeError('raw bytes changed '+rec['record_id'])
   recomputed=frequency_rows(np.load(rec['path']),meta)
  for r in sorted(rr,key=lambda a:(int(a['line']),int(a['window']))):
   line=int(r['line']);offset=float(r['rx_center_offset_hz']);fhat=meta['digital_reference_hz'][line]+offset;detection=float(r['detection_ratio']);valid=np.isfinite(fhat) and np.isfinite(offset) and np.isfinite(detection) and detection>=6
   ledger.append(dict(record_id=r['record_id'],line=line,window=int(r['window']),f_hat_hz=fhat,saved_f_rx_hz=float(r['f_rx_hz']),center_offset_hz=offset,public_LO_hz=meta['digital_reference_hz'][line],carrier_valid=bool(valid),detection_ratio=detection,old_residual_nu_hz=float(r['nu_hz']),old_residual_q_eff_mps=float(r['q_eff_mps']),old_correlation_accepted=r['accepted'],old_all_peaks=r['all_peaks'],time_start_s=float(r['time_start_s']),time_end_s=float(r['time_end_s']),time_center_s=float(r['time_center_s'])))
   if recomputed is not None:
    f,of,ratio=recomputed[(line,int(r['window']))];comp.append(dict(record_id=r['record_id'],line=line,window=int(r['window']),saved_frequency_hz=fhat,recomputed_frequency_hz=f,frequency_difference_hz=f-fhat,offset_difference_hz=of-offset,detection_relative_difference=abs(ratio-detection)/max(detection,1e-300)))
 table(OUT/'RECEIVER_FREQUENCY_LEDGER.csv',ledger);table(OUT/'CARRIER_RECOMPUTATION.csv',comp)
 features=[]
 for rec in index:
  for line in range(3):
   a=sorted([r for r in ledger if r['record_id']==rec['record_id'] and r['line']==line],key=lambda r:r['window']);complete=len(a)==6 and all(r['carrier_valid'] for r in a);fbar=float(np.mean([r['f_hat_hz'] for r in a])) if complete else None
   for j in range(5):
    lo=next((r for r in a if r['window']==j),None);hi=next((r for r in a if r['window']==j+1),None);valid=lo is not None and hi is not None and lo['carrier_valid'] and hi['carrier_valid'];d=float(hi['f_hat_hz']-lo['f_hat_hz']) if valid else None
    features.append(dict(record_id=rec['record_id'],line=line,difference_j=j,group='TURN' if j==2 else 'SAME_LEG',valid=bool(valid),status='AVAILABLE' if valid else 'NOT_AVAILABLE',d_frequency_hz=d,Fbar_hz=fbar,apparent_radial_change_mps=-meta['c_ref_mps']*d/fbar if valid and complete else None,quantity='APPARENT_RADIAL_CHANGE_NOT_ABSOLUTE_Q',from_start_s=lo['time_start_s'] if lo else None,from_end_s=lo['time_end_s'] if lo else None,to_start_s=hi['time_start_s'] if hi else None,to_end_s=hi['time_end_s'] if hi else None))
 table(OUT/'INTERWINDOW_FEATURES.csv',features)
 controls=[]
 for snr in [20,5]:
  rec=next(r for r in index if int(r['SNR_dB'])==snr and r['raw_available']=='True');x=np.load(rec['path']);t=np.arange(len(x))/meta['sample_rate_hz'];delta=.031
  changed=x*np.exp(-2j*np.pi*delta*t[:,None,None]);changed_lo=np.array(meta['digital_reference_hz'])+delta
  # Restore the original common physical receiver coordinate BEFORE fixed beam/filter.
  restored=changed*np.exp(2j*np.pi*(changed_lo-np.array(meta['digital_reference_hz']))[None,:,None]*t[:,None,None]);a=frequency_rows(restored,meta)
  rr=[r for r in ledger if r['record_id']==rec['record_id']];e=max(abs(a[(r['line'],r['window'])][0]-r['f_hat_hz']) for r in rr)
  controls.append(dict(record_id=rec['record_id'],SNR_dB=snr,delta_LO_hz=delta,max_frequency_difference_hz=e,relative_sample_roundtrip=float(np.linalg.norm(restored-x)/np.linalg.norm(x)),limit_hz=1e-8,PASS=e<=1e-8,scope='same physical reference restored before unchanged beam and band; not source identification'))
 table(OUT/'DIGITAL_REFERENCE_CONTROLS.csv',controls)
 elapsed,rss=guard();memory=max(memory,rss)
 names=['RECEIVER_FREQUENCY_LEDGER.csv','INTERWINDOW_FEATURES.csv','CARRIER_RECOMPUTATION.csv','DIGITAL_REFERENCE_CONTROLS.csv']
 dump(OUT/'FEATURE_OUTPUT_MANIFEST.json',{'time_utc':stamp(),'feature_sha256':{n:sha(OUT/n) for n in names},'receiver_process_truth_reads':0,'receiver_process_source_reference_reads':0,'records':len(index),'line_windows':len(ledger),'frequency_differences':len(features),'raw_records_recomputed':len({r['record_id'] for r in comp}),'missing_raw_ids':[r['record_id'] for r in index if r['raw_available']!='True'],'max_frequency_recompute_difference_hz':max(abs(r['frequency_difference_hz']) for r in comp) if comp else None,'max_f_rx_vs_LO_plus_offset_hz':max(abs(r['f_hat_hz']-r['saved_f_rx_hz']) for r in ledger),'frequency_recompute_limit_hz':1e-8,'digital_reference_PASS':all(r['PASS'] for r in controls),'elapsed_s':elapsed,'sampled_RSS_max_bytes':memory})
 print('RECEIVER FROZEN',len(ledger),len(features),'max Hz error',max(abs(r['frequency_difference_hz']) for r in comp),flush=True)
if __name__=='__main__':main()
