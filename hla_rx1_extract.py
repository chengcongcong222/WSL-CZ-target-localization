"""Receive-only entry. It cannot open generator/provider/truth files.
REF permissions are an explicit separate optional phase reference argument.
"""
import json
import numpy as np
from scipy import signal,fft
from hla_rx1_common import read,OUT,table

def combine(x,meta):
 # Fixed global +x endfire steering: hardware geometry/receiver LO, no bearing truth.
 offsets=np.array(meta['offsets_m'])
 w=np.exp(-2j*np.pi*np.array(meta['digital_reference_hz'])[:,None]*offsets[None,:]/meta['c_ref_mps'])/len(offsets)
 return np.einsum('tfm,fm->tf',x,w)

def stft_line(x,fs):
 # Global t=0 complex carrier convention; each local coefficient is at LO,
 # not reset to a segment-local source carrier. Eight-second Hann, one-second hop.
 win=signal.windows.hann(int(8*fs)+1,sym=True);win/=win.sum()
 y=signal.convolve(x,win[:,None],mode='valid',method='direct')
 ix=np.arange(0,len(y),int(fs));return y[ix],(ix+(len(win)-1)/2)/fs

def peak_frequency(y,t):
 # Receiver-only whole analysis window line-center estimate. All phases retained.
 fs=1/(t[1]-t[0]);n=16384
 a=abs(np.fft.fft(y*signal.windows.hann(len(y)),n))**2
 f=np.fft.fftfreq(n,1/fs);use=np.flatnonzero(abs(f)<=.6)
 k=use[np.argmax(a[use])];delta=0.
 if 0<k<n-1:
  v=np.log(np.maximum(a[k-1:k+2],1e-300));den=v[0]-2*v[1]+v[2]
  if den!=0:delta=.5*(v[0]-v[2])/den
 center=f[k]+delta*fs/n
 noise=float(np.median(a[use]));ratio=float(a[k]/max(noise,1e-300))
 return float(center),ratio

def correlation(y,dt,lag):
 n=len(y);L=min(n-1,int(lag/dt))
 a=signal.correlate(y,y,mode='full',method='fft')[n-1:n+L]
 a=a/max(float(a[0].real),1e-300)
 # Biased finite-record ACF: all available pairs; no simulated correlation.
 c=a.real
 # DCT-I after numerical zero extension of positive lags: no added samples.
 s=np.maximum(fft.dct(c,type=1,n=8193),0)
 f=np.fft.rfftfreq(16384,dt);valid=f<=.5
 f=f[valid];s=s[valid]
 peaks=signal.find_peaks(s)[0].tolist()
 if s[0]>=s[1]:peaks=[0]+peaks
 if s[-1]>=s[-2]:peaks+=[len(s)-1]
 peaks=sorted(peaks,key=lambda k:(-s[k],f[k]))
 top=float(max(s))
 kept=[{'nu_hz':float(f[k]),'relative_power':float(s[k]/max(top,1e-300))} for k in peaks if s[k]>=.25*top]
 return c,f,s,kept

def extract(x,meta,reference=None,keep_spectra=False):
 y,t=stft_line(combine(x,meta),meta['sample_rate_hz']);rows=[];spectra={}
 for j in range(6):
  use=(t>=j*200)&(t<(j+1)*200);tt=t[use];yy=y[use]
  for line in range(3):
   offset,detect=peak_frequency(yy[:,line],tt)
   frx=meta['digital_reference_hz'][line]+offset
   for method in ['RX','REF'] if reference is not None else ['RX']:
    if method=='RX':phase=2*np.pi*offset*tt
    else:
     phase=2*np.pi*(reference['source_frequency_hz'][line]-meta['digital_reference_hz'][line])*tt+np.pi*reference['drift_hz_per_s'][line]*tt**2
    residual=yy[:,line]*np.exp(-1j*phase)
    for lag in [200,100]:
     c,f,s,p=correlation(residual,1.,lag)
     nu=p[0]['nu_hz'] if p else None
     ok=detect>=6 and nu is not None and np.isfinite(residual).all()
     q=meta['c_ref_mps']*nu/frx if nu is not None else None
     rows.append(dict(window=j,line=line,method=method,lag_s=lag,f_rx_hz=frx,rx_center_offset_hz=offset,detection_ratio=detect,
      nu_hz=nu,q_eff_mps=q,accepted=bool(ok),reject_reason='' if ok else 'LINE_OR_CORRELATION_NOT_DETECTED',
      all_peaks=json.dumps(p,separators=(',',':')),direction='UNRESOLVED_PLUS_MINUS',sign_branches='[-1,1]',
      time_center_s=float((tt[0]+tt[-1])/2),time_start_s=float(tt[0]),time_end_s=float(tt[-1]),
      weighting='biased ACF: pair endpoints; lag spectrum nonlinear, no certified scalar radial weight',
      reference_permissions='receiver-time source phase only; no delay/mode/state' if method=='REF' else 'receiver samples and public metadata only'))
     if keep_spectra:spectra[f'w{j}_f{line}_{method}_L{lag}']=(tt,residual,c,f,s)
 return rows,spectra

def main():
 import argparse
 p=argparse.ArgumentParser();p.add_argument('samples');p.add_argument('metadata');p.add_argument('output');a=p.parse_args()
 meta=read(a.metadata);x=np.load(a.samples);rows,_=extract(x,meta);table(a.output,rows)
if __name__=='__main__':main()
