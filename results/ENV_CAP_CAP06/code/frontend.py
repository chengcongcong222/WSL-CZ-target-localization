"""Reads receiver samples and legitimate calibration only. No truth/generator imports."""
import numpy as np
from scipy.special import gammaincinv
X=np.arange(-7,8,2,dtype=float);U=np.linspace(-1,1,1001)
def calibration(raw,known_angles_deg):
 # Four independent known-angle calibration transmissions, NOT the target.
 v=raw.mean(axis=-1);rat=v/v[:,0,None];ph=np.unwrap(np.angle(rat),axis=1)
 angles=np.deg2rad(known_angles_deg);dx=X-X[0];k=2*np.pi*200/1500
 # +/- angle pairs remove each channel's fixed complex response phase.
 dif=[];coef=[]
 for pos,neg in [(0,1),(2,3)]:
  dif.extend(np.angle(rat[pos,1:]*rat[neg,1:].conj()))
  coef.extend(2*k*dx[1:]*np.sin(angles[pos]))
 coef=np.array(coef);b=float(np.arcsin(np.clip(coef@dif/(coef@coef),-1,1)))
 pred=np.exp(1j*k*dx[None,:]*np.cos(angles[:,None]-b))
 gr=rat/pred;g=gr.mean(axis=0)
 residual=float(np.sqrt(np.mean(abs(gr-g)**2))/np.sqrt(np.mean(abs(g)**2)))
 return b,g,residual

def extract(samples,background,cal,cal_angles,calibrated=True):
 beta,g,cerr=calibration(cal,cal_angles) if calibrated else (0.,np.ones(8,complex),None)
 n=samples.shape[-1];fs=128.;fft_n=1024;freq=np.fft.fftfreq(fft_n,1/fs)
 sel=np.flatnonzero(abs(freq)<=8)
 Z=np.fft.fft(samples,n=fft_n,axis=-1)[:,sel]/n;ff=freq[sel]
 energy=np.sum(abs(Z)**2,axis=0);nv=float(np.mean(abs(background)**2))
 # Upper variance allowance for finite background. No target noise truth used.
 nbg=background.size;nv_hi=nv*nbg/gammaincinv(nbg,.0001)
 det=nv_hi/n*gammaincinv(8,1-1e-3/len(sel))
 order=np.argsort(-energy);chosen=[]
 for j in order:
  if energy[j]<=det:break
  if all(abs(ff[j]-ff[q])>=1.5 for q in chosen):chosen.append(int(j))
  if len(chosen)==3:break
 products=[]
 for j in chosen:
  f=200+ff[j];k=2*np.pi*f/1500
  steering=g[:,None]*np.exp(1j*k*X[:,None]*U[None,:]);norm=np.sum(abs(steering)**2,axis=0)
  y=Z[:,j];res=np.maximum(0,energy[j]-abs(steering.conj().T@y)**2/norm)
  threshold=nv_hi/n*gammaincinv(7,1-.0001)+(.02**2)*energy[j]
  # |R''| <= 8 k^2 max|x|^2 ||y||^2; retain all uncertain u-cells.
  interpolation_allowance=(8*k*k*49*energy[j])*(U[1]-U[0])**2/8
  mask=res<=threshold+interpolation_allowance
  if not mask.any():continue
  best=int(np.argmin(res));lo=max(-1,float(U[mask].min()-.002));hi=min(1,float(U[mask].max()+.002))
  products.append(dict(f_rx_hz=float(f),u_best=float(U[best]),u_lo=lo,u_hi=hi,axis_angle_deg=float(np.rad2deg(np.arccos(U[best]))),score=float(res[best]),threshold=float(threshold),power=float(energy[j]),mirror_angles_deg=[float(np.rad2deg(beta)+np.rad2deg(np.arccos(U[best]))),float(np.rad2deg(beta)-np.rad2deg(np.arccos(U[best])))],association='UNASSIGNED'))
 return dict(beta_est_rad=beta,calibration_residual=cerr,g_est_real=g.real.tolist(),g_est_imag=g.imag.tolist(),products=products,detected=bool(products),noise_estimate=nv,noise_upper=nv_hi,detection_threshold=float(det),physical_model='PLANE_WAVE_CONTROL',PHYSICAL_CZ_LINK='UNVERIFIED')
