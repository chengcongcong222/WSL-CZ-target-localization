"""Calibration-only truth labels; freeze thresholds before any test inference."""
from hla_sp1_core import *
def main():
 if (OUT/'THRESHOLD_FREEZE.json').exists():raise RuntimeError('no recalibration')
 dest=LOCAL/'calibration';m=read(OUT/'CALIBRATION_DFT_FREEZE.json');assert sha(m['file'])==m['sha256'];y=np.load(m['file'],mmap_mode='r');private=read(dest/'PRIVATE_MANIFEST.json')['records'];pub=rows(dest/'PUBLIC_RECORD_INDEX.csv');provider=Provider(160001);dist=[];taus=[]
 for step in [.10,.05]:
  rs,ths,cat=catalogue(step);states=np.array([[rs[nearest(rs,r['truth_r_km'])],ths[nearest(ths,abs(r['truth_theta_deg']))],r['truth_z_m']] for r in private]);h=provider.field(states,exact=True)
  for method in range(3):
   ds=[]
   for i in range(len(private)):
    d=1-float(score(h[i],y[i],method));ds.append(d);dist.append(dict(record_id=pub[i]['record_id'],geometry=private[i]['geometry'],replicate=private[i]['replicate'],C=pub[i]['C'],reference_SNR_dB=pub[i]['reference_SNR_dB'],grid_step_km=step,method='P'+str(method),D_nearest_joint_label=d))
   for c in ['C0','C1']:
    for snr in [20,5]:
     values=sorted(ds[i] for i,r in enumerate(pub) if r['C']==c and int(r['reference_SNR_dB'])==snr);assert len(values)==256;tau=values[int(np.ceil(.95*256))-1];taus.append(dict(C=c,reference_SNR_dB=snr,grid_step_km=step,method='P'+str(method),count=256,nearest_rank=int(np.ceil(.95*256)),tau=tau))
 table(OUT/'CALIBRATION_DISTRIBUTION.csv',dist);table(OUT/'CALIBRATION_THRESHOLDS.csv',taus);dump(OUT/'THRESHOLD_FREEZE.json',{'time_utc':stamp(),'calibration_DFT_sha256':m['sha256'],'threshold_sha256':sha(OUT/'CALIBRATION_THRESHOLDS.csv'),'distribution_sha256':sha(OUT/'CALIBRATION_DISTRIBUTION.csv'),'test_inference_not_started':not (OUT/'TEST_INFERENCE_STARTED.json').exists(),'rule':'nearest rank95%,256 per method/C/SNR/grid; not theoretical95% confidence'});print('THRESHOLDS FROZEN',len(taus),flush=True)
if __name__=='__main__':main()
