"""Receiver reads public raw waveform whitelist only, never generator truth/source/gain."""
from t2_tilt_g0_core import *
import sys

def main(panel):
 if (OUT/f'{panel.upper()}_RECEIVER_STARTED.json').exists():raise RuntimeError('receiver once guard')
 dump(OUT/f'{panel.upper()}_RECEIVER_STARTED.json',{'utc':utc()});start=time.monotonic();a=rows(OUT/f'{panel.upper()}_PUBLIC_RECORD_INDEX.csv');outdir=LOCAL/panel/'DFT';outdir.mkdir();idx=[]
 for r in a:
  guard(start);fp=Path(r['waveform_path']);assert sha(fp)==r['waveform_SHA256'];x=np.load(fp);assert x.shape==(K,int(r['channels']),N)
  y=dft(x);yp=outdir/f"DFT_{int(r['record_id']):05d}.npy";np.save(yp,y);idx.append({k:r[k] for k in ['record_id','base_id','resource','calibration','reference_SNR_dB']}|{'DFT_path':str(yp),'DFT_SHA256':sha(yp)})
  if len(idx)%1024==0:print('DFT',panel,len(idx),'/',len(a),flush=True)
 table(OUT/f'{panel.upper()}_DFT_INDEX.csv',idx)
 dump(OUT/('RECEIVER_DFT_FREEZE.json' if panel=='test' else 'CALIBRATION_DFT_FREEZE.json'),{'utc':utc(),'panel':panel,'index_SHA256':sha(OUT/f'{panel.upper()}_DFT_INDEX.csv'),'records':len(idx),'only_sampling_and_public_waveforms_accessed':True,'science_wall_s':time.monotonic()-start})
if __name__=='__main__':main(sys.argv[1])