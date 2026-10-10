"""Private static generator. Receiver is a separate restricted-input process."""
from hla_sp1_core import *

def main(kind):
 dest=LOCAL/kind;dest.mkdir(exist_ok=True)
 if (dest/'PRIVATE_MANIFEST.json').exists():raise RuntimeError('no regeneration')
 panel=rows(OUT/'CALIBRATION_PANEL.csv') if kind=='calibration' else rows(ROOT/'results/R4_SINGLE_HLA_SOURCE_PROFILED_HORIZONTAL_PILOT/OFFGRID_TRUTH_EVALUATION_ONLY.csv');states=np.array([[float(r['r_km']),float(r['theta_deg']),float(r['z_label'])] for r in panel]);provider=Provider(160001);h=provider.field(states,exact=True);pref=read(OUT/'PREFLIGHT.json')['reference_power'];root=np.random.SeedSequence(2026101018 if kind=='calibration' else 2026101019);streams=root.spawn(len(panel)*16);n=np.arange(2048);tone=np.exp(2j*np.pi*FREQ[:,None]*n[None,:]/1024);public=[];private=[]
 for g in range(len(panel)):
  for rep in range(16):
   guard();seeds=streams[g*16+rep].spawn(3);rs,rn,rg=[np.random.default_rng(s) for s in seeds];s=10**(rs.uniform(-3,3,(8,3))/20)*np.exp(1j*rs.uniform(-np.pi,np.pi,(8,3)));noise=rn.normal(size=(8,2048,8));gain=10**(rg.uniform(-1,1,8)/20)*np.exp(1j*np.radians(rg.uniform(-30,30,8)));privpath=dest/f'innov_g{g:02d}_r{rep:02d}.npz';np.savez(privpath,source=s,unit_real_noise=noise,gain=gain,unit_source_H=h[g],seed_spawn_keys=np.array([x.spawn_key for x in seeds],dtype=int))
   for c in ['C0','C1']:
    for snr in [20,5]:
     record=f'{kind.upper()}_{len(public):04d}';nu=pref/10**(snr/10);y=s[:,:,None]*h[g][None]*(gain[None,None,:] if c=='C1' else 1);x=2*np.real(np.einsum('kfm,fn->knm',y,tone))+np.sqrt(2048*nu)*noise;path=dest/(record+'.npy');np.save(path,x);public.append(dict(record_id=record,path=str(path),sha256=sha(path),C=c,reference_SNR_dB=snr,fs_hz=1024,N=2048,K=8));actual=np.mean(abs(y)**2,axis=(0,2))/nu;private.append(dict(record_id=record,geometry=g,replicate=rep,private_innovations=str(privpath),private_sha256=sha(privpath),truth_r_km=float(panel[g]['r_km']),truth_theta_deg=float(panel[g]['theta_deg']),truth_z_m=float(panel[g]['z_label']),actual_SNR_by_line_dB=(10*np.log10(actual)).tolist()))
  print('GENERATED',kind,'geometry',g+1,'/',len(panel),flush=True)
 table(dest/'PUBLIC_RECORD_INDEX.csv',public);dump(dest/'PRIVATE_MANIFEST.json',{'records':private,'source_seed':root.entropy,'permission':'generator and retrospective audit only, forbidden receiver and test estimator','source_noise_gain_paired_across_C_SNR':True});dump(dest/'PUBLIC_METADATA.json',{'fs_hz':1024,'N':2048,'K':8,'analysis_frequencies_hz':FREQ.tolist(),'array_x_m':OFF.tolist(),'array_z_m':200,'scope':'BIN_ALIGNED_STATIC_CONTROL','forbidden':'source,H,gain,truth,REF'})
 dump(OUT/(kind.upper()+'_RECORD_MANIFEST.json'),{'public_index':str(dest/'PUBLIC_RECORD_INDEX.csv'),'public_index_sha256':sha(dest/'PUBLIC_RECORD_INDEX.csv'),'private_manifest':str(dest/'PRIVATE_MANIFEST.json'),'private_manifest_sha256':sha(dest/'PRIVATE_MANIFEST.json'),'records':len(public),'local_raw_bytes':sum(Path(r['path']).stat().st_size for r in public),'time_utc':stamp()})
if __name__=='__main__':
 import sys
 main(sys.argv[1])
