"""Generator only: private states never cross into receiver files."""
from common import *
def make_cal(rng,g,beta):
 ang=np.array([30.,-30.,60.,-60.]);x=np.arange(-7,8,2);k=2*np.pi*200/1500
 a=g[None,:]*np.exp(1j*k*x[None,:]*np.cos(np.deg2rad(ang[:,None])-beta))
 phases=np.exp(1j*rng.uniform(-np.pi,np.pi,4))
 return (100*a[:,:,None]*phases[:,None,None]+cn(rng,(4,8,2048))).astype(np.complex64),ang

def waveform(rng,t,target,node,beta,g,snr,stress='SINGLE'):
 n=128;ts=np.arange(n)/128.;x=np.arange(-7,8,2);epoch=len(t)
 trueang=np.arctan2(target[:,1]-node[:,1],target[:,0]-node[:,0]);heading=axis(t)
 offset=rng.uniform(-.2,.2,epoch);ph=rng.uniform(-np.pi,np.pi,epoch)
 samples=cn(rng,(epoch,8,n));bg=cn(rng,(epoch,8,n))
 present=rng.random(epoch)>.05;false=rng.random(epoch)<.10
 if stress!='SINGLE':present[:]=True;false[:]=False
 if stress=='NOISE_ONLY':present[:]=False
 amp=np.sqrt(10**(snr/10));actual_g=np.tile(g,(epoch,1))
 if stress=='BAD_CHANNEL':actual_g[:,3]=0
 if stress=='PHASE_DRIFT':actual_g*=np.exp(1j*np.deg2rad(20)*np.arange(epoch)[:,None]/max(1,epoch-1)*np.linspace(-1,1,8))
 if stress=='WEAK':amp*=.01
 yy=.25*np.sin(np.arange(8)) if stress=='BENT_ARRAY' else np.zeros(8)
 xx=x*np.cos(np.deg2rad(5)) if stress=='BENT_ARRAY' else x
 for k in range(epoch):
  f=200+offset[k];a=actual_g[k]*np.exp(1j*2*np.pi*f/1500*(xx*np.cos(trueang[k]-heading[k]-beta)+yy*np.sin(trueang[k]-heading[k]-beta)))
  tone=np.exp(1j*(2*np.pi*offset[k]*ts+ph[k]));samples[k]+=present[k]*amp*a[:,None]*tone
  if stress=='COHERENT_ECHO':
   b=actual_g[k]*np.exp(1j*2*np.pi*f/1500*x*np.cos(trueang[k]-heading[k]-beta+np.deg2rad(12)))
   samples[k]+=.8*amp*b[:,None]*tone*np.exp(1j*.9)
  if stress=='CLOCK':
   samples[k]*=np.exp(1j*2*np.pi*200*np.linspace(-.0003,.0003,8))[:,None]
  if false[k]:
   fa=rng.uniform(-np.pi,np.pi);b=g*np.exp(1j*2*np.pi*203/1500*x*np.cos(fa));samples[k]+=amp*b[:,None]*np.exp(1j*(2*np.pi*3*ts+rng.uniform(-np.pi,np.pi)))
 return samples.astype(np.complex64),bg.astype(np.complex64),dict(angle=trueang,present=present,false=false,frequency=200+offset)

def main():
 spec=read('CAP06_FROZEN_SPEC.json');panel=read('P0_FROZEN_PANEL.json')
 public=RAW/'p2/receiver_public';private=RAW/'p2/truth_private';public.mkdir(parents=True);private.mkdir(parents=True)
 index=[];records=[]
 for entry in panel['states']:
  reps=spec['repetitions'][entry['split']]
  for rep in range(reps):
   side=1 if rep%2==0 else -1;s=np.array(entry['state'],float);s[1]*=side;s[4]*=side
   epoch=np.array(spec['times_s']);nom=platform(epoch)
   sid=entry['seed_base']+rep*1000;common=np.random.default_rng(sid+11).normal(size=(len(epoch),2))*20
   beta_common=np.random.default_rng(sid+12).uniform(-.05,.05)*np.pi/180
   for nodeid,off in [('MAIN',[0,0]),('B0',[0,side*5000]),('B3',[3000,side*15000])]:
    ni=['MAIN','B0','B3'].index(nodeid);seed=sid+100+ni*100
    rg=np.random.default_rng(seed);b=beta_common+rg.uniform(-.025,.025)*np.pi/180
    g=10**(rg.uniform(-1,1,8)/20)*np.exp(1j*np.deg2rad(rg.uniform(-20,20,8)))
    cal,angles=make_cal(rg,g,b);clock=rg.uniform(-.25,.25);actual_t=epoch+clock
    truepos=platform(actual_t)+off;reported=nom+off+common+rg.normal(size=nom.shape)*15
    target=s[0]*np.array([np.cos(s[1]),np.sin(s[1])])+actual_t[:,None]*s[3]*np.array([np.cos(s[4]),np.sin(s[4])])
    for snr in spec['snr_db']:
     # Same standardized innovations for quality levels, separate from calibration and navigation.
     y,bg,truth=waveform(np.random.default_rng(seed+17),actual_t,target,truepos,b,g,snr)
     token=f"{entry['id']}_{rep:03}_{nodeid}_{snr}"
     np.savez_compressed(public/(token+'.npz'),samples=y,background=bg,calibration=cal,calibration_angles_deg=angles,timestamps=epoch,reported_positions=reported,nominal_heading=axis(epoch),nav_sigma_axis_m=np.array(25.),fs_hz=np.array(128.),LO_hz=np.array(200.))
     np.savez_compressed(private/(token+'.npz'),state=s,true_positions=truepos,true_times=actual_t,beta=np.array(b),gain=g,**truth)
     records.append(dict(token=token,split=entry['split'],state_id=entry['id'],repeat=rep,node=nodeid,snr=snr,stress='MAIN_TRACK',seed=seed+17))
   if rep%8==0:print(entry['id'],rep,flush=True)
 # Predetermined 8 stress families x 8 signed directions x 8 innovations x 3 quality levels.
 angles_test=[-60,-30,-5,-1,1,5,30,60]
 for ci,stress in enumerate(spec['stress_cases']):
  for ai,angle in enumerate(angles_test):
   for rep in range(8):
    seed=2026101950+ci*10000+ai*100+rep
    rg=np.random.default_rng(seed);g=10**(rg.uniform(-1,1,8)/20)*np.exp(1j*np.deg2rad(rg.uniform(-20,20,8)));b=rg.uniform(-.075,.075)*np.pi/180
    cal,ca=make_cal(rg,g,b);t=np.array([0.,300.,600.,900.,1200.]);pos=platform(t);beta=np.deg2rad(angle)+axis(t)
    target=pos+50000*np.column_stack((np.cos(beta),np.sin(beta)))
    for snr in spec['snr_db']:
     y,bg,truth=waveform(np.random.default_rng(seed+1),t,target,pos,b,g,snr,stress)
     token=f'STRESS_{ci}_{ai}_{rep}_{snr}'
     np.savez_compressed(public/(token+'.npz'),samples=y,background=bg,calibration=cal,calibration_angles_deg=ca,timestamps=t,reported_positions=pos,nominal_heading=axis(t),nav_sigma_axis_m=np.array(0.),fs_hz=np.array(128.),LO_hz=np.array(200.))
     np.savez_compressed(private/(token+'.npz'),state=np.zeros(5),true_positions=pos,true_times=t,beta=np.array(b),gain=g,**truth)
     records.append(dict(token=token,split='STRESS_CONTROL',state_id=f'A{ai}',repeat=rep,node='CONTROL',snr=snr,stress=stress,seed=seed+1))
 tab('P2_RECORD_INDEX.csv',records)
 for folder in [public,private]:
  for p in sorted(folder.glob('*.npz')):index.append(dict(path=str(p),bytes=p.stat().st_size,sha256=sha(p),domain=folder.name))
 tab('P2_RAW_DATA_INDEX.csv',index)
 write('P2_GENERATION_SUMMARY.json',dict(records=len(records),unique_trajectories=sum(r['stress']=='MAIN_TRACK' for r in records),stress_records=sum(r['stress']!='MAIN_TRACK' for r in records),model='PLANE_WAVE_CONDITIONAL_NOT_CZ',solver_calls=0))
if __name__=='__main__':main()
