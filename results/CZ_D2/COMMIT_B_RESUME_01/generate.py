"""Frozen SHA panel and archive generator. Private stage; never imported by worker."""
from common import *

def panel():
 verify_code();pc=CFG['panel'];old=[[float(r[k]) if r[k] else None for k in ['r_km','theta_deg','v_mps','psi_deg']] for r in csvread(V3)];accepted=[];attempts=[]
 def reason(h):
  for j,(a,s) in enumerate(zip(pc['grid_origin'],pc['grid_step'])):
   f=(h[j]-a)/s
   if abs(f-round(f))<=pc['grid_fraction_guard']:return 'GRID_'+str(j)
  for family,states in [('HISTORY',old),('NEW',accepted)]:
   for i,state in enumerate(states):
    w=pc['static_position_guard'] if state[2] is None else pc['old_and_new_four_state_neighborhood_half_width']
    for sign in (1,-1):
     if all((abs((h[j]-sign*state[j]+180)%360-180) if j in [1,3] else abs(h[j]-state[j]))<=w[j] for j in range(len(w))):return family+'_'+str(i)+'_'+str(sign)
  return 'ACCEPT'
 for i in range(12):
  for a in range(pc['max_attempts_per_stratum']):
   b=hashlib.sha256((pc['seed_text']+'/'+str(i)+'/'+str(a)).encode('ascii')).digest();u=[((int.from_bytes(b[j:j+8],'big')>>12)+.5)/2**52 for j in range(0,32,8)]
   h=[pc['lower'][0]+(i+u[0])*(pc['upper'][0]-pc['lower'][0])/12]+[pc['lower'][j]+u[j]*(pc['upper'][j]-pc['lower'][j]) for j in range(1,4)]
   why=reason(h);attempts.append(dict(case_index=i,attempt=a,digest=b.hex(),h=h,decision=why))
   if why=='ACCEPT':accepted.append(h);break
  else:raise RuntimeError('PANEL_ATTEMPT_LIMIT')
 dump(DATA/'private/PANEL_AND_REJECTIONS.json',dict(states=accepted,attempts=attempts,registry_sha256=sha(V3)),True)
 # Independent rebuild before any observations. Pure audit function has no generator import.
 from independent_math import rebuild_panel
 other=rebuild_panel(CFG,csvread(V3));assert other==dict(states=accepted,attempts=attempts), 'PANEL_INDEPENDENT_REBUILD'
 dump(O/'PANEL_SEAL.json',dict(path=str(DATA/'private/PANEL_AND_REJECTIONS.json'),sha256=sha(DATA/'private/PANEL_AND_REJECTIONS.json'),n=12,attempts=len(attempts),independent_rebuild=True,code_freeze_sha256=sha(O/'PRODUCTION_CODE_FREEZE.json')),True)
 csvwrite(O/'PANEL_EVALUATION_ONLY.csv',[dict(case_id=f'D{i+1:02d}',**dict(zip(['r_km','theta_deg','v_mps','psi_deg'],h))) for i,h in enumerate(accepted)])
 return accepted

def generate(states):
 gen=module('t3_w1_generate');manifest=[];priv=[]; sides=np.r_[np.full(250,-1),np.full(250,1)]
 for ci,h in enumerate(states):
  rng=np.random.Generator(np.random.PCG64(np.random.SeedSequence([CFG['random']['seed'],ci])));u=rng.uniform(size=(500,2));n=rng.standard_normal((500,121,6));pn=DATA/'private'/f'D{ci+1:02d}_innovations.npz';np.savez(pn,uniforms=u,normals=n,side=sides)
  innovhash=arrayhash(u,n,sides);priv.append(dict(case=ci,path=str(pn),sha256=sha(pn),array_sha256=innovhash))
  tr=dict(zip(['r0_km','theta0_deg','v_mps','psi_deg'],h))
  for q in ['Q1','Q2']:
   mainhash=None
   for b in ['B0','B3']:
    t,P,y,bias=gen.observations(tr,CFG['geometry'][b],CFG['quality'][q],n,u,sides);p=recordpath(q,b,ci);np.savez(p,times=t,reported=P,bearing=y)
    hh=[arrayhash(t,P[j, :,0],y[j,:,0]) for j in range(500)];rh=[arrayhash(t,P[j],y[j]) for j in range(500)]
    if mainhash is None:mainhash=hh
    else:assert mainhash==hh,'B0_B3_MAIN_BYTES'
    for j in range(500):manifest.append(dict(quality=q,geometry=b,case_id=f'D{ci+1:02d}',replicate=j,side=int(sides[j]),path=str(p),record_sha256=rh[j],main_sha256=hh[j],innovation_sha256=innovhash))
  print('GENERATED',ci+1,flush=True)
 csvwrite(O/'PUBLIC_RECORD_IDENTITIES.csv',manifest);dump(O/'PRIVATE_ARCHIVE_INDEX.json',priv,True)
 pub=[dict(path=str(p),sha256=sha(p),bytes=p.stat().st_size) for p in sorted((DATA/'public').glob('*.npz'))];dump(O/'PUBLIC_ARCHIVE_SEAL.json',dict(files=pub,records=len(manifest),B0_B3_main_identical=True,P0_P1_shared=True),True)
