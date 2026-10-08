"""Independent analytic horizontal chain-rule and nonbearing rank audit. No new states/solver."""
import numpy as np
from scipy import linalg
import r4_e2_nine_pair_pilot as p
def main():
 out=[]
 for sc in p.read(p.O/'DESIGN_FREEZE.json')['scenes']:
  x=p.state(sc)
  for grid in [80001,160001]:
   a=np.load(p.O/f"INPUT_{sc['scene_id']}_n{grid}.npz");P=a['pressure'];J=a['J'];B=a['bearing_tangent']
   ms=[p.old.parse_mod(p.paired.path_for(ff,grid)) for ff in p.RAW]
   rr,nodes=p.old.geometry(sc['state'],sc['mirror'])
   target=p.cart(x)[:2]+p.old.TIMES[:,None]*p.cart(x)[2:]
   offsets=np.stack([np.arange(-7.,8.,2.),np.zeros(8)],axis=-1)
   delta=target[:,None,None,:]-nodes[:,:,None,:]-offsets[None,None,:,:]
   er=delta/rr[:,:,:,None]
   grad=[]
   for m in ms:
    zz,phi,k=m;w=phi[list(zz).index(200.)]**2;r=rr.ravel()
    v=w[:,None]*np.sqrt(2*np.pi/(k[:,None]*r))*np.exp(-1j*k[:,None]*r-1j*np.pi/4)
    dr=np.sum(v*(-1j*k[:,None]-1/(2*r)),axis=0).reshape(rr.shape)
    grad.append(dr[:,:,:,None]*er)
   grad=np.stack(grad,axis=-2)/P[:,:,:,:,None]
   geo=np.zeros((3,2,4));geo[:,:,0]=[np.cos(x[1]),np.sin(x[1])]
   geo[:,:,1]=[-x[0]*np.sin(x[1]),x[0]*np.cos(x[1])]
   geo[:,0,2]=p.old.TIMES;geo[:,1,3]=p.old.TIMES
   exact=J.copy();exact[...,:4]=np.einsum('tnmfd,tdc->tnmfc',grad,geo)*p.SCALE[:4]
   v=target[:,None,:]-nodes;bv=np.stack([-v[:,:,1],v[:,:,0]],axis=-1)
   be=np.einsum('tnmfd,tnd->tnmf',grad,bv)
   for rep,nn in [('P1',1),('P2',2)]:
    y,g,b,_=p.system(P,exact,be,rep,nn,'RELATIVE',1);z=p.remove(y,np.column_stack([g,b]))
    yf,gf,bf,_=p.system(P,J,B,rep,nn,'RELATIVE',1);zf=p.remove(yf,np.column_stack([gf,bf]))
    s=linalg.svdvals(z);sf=linalg.svdvals(zf)
    me=p.metrics(z,.01);mf=p.metrics(zf,.01)
    out.append(dict(scene_id=sc['scene_id'],mesh=grid,representation=rep,analytic_horizontal_chain_rank=int(sum(s>s[0]*1e-10)),frozen_fd_rank=int(sum(sf>sf[0]*1e-10)),trace_relative_difference=float(abs(np.sum(z*z)-np.sum(zf*zf))/np.sum(z*z)),range_information_relative_difference=float(abs(me['range_information_depth_profiled']-mf['range_information_depth_profiled'])/max(me['range_information_depth_profiled'],mf['range_information_depth_profiled'],1e-30)),depth_information_relative_difference=float(abs(me['depth_information_horizontal_profiled']-mf['depth_information_horizontal_profiled'])/max(me['depth_information_horizontal_profiled'],mf['depth_information_horizontal_profiled'],1e-30)),weak_singular_ratio_exact=float(s[-1]/s[0]),weak_singular_ratio_frozen=float(sf[-1]/sf[0])))
 p.write(p.O/'ANALYTIC_CHAIN_RULE_AND_RANK_AUDIT.csv',out)
 for rep in ['P1','P2']:
  r=[x for x in out if x['representation']==rep]
  print(rep,'ranks',sorted(set((x['analytic_horizontal_chain_rank'],x['frozen_fd_rank']) for x in r)),'range max',max(x['range_information_relative_difference'] for x in r),'depth max',max(x['depth_information_relative_difference'] for x in r),'trace max',max(x['trace_relative_difference'] for x in r))
if __name__=='__main__':main()
