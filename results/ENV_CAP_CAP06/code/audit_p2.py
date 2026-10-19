"""Independent raw-sample DFT, calibration geometry and saved-output accounting."""
from pathlib import Path
import json,csv,hashlib,numpy as np
from scipy.special import gammaincinv
O=Path(__file__).resolve().parents[1];W=Path(r'D:\ProjectStorage\WSL-CZ\ENV_CAP_CAP06/raw/p2')
def main():
 checks=[]
 def ck(n,b,e=0):checks.append(dict(check_id=n,PASS=bool(b),error=float(e)))
 records=list(csv.DictReader((O/'P2_RECORD_INDEX.csv').open(encoding='utf-8')))
 # Fixed audit panel: first record per main node/quality/split plus first per stress/quality.
 selected={}
 for r in records:selected.setdefault((r['split'],r['node'],r['snr'],r['stress']),r)
 for r in selected.values():
  token=r['token'];a=np.load(W/'receiver_public'/(token+'.npz'));z=json.loads((W/'products'/(token+'_CAL.json')).read_text())
  v=a['calibration'].astype(complex).sum(-1)/a['calibration'].shape[-1];ratio=v/v[:,0,None];dx=np.arange(0,16,2);k=2*np.pi*200/1500
  numerator=denom=0.
  for j,l in [(0,1),(2,3)]:
   for m in range(1,8):
    b=2*k*dx[m]*np.sin(np.deg2rad(a['calibration_angles_deg'][j]));y=np.angle(ratio[j,m]/ratio[l,m]);numerator+=b*y;denom+=b*b
  beta=np.arcsin(np.clip(numerator/denom,-1,1));ck(token+':calibration_beta',abs(beta-z['epochs'][0]['beta_est_rad'])<1e-6,abs(beta-z['epochs'][0]['beta_est_rad']))
  # Direct complex sums at every exported candidate frequency, not production FFT.
  for ep in [0,len(z['epochs'])-1]:
   e=z['epochs'][ep];g=np.array(e['g_est_real'])+1j*np.array(e['g_est_imag']);data=a['samples'][ep].astype(complex)
   for p in e['products']:
    f=p['f_rx_hz'];v=(data*np.exp(-2j*np.pi*(f-200)*np.arange(128)/128)).sum(-1)/128
    u=np.linspace(-1,1,1001);A=g[:,None]*np.exp(2j*np.pi*f/1500*np.arange(-7,8,2)[:,None]*u);score=(abs(v)**2).sum()-abs(np.einsum('ij,i->j',A.conj(),v))**2/(abs(g)**2).sum()
    j=int(score.argmin());ck(token+f':DFT_beam_{ep}_{f}',abs(u[j]-p['u_best'])<1e-12,abs(u[j]-p['u_best']))
    ck(token+f':residual_{ep}_{f}',abs(max(0,score[j])-p['score'])<1e-6,abs(max(0,score[j])-p['score']))
 # Verify every exported row's private-evaluation linkage without calling generator/front end.
 cases=list(csv.DictReader((O/'P2_RECEIVER_BEARING_CASES.csv').open(encoding='utf-8')));by={}
 for r in cases:by.setdefault((r['token'],r['mode']),[]).append(r)
 for (token,mode),rr in by.items():
  truth=np.load(W/'truth_private'/(token+'.npz'));z=json.loads((W/'products'/(token+'_'+mode+'.json')).read_text());ck(token+mode+':denominator',len(rr)==len(truth['present'])==len(z['epochs']))
  err=0;ok=True
  for row,e in zip(rr,z['epochs']):
   j=int(row['epoch']);heading=np.deg2rad(15) if truth['true_times'][j]>600 else 0;target=np.degrees(np.arccos(np.clip(np.cos(truth['angle'][j]-heading-e['beta_est_rad']),-1,1)));match=[p for p in e['products'] if abs(p['f_rx_hz']-truth['frequency'][j])<=.75]
   valid=bool(match) and bool(truth['present'][j]);ok &= valid==(row['target_component_extracted']=='True')
   if valid:err=max(err,abs(min(abs(target-p['axis_angle_deg']) for p in match)-float(row['axis_error_deg'])))
  ck(token+mode+':errors_and_association_evaluation_only',ok and err<1e-7,err)
 # Summary counts from immutable row table, no reuse of scientific summarizer.
 for s in csv.DictReader((O/'P2_EXTRACTION_SUMMARY.csv').open(encoding='utf-8')):
  rr=[r for r in cases if all(r[k]==s[k] for k in ['split','stress','snr','mode'])];n=sum(r['present']=='True' for r in rr);k=sum(r['target_component_extracted']=='True' for r in rr);cov=sum(r['target_u_covered']=='True' for r in rr)
  ck('summary:'+str((s['split'],s['stress'],s['snr'],s['mode'])),len(rr)==int(s['epochs']) and n==int(s['present_epochs']) and k==int(s['extracted']) and cov==int(s['covered']))
 # All raw identities; independent streaming SHA.
 for name in ['P2_RAW_DATA_INDEX.csv','P2_PRODUCT_INDEX.csv']:
  for r in csv.DictReader((O/name).open(encoding='utf-8')):
   h=hashlib.sha256()
   with open(r['path'],'rb') as f:
    for chunk in iter(lambda:f.read(1048576),b''):h.update(chunk)
   ck('sha:'+Path(r['path']).name,h.hexdigest()==r['sha256'])
 mats=np.load(W.parent/'P2_LOCAL_INFORMATION_MATRICES.npz');panel=json.loads((O/'P0_FROZEN_PANEL.json').read_text())['states'];t=np.arange(0,1201,40.);heading=np.where(t>600,np.pi/12,0);u=np.minimum(t,600);v=np.maximum(t-600,0);platform=np.column_stack((2*u+2*v*np.cos(np.pi/12),2*v*np.sin(np.pi/12)));scale=np.array([1000,1000,1,1]);weight=2*128*100*(2*np.pi*200/1500)**2*168
 for en in panel:
  for side in [-1,1]:
   st=np.array(en['state']);st[1]*=side;st[4]*=side;x=np.r_[st[0]*np.array([np.cos(st[1]),np.sin(st[1])]),st[3]*np.array([np.cos(st[4]),np.sin(st[4])])];ff={};free={}
   for node,off in [('MAIN',[0,0]),('B0',[0,side*5000]),('B3',[3000,side*15000])]:
    def response(xx):
     d=xx[:2]+t[:,None]*xx[2:]-platform-off;return np.cos(np.arctan2(d[:,1],d[:,0])-heading)
    J=np.column_stack([(response(x+np.eye(4)[j]*scale[j]*.0002)-response(x-np.eye(4)[j]*scale[j]*.0002))/.0004 for j in range(4)]);d=x[:2]+t[:,None]*x[2:]-platform-off;B=np.sin(np.arctan2(d[:,1],d[:,0])-heading)[:,None];R=J-B@np.linalg.lstsq(B,J,rcond=1e-12)[0];ff[node]=weight*J.T@J;free[node]=weight*R.T@R
   for res,nodes in [('S0',['MAIN']),('S2_B0',['MAIN','B0']),('S2_B3',['MAIN','B3'])]:
    for mode,dd in [('CAL_ORIENTATION_IDEAL',ff),('FREE_NODE_ORIENTATION',free)]:
     F=sum(dd[n] for n in nodes);saved=mats[f"{en['id']}_{side}_{res}_{mode}"];err=np.linalg.norm(F-saved)/max(np.linalg.norm(saved),1e-12);ck('information:'+en['id']+str(side)+res+mode,err<1e-5,err)
 gate=json.loads((O/'P2_GATE.json').read_text());ck('physical_scope',gate['PHYSICAL_CZ_LINK']=='UNVERIFIED' and gate['z']=='NOT_EVALUATED')
 out=dict(PASS=all(c['PASS'] for c in checks),checks=len(checks),FAIL=sum(not c['PASS'] for c in checks),records=checks,scope='independent selected DFT/control reconstruction plus all row/identity accounting; not CZ propagation certificate')
 (O/'P2_INDEPENDENT_COLD_AUDIT.json').write_text(json.dumps(out,indent=2)+'\n');print(out['PASS'],out['checks'])
if __name__=='__main__':main()
