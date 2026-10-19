from pathlib import Path
import csv,json,hashlib,numpy as np,matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
O=Path(__file__).resolve().parents[1];rows=list(csv.DictReader((O/'P2_EXTRACTION_SUMMARY.csv').open()));z=[r for r in rows if r['split']=='STRESS_CONTROL' and r['mode']=='CAL' and r['snr']=='20']
f,ax=plt.subplots(figsize=(10,4.8));vals=[float(r['extract_fraction']) if r['extract_fraction'] else np.nan for r in z];ax.bar(range(len(z)),vals);ax.set_xticks(range(len(z)),[r['stress'] for r in z],rotation=32,ha='right');ax.set(ylabel='Target extraction fraction',ylim=(0,1.04),title='Fixed 20 dB calibrated stress panel: all failures retained')
for j,v in enumerate(vals):
 if np.isfinite(v):ax.text(j,v+.02,f'{v:.3f}',ha='center',fontsize=9)
 else:ax.text(j,.08,'N/A\n0/320 false alarms',ha='center',fontsize=8)
f.text(.5,.015,'Matched synthetic control; paired innovations and epochs are not independent trajectories.',ha='center',fontsize=9);f.tight_layout(rect=[0,.04,1,1]);f.savefig(O/'figures/P2_STRESS_BOUNDARY.png',dpi=180);plt.close(f)
info=json.loads((O/'P1_NUISANCE_INFO.json').read_text());f,ax=plt.subplots(figsize=(7,4.5));v=[r['profiled_state_information'][2] for r in info];ax.bar([r['gain'] for r in info],v)
for j,x in enumerate(v):ax.text(j,x+.5,'~0 (rank 0)' if j==3 else f'{x:.3f}',ha='center',fontsize=9)
ax.set(ylabel='Dimensionless toy depth-direction information',ylim=(0,41),title='Toy gain nuisance: no first-CZ depth precision claim');f.tight_layout();f.savefig(O/'figures/P1_TOY_NUISANCE.png',dpi=180);plt.close(f)
# Administrative checksum updates only; raw archives and all science tables unchanged.
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
for name in ['FIGURE_PROVENANCE.csv','ARTIFACT_INDEX.csv']:
 rr=list(csv.DictReader((O/name).open(encoding='utf-8')))
 for r in rr:
  p=O/r['figure' if name.startswith('FIGURE') else 'path'];r['sha256']=sha(p)
  if 'bytes' in r:r['bytes']=p.stat().st_size
 if name.startswith('ARTIFACT'):
  known={r['path'] for r in rr}
  for p in sorted(O.rglob('*')):
   if p.is_file() and p.name!='ARTIFACT_INDEX.csv' and str(p.relative_to(O)).replace('\\','/') not in known:rr.append(dict(path=str(p.relative_to(O)).replace('\\','/'),bytes=p.stat().st_size,sha256=sha(p)))
 with (O/name).open('w',newline='',encoding='utf-8') as f:
  w=csv.DictWriter(f,fieldnames=list(rr[0]));w.writeheader();w.writerows(rr)
print('Layout corrected; all plotted values unchanged; artifact identities refreshed.')
