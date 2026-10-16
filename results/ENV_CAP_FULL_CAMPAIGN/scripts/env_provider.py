"""ENV-CAP new ray/TL provider. No old experiment is imported or executed."""
from pathlib import Path
import numpy as np,json,csv,hashlib,subprocess,time,sys,os
from scipy.ndimage import uniform_filter1d
from scipy.signal import find_peaks
R=Path(__file__).resolve().parents[3];O=R/'results/ENV_CAP_FULL_CAMPAIGN';W=Path('D:/ProjectStorage/WSL-CZ/ENV_CAP')
EXE=R/'tools/acoustics_toolbox/atWin10/at/bin/bellhop.exe'
def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def js(p,x):Path(p).write_text(json.dumps(x,ensure_ascii=False,indent=2,default=lambda v:v.item() if hasattr(v,'item') else str(v))+'\n',encoding='utf-8',newline='\n')
def table(p,rows,cols=None):
 if cols is None:cols=list(dict.fromkeys(k for row in rows for k in row))
 with Path(p).open('w',encoding='utf-8',newline='') as f:
  a=csv.DictWriter(f,cols);a.writeheader();a.writerows(rows)
def loadprofile(path):
 a=list(csv.DictReader(Path(path).open(encoding='utf-8-sig')));return np.array([float(x['depth_m']) for x in a]),np.array([float(x['sound_speed_mps']) for x in a])
def cut(z,c,H):
 assert z[0]==0 and z[-1]>=H and np.all(np.diff(z)>0)
 return np.r_[z[z<H],H],np.r_[c[z<H],np.interp(H,z,c)]
def integral(z,c,p):
 """Exact dx for C-linear segments; singular turning endpoint integrated analytically."""
 a=p*c[:-1];b=p*c[1:];g=np.diff(c)/np.diff(z)
 out=np.empty(len(g));m=np.abs(g)>1e-12
 out[m]=(np.sqrt(np.maximum(0,1-a[m]**2))-np.sqrt(np.maximum(0,1-b[m]**2)))/(p*g[m])
 out[~m]=np.diff(z)[~m]*a[~m]/np.sqrt(1-a[~m]**2)
 return np.sum(out)
def return_curve(z,c,zs,angles):
 rows=[];cs=np.interp(zs,z,c)
 for a in angles:
  p=np.cos(np.deg2rad(a))/cs;cp=1/p
  deep=np.flatnonzero((z[1:]>zs)&(c[:-1]<cp)&(c[1:]>=cp)&(np.diff(c)>0))
  if not len(deep):continue
  j=deep[0];turn=z[j]+(cp-c[j])/(c[j+1]-c[j])*(z[j+1]-z[j]);zz=np.r_[zs,z[(z>zs)&(z<turn)],turn];cc=np.interp(zz,z,c);cc[-1]=cp
  rr=2*integral(zz,cc,p)/1000
  rows.append({'angle_deg':float(a),'first_return_km':float(rr),'turn_z_m':float(turn),'c_phase_mps':cp,'surface_bounces':0,'bottom_bounces':0})
 return rows

def write_env(folder,tag,z,c,H,zs,kind='I',nbeam=2001,step=50,angles=(-25,25),rmax=120,depths=None):
 folder=Path(folder);folder.mkdir(parents=True,exist_ok=True);p=folder/(tag+'.env')
 if depths is None:depths=np.arange(100,300.0001,2.5)
 ranges=np.arange(10,rmax+.001,.125)
 lines=[f"'ENV_CAP_{tag}'",'200.0','1',"'CVN'",f'0 0.0 {H:.10f}']
 lines.extend(f'{zi:.10f} {ci:.12f} /' for zi,ci in zip(z,c))
 # Nominal rigid bottom reproduces no real sediment contract; refractive ray identity is separate.
 lines.extend(["'R' 0.0",'1',f'{zs:.9f} /',str(len(depths)),' '.join(f'{v:.5f}' for v in depths)+' /',str(len(ranges)),' '.join(f'{v:.5f}' for v in ranges)+' /',f"'{kind}'",str(nbeam),f'{angles[0]:.9f} {angles[1]:.9f} /',f'{step:.9f} {H+10:.9f} {rmax+1:.9f}'])
 p.write_text('\n'.join(lines)+'\n',encoding='ascii',newline='\n');return p

def readshd(p):
 b=Path(p).read_bytes();rec=4*int(np.frombuffer(b[:4],'<i4')[0]);dims=np.frombuffer(b[2*rec:2*rec+28],'<i4');nf,nt,nx,ny,ns,nz,nr=map(int,dims)
 assert (nf,nt,nx,ny,ns)==(1,1,1,1,1),dims
 rr=np.frombuffer(b[9*rec:9*rec+4*nr],'<f4').astype(float)/1000;zz=np.frombuffer(b[8*rec:8*rec+4*nz],'<f4').astype(float)
 assert len(b)>=(10+nz)*rec
 pp=np.stack([np.frombuffer(b[(10+i)*rec:(10+i)*rec+nr*8],'<c8').astype(complex) for i in range(nz)])
 assert np.all(np.isfinite(pp)) and np.all(np.diff(rr)>0)
 return rr,zz,pp

def readrays(p):
 lines=Path(p).read_text().splitlines();assert 'rz' in lines[6].lower() or 'xyz' in lines[6].lower(),lines[:8]
 dim=2 if 'rz' in lines[6].lower() else 3;out=[];i=7
 while i<len(lines):
  if not lines[i].strip():i+=1;continue
  angle=float(lines[i]);n,ns,nb=map(int,lines[i+1].split());i+=2
  ar=np.array([[float(v) for v in s.split()] for s in lines[i:i+n]]);i+=n
  assert ar.shape==(n,dim)
  out.append((angle,ns,nb,ar))
 return out

def ray_first_returns(rays,zs,H):
 rows=[]
 for a,ns,nb,x in rays:
  zz=x[:,-1];rr=x[:,0]/1000
  # Registered downward launches only. First upward crossing after depth exceeded source.
  hits=np.flatnonzero((zz[:-1]>zs)&(zz[1:]<=zs))
  if not len(hits):continue
  j=hits[0];seg=zz[:j+2]
  if np.min(seg)<1 or np.max(seg)>=H-1:continue
  q=(zs-zz[j])/(zz[j+1]-zz[j]);r=rr[j]+q*(rr[j+1]-rr[j]);rows.append({'angle_deg':a,'first_return_km':r,'turn_z_m':np.max(seg),'surface_bounces':0,'bottom_bounces':0})
 return rows

def execute(p,unit,timeout=120):
 start=time.perf_counter();rc=None;error=None
 try:
  v=subprocess.run([str(EXE),p.stem],cwd=p.parent,capture_output=True,timeout=timeout)
  (p.parent/(p.stem+'.stdout')).write_bytes(v.stdout+b'\n'+v.stderr);rc=v.returncode
 except subprocess.TimeoutExpired as e:error='TIMEOUT';(p.parent/(p.stem+'.stdout')).write_bytes((e.stdout or b'')+(e.stderr or b''))
 row={'unit':unit,'tag':p.stem,'input_sha256':sha(p),'solver_sha256':sha(EXE),'elapsed_s':time.perf_counter()-start,'returncode':rc,'error':error}
 with (W/'SOLVER_CALLS.jsonl').open('a',encoding='utf-8') as f:f.write(json.dumps(row)+'\n')
 if rc!=0:raise RuntimeError(row)
 # Bellhop can return zero on Fortran STOP; reject diagnostic failure strings.
 prt=p.with_suffix('.prt').read_text(errors='replace')
 if any(s in prt.lower() for s in ['fatal','insufficient storage','error in']):raise RuntimeError(prt[-1500:])
 return row

def energycurve(rr,zz,p,zs):
 band=np.abs(zz-zs)<=25
 e=np.mean(np.abs(p[band])**2,axis=0)
 # 1.125 km moving mean ENERGY, never dB; fixed 9 samples at 125m.
 e=uniform_filter1d(e,9,mode='nearest');return -10*np.log10(np.maximum(e,1e-30))

def focus(rr,tl,expected):
 valid=(rr>=expected-5)&(rr<=expected+5)
 if not valid.any():return None
 idx=np.flatnonzero(valid)[np.argmin(tl[valid])];r=rr[idx]
 # fixed 5..10km flanking background excludes proposed focus.
 flank=(np.abs(rr-r)>=5)&(np.abs(rr-r)<=10)
 prom=float(np.median(tl[flank])-tl[idx]) if flank.any() else -np.inf
 half=tl[idx]+max(prom,0)/2;l=u=idx
 while l>0 and tl[l]<=half:l-=1
 while u<len(rr)-1 and tl[u]<=half:u+=1
 return {'peak_km':r,'peak_TL_dB':tl[idx],'prominence_dB':prom,'width_km':rr[u]-rr[l],'boundary_selected':idx in [np.flatnonzero(valid)[0],np.flatnonzero(valid)[-1]]}
