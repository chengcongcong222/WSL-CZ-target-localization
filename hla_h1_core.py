"""HLA-H1 pure model and likelihood; no historical module imports or writes."""
from pathlib import Path
import hashlib, json, time
import numpy as np
from scipy.interpolate import CubicSpline
from scipy.linalg import helmert
from scipy.stats import chi2
ROOT=Path(__file__).resolve().parent
OUT=ROOT/'results/R4_SINGLE_HLA_SOURCE_PROFILED_HORIZONTAL_PILOT'
LOCAL=Path('D:/ProjectStorage/WSL-CZ/HLA_H1')
MODES=ROOT/'results/R3_C2_Yang_SA_depth/R3_C2_1/_kraken_zgrid'
AXES=('r_km','theta_deg','v_mps','psi_deg')
SHAPE=(16,21,11,31)
ORIGIN=np.array([45.,-5.,1.,-15.])
STEP=np.array([1.,.5,.2,1.])
LOW=ORIGIN.copy()
HIGH=np.array([60.,5.,3.,15.])
TIMES=np.arange(0.,1201.,10.)
PROFILE=np.arange(150.,251.,5.)
FREQS=(201,235,283)
GEOMS=np.array([[50,0,2,5],[47,-1,1.4,-7],[54,1,2.6,9],[59,-.5,1.8,-11],[56,-1.5,2.2,3],[48,1.5,2.8,-5]],float)
ZS=np.array([200,180,220,200,180,220])
METHODS=('BEARING','M0','M1','ORACLE')
DF={'BEARING':121,'M0':478,'M1':359,'ORACLE':478}
CUT={k:float(chi2.ppf(.95,v)) for k,v in DF.items()}
UF=helmert(3,full=False).T
UT=[helmert(n,full=False).T for n in (61,60)]
SECTIONS=(slice(0,61),slice(61,121))
def sha(p):
    h=hashlib.sha256()
    with Path(p).open('rb') as f:
        for b in iter(lambda:f.read(1048576),b''):h.update(b)
    return h.hexdigest()
def write_json(p,obj):
    Path(p).write_text(json.dumps(obj,indent=2,ensure_ascii=False,allow_nan=False)+'\n',encoding='utf-8')
def grid():
    return np.round(ORIGIN+np.indices(SHAPE).reshape(4,-1).T*STEP,10)
def wrap(x):return (np.asarray(x)+np.pi)%(2*np.pi)-np.pi
def geometry(states):
    s=np.atleast_2d(states);r,th,v,ps=s.T;th=np.radians(th);ps=np.radians(ps)
    dt=np.maximum(TIMES-600,0)
    px=2*np.minimum(TIMES,600)+2*dt*np.cos(np.radians(15))
    py=2*dt*np.sin(np.radians(15))
    x=1000*r[:,None]*np.cos(th[:,None])+v[:,None]*TIMES*np.cos(ps[:,None])-px
    y=1000*r[:,None]*np.sin(th[:,None])+v[:,None]*TIMES*np.sin(ps[:,None])-py
    return np.arctan2(y,x),np.hypot(x,y)
def parse_mod(path):
    b=Path(path).read_bytes();rec=4*int(np.frombuffer(b[:4],'<i4')[0])
    hdr=np.frombuffer(b[84:108],'<i4');nt,nm=map(int,hdr[2:4])
    depths=np.frombuffer(b[4*rec:5*rec],'<f4')[:nt].astype(float)
    count=int(np.frombuffer(b[5*rec:5*rec+4],'<i4')[0])
    phi=np.array([np.frombuffer(b[(7+i)*rec:(8+i)*rec],'<c8')[:nm] for i in range(count)]).T.astype(complex)
    k=np.frombuffer(b[(7+count)*rec:(7+count)*rec+8*count],'<c8').copy()
    assert nt==nm and np.isfinite(phi).all() and np.all(np.diff(depths)>0)
    return dict(depths=depths,phi=phi,k=k)
def project(levels,method):
    chunks=[]
    for sl,u in zip(SECTIONS,UT):
        q=np.asarray(levels)[...,sl]@u
        if method=='M1':q=np.einsum('af,...ft->...at',UF.T,q)
        chunks.append(q.reshape(*q.shape[:-2],-1))
    return np.concatenate(chunks,axis=-1)
def norm_rms(x):
    x=np.asarray(x)-np.mean(x)
    return x/np.sqrt(np.mean(x*x))
def sources():
    ramp=norm_rms(TIMES/1200-.5);sine=norm_rms(np.sin(2*np.pi*TIMES/1200))
    return np.array([np.zeros((3,121)),np.tile(ramp,(3,1)),np.tile(sine,(3,1)),
                     np.tile(sine,(3,1))+.25*np.array([1,-1,1])[:,None]*sine])
class Model:
    def __init__(self):
        self.mods={f:parse_mod(MODES/f'zgrid_f{f}.mod') for f in FREQS}
        self.splines={};self.weights={};self.carrier={};self.mapping=[]
        x=np.arange(39000.,66001.,1.)
        for f,m in self.mods.items():
            iz=np.abs(m['depths'][:,None]-PROFILE).argmin(axis=0)
            ir=int(np.abs(m['depths']-200).argmin())
            self.weights[f]=m['phi'][iz]*m['phi'][ir]
            self.carrier[f]=float((m['k'].real.min()+m['k'].real.max())/2)
            self.mapping.append(dict(frequency=f,labels=PROFILE.tolist(),actual=m['depths'][iz].tolist(),receiver_actual=float(m['depths'][ir])))
            p=np.empty((len(x),21),complex)
            for start in range(0,len(x),1024):p[start:start+1024]=self.pressure(f,x[start:start+1024],demod=True).T
            self.splines[f]=CubicSpline(x,p,axis=0,extrapolate=False)
    def pressure(self,f,ranges,demod=False):
        r=np.asarray(ranges).ravel()[None,:];m=self.mods[f]
        k=m['k'].real.astype(float)[:,None];a=-m['k'].imag.astype(float)[:,None]
        phase=k-(self.carrier[f] if demod else 0)
        return self.weights[f]@(np.sqrt(2*np.pi/(k*r))*np.exp(-1j*phase*r-a*r-1j*np.pi/4))
    def levels(self,states,exact=False):
        ranges=geometry(states)[1];n,nt=ranges.shape
        out=np.empty((n,21,3,nt))
        for j,f in enumerate(FREQS):
            if exact:
                p=self.pressure(f,ranges).reshape(21,n,nt).transpose(1,0,2)
            else:
                p=self.splines[f](ranges).transpose(0,2,1)
            out[:,:,j]=20*np.log10(np.maximum(np.abs(p),1e-30))
        if not np.isfinite(out).all():raise ValueError('Nonfinite inherited level field')
        return out
def bearing_cost(pred,obs):
    return np.square(wrap(pred-np.asarray(obs))/np.radians(.1)).sum(axis=-1)
def errors(s,true):
    e=np.abs(np.asarray(s)-true)
    e[...,1]=np.abs(wrap(np.radians(np.asarray(s)[...,1]-true[1])))*180/np.pi
    e[...,3]=np.abs(wrap(np.radians(np.asarray(s)[...,3]-true[3])))*180/np.pi
    return e
