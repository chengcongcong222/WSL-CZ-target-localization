"""Outward rounded modal level and Q2 enclosures; no archived observation loader."""
from __future__ import annotations
from contextlib import contextmanager
from pathlib import Path
import numpy as np
from flint import arb, acb, fmpq, ctx
from r4_cz_envelope_representation import FREQUENCIES, TIMES, WINDOWS
from r4_cz_envelope_support import exact, Cell, Bounds, Certificate, DEPTH_LABELS

PRECISION = 160
MODEL_DIRECTORY = Path(__file__).resolve().parent/'results/R3_C2_Yang_SA_depth/R3_C2_1/_kraken_zgrid'

@contextmanager
def precision(bits=PRECISION):
    old=ctx.prec;ctx.prec=bits
    try:yield
    finally:ctx.prec=old

def platform_xy():
    dt=np.maximum(TIMES-600.,0.)
    return (2*np.minimum(TIMES,600.)+2*dt*np.cos(np.deg2rad(15.)),2*dt*np.sin(np.deg2rad(15.)))

def geometry(cell):
    r,theta,v,psi=[exact(a).union(exact(b)) for a,b in zip(cell.low,cell.high)]
    theta=theta*arb.pi()/180;psi=psi*arb.pi()/180
    xp,yp=platform_xy(); ranges=[];angles=[]
    for t,px,py in zip(TIMES,xp,yp):
        dx=r*1000*theta.cos()+v*exact(t)*psi.cos()-exact(px)
        dy=r*1000*theta.sin()+v*exact(t)*psi.sin()-exact(py)
        if not bool(dx>0):raise ArithmeticError('Cannot certify positive horizontal x')
        rr=(dx**2+dy**2).sqrt()
        if not rr.is_finite() or not bool(rr>0):raise ArithmeticError('Cannot certify positive range')
        ranges.append(rr);angles.append((dy/dx).atan())
    return angles,ranges

def point_geometry(state):
    r,theta,v,psi=map(float,state);theta,psi=np.deg2rad([theta,psi]);xp,yp=platform_xy()
    dx=r*1000*np.cos(theta)+v*TIMES*np.cos(psi)-xp
    dy=r*1000*np.sin(theta)+v*TIMES*np.sin(psi)-yp
    return np.arctan2(dy,dx),np.hypot(dx,dy)

def certified_basis(n):
    x=[fmpq(2*i-(n-1)) for i in range(n)]
    den1=arb(sum(v*v for v in x)).sqrt()
    avg=sum(v*v for v in x)/n
    q=[v*v-avg for v in x];den2=arb(sum(v*v for v in q)).sqrt()
    return [(arb(v)/den1,arb(w)/den2) for v,w in zip(x,q)]

def project_balls(levels):
    out=[]
    for line in levels:
        if len(line)!=121:raise ValueError('121 level balls per line required')
        for section in WINDOWS:
            values=line[section];mean=sum(values,arb(0))/len(values)
            b=certified_basis(len(values))
            out.extend(sum(((val-mean)*row[k] for val,row in zip(values,b)),arb(0)) for k in (0,1))
    return tuple(out)

def distance_bounds(predicted, coefficients):
    lo=arb(0);hi=arb(0)
    for p,o in zip(predicted,coefficients):
        delta=p-exact(o);lo+=delta.abs_lower()**2;hi+=delta.abs_upper()**2
    return Bounds((lo/6).sqrt().lower().max(0).lower(),(hi/6).sqrt().upper())

def lower_bound_bearing_cost(cell, observation, budget):
    budget.charge(n_bound_evaluations=1)
    angles,_=geometry(cell)
    return _bearing_bounds(angles,observation)

def _bearing_bounds(angles,observation):
    lo=arb(0);hi=arb(0);pi=arb.pi()
    for p,o in zip(angles,observation.bearing_rad):
        delta=p-exact(o)
        if bool(delta>=-pi) and bool(delta<=pi):
            lo+=delta.abs_lower()**2;hi+=delta.abs_upper()**2
        else:
            # Uncertain principal wrap: universal nonnegative SSE enclosure.
            hi+=pi.upper()**2
    return Bounds(lo.lower().max(0).lower(),hi.upper())

class ModalModel:
    def __init__(self,directory=MODEL_DIRECTORY):
        self.lines={};self.mapping=[]
        for f in FREQUENCIES:
            data=(Path(directory)/f'zgrid_f{f}.mod').read_bytes()
            recl=4*int(np.frombuffer(data[:4],dtype='<i4')[0]);hdr=np.frombuffer(data[84:108],dtype='<i4')
            ntot,nmat=int(hdr[2]),int(hdr[3])
            depths=np.frombuffer(data[4*recl:5*recl],dtype='<f4')[:ntot].astype(float)
            count=int(np.frombuffer(data[5*recl:5*recl+4],dtype='<i4')[0])
            phi=np.empty((nmat,count),dtype=complex)
            for im in range(count):phi[:,im]=np.frombuffer(data[(7+im)*recl:(8+im)*recl],dtype='<c8')[:nmat]
            k=np.frombuffer(data[(7+count)*recl:(7+count)*recl+count*8],dtype='<c8').copy()
            if (ntot!=nmat or len(depths)!=nmat or count<=0 or not np.all(np.diff(depths)>0) or
                not np.isfinite(phi).all() or not np.isfinite(k).all() or not np.all(k.real>0)):
                raise ValueError('Invalid inherited modal layout')
            izr=int(np.argmin(abs(depths-200.)))
            self.lines[f]={'depths':depths,'phi':phi,'k':k,'izr':izr,'M':count}
            for label in DEPTH_LABELS:
                iz=int(np.argmin(abs(depths-label)))
                self.mapping.append({'frequency_hz':f,'profile_label_m':label,'effective_source_depth_m':float(depths[iz]),'receiver_depth_m':float(depths[izr])})
    def levels_numpy(self,state,label,budget):
        budget.charge(n_forward_state_evaluations=1,n_profile_evaluations=1)
        _,ranges=point_geometry(state);out=[]
        for f in FREQUENCIES:
            m=self.lines[f];iz=int(np.argmin(abs(m['depths']-label)));weights=m['phi'][iz]*m['phi'][m['izr']]
            kre=m['k'].real[:,None];alpha=-m['k'].imag[:,None];r=ranges[None,:]
            factors=np.sqrt(2*np.pi/(kre*r))*np.exp(-1j*kre*r-alpha*r-1j*np.pi/4)
            p=weights@factors;out.append(20*np.log10(np.maximum(abs(p),1e-30)))
        return np.asarray(out)
    def levels_arb(self,cell,label,budget,*,ranges=None):
        return self.levels_arb_batch(cell,(label,),budget,ranges=ranges)[0]
    def levels_arb_batch(self,cell,labels,budget,*,ranges=None):
        if not labels or any(label not in DEPTH_LABELS for label in labels):raise ValueError('Only inherited finite labels')
        budget.charge(n_forward_state_evaluations=len(labels),n_profile_evaluations=len(labels))
        return self._levels_arb_batch(cell,labels,ranges=ranges)
    def _levels_arb_batch(self,cell,labels,*,ranges=None):
        if ranges is None:_,ranges=geometry(cell)
        out=[[] for _ in labels];pi=arb.pi();floor=exact(1e-30)
        for f in FREQUENCIES:
            m=self.lines[f];weights=[]
            params=[(exact(k.real),-exact(k.imag)) for k in m['k']]
            for label in labels:
                iz=int(np.argmin(abs(m['depths']-label)));row=[]
                for a,b in zip(m['phi'][iz],m['phi'][m['izr']]):
                    row.append(acb(exact(a.real),exact(a.imag))*acb(exact(b.real),exact(b.imag)))
                weights.append(row)
            lines=[[] for _ in labels]
            for r in ranges:
                # Common expression elimination inside one batch, not a persistent cache.
                factors=[acb(-alpha*r,-k*r-pi/4).exp()*(2*pi/(k*r)).sqrt() for k,alpha in params]
                for line,row in zip(lines,weights):
                    pressure=sum((w*z for w,z in zip(row,factors)),acb(0))
                    plo=pressure.abs_lower();phi=pressure.abs_upper()
                    # Take log of positive endpoints BEFORE ball union. A wide amplitude
                    # ball crossing zero must never cause a spurious log failure.
                    positive_lo=plo if bool(plo>floor) else floor
                    positive_hi=phi.max(floor)
                    lower=20*positive_lo.log()/arb(10).log()
                    upper=20*positive_hi.log()/arb(10).log()
                    value=lower.union(upper)
                    if not value.is_finite():raise ArithmeticError('Nonfinite level enclosure')
                    line.append(value)
            for result,line in zip(out,lines):result.append(line)
        return out
    def lower_bound_feature_mismatch(self,cell,label,observation,budget,*,ranges=None):
        if label not in DEPTH_LABELS:raise ValueError('Only inherited finite labels')
        budget.charge(n_bound_evaluations=1,n_forward_state_evaluations=1,n_profile_evaluations=1)
        levels=self._levels_arb_batch(cell,(label,),ranges=ranges)[0]
        return distance_bounds(project_balls(levels),observation.coefficients)
    def certificate(self,cell,observation,budget,baseline=False):
        with precision():
            budget.charge(n_bound_evaluations=1)
            angles,ranges=geometry(cell)
            bearing=_bearing_bounds(angles,observation)
            if baseline or bool(bearing.lower>exact(observation.bearing_cutoff)):
                return Certificate(cell,bearing,None,(),'ARB_WHOLE_CELL_V1')
            budget.charge(n_bound_evaluations=21,n_forward_state_evaluations=21,n_profile_evaluations=21)
            levels=self._levels_arb_batch(cell,DEPTH_LABELS,ranges=ranges)
            bounds=[distance_bounds(project_balls(row),observation.coefficients) for row in levels]
            lower=bounds[0].lower;upper=bounds[0].upper
            for b in bounds[1:]:lower=lower.min(b.lower).lower();upper=upper.min(b.upper).upper()
            return Certificate(cell,bearing,Bounds(lower.max(0).lower(),upper),DEPTH_LABELS,'ARB_WHOLE_CELL_V1')
