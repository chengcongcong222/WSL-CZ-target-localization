#!/usr/bin/env python3
"""Independent analytic homogeneous ideal-waveguide benchmark (no KRAKEN imports).

Two *mathematical* boundary contracts, not site-specific seabed representations:
PR_R: p(0)=0, dp/dz(H)=0, mu_n=(n+1/2)*pi/H, n>=0
PR_PR: p(0)=p(H)=0, mu_n=n*pi/H, n>=1
For propagating modes kr_n=sqrt((2*pi*f/c)^2-mu_n^2),
phi_n(z)=sqrt(2/H)*sin(mu_n*z).
Exact radial Green factor is (-i*pi)*H_0^(2)(kr*r). Its large-kr asymptote
is sqrt(2*pi/(kr*r))*exp(-i*kr*r-i*pi/4), matching the old project's
far-field convention up to a common coefficient. Finite windows select
phase speed cp=omega/kr strictly below the registered CHIGH (m/s).
Only outgoing propagating modes are compared; evanescent modes are omitted.
"""
from __future__ import annotations
import argparse,csv,math
from pathlib import Path
import numpy as np
from scipy.special import hankel2
ROOT=Path(__file__).resolve().parent
H=5000.0
C=1500.0
FREQS=(201.,283.)
WINDOWS=(1800.,2400.,4000.,float('inf'))
PROBES=(('P01',45.5,0.,180.),('P02',50.,0.,200.),('P03',55.,2.5,200.),('P04',58.,4.5,210.),('P05',59.5,5.,220.))
RESOURCES=('A0_0','A_P4M','B1_V4','B2_DUP')

def modal_spectrum(freq:float, boundary:str):
    # Explicit analytic orthogonal eigenfunctions independent of numerical solver.
    w=2*np.pi*freq;k0=w/C
    # Strict propagation condition is written in an exact-cutoff-stable
    # dimensionless form. Comparing independently rounded mu and k0
    # would incorrectly admit the PR_PR 201 Hz mode exactly at kr=0.
    cutoff_ratio=2*freq*H/C
    n=np.arange(int(np.ceil(cutoff_ratio))+2,dtype=float)
    mode_order=n+.5 if boundary=='PR_R' else n+1.
    mu=(mode_order[mode_order<cutoff_ratio])*np.pi/H
    kr=np.sqrt(k0*k0-mu*mu)
    cp=w/kr
    return mu,kr,cp

def geom_arrays():
    p=ROOT/'T2_TILT_G0_ARRAY_GEOMETRY_20261013_v1.csv'
    by={k:[] for k in RESOURCES}
    with p.open(newline='',encoding='utf-8-sig') as f:
        for x in csv.DictReader(f):
            if x['resource'] in by:
                by[x['resource']].append((int(x['channel']),float(x['x_m']),float(x['y_m']),float(x['z_m'])))
    out={k:np.array([a[1:] for a in sorted(v)],float) for k,v in by.items()}
    assert [len(out[r]) for r in RESOURCES]==[8,8,12,12]
    return out

def pressure(freq,boundary,cp_limit,r_km,theta_deg,z_src,coords,use_exact=True):
    mu,kr,cp=modal_spectrum(freq,boundary)
    selected=(cp<cp_limit)
    mu=mu[selected];kr=kr[selected]
    xy=1000*r_km*np.array([np.cos(np.deg2rad(theta_deg)),np.sin(np.deg2rad(theta_deg))])
    dist=np.linalg.norm(xy[None,:]-coords[:,:2],axis=1)
    phi_src=np.sqrt(2/H)*np.sin(mu*z_src)
    phi_recv=np.sqrt(2/H)*np.sin(mu[:,None]*coords[None,:,2])
    x=kr[:,None]*dist[None,:]
    radial=(-1j*np.pi)*hankel2(0,x) if use_exact else np.sqrt(2*np.pi/x)*np.exp(-1j*x-1j*np.pi/4)
    p=np.sum((phi_src[:,None]*phi_recv)*radial,axis=0,dtype=np.complex128)
    return p,int(selected.sum()),len(selected)

def shape(a,b):
    x=a/np.linalg.norm(a);y=b/np.linalg.norm(b)
    return float(np.linalg.norm(x-y*np.vdot(y,x)))

def main(outdir:Path):
    outdir.mkdir(parents=True,exist_ok=True)
    resources=geom_arrays()
    count=[];pairs=[];exact_vs_far=[]
    for bc in ('PR_R','PR_PR'):
      for f in FREQS:
        mu,kr,cp=modal_spectrum(f,bc)
        for limit in WINDOWS:
          count.append({'boundary':bc,'frequency_hz':int(f),'chigh_mps':'ALL_PROPAGATING' if np.isinf(limit) else int(limit),'analytic_modes':int(np.sum(cp<limit)),'total_propagating_modes':len(cp),'max_cp_full_mps':float(np.max(cp)),'source':'closed_form_kz_eigenvalues_no_solver'})
        for id,r,t,z in PROBES:
          for name,coords in resources.items():
            pressures={}
            for lim in WINDOWS:
              p,m,tot=pressure(f,bc,lim,r,t,z,coords)
              pressures[lim]=p
            far,m,_=pressure(f,bc,4000.,r,t,z,coords,use_exact=False)
            exact_vs_far.append({'boundary':bc,'frequency_hz':int(f),'probe':id,'resource':name,'window_mps':4000,'exact_Hankel_vs_farfield_shape':shape(pressures[4000.],far),'exact_Hankel_vs_farfield_raw_globally_profiled':float(np.linalg.norm(far-pressures[4000.]*np.vdot(pressures[4000.],far)/np.vdot(pressures[4000.],pressures[4000.]))/np.linalg.norm(far))})
            for a,b in zip(WINDOWS,WINDOWS[1:]):
              pairs.append({'boundary':bc,'frequency_hz':int(f),'probe':id,'resource':name,'left_chigh_mps':int(a),'right_chigh_mps':'ALL_PROPAGATING' if np.isinf(b) else int(b),'left_modes':int(np.sum(cp<a)),'right_modes':int(np.sum(cp<b)),'shape_difference':shape(pressures[a],pressures[b]),'model_scope':'constant_c_ideal_boundary_propagating_Hankel_only'})
    def write(name,data):
      with (outdir/name).open('w',encoding='utf-8',newline='') as handle:
        w=csv.DictWriter(handle,fieldnames=list(data[0]));w.writeheader();w.writerows(data)
    write('ANALYTIC_WAVEGUIDE_MODE_COUNTS.csv',count)
    write('ANALYTIC_WAVEGUIDE_WINDOW_SHAPES.csv',pairs)
    write('ANALYTIC_EXACT_RADIAL_VS_FARFIELD.csv',exact_vs_far)
    print('ANALYTIC modes',len(count),'window cases',len(pairs),'Hankel checks',len(exact_vs_far))
    for bc in ('PR_R','PR_PR'):
      for f in FREQS:
        q=[r for r in count if r['boundary']==bc and r['frequency_hz']==int(f)]
        a=[p['shape_difference'] for p in pairs if p['boundary']==bc and p['frequency_hz']==int(f) and p['right_chigh_mps']==2400]
        b=[p['shape_difference'] for p in pairs if p['boundary']==bc and p['frequency_hz']==int(f) and p['right_chigh_mps']=='ALL_PROPAGATING']
        print(bc,int(f),'N1800/2400/4000/FULL',[r['analytic_modes'] for r in q],'shape1800to2400 max',round(max(a),6),'shape4000tofull max',round(max(b),6))
    print('farfield shape max',max(z['exact_Hankel_vs_farfield_shape'] for z in exact_vs_far))

if __name__=='__main__':
    ar=argparse.ArgumentParser();ar.add_argument('--outdir',type=Path,default=ROOT/'analytic_preview');q=ar.parse_args();main(q.outdir)
