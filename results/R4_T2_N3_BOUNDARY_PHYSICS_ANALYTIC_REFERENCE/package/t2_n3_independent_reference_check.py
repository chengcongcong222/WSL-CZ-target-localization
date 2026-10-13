#!/usr/bin/env python3
"""Independent, solver-free check of archived ideal-waveguide analytic preview.

Never imports t2_n3_analytic_reference, KRAKEN, or archived modes.
Computes mode counts by strict half-integer inequalities, and recomputes
all archived shape differences through an independent candidate loop.
"""
from pathlib import Path
import csv,json,math
import numpy as np
from scipy.special import hankel2
P=Path(__file__).resolve().parent
H=5000.;c=1500.; k2=2*math.pi

def read_rows(name):
    with (P/name).open(newline='',encoding='utf-8-sig') as f:return list(csv.DictReader(f))

resources={}
for r in read_rows('T2_TILT_G0_ARRAY_GEOMETRY_20261013_v1.csv'):
    if r['resource'] in ('A0_0','A_P4M','B1_V4','B2_DUP'):
        resources.setdefault(r['resource'],[]).append((int(r['channel']),float(r['x_m']),float(r['y_m']),float(r['z_m'])))
for x in list(resources):resources[x]=np.array([a[1:] for a in sorted(resources[x])],float)
probes={'P01':(45.5,0.,180.),'P02':(50.,0.,200.),'P03':(55.,2.5,200.),'P04':(58.,4.5,210.),'P05':(59.5,5.,220.)}

def count(f,BC,limit):
    s=2*f*H/c
    if limit!='ALL_PROPAGATING':s*=math.sqrt(1-(c/float(limit))**2)
    # n+half < s or n+1 <s ; don't rely on numerical mode root enumeration.
    return int(max(0,math.ceil(s-.5))) if BC=='PR_R' else int(max(0,math.ceil(s-1)))

def psi(f,BC,limit,position,receiver):
    r,th,z=position;xy=1000*r*np.array((math.cos(math.radians(th)),math.sin(math.radians(th))))
    ranges=np.sqrt(np.sum((xy[None,:]-receiver[:,:2])**2,axis=1));omega=k2*f
    nmax=count(f,BC,limit)
    mm=np.arange(nmax,dtype=float)+(.5 if BC=='PR_R' else 1.)
    mu=mm*math.pi/H;kr=np.sqrt((omega/c)**2-mu**2)
    source_sin=np.sin(mu*z)*math.sqrt(2/H)
    recv_sin=np.sin(mu[:,None]*receiver[None,:,2])*math.sqrt(2/H)
    return np.sum(source_sin[:,None]*recv_sin*(-1j*math.pi*hankel2(0,kr[:,None]*ranges[None,:])),axis=0)

def shape(v,w):
    v=v/np.linalg.norm(v);w=w/np.linalg.norm(w)
    return float(np.linalg.norm(v-np.vdot(w,v)*w))

check=0;max_diff=0.
for row in read_rows('ANALYTIC_WAVEGUIDE_MODE_COUNTS.csv'):
    n=count(float(row['frequency_hz']),row['boundary'],row['chigh_mps'])
    assert n==int(row['analytic_modes']),(row,n)
    check+=1

for row in read_rows('ANALYTIC_WAVEGUIDE_WINDOW_SHAPES.csv'):
    BC=row['boundary'];f=float(row['frequency_hz']);pos=probes[row['probe']];rx=resources[row['resource']]
    a=psi(f,BC,row['left_chigh_mps'],pos,rx)
    b=psi(f,BC,row['right_chigh_mps'],pos,rx)
    s=shape(a,b);delta=abs(s-float(row['shape_difference']));max_diff=max(max_diff,delta)
    assert delta<=1e-10,(row,delta)
    check+=1

from scipy.special import hankel2
results={'PASS':True,'mode_count_records':16,'window_shape_records':240,'checks':check,'max_independent_shape_difference_discrepancy':max_diff,'solver_calls':0,'scope':'closed_form_ideal_constant_speed_waveguide_only'}
(P/'ANALYTIC_REFERENCE_INDEPENDENT_CHECK.json').write_text(json.dumps(results,ensure_ascii=False,indent=2)+'\n',encoding='utf-8')
print(json.dumps(results,ensure_ascii=False))
