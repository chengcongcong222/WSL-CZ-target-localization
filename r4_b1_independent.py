"""Independent B1 scalar geometry, struct modal parser and ordered modal sums."""
from pathlib import Path
import struct
import math
import numpy as np
ROOT=Path(__file__).resolve().parent

def independent_models():
    models={}
    for f in (201,235,283,338):
        buf=(ROOT/'results/R3_C2_Yang_SA_depth/R3_C2_1/_kraken_zgrid'/f'zgrid_f{f}.mod').read_bytes()
        record=4*struct.unpack_from('<i',buf,0)[0]
        header=struct.unpack_from('<6i',buf,84); ntot,nmat=header[2],header[3]
        depths=np.array(struct.unpack_from(f'<{ntot}f',buf,4*record),float)
        count=struct.unpack_from('<i',buf,5*record)[0]
        columns=[]
        for m in range(count):
            packed=struct.unpack_from(f'<{2*nmat}f',buf,(7+m)*record)
            columns.append([complex(packed[2*i],packed[2*i+1]) for i in range(nmat)])
        packed=struct.unpack_from(f'<{2*count}f',buf,(7+count)*record)
        ks=np.array([complex(packed[2*m],packed[2*m+1]) for m in range(count)])
        models[f]={'depths':depths,'phi':np.array(columns).T,'k':ks,'M':count}
    return models

def independent_ranges(horizontal):
    r,theta,v,psi=horizontal; th=math.radians(theta); ps=math.radians(psi)
    turn=math.radians(15.)
    out=[]
    for t in range(0,1201,10):
        xp=2*t if t<=600 else 1200+2*(t-600)*math.cos(turn)
        yp=0. if t<=600 else 2*(t-600)*math.sin(turn)
        dx=r*1000*math.cos(th)+v*t*math.cos(ps)-xp
        dy=r*1000*math.sin(th)+v*t*math.sin(ps)-yp
        out.append(math.hypot(dx,dy))
    return np.array(out)

def independent_levels(horizontal,models,labels):
    r=independent_ranges(horizontal); results={}
    for f,mod in models.items():
        indices=[min(range(len(mod['depths'])),key=lambda i:abs(float(mod['depths'][i])-z)) for z in labels]
        receiver=min(range(len(mod['depths'])),key=lambda i:abs(float(mod['depths'][i])-200.))
        total=np.zeros((len(labels),len(r)),complex)
        compensation=np.zeros_like(total)
        for m,k in enumerate(mod['k']):
            factor=np.sqrt((2*math.pi)/(float(k.real)*r))*np.exp(-1j*float(k.real)*r+float(k.imag)*r-1j*math.pi/4)
            weight=mod['phi'][indices,m]*mod['phi'][receiver,m]
            term=weight[:,None]*factor
            increment=term-compensation; updated=total+increment
            compensation=(updated-total)-increment; total=updated
        results[f]=20*np.log10(np.maximum(np.abs(total),1e-30))
    return results

def independent_center(raw):
    centered={}
    for f,levels in raw.items():
        out=levels.copy()
        for row in out:
            for start,end in ((0,61),(61,121)):
                avg=math.fsum(float(x) for x in row[start:end])/(end-start)
                row[start:end]-=avg
        centered[f]=out
    return centered

def independent_score(horizontal,observed,models,labels,frequencies):
    levels=independent_center(independent_levels(horizontal,models,labels))
    return np.array([math.sqrt(math.fsum((float(levels[f][i,t])-float(observed[f][t]))**2 for f in frequencies for t in range(121))/(121*len(frequencies))) for i in range(len(labels))])
