#!/usr/bin/env python3
"""Non-scientific synthetic-only tests; never substitute for archive cold audit."""
import math, tempfile, struct
from pathlib import Path
import numpy as np
from t2_prec2_readonly import parse_mod, parse_prt,env_essentials,lossless_lambda_tridiagonal,richardson
from t2_prec2_cold_audit import prufer_root

def synthetic_files(base):
    rec=160;M=41;nt=20;total=(7+M)*rec+M*8
    b=bytearray(total)
    struct.pack_into('<i',b,0,rec//4)
    struct.pack_into('<4i',b,84,1,1,nt,nt)
    struct.pack_into('<i',b,5*rec,M)
    s=[]
    for j in range(M):
        k=0.821-0.00022*j
        c=np.complex64(complex(k,-1e-6))
        struct.pack_into('<ff',b,(7+M)*rec+8*j,float(c.real),float(c.imag))
        s.append((j+1,k))
    mod=base/'test.mod';mod.write_bytes(b)
    prt=base/'test.prt'
    with prt.open('w') as f:
        f.write('    I    k (1/m)            alpha (1/m)   Phase Speed (m/s) Group Speed (m/s)\n')
        for j,k in s[:30]:f.write(f'{j:5d}  {k:.10f}    -1.000000000E-06     1500.0     1500.0\n')
    env=base/'test.env'
    with env.open('w') as f:
        f.write("'Synthetic test'\n201.0\n1\n'CVWT'\n20000 0 5000\n")
        for d in range(0,5001,20):f.write(f'{d:.3f} 1500.0 0 1 0 0\n')
        f.write("'R' 0.0\n1500 2400\n60\n1\n200.\n1\n200.\nR\n")
    return mod,prt,env

def main():
    checks=0
    with tempfile.TemporaryDirectory() as td:
        mod,prt,env=synthetic_files(Path(td))
        k,meta=parse_mod(mod);r=parse_prt(prt);z,c,opt=env_essentials(env)
        assert meta['nmodes']==41 and len(r)==30 and len(z)==251; checks+=3
        assert all(np.float32(rr['k_re']).view(np.uint32)==np.float32(k[rr['index_1based']-1].real).view(np.uint32) for rr in r);checks+=1
        assert opt['options']=="'CVWT'" and c.min()==1500;checks+=2
        f=201; H=5000;k0=2*np.pi*f/1500
        for mode in [1,33,321,641,961]:
            lam=[lossless_lambda_tridiagonal(f,z,c,N,mode) for N in (20000,40000,80000,160000)]
            ktrue=math.sqrt(k0*k0-((mode-.5)*np.pi/H)**2)
            x=richardson(lam)
            assert abs(x['kR6']-ktrue)<2e-9,(mode,x['kR6']-ktrue)
            checks+=1
        kr=richardson([lossless_lambda_tridiagonal(f,z,c,N,321) for N in (20000,40000,80000,160000)])['kR6']
        kt,n=prufer_root(f,z,c,321,kr,.5,1e-12,120,__import__('time').monotonic())
        assert abs(kt-kr)<2e-9,(kt,kr);checks+=1
    print('SYNTHETIC ONLY PASS',checks,'checks')
if __name__=='__main__':main()
