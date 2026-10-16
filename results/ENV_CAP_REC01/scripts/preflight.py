"""Small synthetic implementation controls only, run before design commit."""
import math,json,py_compile
import numpy as np
import rec01_analysis as m
from scipy.ndimage import uniform_filter1d
from decimal import Decimal,localcontext
checks=[]
def ck(n,b):checks.append({'check':n,'PASS':bool(b)})
for f in ['rec01_analysis.py','rec01_cold.py','run_once.py']:py_compile.compile(str(m.O/'scripts'/f),cfile=str(m.W/(f+'.pyc')),doraise=True)
z=np.arange(100,300.01,2.5);p=np.ones((len(z),881),complex);e,w=m.positive_energy(p,z,200);ck('constant_nonnegative_mean',np.array_equal(w,np.ones(881)))
p[:]=0;e,w=m.positive_energy(p,z,200);ck('exact_zero_no_floor_energy',np.array_equal(w,np.zeros(881)))
p[:,400]=1+1j;e,w=m.positive_energy(p,z,200);ck('isolated_impulse_9point',np.all(w[396:405]==2/9) and np.count_nonzero(w)==9)
p[:]=complex(1e-20,1e-20);p[:,10]=1;e,w=m.positive_energy(p,z,200);ck('weak_positive_after_bright',w[500]>0 and w[500]<1e-35)
# Analytic constant speed segment integral and reversal symmetry.
z=np.array([100.,600.]);c=np.array([1500.,1500.]);a=.2;p=math.cos(a)/1500;ck('flat_C_linear_integral',abs(m.segment(z,c,p)-500/math.tan(a))<1e-9)
z=np.array([0.,500.,1000.]);c=np.array([1500.,1500.,1500.]);q=m.path(z,c,200,1000,10);ck('constant_speed_no_turn',q['deep_bottom_bounce'] and q['upper_surface_bounce'])
# Exact linear-gradient turning path: independent elementary primitive.
z=np.array([0.,1000.]);c=np.array([1500.,1600.]);a=10.;q=m.path(z,c,100,1000,a);p=math.cos(math.radians(a))/1510.;cp=1/p;ref=2*math.sqrt(1-(p*1510)**2)/(.1*p)/1000;ck('gradient_turn_integral',not q['deep_bottom_bounce'] and abs(q['deep_return_km']-ref)<1e-10)
for k,ar in m.arrays().items():
 x=np.array([50000,5000,215],dtype=np.int32);y=x.copy();y[1]*=-1;dd=m.mirror_geometry(x,ar);ck(k+'_float64_int32_input',dd.dtype==np.float64 and np.isfinite(dd).all() and np.array_equal(dd,m.mirror_geometry(y,ar)))
for s in [[float('nan'),0,0],[float('inf'),0,0],[1e10,0,0]]:
 try:m.mirror_geometry(s,m.arrays()['HLA8']);ok=False
 except ValueError:ok=True
 ck('invalid_coordinate_guard',ok)
# Independent Decimal scalar mirror distance, unlike the old int32 square.
with localcontext() as ctx:
 ctx.prec=60
 r=float((Decimal(50007)**2+Decimal(5000)**2+Decimal(15)**2).sqrt());dd=m.mirror_geometry([50000,5000,215],m.arrays()['HLA8'])[0];ck('decimal_reference',abs(dd-r)<=4*math.ulp(r))
res={'PASS':all(a['PASS'] for a in checks),'checks':checks,'actual_archived_case_scientific_evaluation':False};m.js(m.O/'PREFLIGHT_IMPLEMENTATION_CONTROLS.json',res);print(json.dumps(res,indent=2));raise SystemExit(0 if res['PASS'] else 2)
