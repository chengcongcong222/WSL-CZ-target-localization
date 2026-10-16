"""Synthetic controls before Commit A; never reads the nine scientific profiles."""
from physics import *
from cold import IndependentProfile
import py_compile
checks=[]
def ck(n,b,d=None):checks.append({'check':n,'PASS':bool(b),'detail':d})
for f in ['physics.py','analysis.py','cold.py','preflight.py','run_once.py']:py_compile.compile(str(O/'scripts'/f),cfile=str(W/(f+'.pyc')),doraise=True)
z=np.array([0.,1000.,2000.]);c=np.array([1500.,1600.,1700.]);zs=100.;a=10.;cp=1510/math.cos(math.radians(a));expected=2*cp*math.sin(math.radians(a))/.1/1000
for method in ['C_LINEAR','PCHIP']:
 p=Profile(z,c,method);v=p.evaluate(zs,a,64);ck(method+'_linear_gradient_exact_range',abs(v['deep_return_km']-expected)<1e-8)
 ref=IndependentProfile(z,c,method).calc(zs,a);ck(method+'_independent_adaptive_gradient',abs(v['deep_return_km']-ref[0])<1e-7)
 met,rx,yy=p.ode(zs,a,100);ck(method+'_ODE_gradient',met['success'] and met['returned'] and abs(met['range_km']-expected)<.005,met['range_km'])
 flat=Profile([0.,1000.,2000.],[1500.,1500.,1500.],method).evaluate(100,10,64);ck(method+'_flat_no_refracted_return',flat['deep_bottom_bounce'] and flat['upper_surface_bounce'])
z=np.array([0.,400.,1200.,2000.]);c=np.array([1500.,1490.,1510.,1550.]);grid=np.linspace(0,2000,1001)
for method in ['C_LINEAR','PCHIP']:
 p=Profile(z,c,method);ip=IndependentProfile(z,c,method);err=max(abs(float(p.v(x))-ip.val(float(x))) for x in grid);ck(method+'_own_interpolant_dense',err<1e-10,err)
 for aa in [3.,5.,8.]:
  v=p.evaluate(100,aa,64);ref=ip.calc(100,aa);ck(method+f'_turn_and_quad_a{aa}',abs(v['deep_return_km']-ref[0])<.001 and abs(v['upper_return_km']-ref[1])<.001)
 p32=p.evaluate(100,5,32);p64=p.evaluate(100,5,64);ck(method+'_order_comparison',abs(p32['deep_return_km']-p64['deep_return_km'])<.005)
 if method=='PCHIP':
  for i in range(len(z)-1):
   vals=p.v(np.linspace(z[i],z[i+1],101));ck('PCHIP_no_overshoot_segment'+str(i),vals.min()>=min(c[i:i+2])-1e-10 and vals.max()<=max(c[i:i+2])+1e-10)
  ck('PCHIP_C1_jumps',max(abs(float(p.df(x-1e-6))-float(p.df(x+1e-6))) for x in z[1:-1])<1e-7)
# Independent ray equation through a derivative discontinuity; synthetic, not a WOA case.
p=Profile([0.,500.,1000.,2000.],[1500.,1450.,1495.,1550.],'C_LINEAR');v=p.evaluate(200,5,64);met,rx,yy=p.ode(200,5,50);ck('ODE_piecewise_gradient',met['returned'] and abs(v['deep_return_km']-met['range_km'])<.01,met)
# Archive parser must label prefix contacts rather than whole-trace counters.
ray=W/'synthetic_prefix.ray';ray.write_text("test\n0\n0\n0\n0\n0\nrz\n5\n5 1 0\n0 100\n10000 500\n20000 100\n30000 0\n40000 100\n",encoding='ascii');rows=path_archive(ray,100,1000,10,50);ck('ray_prefix_not_whole_trace',len(rows)==2 and rows[0]['identity']=='ZERO_BOUNDARY' and rows[1]['identity']=='SURFACE_ONLY')
ck('cold_launcher_contract',"['completed']" in (O/'scripts/run_once.py').read_text(encoding='utf-8') and 'PRIMARY_OUTPUT_FREEZE.json' in (O/'scripts/run_once.py').read_text(encoding='utf-8'))
res={'PASS':all(a['PASS'] for a in checks),'controls':checks,'actual_nine_case_science_executed':False};js(O/'PREFLIGHT_IMPLEMENTATION_CONTROLS.json',res);print(json.dumps(res,ensure_ascii=False,indent=2));raise SystemExit(0 if res['PASS'] else 2)
