"""Final synthetic harness controls; imports routines but never executes actual-case main()."""
from integrator import *
from independent_reference import mp_reference
from cold import golden
from analysis import refine
from decimal import Decimal,localcontext
import py_compile,time
D=O/'development/v004';D.mkdir();start=time.time();checks=[]
for p in (O/'scripts').glob('*.py'):
 (D/p.name).write_bytes(p.read_bytes());py_compile.compile(str(p),cfile=str(W/(p.name+'.pyc')),doraise=True)
def ck(n,b,d=None):checks.append({'check':n,'PASS':bool(b),'detail':d})
class Fake:
 kind='PCHIP'
 def calc(self,zs,a):return {'range_m':50000+1000*(a-2.3456789)**2,'turn_z_m':2000+10*a,'angle_deg':a}
f=Fake();r=golden(f,200,2.,3.);ck('independent_golden_quadratic',abs(r['angle_deg']-2.3456789)<1e-6)
rows=[dict(f.calc(200,a),angle_deg=a) for a in [2.,2.35,3.]];r=refine(f,200,rows,1);ck('primary_bounded_refinement_quadratic',r['strict_interior_minimum'] and abs(r['angle_deg']-2.3456789)<1e-6)
for kind in ['C_LINEAR','PCHIP']:
 z=[0.,600.,1600.,2600.,3600.];c=[1450.,1453.,1462.,1479.,1500.];p=Profile(z,c,kind);cp=float(p.f(1600.0038));angle=math.degrees(math.acos(float(p.f(200))/cp));a=mp_reference(z,c,200,angle,kind,60);b=mp_reference(z,c,200,angle,kind,80)
 with localcontext() as ctx:ctx.prec=100;delta=float(abs(Decimal(a['range_decimal'])-Decimal(b['range_decimal'])))
 ck(kind+'_decimal_MP60_MP80',delta<1e-12,delta)
ck('cold_imports_no_main_integrator',not any(s in (O/'scripts/cold.py').read_text(encoding='utf-8') for s in ['from integrator import','import integrator','from analysis import','import analysis']))
ck('guard_seals_output_and_requires_completed',"['completed']" in (O/'scripts/run_once.py').read_text(encoding='utf-8-sig') and 'PRIMARY_OUTPUT_FREEZE.json' in (O/'scripts/run_once.py').read_text(encoding='utf-8-sig'))
r={'PASS':all(a['PASS'] for a in checks),'version':'v004','reason':'Final synthetic harness validation and Decimal cross-precision comparison; numerical integrator unchanged from v001','wall_s':time.time()-start,'actual_environment_science':False,'checks':len(checks),'FAIL_count':sum(not a['PASS'] for a in checks),'records':checks};js(D/'RESULTS.json',r)
with (O/'DEVELOPMENT_LOG.jsonl').open('a',encoding='utf-8') as f:f.write(json.dumps({k:v for k,v in r.items() if k!='records'})+'\n')
print(json.dumps(r,ensure_ascii=False,indent=2));raise SystemExit(0 if r['PASS'] else 2)
