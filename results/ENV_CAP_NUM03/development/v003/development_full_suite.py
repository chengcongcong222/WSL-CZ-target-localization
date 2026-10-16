"""Separate full synthetic-domain development suite, no archived data import."""
from integrator import *
from independent_reference import Reference,mp_reference
import time
D=O/'development/v003'
if D.exists():raise RuntimeError('Unique development attempt required')
D.mkdir();start=time.time();checks=[]
for name in ['integrator.py','independent_reference.py','development_full_suite.py']:(D/name).write_bytes((O/'scripts'/name).read_bytes())
js(D/'START.json',{'reason':'Fix the test driver only after v002 parse failure; integrator/reference unchanged. Full synthetic angular, MP60/80 and ODE tests.','time_unix':start,'actual_environment_science':False})
def test(name,fn):
 try:
  info=fn();checks.append(dict(check=name,PASS=bool(info.pop('PASS')),**info))
 except Exception as e:checks.append({'check':name,'PASS':False,'error':repr(e)})
profiles=[([0.,500.,1500.,2500.,3500.,4500.],[1515.,1495.,1500.,1510.,1525.,1540.],200.,2500.),([0.,600.,1600.,2600.,3600.],[1450.,1453.,1462.,1479.,1500.],200.,1600.)]
for pi,(z,c,zs,node) in enumerate(profiles):
 for kind in ['C_LINEAR','PCHIP']:
  p=Profile(z,c,kind);ref=Reference(z,c,kind)
  for ai,angle in enumerate(np.linspace(.1,20,4001)):
   def run(p=p,ref=ref,angle=float(angle),zs=zs):
    a=p.calc(zs,angle);b=ref.calc(zs,angle);d=abs(a['range_m']-b['range_m']);return {'PASS':d<=.001 and abs(a['turn_z_m']-b['turn_z_m'])<1e-5 and a['bottom_encountered']==b['bottom_encountered'],'delta_m':d,'algorithm':a['algorithm']}
   test(f'synthetic_full_domain_{pi}_{kind}_{ai}',run)
  for delta in [.0038,.0159]:
   cp=float(p.f(node+delta));angle=math.degrees(math.acos(float(p.f(zs))/cp))
   def run(p=p,z=z,c=c,zs=zs,angle=angle,kind=kind):
    a=mp_reference(z,c,zs,angle,kind,60);b=mp_reference(z,c,zs,angle,kind,80);v=p.calc(zs,angle);return {'PASS':abs(a['range_m']-b['range_m'])<1e-8 and abs(v['range_m']-b['range_m'])<.001,'MP60_MP80_delta_m':abs(a['range_m']-b['range_m']),'double_MP80_delta_m':abs(v['range_m']-b['range_m'])}
   test(f'cross_precision_{pi}_{kind}_{delta}',run)
  print('COMPLETE_SYNTHETIC',pi,kind,flush=True)
for kind in ['C_LINEAR','PCHIP']:
 p=Profile([0.,500.,1000.,2000.],[1500.,1450.,1495.,1550.],kind)
 def run(p=p):
  v=p.calc(200,5);a,rr,yy=p.ode(200,5,50);b,rr,yy=p.ode(200,5,25);return {'PASS':a['returned'] and b['returned'] and abs(a['range_m']-v['range_m'])<.1 and abs(a['range_m']-b['range_m'])<.1 and a['invariant_error']<1e-8,'ODE_Snell_delta_m':abs(a['range_m']-v['range_m']),'ODE_step_delta_m':abs(a['range_m']-b['range_m'])}
 test(kind+'_differential_ray_control',run)
res={'PASS':all(a['PASS'] for a in checks),'version':'v003','reason':'Test-driver correction only; integrator and independent reference unchanged since v001','wall_s':time.time()-start,'actual_environment_science':False,'checks':len(checks),'FAIL_count':sum(not a['PASS'] for a in checks),'records':checks};js(D/'RESULTS.json',res)
with (O/'DEVELOPMENT_LOG.jsonl').open('a',encoding='utf-8') as f:f.write(json.dumps({k:v for k,v in res.items() if k!='records'})+'\n')
print(json.dumps({k:v for k,v in res.items() if k!='records'}));print(json.dumps([a for a in checks if not a['PASS']],indent=2));raise SystemExit(0 if res['PASS'] else 2)
