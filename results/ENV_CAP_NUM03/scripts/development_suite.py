"""NUM03 synthetic-only development suite. Does not read archived environmental SSPs."""
from integrator import *
from independent_reference import Reference,mp_reference
from scipy.interpolate import PchipInterpolator
import time,sys,shutil
version=sys.argv[1];reason=sys.argv[2];D=O/'development'/version
if D.exists():raise RuntimeError('Development attempt ID must be unique')
D.mkdir();start=time.time();records=[]
for name in ['integrator.py','independent_reference.py','development_suite.py']:(D/name).write_bytes((O/'scripts'/name).read_bytes())
js(D/'START.json',{'reason':reason,'time_unix':start,'actual_SCS_profiles_used':False,'code':{n:sha(D/n) for n in ['integrator.py','independent_reference.py','development_suite.py']}})
def test(name,fn):
 try:
  info=fn();records.append(dict(check=name,PASS=bool(info.pop('PASS')),**info))
 except Exception as e:records.append({'check':name,'PASS':False,'error':repr(e)})
profiles=[([0.,500.,1500.,2500.,3500.,4500.],[1515.,1495.,1500.,1510.,1525.,1540.],200.,2500.),([0.,600.,1600.,2600.,3600.],[1450.,1453.,1462.,1479.,1500.],200.,1600.)]
offsets=[-25.,-1.,-.1,-.01,-.001,-.000001,0.,.000001,.001,.0038,.0083,.0159,.01,.1,1.,25.]
for pi,(z,c,zs,node) in enumerate(profiles):
 for kind in ['C_LINEAR','PCHIP']:
  p=Profile(z,c,kind);ref=Reference(z,c,kind)
  for delta in offsets:
   target=node+delta;cp=float(p.f(target));angle=math.degrees(math.acos(float(p.f(zs))/cp))
   def run(p=p,ref=ref,angle=angle,target=target,delta=delta,pi=pi,kind=kind):
    v=p.calc(zs,angle);q=ref.calc(zs,angle);difference=abs(v['range_m']-q['range_m']);res={'PASS':difference<=.001 and abs(v['turn_z_m']-q['turn_z_m'])<=1e-5 and v['bottom_encountered']==q['bottom_encountered'],'delta_m':difference,'turn_delta_m':abs(v['turn_z_m']-q['turn_z_m']),'nominal_turn_node_offset_m':delta,'algorithm':v['algorithm'],'Gauss_difference_m':v['G32_G64_difference_m'],'estimated_error_m':v['error_estimate_m']}
    if delta in [-.000001,.0038,.0083,.0159,25.]:
     mpv=mp_reference(z,c,zs,angle,kind,60);res['mp_delta_m']=abs(v['range_m']-mpv['range_m']);res['mp_error_estimate_m']=mpv['quad_error_estimate_m'];res['PASS']=res['PASS'] and res['mp_delta_m']<=.001 and mpv['quad_error_estimate_m']<=1e-6
    return res
   test(f'profile{pi}_{kind}_node_offset_{delta}',run)
# Near-bottom but unturned branches, and nearly horizontal launches, are explicitly tested too.
z,c,zs,node=profiles[1]
for kind in ['C_LINEAR','PCHIP']:
 p=Profile(z,c,kind);ref=Reference(z,c,kind)
 for gap in [1e-8,1e-6,1e-4,.01,1.]:
  cp=c[-1]+gap;angle=math.degrees(math.acos(float(p.f(zs))/cp))
  def run(angle=angle,p=p,ref=ref):
   v=p.calc(zs,angle);q=ref.calc(zs,angle);mpv=mp_reference(z,c,zs,angle,kind,60);d=abs(v['range_m']-q['range_m']);dm=abs(v['range_m']-mpv['range_m']);return {'PASS':d<.001 and dm<.001 and v['bottom_encountered'],'delta_m':d,'mp_delta_m':dm,'algorithm':v['algorithm'],'Gauss_difference_m':v['G32_G64_difference_m']}
  test(kind+f'_near_bottom_cp_gap{gap}',run)
 for a in [.1,.4,1.,3.,5.,10.,20.]:
  def run(a=a,p=p,ref=ref):
   v=p.calc(zs,a);q=ref.calc(zs,a);d=abs(v['range_m']-q['range_m']);return {'PASS':d<.001,'delta_m':d,'algorithm':v['algorithm']}
  test(kind+f'_launch_{a}',run)
for kind in ['C_LINEAR','PCHIP']:
 def run(kind=kind):
  p=Profile([0.,1000.,2000.],[1500.,1600.,1700.],kind);v=p.calc(100,10);cp=1510/math.cos(math.radians(10));exact=2*cp*math.sin(math.radians(10))/.1;return {'PASS':abs(v['range_m']-exact)<1e-5,'delta_m':abs(v['range_m']-exact)}
 test(kind+'_exact_linear_gradient',run)
 def flat(kind=kind):
  p=Profile([0.,1000.,2000.],[1500.,1500.,1500.],kind);v=p.calc(100,10);exact=2*1900/math.tan(math.radians(10));return {'PASS':abs(v['range_m']-exact)<1e-5 and v['bottom_encountered'],'delta_m':abs(v['range_m']-exact)}
 test(kind+'_constant_speed_boundary',flat)
def critical():
 p=Profile([0.,1000.,2000.,3000.],[1500.,1520.,1480.,1550.],'PCHIP')
 try:p.integrate(200,1000,1520.,True);ok=False
 except ArithmeticError:ok=True
 return {'PASS':ok}
test('PCHIP_exact_nonsimple_turn_rejected',critical)
res={'PASS':all(a['PASS'] for a in records),'version':version,'reason':reason,'wall_s':time.time()-start,'actual_environment_science':False,'checks':len(records),'FAIL_count':sum(not a['PASS'] for a in records),'records':records};js(D/'RESULTS.json',res)
with (O/'DEVELOPMENT_LOG.jsonl').open('a',encoding='utf-8') as f:f.write(json.dumps({k:v for k,v in res.items() if k!='records'},ensure_ascii=False)+'\n')
print(json.dumps({k:v for k,v in res.items() if k!='records'},ensure_ascii=False));print(json.dumps([a for a in records if not a['PASS']],ensure_ascii=False,indent=2));raise SystemExit(0 if res['PASS'] else 2)
