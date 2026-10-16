from pathlib import Path
import sys,ast,json,importlib.util,numpy as np
O=Path.cwd()/'results/ENV_CAP_CZO04';sys.path.insert(0,str(O/'scripts'))
import model
sp=importlib.util.spec_from_file_location('indep',O/'scripts/cold.py');c=importlib.util.module_from_spec(sp);sp.loader.exec_module(c)
ck=[]
def check(n,b):ck.append(dict(check=n,PASS=bool(b)))
for p in (O/'scripts').glob('*.py'):ast.parse(p.read_text(encoding='utf-8-sig'));check('syntax/'+p.name,True)
for rho in [0,.9]:
 wh=model.weight(rho);a=np.arange(3*9*121).reshape(3,9,121)/100
 for cond in ['C0','C1','C2']:
  check(f'project-crosscheck/{rho}/{cond}',np.max(abs(model.project(a,wh,cond)-c.projected(a,wh,cond)))<1e-10)
for s in [[50,0,2,5],[47,-1,1.4,-7],[59,-.5,1.8,-11]]:
 a,b=model.geom(s);x,y=c.geometry(s);check('geometry/'+str(s),np.max(abs(a-x))<1e-12 and np.max(abs(b-y))<1e-8)
# No archive read here: toy nine-label tensor reduction and empty candidate behavior.
for n in [0,3]:
 a=np.zeros((n,9,121));tx=np.ones((9,121));px=a.reshape(-1,121);d=np.maximum((px*px).sum(1)[:,None]+(tx*tx).sum(1)[None,:]-2*px@tx.T,0).reshape(n,9,9)
 check('tensor-'+str(n),d.min(axis=1).T.shape==(9,n))
check('scenario-count',len(model.states())==36)
model.dump(O/'IMPLEMENTATION_PREFLIGHT_CONTROLS.json',dict(PASS=all(x['PASS'] for x in ck),checks=ck,scope='synthetic fields only; historical scenario list is metadata'))
print(json.dumps({'checks':len(ck),'PASS':all(x['PASS'] for x in ck)}));assert all(x['PASS'] for x in ck)
