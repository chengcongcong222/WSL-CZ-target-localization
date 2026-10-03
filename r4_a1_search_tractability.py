"""Pre-run deterministic search components; no scientific case driver is enabled."""
from __future__ import annotations
from dataclasses import dataclass
from pathlib import Path
import hashlib
import json
import time
import numpy as np
import pandas as pd
from scipy.optimize import shgo,direct
from scipy.spatial import cKDTree
import r4_a1_fix2_coverage as inherited

L=inherited.L
ROOT=Path(__file__).resolve().parent
OUT=ROOT/'results/R4_A1_SEARCH_TRACTABILITY_AUDIT'
BASELINE='7ab24845e6e1551b75287fefb1ab662e92b395b8'
R_BOUNDS=(45000.,60000.)
U_BOUNDS=(.94,3.)
STATE_AGREEMENT=(1e-5,.001,1e-4,.01)
DEDUP=(.005,.001,.005,.05)


class ObjectiveBudgetExceeded(RuntimeError):
    pass


@dataclass(frozen=True)
class Prepared:
    state: np.ndarray
    ranges: np.ndarray
    bearing_cost: float
    feasible: bool


@dataclass(frozen=True)
class Evaluation:
    value: float
    J_exact: float | None
    state: np.ndarray
    bearing_cost: float
    depth_label_m: float
    feasible: bool
    exact_forward_calls: int


class ExactObjective:
    """Observation-conditioned two-coordinate scalar objective for one branch."""
    def __init__(self,models,bearing_rad,relative_tl,cutoff,depth_label_m):
        if depth_label_m not in L.PROFILE:raise ValueError('Unknown inherited depth label')
        self.models=models
        self.bearing=np.array(bearing_rad,dtype=float,copy=True)
        self.relative_tl=np.array(relative_tl,dtype=float,copy=True)
        if self.bearing.shape!=(121,) or self.relative_tl.shape!=(3,121):raise ValueError('Observation shapes')
        if not np.isfinite(self.bearing).all() or not np.isfinite(self.relative_tl).all():raise ValueError('Nonfinite observation')
        if not np.isfinite(cutoff) or cutoff<0:raise ValueError('Invalid cutoff')
        self.bearing.setflags(write=False);self.relative_tl.setflags(write=False)
        self.cutoff=float(cutoff);self.depth_label_m=float(depth_label_m)

    def prepare(self,q):
        q=np.asarray(q,dtype=float)
        if q.shape!=(2,) or not np.isfinite(q).all():raise ValueError('Coordinate shape/finite check')
        r=R_BOUNDS[0]+(R_BOUNDS[1]-R_BOUNDS[0])*q[0]
        u=U_BOUNDS[0]+(U_BOUNDS[1]-U_BOUNDS[0])*q[1]
        s=inherited.radial_profile([r],[u],self.bearing)[0]
        bearing,ranges=L.geometry(s)
        residual=(bearing[0]-self.bearing+np.pi)%(2*np.pi)-np.pi
        cost=float(residual@residual)
        inside=np.all(q>=0) and np.all(q<=1) and np.all(s>=inherited.prior.LOW) and np.all(s<=inherited.prior.HIGH)
        return Prepared(s,ranges,float(cost),bool(inside and cost<=self.cutoff+1e-14))

    def key(self,p):return (self.depth_label_m,p.state.astype('<f8').tobytes())

    def evaluate(self,p):
        if not p.feasible:
            # A deterministic feasibility penalty is never a Gate-bearing score.
            return Evaluation(1e6+p.bearing_cost,None,p.state,p.bearing_cost,self.depth_label_m,False,0)
        features=L.direct_features(self.models,p.ranges,[self.depth_label_m])[0,0]
        score=float(np.sqrt(np.mean((features-self.relative_tl)**2)))
        return Evaluation(score,score,p.state,p.bearing_cost,self.depth_label_m,True,1)


class BudgetObjective:
    """Admitted requests are capped before preparation, cache and forward work."""
    def __init__(self,objective,limit,cache_enabled=True):
        if isinstance(limit,bool) or not isinstance(limit,int) or limit<1:raise ValueError('Positive integer hard cap required')
        self.objective=objective;self.limit=limit;self.cache_enabled=bool(cache_enabled)
        self.n_objective_requests=0;self.n_objective_attempts=0;self.n_blocked_requests=0
        self.n_exact_forward_evaluations=0;self.n_cache_hits=0
        self.hard_budget_triggered=False;self.cache={};self.seen=set();self.ledger=[]

    @property
    def n_unique_states_evaluated(self):return len(self.seen)

    def __call__(self,q):
        self.n_objective_attempts+=1
        if self.n_objective_requests>=self.limit:
            self.hard_budget_triggered=True;self.n_blocked_requests+=1
            self.ledger.append(dict(attempt=self.n_objective_attempts,admitted=False,request=self.n_objective_requests,cache_hit=False,exact_forward_calls=0))
            raise ObjectiveBudgetExceeded(f'admitted request cap {self.limit} reached')
        self.n_objective_requests+=1
        p=self.objective.prepare(q);key=self.objective.key(p)
        hit=self.cache_enabled and key in self.cache
        if hit:
            self.n_cache_hits+=1;e=self.cache[key]
        else:
            e=self.objective.evaluate(p);self.n_exact_forward_evaluations+=e.exact_forward_calls
            self.seen.add(key);self.cache[key]=e
        self.ledger.append(dict(attempt=self.n_objective_attempts,admitted=True,request=self.n_objective_requests,cache_hit=bool(hit),
            exact_forward_calls=0 if hit else e.exact_forward_calls,state_sha256=hashlib.sha256(key[1]).hexdigest(),
            q0=float(q[0]),q1=float(q[1]),J_exact=e.J_exact,feasible=e.feasible,depth_label_m=e.depth_label_m,
            **dict(zip(L.AXES,map(float,e.state))),bearing_cost=e.bearing_cost))
        return e.value

    def counters(self):
        return dict(n_objective_requests=self.n_objective_requests,n_objective_attempts=self.n_objective_attempts,
            n_blocked_requests=self.n_blocked_requests,n_exact_forward_evaluations=self.n_exact_forward_evaluations,
            n_unique_states_evaluated=self.n_unique_states_evaluated,n_cache_hits=self.n_cache_hits,
            hard_budget_limit=self.limit,hard_budget_triggered=self.hard_budget_triggered)


def exception_chain(error):
    result=[];seen=set()
    while error is not None and id(error) not in seen:
        seen.add(id(error));result.append(type(error).__name__)
        error=error.__cause__ if error.__cause__ is not None else error.__context__
    return result


def run_solver(name,controlled,options):
    """Fresh optimizer state per call; preserves native failures and causal chain."""
    started=time.perf_counter();native=None;chain=[]
    try:
        if name=='SHGO':
            native=shgo(controlled,[(0.,1.),(0.,1.)],**options)
        elif name=='DIRECT':
            native=direct(controlled,[(0.,1.),(0.,1.)],**options)
        else:raise ValueError('Only SHGO and DIRECT are authorized')
        reason='SOLVER_RETURN';message=str(native.message);success=bool(native.success)
    except (ObjectiveBudgetExceeded,SystemError) as error:
        chain=exception_chain(error)
        # SciPy DIRECT's C boundary retains the original exception as its cause.
        # This native chain is disclosed; an unrelated SystemError is re-raised.
        if 'ObjectiveBudgetExceeded' not in chain or not controlled.hard_budget_triggered:raise
        reason='HARD_BUDGET_EXCEEDED';message=repr(error);success=False
    supported=not controlled.hard_budget_triggered or (controlled.n_objective_requests==controlled.limit and controlled.n_blocked_requests==1 and 'ObjectiveBudgetExceeded' in chain)
    if not supported:reason='HARD_BUDGET_ENFORCEMENT_FAILED';success=False
    return native,dict(solver=name,termination_reason=reason,solver_success=success,solver_message=message,
        exception_chain=json.dumps(chain),hard_budget_enforcement_supported=supported,elapsed_seconds=time.perf_counter()-started,**controlled.counters())


def candidate_catalog(controlled):
    entries=[e for e in controlled.cache.values() if e.feasible and e.J_exact is not None]
    entries.sort(key=lambda e:e.J_exact);selected=[]
    if entries:
        states=np.array([e.state for e in entries]);scale=np.array(DEDUP)
        tree=cKDTree(states/scale);removed=np.zeros(len(entries),bool)
        for i,e in enumerate(entries):
            if removed[i]:continue
            selected.append(e)
            neighbors=np.array(tree.query_ball_point(states[i]/scale,1+1e-9,p=np.inf),dtype=int)
            close=np.all(abs(states[neighbors]-states[i])<=scale,axis=1)
            removed[neighbors[close]]=True
    rows=[dict(candidate_kind='EVALUATED_WITNESS_NOT_CERTIFIED_LOCAL_MINIMUM',J_exact=e.J_exact,
        z_star_label_m=e.depth_label_m,bearing_cost=e.bearing_cost,**dict(zip(L.AXES,map(float,e.state)))) for e in selected]
    return pd.DataFrame(rows,columns=['candidate_kind','J_exact','z_star_label_m','bearing_cost',*L.AXES])


def verify_catalog(models,bearing_rad,relative_tl,cutoff,catalog):
    rows=[]
    for row in catalog.itertuples():
        s=np.array([getattr(row,a) for a in L.AXES]);b,r=L.geometry(s)
        residual=(b[0]-bearing_rad+np.pi)%(2*np.pi)-np.pi;cost=float(residual@residual)
        features=L.direct_features(models,r,[row.z_star_label_m])[0,0]
        j=float(np.sqrt(np.mean((features-relative_tl)**2)))
        if abs(j-row.J_exact)>1e-7 or abs(cost-row.bearing_cost)>1e-12 or cost>cutoff+1e-14:raise AssertionError('Candidate verification mismatch')
        rows.append(dict(recomputed_J_exact=j,recomputed_bearing_cost=cost,pass_check=True))
    return pd.DataFrame(rows)


def save(name,frame):
    target=OUT/name;target.parent.mkdir(parents=True,exist_ok=True)
    frame.to_csv(target,index=False,float_format='%.17g',lineterminator='\n')


def json_write(name,value):
    OUT.mkdir(exist_ok=True);(OUT/name).write_bytes((json.dumps(value,indent=2,allow_nan=False)+'\n').encode())


if __name__=='__main__':
    raise SystemExit('PRE_RUN_CHECKPOINT_ONLY: no development driver enabled; research-lead approval required')
