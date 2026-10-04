"""Conservative continuous-cell support, immutable observations, strict accounting."""
from __future__ import annotations
from collections import deque
from dataclasses import dataclass, field
import math
from flint import arb, fmpq
import numpy as np
from r4_cz_envelope_guard import estimator_scope

DOMAIN_LOW = (45., -5., 1., -15.)
DOMAIN_HIGH = (60., 5., 3., 15.)
DEPTH_LABELS = tuple(range(150, 251, 5))
COUNTERS = ('n_cell_requests','n_forward_state_evaluations','n_profile_evaluations',
            'n_bound_evaluations','n_cache_hits','n_validation_evaluations')
COARSE_WIDTHS = (.25,.625,.125,1.875)
FINE_WIDTHS = tuple(x/2 for x in COARSE_WIDTHS)

def exact(value):
    if isinstance(value, arb):
        return value
    if isinstance(value, (int, np.integer)):
        return arb(int(value))
    n, d = float(value).as_integer_ratio()
    return arb(fmpq(n, d))

@dataclass(frozen=True)
class Cell:
    low: tuple
    high: tuple
    def __post_init__(self):
        lo, hi = tuple(map(float,self.low)), tuple(map(float,self.high))
        if len(lo)!=4 or len(hi)!=4 or not all(math.isfinite(a) and math.isfinite(b) and a<=b for a,b in zip(lo,hi)):
            raise ValueError('Invalid 4D cell')
        if not all(a<=x<=y<=b for a,x,y,b in zip(DOMAIN_LOW,lo,hi,DOMAIN_HIGH)):
            raise ValueError('Cell outside inherited domain')
        object.__setattr__(self,'low',lo);object.__setattr__(self,'high',hi)
    @property
    def widths(self): return tuple(b-a for a,b in zip(self.low,self.high))
    @property
    def boundary_censored(self):
        return any(a==x or b==y for a,b,x,y in zip(self.low,self.high,DOMAIN_LOW,DOMAIN_HIGH))
    def split(self, terminal_widths):
        ratios=[w/t for w,t in zip(self.widths,terminal_widths)]
        axis=max(range(4),key=lambda i:(ratios[i],-i))
        mid=(self.low[axis]+self.high[axis])/2
        if mid in (self.low[axis],self.high[axis]):
            raise ValueError('Cannot split at binary64 precision')
        hi=list(self.high);hi[axis]=mid
        lo=list(self.low);lo[axis]=mid
        return Cell(self.low,tuple(hi)),Cell(tuple(lo),self.high)

@dataclass(frozen=True)
class Observation:
    bearing_rad: np.ndarray
    coefficients: np.ndarray
    bearing_cutoff: float
    def __post_init__(self):
        for name,shape in [('bearing_rad',(121,)),('coefficients',(12,))]:
            a=np.asarray(getattr(self,name))
            if np.iscomplexobj(a) or a.dtype.kind not in 'fiu' or a.shape!=shape or not np.isfinite(a).all():
                raise ValueError('Invalid '+name)
            # immutable bytes backing prevents even setflags(write=True)
            a=np.frombuffer(a.astype(float).tobytes(),dtype=float).reshape(shape)
            object.__setattr__(self,name,a)
        if not math.isfinite(self.bearing_cutoff) or self.bearing_cutoff<0:
            raise ValueError('Invalid observation-side cutoff')

class BudgetExceeded(RuntimeError): pass
@dataclass
class HardBudget:
    limits: dict
    hard_budget_limit: int
    counts: dict = field(default_factory=lambda:dict.fromkeys(COUNTERS,0))
    attempts: dict = field(default_factory=lambda:dict.fromkeys(COUNTERS,0))
    blocked_transactions: int = 0
    def __post_init__(self):
        if set(self.limits)!=set(COUNTERS) or any(type(v)!=int or v<0 for v in self.limits.values()) or type(self.hard_budget_limit)!=int or self.hard_budget_limit<0:
            raise ValueError('All hard caps must be nonnegative integers')
        self.limits=dict(self.limits)
    def charge(self, **units):
        if any(k not in COUNTERS or type(v)!=int or v<0 for k,v in units.items()):
            raise ValueError('Invalid accounting transaction')
        for k,v in units.items():self.attempts[k]+=v
        if (sum(self.counts.values())+sum(units.values())>self.hard_budget_limit or
            any(self.counts[k]+v>self.limits[k] for k,v in units.items())):
            self.blocked_transactions+=1
            raise BudgetExceeded('HARD_CAP; transaction rejected before calculation')
        for k,v in units.items():self.counts[k]+=v
    def report(self):
        return {**self.counts,'hard_budget_limit':self.hard_budget_limit,'charged_total':sum(self.counts.values()),
                'attempts':dict(self.attempts),'blocked_transactions':self.blocked_transactions,'limits':dict(self.limits)}

def search_budget(cell_cap, baseline=False):
    limits=dict.fromkeys(COUNTERS,0)
    limits.update(n_cell_requests=cell_cap,n_cache_hits=cell_cap,n_bound_evaluations=cell_cap*(1 if baseline else 22))
    if not baseline:limits.update(n_forward_state_evaluations=21*cell_cap,n_profile_evaluations=21*cell_cap)
    return HardBudget(limits,cell_cap*(3 if baseline else 66))

@dataclass(frozen=True)
class Bounds:
    lower: arb
    upper: arb
    def valid(self):
        return self.lower.is_finite() and self.upper.is_finite() and bool(self.lower>=0) and bool(self.upper>=self.lower)

@dataclass(frozen=True)
class Certificate:
    cell: Cell
    bearing: Bounds
    feature: Bounds|None
    depth_labels: tuple
    method: str
    def verdict(self, observation, tau, baseline):
        # A missing/uncertain certificate never permits whole-cell rejection.
        if self.method not in ('ARB_WHOLE_CELL_V1','ANALYTIC_TOY_INTERVAL_V1') or not self.bearing.valid():
            return 'RETAINED_UNRESOLVED'
        bc=exact(observation.bearing_cutoff)
        if bool(self.bearing.lower>bc):return 'REJECTED_BY_CERTIFIED_BOUND'
        if baseline:
            return 'RETAINED_COMPATIBLE' if bool(self.bearing.upper<=bc) else 'RETAINED_UNRESOLVED'
        if self.depth_labels!=DEPTH_LABELS or self.feature is None or not self.feature.valid():
            return 'RETAINED_UNRESOLVED'
        if bool(self.feature.lower>exact(tau)):return 'REJECTED_BY_CERTIFIED_BOUND'
        return 'RETAINED_COMPATIBLE' if bool(self.bearing.upper<=bc) and bool(self.feature.upper<=exact(tau)) else 'RETAINED_UNRESOLVED'

@dataclass(frozen=True)
class CellRecord:
    cell: Cell
    status: str
    certificate: Certificate|None
    reason: str
    @property
    def flags(self):return ('BOUNDARY_CENSORED',) if self.cell.boundary_censored else ()

def adjacent(a,b):
    # Closed-face adjacency; overlap of positive measure on other axes.
    for axis in range(4):
        if a.high[axis]==b.low[axis] or b.high[axis]==a.low[axis]:
            if all(min(a.high[i],b.high[i])>max(a.low[i],b.low[i]) for i in range(4) if i!=axis):return True
    return False

def components(cells):
    cells=sorted(cells,key=lambda c:(c.low,c.high)); unseen=set(range(len(cells)));n=0
    while unseen:
        n+=1; queue=[min(unseen)];unseen.remove(queue[0])
        while queue:
            i=queue.pop()
            connected=[j for j in sorted(unseen) if adjacent(cells[i],cells[j])]
            unseen.difference_update(connected);queue.extend(connected)
    return n

def projection(cells):
    merged=[]
    for c in sorted(cells,key=lambda c:(c.low[0],c.high[0])):
        lo,hi=c.low[0],c.high[0]
        if merged and lo<=merged[-1][1]:merged[-1]=(merged[-1][0],max(hi,merged[-1][1]))
        else:merged.append((lo,hi))
    return tuple(merged)

@dataclass(frozen=True)
class SupportResult:
    records: tuple
    accounting: dict
    coverage_closed: bool
    termination_reason: str
    execution_valid: bool
    initial_domain: Cell
    @property
    def retained(self):return tuple(r.cell for r in self.records if r.status!='REJECTED_BY_CERTIFIED_BOUND')
    @property
    def range_intervals(self):return projection(self.retained)
    def metrics(self):
        iv=self.range_intervals
        return {'range_hull_km':[iv[0][0],iv[-1][1]] if iv else None,
                'hull_width_km':iv[-1][1]-iv[0][0] if iv else None,
                'union_length_km':sum(b-a for a,b in iv),'range_component_count':len(iv),
                'joint_broad_component_count':components(self.retained),'empty_support':not bool(iv),
                'coverage_closed':self.coverage_closed,'execution_valid':self.execution_valid,
                'unresolved_cells':sum(r.status=='RETAINED_UNRESOLVED' for r in self.records),
                'boundary_censored':any(r.flags for r in self.records if r.status!='REJECTED_BY_CERTIFIED_BOUND')}

class SupportEngine:
    def __init__(self, provider, observation, tau, terminal_widths, budget, *, baseline=False):
        if len(terminal_widths)!=4 or any(not math.isfinite(w) or w<=0 for w in terminal_widths) or not math.isfinite(tau) or tau<0:
            raise ValueError('Invalid frozen resolution/tolerance')
        self.provider,self.observation,self.tau=provider,observation,tau
        self.terminal_widths=tuple(terminal_widths);self.budget=budget;self.baseline=baseline
        self.cache={};self.queue=deque();self._ran=False
    def certificate(self, cell):
        with estimator_scope():
            self.budget.charge(n_cell_requests=1)
            if cell in self.cache:
                self.budget.charge(n_cache_hits=1)
                return self.cache[cell]
            cert=self.provider(cell,self.observation,self.budget,self.baseline)
            if not isinstance(cert,Certificate) or cert.cell!=cell:raise ValueError('Whole-cell certificate required')
            self.cache[cell]=cert
            return cert
    def run(self):
        with estimator_scope():return self._run()
    def _run(self):
        if self._ran or self.queue or self.cache:raise RuntimeError('Raw run requires fresh engine/cache/queue/counters')
        if any(self.budget.counts.values()):raise RuntimeError('Fresh accounting required')
        self._ran=True;root=Cell(DOMAIN_LOW,DOMAIN_HIGH);self.queue.append(root);records=[];reason='QUEUE_EXHAUSTED';valid=True
        while self.queue:
            cell=self.queue.popleft();cert=None
            try:
                cert=self.certificate(cell)
                status=cert.verdict(self.observation,self.tau,self.baseline)
            except BudgetExceeded:
                status='RETAINED_UNRESOLVED';reason='HARD_CAP'
                records.append(CellRecord(cell,status,None,reason))
                records.extend(CellRecord(c,status,None,reason) for c in self.queue);self.queue.clear();break
            except (ArithmeticError,ValueError):
                valid=False;status='RETAINED_UNRESOLVED';reason='INVALID_BOUND'
                records.append(CellRecord(cell,status,None,reason))
                records.extend(CellRecord(c,status,None,reason) for c in self.queue);self.queue.clear();break
            if status=='RETAINED_UNRESOLVED' and any(w>t for w,t in zip(cell.widths,self.terminal_widths)):
                self.queue.extend(cell.split(self.terminal_widths))
            else:records.append(CellRecord(cell,status,cert,'CERTIFIED' if status!='RETAINED_UNRESOLVED' else 'LOOSE_BOUND_AT_TERMINAL'))
        closed=valid and not any(r.status=='RETAINED_UNRESOLVED' for r in records)
        if reason=='QUEUE_EXHAUSTED' and not closed:reason='TERMINAL_UNRESOLVED'
        return SupportResult(tuple(records),self.budget.report(),closed,reason,valid,root)

def raw_pair(provider, observation, tau, budgets, *, baseline=False):
    # Each call constructs two independent full-domain engines; no inherited support.
    return tuple(SupportEngine(provider,observation,tau,widths,search_budget(budgets[label],baseline),baseline=baseline).run()
                 for label,widths in [('COARSE',COARSE_WIDTHS),('FINE',FINE_WIDTHS)])
