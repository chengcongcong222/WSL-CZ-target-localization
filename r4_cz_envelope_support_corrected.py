"""Corrected conservative outer partition; original ball certificates unchanged."""
from __future__ import annotations
from dataclasses import dataclass
import math
from r4_cz_envelope_support import (SupportEngine as HistoricalEngine,Cell,Certificate,BudgetExceeded,
    DOMAIN_LOW,DOMAIN_HIGH,DEPTH_LABELS,COARSE_WIDTHS,FINE_WIDTHS,search_budget,exact,projection,components)

RETAINED_STATUSES=('RETAINED_CERTIFIED_COMPATIBLE','RETAINED_POSSIBLE_AT_TERMINAL')
ALL_STATUSES=('REJECTED_BY_CERTIFIED_BOUND',*RETAINED_STATUSES,'RETAINED_UNRESOLVED_BUDGET','RETAINED_INVALID_BOUND')

def structural_analysis(terminal_widths):
    if len(terminal_widths)!=4 or any(not math.isfinite(w) or w<=0 for w in terminal_widths):raise ValueError('Invalid terminal widths')
    splits=[]
    for lo,hi,w in zip(DOMAIN_LOW,DOMAIN_HIGH,terminal_widths):
        n=0
        while (hi-lo)/(2**n)>w:n+=1
        splits.append(n)
    depth=sum(splits)
    return {'axis_splits':tuple(splits),'path_split_depth':depth,'structural_minimum_cell_requests':2*depth+1}

@dataclass(frozen=True)
class OuterRecord:
    cell: Cell
    status: str
    certificate: Certificate|None
    depth: int
    reason: str
    @property
    def compatibility_certified(self):return self.status=='RETAINED_CERTIFIED_COMPATIBLE'
    @property
    def flags(self):return ('BOUNDARY_CENSORED',) if self.cell.boundary_censored else ()

@dataclass(frozen=True)
class OuterResult:
    records: tuple
    accounting: dict
    domain_partition_closed: bool
    compatibility_certified: bool
    search_budget_closed: bool
    termination_reason: str
    execution_valid: bool
    initial_domain: Cell
    terminal_widths: tuple
    maximum_depth: int
    queue_peak: int
    estimated_storage_peak_bytes: int
    @property
    def retained(self):
        # Safe conservative diagnostics include pending/invalid regions too.
        return tuple(r.cell for r in self.records if r.status!='REJECTED_BY_CERTIFIED_BOUND')
    @property
    def range_intervals(self):return projection(self.retained)
    @property
    def scientific_cells(self):
        if not self.domain_partition_closed:raise RuntimeError('Unclosed partition cannot enter scientific Gate')
        return tuple(r.cell for r in self.records if r.status in RETAINED_STATUSES)
    def diagnostics(self):
        iv=self.range_intervals
        return {**self.accounting,'n_rejected':sum(r.status=='REJECTED_BY_CERTIFIED_BOUND' for r in self.records),
            'n_terminal_possible':sum(r.status=='RETAINED_POSSIBLE_AT_TERMINAL' for r in self.records),
            'n_certified_compatible':sum(r.status=='RETAINED_CERTIFIED_COMPATIBLE' for r in self.records),
            'n_budget_unresolved':sum(r.status=='RETAINED_UNRESOLVED_BUDGET' for r in self.records),
            'n_invalid_bound':sum(r.status=='RETAINED_INVALID_BOUND' for r in self.records),
            'maximum_depth':self.maximum_depth,'queue_peak':self.queue_peak,'estimated_storage_peak_bytes':self.estimated_storage_peak_bytes,
            'retained_range_hull':(iv[0][0],iv[-1][1]) if iv else None,
            'domain_partition_closed':self.domain_partition_closed,'search_budget_closed':self.search_budget_closed,
            'compatibility_certified':self.compatibility_certified,'execution_valid':self.execution_valid,'termination_reason':self.termination_reason}
    def scientific_metrics(self):
        cells=self.scientific_cells;iv=projection(cells)
        return {'range_hull':(iv[0][0],iv[-1][1]) if iv else None,'hull_width_km':iv[-1][1]-iv[0][0] if iv else None,
            'union_length_km':sum(b-a for a,b in iv),'range_component_count':len(iv),'joint_broad_component_count':components(cells),
            'empty_support':not bool(cells),'domain_partition_closed':True,'search_budget_closed':self.search_budget_closed}

class CorrectedEngine(HistoricalEngine):
    def __init__(self,*args,progress=None,storage_envelope_bytes=32*1024*1024,**kwargs):
        super().__init__(*args,**kwargs)
        self.progress=progress;self.storage_envelope_bytes=storage_envelope_bytes
        if type(storage_envelope_bytes)!=int or storage_envelope_bytes<8192:raise ValueError('Invalid storage cap')
    def _validate(self,cert):
        if cert.method not in ('ARB_WHOLE_CELL_V1','ANALYTIC_TOY_INTERVAL_V1') or not cert.bearing.valid():raise ValueError('Invalid whole-cell bearing certificate')
        if not self.baseline and not bool(cert.bearing.lower>exact(self.observation.bearing_cutoff)):
            if cert.feature is None or cert.depth_labels!=DEPTH_LABELS or not cert.feature.valid():raise ValueError('Incomplete or invalid feature certificate')
    def _run(self):
        if self._ran or self.queue or self.cache or any(self.budget.counts.values()):raise RuntimeError('Raw run requires fresh full-domain state')
        self._ran=True;root=Cell(DOMAIN_LOW,DOMAIN_HIGH);self.queue.append((root,0));records=[]
        reason='QUEUE_EXHAUSTED';valid=True;natural_exhaustion=False;maximum_depth=0;peak=1;storage_peak=8192
        while self.queue:
            storage=(len(self.cache)+len(records)+len(self.queue))*8192;storage_peak=max(storage_peak,storage)
            if storage>self.storage_envelope_bytes:
                reason='STORAGE_CAP'
                records.extend(OuterRecord(c,'RETAINED_UNRESOLVED_BUDGET',None,d,reason) for c,d in self.queue)
                self.queue.clear();break
            cell,depth=self.queue.popleft();maximum_depth=max(maximum_depth,depth)
            try:
                certificate=self.certificate(cell);self._validate(certificate)
                verdict=certificate.verdict(self.observation,self.tau,self.baseline)
            except BudgetExceeded:
                reason='HARD_CAP'
                records.append(OuterRecord(cell,'RETAINED_UNRESOLVED_BUDGET',None,depth,reason))
                records.extend(OuterRecord(c,'RETAINED_UNRESOLVED_BUDGET',None,d,reason) for c,d in self.queue);self.queue.clear();break
            except (ArithmeticError,ValueError):
                reason='INVALID_BOUND';valid=False
                records.append(OuterRecord(cell,'RETAINED_INVALID_BOUND',None,depth,reason))
                records.extend(OuterRecord(c,'RETAINED_INVALID_BOUND',None,d,reason) for c,d in self.queue);self.queue.clear();break
            terminal=all(w<=t for w,t in zip(cell.widths,self.terminal_widths))
            if verdict=='REJECTED_BY_CERTIFIED_BOUND':records.append(OuterRecord(cell,verdict,certificate,depth,'CERTIFIED_LOWER_BOUND'))
            elif not terminal:
                # All retained regions reach the declared resolution, even if upper-certified.
                self.queue.extend((c,depth+1) for c in cell.split(self.terminal_widths))
                peak=max(peak,len(self.queue))
            else:
                status='RETAINED_CERTIFIED_COMPATIBLE' if verdict=='RETAINED_COMPATIBLE' else 'RETAINED_POSSIBLE_AT_TERMINAL'
                records.append(OuterRecord(cell,status,certificate,depth,'TERMINAL_OUTER_SUPPORT'))
            if self.progress is not None:self.progress(self.budget.counts['n_cell_requests'],depth,len(self.queue),len(records))
        else:natural_exhaustion=True
        closed=valid and natural_exhaustion and all(r.status in ('REJECTED_BY_CERTIFIED_BOUND',*RETAINED_STATUSES) for r in records)
        retained=[r for r in records if r.status!='REJECTED_BY_CERTIFIED_BOUND']
        uniformly=bool(retained) and all(r.compatibility_certified for r in retained)
        return OuterResult(tuple(records),self.budget.report(),closed,uniformly,natural_exhaustion,reason,valid,root,
            self.terminal_widths,maximum_depth,peak,storage_peak)

def raw_pair(provider,observation,tau,caps,*,baseline=False):
    return tuple(CorrectedEngine(provider,observation,tau,widths,search_budget(caps[label],baseline),baseline=baseline).run()
        for label,widths in [('COARSE',COARSE_WIDTHS),('FINE',FINE_WIDTHS)])
