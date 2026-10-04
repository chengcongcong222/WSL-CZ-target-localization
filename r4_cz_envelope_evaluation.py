"""Reference retention is evaluation only, after support is sealed."""
from r4_cz_envelope_guard import evaluation_only
@evaluation_only
def retention(result, reference_state):
    x = tuple(reference_state)
    return {'joint_retained': any(all(a <= v <= b for a,v,b in zip(c.low,x,c.high)) for c in result.retained),
            'range_retained': any(a <= x[0] <= b for a,b in result.range_intervals)}
