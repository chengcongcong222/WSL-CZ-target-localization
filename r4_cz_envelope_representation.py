"""Fixed M0 WINDOW_SHAPE_Q2; real level input only."""
from __future__ import annotations
import numpy as np

FREQUENCIES = (201, 235, 283)
WINDOWS = (slice(0, 61), slice(61, 121))
TIMES = np.arange(121, dtype=float) * 10

def basis(n):
    if n not in (61, 60):
        raise ValueError('Only inherited windows are admitted')
    x = 2 * np.arange(n, dtype=float) - (n - 1)
    b1 = x / np.linalg.norm(x)
    q = x*x - np.mean(x*x)
    b2 = q / np.linalg.norm(q)
    return np.column_stack((b1, b2))

def feature(relative_level, *, frequencies=FREQUENCIES, times=TIMES):
    a = np.asarray(relative_level)
    if np.iscomplexobj(a) or a.dtype.kind not in 'fiu' or a.shape != (3, 121):
        raise ValueError('M0 requires finite real 3 x 121 levels')
    if tuple(frequencies) != FREQUENCIES or not np.array_equal(times, TIMES):
        raise ValueError('Frozen frequency/time metadata required')
    if not np.isfinite(a).all():
        raise ValueError('Nonfinite levels')
    a = a.astype(float)
    result = []
    for line in a:
        for section in WINDOWS:
            v = line[section]
            result.extend((v - v.mean()) @ basis(len(v)))
    return np.asarray(result, dtype=float)

def score(a, b):
    a, b = np.asarray(a), np.asarray(b)
    if a.shape != (12,) or b.shape != (12,) or not np.isfinite(a).all() or not np.isfinite(b).all():
        raise ValueError('12 finite coefficients required')
    return float(np.sqrt(np.sum((a-b)**2) / 6))
