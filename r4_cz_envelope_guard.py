"""Runtime isolation shared by the estimator and the separate evaluator."""
from contextlib import contextmanager
from contextvars import ContextVar
from functools import wraps
_ACTIVE = ContextVar('cz_estimator_active', default=False)
@contextmanager
def estimator_scope():
    token = _ACTIVE.set(True)
    try:
        yield
    finally:
        _ACTIVE.reset(token)

def evaluation_only(fn):
    @wraps(fn)
    def guarded(*args, **kwargs):
        if _ACTIVE.get():
            raise RuntimeError('Evaluation helper called inside estimator')
        return fn(*args, **kwargs)
    return guarded
