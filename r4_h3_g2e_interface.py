"""Truth-free H3-G2E statistical interface. No modal field, source or noise truth input."""
import numpy as np
from scipy.linalg import cho_factor, cho_solve

def covariance(samples):
    """Rows are proper complex sample vectors; R_ij=E[y_i conj(y_j)]."""
    return samples.T @ samples.conj() / len(samples)

def herm(a):
    return (a + a.conj().T) / 2

def joint_model(noise, response, powers):
    # response: channels x frequencies, powers: one positive power per frequency.
    n, f = response.shape
    a = np.zeros((n*f, f), complex)
    for j in range(f):
        a[np.arange(n)*f+j, j] = response[:, j]
    return herm(noise + (a * powers[None, :]) @ a.conj().T)

def likelihoods(samples, model):
    # Never invert the sample CSD. Only the positive-definite MODEL is factored.
    cf = cho_factor(model, lower=True, check_finite=True)
    ld = 2*np.log(np.diag(cf[0]).real).sum()
    r = covariance(samples)
    raw = len(samples)*ld + np.sum(samples.conj().T * cho_solve(cf, samples.T)).real
    csd = len(samples)*(ld + np.trace(cho_solve(cf, r)).real)
    return float(raw), float(csd)

def estimate_noise(background):
    r = herm(covariance(background))
    np.linalg.cholesky(r)  # no ridge, no eigenvalue clipping, no diagonalization
    return r

def response_block(samples, background, supplied_oracle=None, alpha=.01):
    """Returns a diagnostic direction only after a frozen conservative energy gate.
    Unknown q and fixed gains are not supplied. Output is gain-distorted response,
    not an estimate of the uncalibrated propagation vector.
    """
    n = samples.shape[1]
    r = herm(covariance(samples))
    rn = estimate_noise(background) if supplied_oracle is None else supplied_oracle
    cf = cho_factor(rn, lower=True)
    energy = float(np.trace(cho_solve(cf, r)).real)
    # Complex Wishart E[Rn_hat^-1]=B/(B-n) Rn^-1, B independent proper samples.
    # Markov's inequality gives conservative per-block null false-alarm <=alpha.
    factor = len(background)/(len(background)-n) if supplied_oracle is None else 1
    threshold = n*factor/alpha
    ev0, v0 = np.linalg.eigh(r)
    ev, v = np.linalg.eigh(herm(r-rn))
    detected = bool(energy > threshold and ev[-1] > 0)
    u0 = v0[:, -1]
    u = v[:, -1] if detected else None
    return {"before": u0, "after": u, "energy": energy, "threshold": threshold,
            "detected": detected, "largest_signal_eigenvalue": float(ev[-1]),
            "rn": rn, "sample_rank": int(np.linalg.matrix_rank(r, tol=np.linalg.norm(r,2)*1e-10))}

def projector(u):
    return np.outer(u, u.conj()) / np.vdot(u,u).real

def response_error(u, nominal):
    if u is None:
        return None
    return float(np.linalg.norm(projector(u)-projector(nominal), "fro"))
