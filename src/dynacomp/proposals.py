"""Snapshot-based numerical proposals, rationalized before algebraic use."""
from fractions import Fraction
from .rational import matrix, transpose, mul, inverse, identity, encoded


def reduced_candidate(model, V):
    """Any independent rational lift gets an exact left inverse and residuals.

    This constructs a candidate, not a claim of invariant closure.
    """
    v = matrix(V)
    if len(v) != len(model.A):
        raise ValueError("Lift must have one row per full-state coordinate")
    vt = transpose(v)
    e = mul(inverse(mul(vt, v)), vt)
    assert mul(e, v) == identity(len(v[0]))
    return {"V": encoded(v), "E": encoded(e), "G": encoded(mul(mul(e, model.A), v)),
            "H": encoded(mul(e, model.B)), "D": encoded(mul(model.C, v)),
            "status": "candidate_requires_error_bound"}


def propose_from_snapshots(model, snapshots, rank, max_denominator=10000):
    """Optional POD proposal. Rationalization can alter rank and is checked."""
    import numpy as np
    x = np.asarray(snapshots, dtype=float)
    if x.ndim != 2 or x.shape[0] != len(model.A) or not np.isfinite(x).all():
        raise ValueError("Expected finite snapshots[n,samples]")
    if type(rank) is not int or not 1 <= rank <= min(x.shape) or max_denominator < 1:
        raise ValueError("Invalid rank or denominator budget")
    u, singular, _ = np.linalg.svd(x, full_matrices=False)
    v = [[Fraction(float(value)).limit_denominator(max_denominator) for value in row]
         for row in u[:, :rank]]
    result = reduced_candidate(model, v)
    result["proposal_method"] = "POD with rationalized basis"
    result["snapshot_singular_values"] = [float(value) for value in singular]
    return result
