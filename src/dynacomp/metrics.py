"""Exact weighted infinity-norm margins; heuristic discovery, exact acceptance."""
from .rational import matrix, q


def metric_margin(A, weights):
    a = matrix(A)
    weights = tuple(q(v) for v in weights)
    if len(a) != len(a[0]) or len(weights) != len(a) or any(w <= 0 for w in weights):
        raise ValueError("Need square A and strictly positive weights")
    # ||x||_w = max_i |x_i| / w_i.
    rho = max(sum(abs(value)*weights[j] for j, value in enumerate(row))/weights[i]
              for i, row in enumerate(a))
    return {"arithmetic": "rational", "weights": [str(w) for w in weights],
            "factor": str(rho), "margin": str(1-rho), "contractive": rho < 1}


def discover_metric(A, rounds=32):
    """Search positive diagonal metrics. Failure is not proof of instability."""
    a = matrix(A)
    if len(a) != len(a[0]) or rounds < 0:
        raise ValueError("Invalid square operator or search budget")
    weights = [q(1)]*len(a)
    best = metric_margin(a, weights)
    for _ in range(rounds):
        # Positive resolvent-like iteration on |A|, with an additive floor.
        weights = [1+sum(abs(v)*weights[j] for j, v in enumerate(row)) for row in a]
        candidate = metric_margin(a, weights)
        if q(candidate["factor"]) < q(best["factor"]):
            best = candidate
    return best
