"""Exact rational output-error envelopes for approximate lifted dynamics.

Bounds are algebraic finite-horizon inequalities, not floating-point estimates.
They need no asymptotic stability assumption and include hidden initial state.
"""
from .rational import matrix, mul, sub, mv, identity, norm_inf, vector_norm, q, addv, encoded


def finite_horizon_bound(model, V, G, H, x0, inputs):
    v, g, h = matrix(V), matrix(G), matrix(H)
    n, r = len(model.A), len(g)
    if len(v) != n or len(v[0]) != r or len(g[0]) != r:
        raise ValueError("Expected V[n,r], G[r,r]")
    if len(h) != r or len(h[0]) != len(model.B[0]):
        raise ValueError("Expected H[r,m]")
    # Initial reduced state is supplied separately via projection by caller.
    if not isinstance(x0, dict) or set(x0) != {"full", "reduced"}:
        raise ValueError("x0 must contain full and reduced initial states")
    full0 = tuple(q(x) for x in x0["full"])
    z = tuple(q(x) for x in x0["reduced"])
    e0 = tuple(a-b for a, b in zip(full0, mv(v, z)))
    if len(full0) != n:
        raise ValueError("Wrong full initial dimension")
    residual = sub(mul(model.A, v), mul(v, g))
    forcing = sub(model.B, mul(v, h))
    us = [tuple(q(x) for x in u) for u in inputs]
    zs = [z]
    for u in us:
        zs.append(addv(mv(g, zs[-1]), mv(h, u)))
    powers = [identity(n)]
    for _ in us:
        powers.append(mul(powers[-1], model.A))
    # e_t=A^t e0 + sum A^(t-1-j)(R z_j + S u_j).
    envelopes = []
    for t in range(len(us)+1):
        bound = vector_norm(mv(mul(model.C, powers[t]), e0))
        for j in range(t):
            ca = mul(model.C, powers[t-1-j])
            bound += norm_inf(mul(ca, residual))*vector_norm(zs[j])
            bound += norm_inf(mul(ca, forcing))*vector_norm(us[j])
        envelopes.append(bound)
    return {"schema": "dynacomp.bound/1", "arithmetic": "rational",
            "bound": [str(x) for x in envelopes],
            "reduced_states": [[str(x) for x in z] for z in zs],
            "transition_residual": encoded(residual), "input_residual": encoded(forcing)}


def usable_horizon(bounds, tolerance):
    """Largest prefix whose every bound meets tolerance; -1 if t=0 fails."""
    tol = q(tolerance)
    if tol < 0:
        raise ValueError("Tolerance must be nonnegative")
    for t, b in enumerate(bounds):
        if q(b) > tol:
            return t-1
    return len(bounds)-1
