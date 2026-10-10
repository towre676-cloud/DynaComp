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
    targets = [model.C]
    for _ in us:
        targets.append(mul(targets[-1], model.A))
    transition_norms = [norm_inf(mul(ca, residual)) for ca in targets]
    input_norms = [norm_inf(mul(ca, forcing)) for ca in targets]
    # e_t=A^t e0 + sum A^(t-1-j)(R z_j + S u_j).
    envelopes = []
    for t in range(len(us)+1):
        bound = vector_norm(mv(targets[t], e0))
        for j in range(t):
            lag = t-1-j
            bound += transition_norms[lag]*vector_norm(zs[j])
            bound += input_norms[lag]*vector_norm(us[j])
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


def absolute_apply(a, radii):
    """Componentwise envelope |a x| <= |a| radii, if |x| <= radii."""
    radii = tuple(q(v) for v in radii)
    if len(a[0]) != len(radii) or any(v < 0 for v in radii):
        raise ValueError("Nonnegative radius vector with matching dimension required")
    return tuple(sum(abs(v)*r for v, r in zip(row, radii)) for row in a)


def box_horizon_bound(model, V, E, G, H, initial_center, initial_radius, input_radii):
    """Uniform output bound for a full initial-state box and arbitrary bounded inputs.

    Reduced initial state is E x0. Every input sequence satisfying the given
    componentwise bounds is covered, not only retained simulation traces.
    """
    v, e, g, h = (matrix(a) for a in (V, E, G, H))
    n, r = len(model.A), len(g)
    if len(v) != n or len(v[0]) != r or len(e) != r or len(e[0]) != n or len(g[0]) != r:
        raise ValueError("Expected V[n,r], E[r,n], G[r,r]")
    if len(h) != r or len(h[0]) != len(model.B[0]):
        raise ValueError("Expected H[r,m]")
    if mul(e, v) != identity(r):
        raise ValueError("Candidate must satisfy EV=I")
    center = tuple(q(x) for x in initial_center)
    radius = tuple(q(x) for x in initial_radius)
    if len(center) != n or len(radius) != n or any(x < 0 for x in radius):
        raise ValueError("Initial box must have n coordinates and nonnegative radii")
    us = [tuple(q(x) for x in u) for u in input_radii]
    if any(len(u) != len(model.B[0]) or any(x < 0 for x in u) for u in us):
        raise ValueError("Invalid input-radius vectors")
    residual = sub(mul(model.A, v), mul(v, g))
    forcing = sub(model.B, mul(v, h))
    initial_map = sub(identity(n), mul(v, e))
    initial_center_error = mv(initial_map, center)
    initial_radius_error = absolute_apply(initial_map, radius)
    zcenter = mv(e, center)
    zradius = absolute_apply(e, radius)
    zs = [tuple(abs(c)+a for c, a in zip(zcenter, zradius))]
    for u in us:
        zs.append(addv(absolute_apply(g, zs[-1]), absolute_apply(h, u)))
    # Store C A^k directly: output space can be much smaller than state space.
    ca = [model.C]
    for _ in us:
        ca.append(mul(ca[-1], model.A))
    envelopes = []
    for t in range(len(us)+1):
        values = addv(tuple(abs(x) for x in mv(ca[t], initial_center_error)),
                      absolute_apply(ca[t], initial_radius_error))
        for j in range(t):
            values = addv(values, absolute_apply(mul(ca[t-1-j], residual), zs[j]))
            values = addv(values, absolute_apply(mul(ca[t-1-j], forcing), us[j]))
        envelopes.append(values)
    return {"schema": "dynacomp.box-bound/1", "arithmetic": "rational",
            "coverage": "all initial states in the box and all inputs within componentwise radii",
            "bound": [str(vector_norm(v)) for v in envelopes],
            "output_bounds": [[str(x) for x in v] for v in envelopes],
            "reduced_state_radii": [[str(x) for x in v] for v in zs],
            "transition_residual": encoded(residual), "input_residual": encoded(forcing)}
