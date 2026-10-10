"""Certified greedy enrichment of approximate linear reductions.

Start from the target directions, project the dynamics, and certify every
scenario with the exact rational envelope of ``bounds.finite_horizon_bound``.
While a certified bound exceeds the tolerance, add the full-state error
snapshot whose future is most visible to the target. Each accepted step is a
proof for the supplied scenarios, not a training-set fit. Because the lift
can grow to the full space, where all residuals vanish, the loop always ends.
"""
from .bounds import finite_horizon_bound
from .krylov import Echelon
from .rational import matrix, sub, mul, transpose, inverse, mv, addv, q, vector_norm, encoded, row_basis


def _project(model, columns):
    """Orthogonal projector over Q: W = (V^T V)^(-1) V^T, G = WAV, H = WB."""
    v = transpose(columns)
    w = mul(inverse(mul(columns, v)), columns)
    return v, w, mul(mul(w, model.A), v), mul(w, model.B)


def _trajectory(model, x0, inputs):
    xs = [tuple(q(x) for x in x0)]
    for u in inputs:
        xs.append(addv(mv(model.A, xs[-1]), mv(model.B, u)))
    return xs


def _visibility(model, e, horizon):
    """Target-weighted size of an error direction: sum_k |C A^k e| for k <= horizon."""
    total, x = q(0), e
    for _ in range(horizon+1):
        total += vector_norm(mv(model.C, x))
        x = mv(model.A, x)
    return total


def simplify_direction(e, denominator):
    """Scale to unit max-norm and round to small denominators.

    Soundness does not depend on which rational lift is used: the envelope is
    recomputed exactly for whatever V results. Rounding only controls growth.
    """
    if not denominator:
        return e
    scale = vector_norm(e)
    return tuple((x/scale).limit_denominator(denominator) for x in e)


def greedy_enrichment(model, scenarios, tolerance, max_dimension=None, lookahead=None, denominator=64):
    """Grow a lift until every scenario's certified output envelope meets tolerance.

    scenarios: [(x0, inputs)]. ``denominator`` bounds the entries of each added
    direction (None keeps exact snapshots). Returns the final reduced model, the per-round
    worst certified bounds, and status "certified" or "budget_exhausted".
    """
    tol = q(tolerance)
    n = len(model.A)
    budget = n if max_dimension is None else min(n, max_dimension)
    horizon = lookahead if lookahead is not None else min(n, max(len(u) for _, u in scenarios))
    columns = list(row_basis(model.C))
    if len(columns) > budget:
        raise ValueError("Budget is smaller than the number of independent targets")
    trajectories = [_trajectory(model, x0, inputs) for x0, inputs in scenarios]
    history = []
    while True:
        v, w, g, h = _project(model, tuple(columns))
        worst, certificates, candidates = q(0), [], []
        for (x0, inputs), xs in zip(scenarios, trajectories):
            z0 = mv(w, xs[0])
            cert = finite_horizon_bound(model, v, g, h, {"full": xs[0], "reduced": z0}, inputs)
            certificates.append(cert)
            worst = max([worst] + [q(b) for b in cert["bound"]])
            z = z0
            for t, x in enumerate(xs):
                error = tuple(a-b for a, b in zip(x, mv(v, z)))
                if any(error):
                    candidates.append(error)
                if t < len(inputs):
                    z = addv(mv(g, z), mv(h, inputs[t]))
        history.append({"dimension": len(columns), "worst_certified_bound": str(worst),
                        "approximately": f"{float(worst):.6g}"})
        if worst <= tol:
            status = "certified"
            break
        if len(columns) >= budget:
            status = "budget_exhausted"
            break
        span = Echelon()
        for c in columns:
            span.add(c)
        # Residual directions cover bounds driven by A V - V G when snapshots agree.
        candidates += list(transpose(sub(mul(model.A, v), mul(v, g))))
        ranked = sorted(candidates, key=lambda e: (_visibility(model, e, horizon), vector_norm(e)), reverse=True)
        chosen = next((s for s in (simplify_direction(e, denominator) for e in ranked)
                       if any(span.reduce(s))), None)
        if chosen is None:
            chosen = next(tuple(int(i == j) for i in range(n)) for j in range(n)
                          if any(span.reduce(tuple(int(i == j) for i in range(n)))))
        columns.append(chosen)
    return {"schema": "dynacomp.enrichment/1", "status": status, "tolerance": str(tol),
            "full_dimension": n, "reduced_dimension": len(columns),
            "V": encoded(v), "W": encoded(w), "G": encoded(g), "H": encoded(h),
            "A": encoded(model.A), "B": encoded(model.B), "C": encoded(model.C),
            "scenarios": [{"x0": [str(x) for x in x0], "inputs": encoded(inputs) if inputs else []}
                          for x0, inputs in scenarios],
            "history": history, "certificates": [c["bound"] for c in certificates],
            "claim": "rational output-error envelopes for the supplied scenarios only"}


def replay_enrichment(packet):
    """Recompute every envelope from the stored lift; no greedy search is rerun."""
    from .linear import LinearModel
    from .rational import identity
    if packet.get("schema") != "dynacomp.enrichment/1":
        raise ValueError("Unknown schema")
    model = LinearModel.from_rows(packet["A"], packet["B"], packet["C"])
    v, w, g, h = (matrix(packet[k]) for k in ("V", "W", "G", "H"))
    if mul(w, v) != identity(len(g)):
        raise ValueError("W is not a left inverse of V")
    if mul(mul(w, model.A), v) != g or mul(w, model.B) != h:
        raise ValueError("Reduced matrices are not the stated projection")
    tol, worst = q(packet["tolerance"]), q(0)
    if len(packet["certificates"]) != len(packet["scenarios"]):
        raise ValueError("One certificate per scenario required")
    for scenario, stored in zip(packet["scenarios"], packet["certificates"]):
        x0 = tuple(q(x) for x in scenario["x0"])
        cert = finite_horizon_bound(model, v, g, h, {"full": x0, "reduced": mv(w, x0)}, scenario["inputs"])
        if cert["bound"] != stored:
            raise ValueError("Stored envelope does not match recomputation")
        worst = max([worst] + [q(b) for b in stored])
    if packet["status"] == "certified" and worst > tol:
        raise ValueError("Certified status but an envelope exceeds tolerance")
    return True
