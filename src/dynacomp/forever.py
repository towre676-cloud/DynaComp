"""All-time uniform output-error brackets for persistent input boxes.

Question: does a fixed reduced model keep |y_t - y_hat_t| <= tol for EVERY
t >= 0, every initial state in a box and every input sequence that stays in a
constant box forever? Finite-horizon certificates cannot answer this: a model
accepted for twenty steps can drift away at steady state.

Work on the exact observable quotient (Eo, Go, Ho, Do), whose outputs equal
the full model's, and stack w = (Eo x, z):
    Aw = diag(Go, G), Bw = [Ho; H], Cw = [Do, -D],  e_t = Cw w_t.
With rows R_t = Cw Aw^t, steady state w* = (I-Aw)^(-1) Bw a, e* = Cw w*,
d0 = lift(c) - w*, the exact supremum of |e_t,i| over the domain at time t is
    S_t,i = |e*_i + R_t,i d0| + sum_k |R_t,i lift_k| r_k + sum_(s<t) sum_k |R_s,i Bw_k| b_k,
attained by a sign-pattern witness. DynaComp encloses S_t from both sides for
t = 0..T with dyadic arithmetic whose rounding error is tracked in the dual
norm ||.||_(P^-1) of a certified quadratic metric P (rows contract there, so
rounding errors stay bounded). For t > T every term is bounded with
||R_t v|| <= kappa rho^t ||v||_P, giving
    S_t,i <= |e*_i| + sum_(s<T) |...| + kappa_i rho^T (||d0||_P + nu_x + beta/(1-rho)).
The result is a rigorous bracket [lower, upper] for sup_t sup_domain |e_t|.
If lower > tol the tolerance is REFUTED: some admissible scenario violates it.
"""
from fractions import Fraction
from .rational import matrix, mul, identity, sub, mv, q, encoded, inverse, transpose
from .sparse import compile_sparse
from .stability import quadratic_certificate, contraction_holds, positive_definite, sqrt_up
from .uniform import as_sparse, sharp_box_bound, SparseLinearModel


def _blockdiag(a, b):
    n, m = len(a), len(b)
    return tuple(tuple(a[i][j] if j < n else Fraction(0) for j in range(n+m)) for i in range(n)) + \
        tuple(tuple(Fraction(0) if j < n else b[i][j-n] for j in range(n+m)) for i in range(m))


def _quad_abs(p, m):
    """Upper bound of sqrt(v^T P v) over |v_i| <= m_i."""
    k = len(m)
    return sqrt_up(sum(abs(p[i][j])*m[i]*m[j] for i in range(k) if m[i] for j in range(k) if m[j]))


def _pnorm(p, v):
    return sqrt_up(sum(v[i]*p[i][j]*v[j] for i in range(len(v)) if v[i] for j in range(len(v)) if v[j]))


def _round(x, bits):
    """Nearest integer to x * 2^bits."""
    y = q(x)*(1 << bits)
    return (y.numerator*2 + y.denominator)//(2*y.denominator)


_QUOTIENTS = {}


def _quotient(model):
    """Exact observable quotient and its certified metric, cached per model."""
    key = (model.operator.n, model.operator.entries, model.B, model.C)
    if key not in _QUOTIENTS:
        op = model.operator
        packet = compile_sparse(op, model.B, model.C, max_rank=op.n)
        parts = tuple(matrix(packet[k]) for k in ("E", "G", "H", "D"))
        _QUOTIENTS[key] = (parts, quadratic_certificate(parts[1]))
    return _QUOTIENTS[key]


def stacked_error_system(model, E, G, H, D):
    (eo, go, ho, do), _ = _quotient(model)
    return {"Go": go, "Eo": eo, "Ho": ho, "Do": do, "Aw": _blockdiag(go, G),
            "Bw": tuple(ho) + tuple(H), "Cw": tuple(tuple(r) + tuple(-x for x in s) for r, s in zip(do, D)),
            "lift": tuple(eo) + tuple(E)}


def all_time_bound(model, V, E, G, H, initial_center, initial_radius, input_center, input_radius,
                   tolerance=None, slack="1/1000000", max_steps=20000, bits=64):
    model = as_sparse(model)
    v, e, g, h = map(matrix, (V, E, G, H))
    if mul(e, v) != identity(len(g)):
        raise ValueError("Need EV = I")
    d = mul(model.C, v)
    c0, r0 = tuple(q(x) for x in initial_center), tuple(q(x) for x in initial_radius)
    a, b = tuple(q(x) for x in input_center), tuple(q(x) for x in input_radius)
    n, m = model.operator.n, len(model.B[0])
    if len(c0) != n or len(r0) != n or len(a) != m or len(b) != m or any(x < 0 for x in r0+b):
        raise ValueError("Boxes need matching dimensions and nonnegative radii")
    tol = None if tolerance is None else q(tolerance)
    slack = q(slack)
    s = stacked_error_system(model, e, g, h, d)
    aw, bw, cw, lift = s["Aw"], s["Bw"], s["Cw"], s["lift"]
    k, p_out = len(aw), len(cw)
    base = {"schema": "dynacomp.all-time-bound/1", "model": model.spec(),
            "V": encoded(v), "E": encoded(e), "G": encoded(g), "H": encoded(h),
            "initial_center": [str(x) for x in c0], "initial_radius": [str(x) for x in r0],
            "input_center": [str(x) for x in a], "input_radius": [str(x) for x in b],
            "tolerance": None if tol is None else str(tol), "slack": str(slack), "bits": bits,
            "quotient_order": len(s["Go"]), "reduced_order": len(g)}
    metrics = [_quotient(model)[1], quadratic_certificate(g)]
    if any(x["status"] != "certified" for x in metrics):
        return {**base, "status": "no_metric", "metrics": [x["status"] for x in metrics],
                "claim": "no all-time statement; finite-horizon certificates still apply"}
    P = _blockdiag(matrix(metrics[0]["P"]), matrix(metrics[1]["P"]))
    gamma = max(q(x["gamma"]) for x in metrics)
    rho = sqrt_up(gamma)
    if rho >= 1:
        return {**base, "status": "no_metric", "claim": "rounded contraction factor reached 1"}
    Pinv = inverse(P)
    w_star = mv(inverse(sub(identity(k), aw)), mv(bw, a))
    steady = mv(cw, w_star)
    d0 = tuple(x - y for x, y in zip(mv(lift, c0), w_star))
    kappa = [sqrt_up(sum(c[i]*Pinv[i][j]*c[j] for i in range(k) for j in range(k))) for c in cw]
    mu = _pnorm(P, d0)
    nu_x = _quad_abs(P, [sum(abs(x)*r for x, r in zip(row, r0)) for row in lift])
    beta = _quad_abs(P, [sum(abs(x)*r for x, r in zip(row, b)) for row in bw])
    lam = sqrt_up(sum(abs(x) for row in Pinv for x in row))  # ||eps||_(P^-1) <= lam * max|eps_k|
    tail_factor = mu + nu_x + beta/(1-rho)

    ulp = Fraction(1, 1 << bits)
    lift_cols = transpose(lift)
    bw_cols = transpose(bw)
    awi = [[_round(x, bits) for x in row] for row in aw]
    d0i = [_round(x, bits) for x in d0]
    lifti = [[_round(x, bits) for x in col] for col in lift_cols]
    bwi = [[_round(x, bits) for x in col] for col in bw_cols]
    lift_norm = [_pnorm(P, col) for col in lift_cols]
    bw_norm = [_pnorm(P, col) for col in bw_cols]
    d0_norm = mu
    active_x = [j for j, r in enumerate(r0) if r]
    active_u = [j for j, r in enumerate(b) if r]

    rows = [[_round(x, bits) for x in row] for row in cw]
    eta = [lam*ulp/2 for _ in cw]  # ||R_t - rows/2^bits||_(P^-1)
    rho_t = Fraction(1)
    s_hi, s_lo = [Fraction(0)]*p_out, [Fraction(0)]*p_out
    upper_so_far, lower_best = Fraction(0), (Fraction(0), 0, 0)
    two = 1 << (2*bits)
    history = []
    status, steps = None, 0
    use_state_metric = [False]*p_out
    for t in range(max_steps+1):
        steps = t
        tails = []
        for i, row in enumerate(rows):
            one_norm = Fraction(sum(abs(x) for x in row), 1 << bits)
            def dot(col):
                return Fraction(sum(x*y for x, y in zip(row, col) if x and y), two)
            def err(norm):
                # stored-row error (eta * ||v||_P) plus rounding of v itself
                return eta[i]*norm + one_norm*ulp/2
            mean = steady[i] + dot(d0i)
            me = err(d0_norm)
            mean_hi, mean_lo = abs(mean) + me, max(abs(mean) - me, 0)
            if use_state_metric[i]:
                state_hi, state_lo = kappa[i]*rho_t*nu_x, Fraction(0)
            else:
                state_hi = state_lo = Fraction(0)
                for j in active_x:
                    val, ev = abs(dot(lifti[j])), err(lift_norm[j])
                    state_hi += (val + ev)*r0[j]
                    state_lo += max(val - ev, 0)*r0[j]
                if kappa[i]*rho_t*nu_x <= slack/4:
                    use_state_metric[i] = True
            hi = mean_hi + state_hi + s_hi[i]
            lo = mean_lo + state_lo + s_lo[i]
            upper_so_far = max(upper_so_far, hi)
            if lo > lower_best[0]:
                lower_best = (lo, t, i)
            for j in active_u:
                val, ev = abs(dot(bwi[j])), err(bw_norm[j])
                s_hi[i] += (val + ev)*b[j]
                s_lo[i] += max(val - ev, 0)*b[j]
            tails.append(abs(steady[i]) + s_hi[i] + kappa[i]*rho_t*rho*tail_factor)
        upper = max(upper_so_far, max(tails))
        if t % 64 == 0:
            history.append({"t": t, "lower": f"{float(lower_best[0]):.9g}", "upper": f"{float(upper):.9g}"})
        if tol is not None and lower_best[0] > tol:
            status = "refuted"
            break
        if tol is not None and upper <= tol:
            status = "certified"
            break
        if max(kappa)*rho_t*rho*tail_factor <= slack and all(use_state_metric):
            status = "converged" if tol is None else "undecided_within_slack"
            break
        # Advance rows: R_(t+1) = R_t Aw, tracking rounding in the dual metric.
        new_rows = []
        for i, row in enumerate(rows):
            nxt = [sum(x*y[c] for x, y in zip(row, awi) if x) for c in range(k)]
            new_rows.append([(z + (1 << (bits-1))) >> bits for z in nxt])
            one_norm = Fraction(sum(abs(x) for x in row), 1 << bits)
            eta[i] = rho*eta[i] + lam*(ulp/2 + one_norm*ulp/2)
            eta[i] = Fraction(_round(eta[i], bits) + 1, 1 << bits)  # keep sizes bounded, round up
        rows = new_rows
        rho_t = Fraction(_round(rho_t*rho, bits) + 1, 1 << bits)
    else:
        status = "step_limit"
    lower, t_low, i_low = lower_best
    result = {**base, "status": status, "steps_examined": steps,
              "lower": str(lower), "upper": str(upper),
              "lower_approximately": f"{float(lower):.9g}", "upper_approximately": f"{float(upper):.9g}",
              "worst_time_lower_bound": t_low, "worst_output": i_low,
              "steady_state_error": [str(x) for x in steady],
              "steady_state_error_approximately": [f"{float(x):.9g}" for x in steady],
              "P": encoded(P), "gamma": str(gamma), "rho": str(rho),
              "history": history,
              "claim": ("rigorous bracket for sup over all t >= 0, all initial states in the box and all input "
                        "sequences remaining in the constant input box; the lower value is attained by a "
                        "sign-pattern scenario at the recorded time")}
    if status == "refuted":
        result["refutation"] = {"time": t_low, "output": i_low, "error_at_least": str(lower)}
    return result


def replay_all_time(packet):
    """Check the metric exactly, then recompute the deterministic bracket."""
    if packet.get("schema") != "dynacomp.all-time-bound/1" or packet.get("status") == "no_metric":
        raise ValueError("Need an all-time bracket packet")
    model = SparseLinearModel.from_spec(packet["model"])
    v, e, g, h = (matrix(packet[k]) for k in ("V", "E", "G", "H"))
    p, gamma = matrix(packet["P"]), q(packet["gamma"])
    if not positive_definite(p) or not 0 <= gamma < 1:
        raise ValueError("Invalid metric")
    aw = stacked_error_system(model, e, g, h, mul(model.C, v))["Aw"]
    if not contraction_holds(aw, p, gamma):
        raise ValueError("Stacked error dynamics do not contract in the stored metric")
    fresh = all_time_bound(model, v, e, g, h, packet["initial_center"], packet["initial_radius"],
                           packet["input_center"], packet["input_radius"], packet["tolerance"],
                           packet["slack"], bits=packet["bits"])
    for key in ("status", "lower", "upper", "P", "gamma", "steps_examined"):
        if fresh[key] != packet[key]:
            raise ValueError(f"Stored {key} does not match recomputation")
    return True


def moment_rows(model, count):
    """Rows C (I - A)^(-j), j = 1..count: left Krylov directions at z = 1.

    Including them in the encoder matches steady-state (DC) gains, which a
    finite-horizon enrichment does not see. Computed on the exact quotient.
    """
    model = as_sparse(model)
    (eo, go, _, do), _ = _quotient(model)
    resolvent = inverse(sub(identity(len(go)), go))
    rows, block = [], do
    for _ in range(count):
        block = mul(block, resolvent)
        rows.extend(mul(block, eo))
    return rows


def all_time_enrichment(model, initial_center, initial_radius, input_center, input_radius, tolerance,
                        max_rank=None, slack="1/1000000"):
    """Grow directions until the all-time bracket certifies the tolerance.

    Two candidate families compete at each round: the next steady-state moment
    row C (I - A)^(-j) and the next target pullback C A^k. The candidate with
    the smaller rigorous lower bound is kept (refutations are found early, so
    rejected candidates are cheap). This is a heuristic; the final packet is
    an ordinary all-time bracket, replayed independently of the search.
    """
    from .uniform import _project
    from .krylov import Echelon
    model = as_sparse(model)
    n = model.operator.n
    budget = n if max_rank is None else max_rank
    moments = moment_rows(model, budget)
    pulls, row = [], model.C
    for _ in range(budget):
        row = model.operator.pullback(row)
        pulls.extend(row)
    span, chosen = Echelon(), []
    for r in model.C:
        if span.add(r):
            chosen.append(r)
    queues = {"moment": list(moments), "pullback": pulls}
    history = []

    def evaluate(rows):
        v, e, g, h, d = _project(model, rows)
        return all_time_bound(model, v, e, g, h, initial_center, initial_radius, input_center,
                              input_radius, tolerance, slack)

    current = evaluate(chosen)
    history.append({"rank": len(chosen), "status": current["status"], "lower": current["lower_approximately"],
                    "upper": current["upper_approximately"]})
    while current["status"] != "certified" and len(chosen) < budget:
        trials = []
        for kind, queue in queues.items():
            while queue and not any(span.reduce(queue[0])):
                queue.pop(0)
            if queue:
                packet = evaluate(chosen + [queue[0]])
                trials.append((packet["status"] != "certified", Fraction(packet["lower"]), kind, packet))
        if not trials:
            break
        _, _, kind, current = min(trials, key=lambda x: (x[0], x[1]))
        added = queues[kind].pop(0)
        span.add(added)
        chosen.append(added)
        history.append({"rank": len(chosen), "added": kind, "status": current["status"],
                        "lower": current["lower_approximately"], "upper": current["upper_approximately"]})
    return {**current, "enrichment_history": history,
            "enrichment_claim": "heuristic choice between steady-state moments and target pullbacks; acceptance is the all-time bracket"}
