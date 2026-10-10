"""Sharp output-error bounds over correlated (zonotopic) scenario domains.

Box domains treat every coordinate and every time independently. Physical
uncertainty is usually correlated: one unknown ambient offset shifts every
cell, a heater runs at one unknown but constant level. A zonotope expresses
this: the scenario (x0, u_0, ..., u_(T-1)) equals a centre plus sum_g xi_g w_g
with |xi_g| <= 1, where a generator w_g may touch any coordinates and times.

The error at time t, output i, is a linear functional L of the scenario, so
its exact supremum is |L(centre)| + sum_g |L(w_g)|, attained at
xi_g = orientation * sign(L(w_g)). Boxes are the special case of one
generator per coordinate and time. Every zonotope inside a box has a bound no
larger than the box's, and usually much smaller.
"""
from .rational import q, matrix, mv, addv, encoded, identity, mul
from .uniform import as_sparse, _kernels, _pullbacks, _project, SparseLinearModel


def _gen(raw, n, m, horizon):
    x = {int(k): q(v) for k, v in (raw.get("x0") or {}).items()}
    u = {}
    for t, vec in (raw.get("inputs") or {}).items():
        t = int(t)
        if not 0 <= t < horizon or len(vec) != m:
            raise ValueError("Generator input entries need a valid time and m coordinates")
        u[t] = tuple(q(v) for v in vec)
    if any(not 0 <= k < n for k in x):
        raise ValueError("Generator state index out of range")
    return x, u


def box_generators(n, radius=None, input_radii=()):
    """One generator per coordinate and time: the box as a zonotope."""
    gens = []
    if radius is not None:
        gens += [{"x0": {k: r}} for k, r in enumerate(radius) if q(r)]
    for t, row in enumerate(input_radii):
        for k, r in enumerate(row):
            if q(r):
                vec = [0]*len(row)
                vec[k] = r
                gens.append({"inputs": {t: vec}})
    return gens


def _serialize(gens):
    return [{"x0": {str(k): str(v) for k, v in x.items()}, "inputs": {str(t): [str(a) for a in vec] for t, vec in u.items()}}
            for x, u in gens]


def sharp_zonotope_bound(model, V, E, G, H, x0_center, input_centers, generators):
    model = as_sparse(model)
    v, e, g, h = map(matrix, (V, E, G, H))
    n, m = model.operator.n, len(model.B[0])
    if mul(e, v) != identity(len(g)):
        raise ValueError("Need EV = I")
    c = tuple(q(x) for x in x0_center)
    a = [tuple(q(x) for x in row) for row in input_centers]
    horizon = len(a)
    if len(c) != n or any(len(row) != m for row in a):
        raise ValueError("Centre dimensions do not match the model")
    gens = [_gen(raw, n, m, horizon) for raw in generators]
    d = mul(model.C, v)
    state, forcing = _kernels(model, e, g, h, d, _pullbacks(model, horizon))
    table, best = [], None
    for t, kernel in enumerate(state):
        row_bounds = []
        for i, krow in enumerate(kernel):
            mean = sum(x*y for x, y in zip(krow, c))
            mean += sum(sum(x*y for x, y in zip(forcing[t-1-j][i], a[j])) for j in range(t))
            coeffs = []
            for gx, gu in gens:
                value = sum(krow[k]*w for k, w in gx.items())
                value += sum(sum(x*y for x, y in zip(forcing[t-1-j][i], vec)) for j, vec in gu.items() if j < t)
                coeffs.append(value)
            bound = abs(mean) + sum(abs(x) for x in coeffs)
            row_bounds.append(bound)
            if best is None or bound > best[0]:
                best = (bound, t, i, mean, coeffs)
        table.append(row_bounds)
    bound, t, i, mean, coeffs = best
    orientation = -1 if mean < 0 else 1
    signs = [orientation*((x > 0) - (x < 0)) for x in coeffs]
    x0 = list(c)
    word = [list(row) for row in a]
    for s, (gx, gu) in zip(signs, gens):
        for k, w in gx.items():
            x0[k] += s*w
        for j, vec in gu.items():
            word[j] = [p + s*w for p, w in zip(word[j], vec)]
    return {"schema": "dynacomp.zonotope-bound/1", "arithmetic": "rational",
            "model": model.spec(), "V": encoded(v), "E": encoded(e), "G": encoded(g), "H": encoded(h),
            "x0_center": [str(x) for x in c], "input_centers": encoded(a) if a else [],
            "generators": _serialize(gens), "generator_count": len(gens),
            "bound": [str(max(row)) for row in table], "worst_bound": str(bound),
            "worst_bound_approximately": f"{float(bound):.9g}",
            "attaining_witness": {"time": t, "output": i, "x0": [str(x) for x in x0],
                                  "inputs": encoded(word) if word else [],
                                  "generator_signs": signs,
                                  "signed_error": str(orientation*bound)},
            "claim": "sharp per-time supremum over the zonotope of initial states and input words, with an attaining scenario"}


def replay_zonotope(packet):
    """Recompute the bound and simulate the witness through both models exactly."""
    if packet.get("schema") != "dynacomp.zonotope-bound/1":
        raise ValueError("Unknown schema")
    model = SparseLinearModel.from_spec(packet["model"])
    v, e, g, h = (matrix(packet[k]) for k in ("V", "E", "G", "H"))
    fresh = sharp_zonotope_bound(model, v, e, g, h, packet["x0_center"], packet["input_centers"], packet["generators"])
    if fresh != packet:
        raise ValueError("Stored bound or witness does not match recomputation")
    w = packet["attaining_witness"]
    x = tuple(q(a) for a in w["x0"])
    z = mv(e, x)
    for u in w["inputs"][:w["time"]]:
        x = addv(model.operator.step(x), mv(model.B, u))
        z = addv(mv(g, z), mv(h, u))
    d = mul(model.C, v)
    actual = mv(model.C, x)[w["output"]] - mv(d, z)[w["output"]]
    if actual != q(w["signed_error"]) or abs(actual) != q(packet["worst_bound"]):
        raise ValueError("Witness does not attain the bound")
    return True
