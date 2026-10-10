"""Exact autonomous polynomial observable closure, with explicit resource limits."""
from .rational import q, row_basis, coordinates, encoded, matrix


def clean(p):
    return {tuple(k): q(v) for k, v in p.items() if q(v)}


def validate(p, n):
    p = clean(p)
    if any(len(k) != n or any(type(e) is not int or e < 0 for e in k) for k in p):
        raise ValueError("Invalid monomial exponents")
    return p


def add(p, r):
    out = dict(p)
    for k, v in r.items():
        out[k] = out.get(k, q(0)) + v
    return clean(out)


def scale(p, c):
    return clean({k: v*c for k, v in p.items()})


def multiply(p, r):
    out = {}
    for k, v in p.items():
        for l, w in r.items():
            m = tuple(a+b for a, b in zip(k, l))
            out[m] = out.get(m, q(0)) + v*w
    return clean(out)


def compose(p, update, max_degree=32, max_terms=10000):
    n = len(update)
    one = {(0,)*n: q(1)}
    result = {}
    for exponent, coefficient in p.items():
        term = scale(one, coefficient)
        for i, power in enumerate(exponent):
            for _ in range(power):
                term = multiply(term, update[i])
                if len(term) > max_terms or any(sum(k) > max_degree for k in term):
                    raise OverflowError("Polynomial expansion exceeded resource limit")
        result = add(result, term)
        if len(result) > max_terms:
            raise OverflowError("Polynomial expansion exceeded term limit")
    return result


def evaluate(p, x):
    total = q(0)
    for k, coefficient in p.items():
        term = coefficient
        for v, e in zip(x, k):
            term *= q(v)**e
        total += term
    return total


def as_basis(polynomials):
    keys = sorted(set().union(*(p.keys() for p in polynomials)))
    if not keys:
        return []
    rows = row_basis([[p.get(k, q(0)) for k in keys] for p in polynomials])
    return [clean(dict(zip(keys, row))) for row in rows]


def coeffs(p, basis):
    keys = sorted(set(p).union(*(b.keys() for b in basis)))
    return coordinates([p.get(k, q(0)) for k in keys],
                       matrix([[b.get(k, q(0)) for k in keys] for b in basis]))


def pack(p):
    return [{"exponent": list(k), "coefficient": str(v)} for k, v in sorted(p.items())]


def unpack(terms, n):
    out = {}
    for term in terms:
        k = tuple(term["exponent"])
        out[k] = out.get(k, q(0)) + q(term["coefficient"])
    return validate(out, n)


def compile_polynomial(update, targets, max_degree=12, max_features=64, max_rounds=32):
    n = len(update)
    if not n or max_features < 1 or max_rounds < 0 or max_degree < 0:
        raise ValueError("Invalid dimensions or limits")
    update = [validate(p, n) for p in update]
    targets = [validate(p, n) for p in targets]
    basis = as_basis(targets)
    if not basis:
        raise ValueError("Need a nonzero target")
    dimensions = [len(basis)]
    for _ in range(max_rounds+1):
        if len(basis) > max_features or any(sum(k) > max_degree for p in basis for k in p):
            return {"status": "resource_limit", "dimensions": dimensions, "reason": "feature/degree limit"}
        try:
            pulled = [compose(p, update, max_degree=max_degree) for p in basis]
        except OverflowError as exc:
            return {"status": "resource_limit", "dimensions": dimensions, "reason": str(exc)}
        enlarged = as_basis(basis+pulled)
        if len(enlarged) == len(basis):
            return {"schema": "dynacomp.polynomial/1", "status": "exact", "variables": n,
                    "update": [pack(p) for p in update], "targets": [pack(p) for p in targets],
                    "features": [pack(p) for p in basis],
                    "G": encoded([coeffs(p, basis) for p in pulled]),
                    "D": encoded([coeffs(p, basis) for p in targets]), "dimensions": dimensions}
        basis = enlarged
        dimensions.append(len(basis))
    return {"status": "resource_limit", "dimensions": dimensions, "reason": "round limit"}


def replay_polynomial(packet):
    if packet.get("schema") != "dynacomp.polynomial/1" or packet.get("status") != "exact":
        raise ValueError("Need an exact polynomial packet")
    n = packet["variables"]
    update = [unpack(p, n) for p in packet["update"]]
    targets = [unpack(p, n) for p in packet["targets"]]
    basis = [unpack(p, n) for p in packet["features"]]
    if len(update) != n or not basis or len(as_basis(basis)) != len(basis):
        raise ValueError("Invalid update or dependent features")
    g, d = matrix(packet["G"]), matrix(packet["D"])
    if len(g) != len(basis) or len(g[0]) != len(basis) or len(d) != len(targets) or len(d[0]) != len(basis):
        raise ValueError("Coefficient matrix dimensions do not match")
    # Replay arbitrary accepted degrees with a limit informed by actual packet data.
    degree = max([sum(k) for p in update+basis+targets for k in p] or [0])
    for feature, row in zip(basis, g):
        rhs = {}
        for b, c in zip(basis, row):
            rhs = add(rhs, scale(b, c))
        if compose(feature, update, max_degree=max(32, degree*degree)) != rhs:
            raise ValueError("Polynomial transition identity fails")
    for target, row in zip(targets, d):
        rhs = {}
        for b, c in zip(basis, row):
            rhs = add(rhs, scale(b, c))
        if target != rhs:
            raise ValueError("Polynomial target identity fails")
    return True
