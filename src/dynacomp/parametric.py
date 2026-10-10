"""Structural versus accidental closure for one real parameter.

Entries of A(t), B(t), C(t) are polynomials in t. A reduction compiled at one
parameter value may be smaller than at nearby values: closure there is an
accident of that coefficient choice. For a real polynomial matrix M(t), let
e_k(t) be the k-th elementary symmetric function of the eigenvalues of
M(t)^T M(t). By Cauchy-Binet, e_k is the sum of squares of all k x k minors of
M, so for real t, rank M(t) >= k exactly when e_k(t) != 0. The generic rank r
is the largest k with e_k not identically zero, and the exceptional real
parameters are precisely the real roots of e_r. Sturm sequences count and
isolate them exactly; rational exceptional values are found exactly.
"""
from fractions import Fraction
from . import upoly as P
from .rational import matrix, rref


def poly_matrix(rows):
    out = tuple(tuple(P.parse(x) for x in row) for row in rows)
    if not out or not out[0] or any(len(r) != len(out[0]) for r in out):
        raise ValueError("Expected a nonempty rectangular parametric matrix")
    return out


def pm_mul(a, b):
    if len(a[0]) != len(b):
        raise ValueError("Matrix dimensions do not match")
    cols = tuple(zip(*b))
    out = []
    for row in a:
        new = []
        for col in cols:
            s = ()
            for x, y in zip(row, col):
                if x and y:
                    s = P.add(s, P.mul(x, y))
            new.append(s)
        out.append(tuple(new))
    return tuple(out)


def pm_transpose(a):
    return tuple(zip(*a))


def pm_eval(a, t):
    return matrix([[P.evaluate(x, t) for x in row] for row in a])


def elementary_symmetric(m):
    """[e_0, ..., e_n] of the eigenvalues of square polynomial matrix m (Faddeev-LeVerrier)."""
    n = len(m)
    coefficient = [P.const(1)]
    mk = tuple(tuple(() for _ in range(n)) for _ in range(n))
    for k in range(1, n+1):
        # M_k = M (M_(k-1) + c_(k-1) I), c_k = -tr(M_k)/k.
        shifted = tuple(tuple(P.add(mk[i][j], coefficient[-1]) if i == j else mk[i][j]
                              for j in range(n)) for i in range(n))
        mk = pm_mul(m, shifted)
        trace = ()
        for i in range(n):
            trace = P.add(trace, mk[i][i])
        coefficient.append(P.scale(trace, Fraction(-1, k)))
    # coefficient[k] multiplies lambda^(n-k); e_k = (-1)^k coefficient[k].
    return [P.scale(c, (-1)**k) for k, c in enumerate(coefficient)]


def rank_profile(m):
    """Generic rank, drop polynomial and exceptional real parameters of m(t)."""
    gram = pm_mul(pm_transpose(m), m)
    e = elementary_symmetric(gram)
    r = max(k for k, p in enumerate(e) if p)
    if r == 0:
        return {"generic_rank": 0, "drop_polynomial": ["1"], "rational_exceptions": [],
                "irrational_exceptions": [], "real_exception_count": 0, "witness": "0"}
    drop = P.squarefree(e[r])
    rational, irrational = P.real_roots(drop)
    witness = next(Fraction(k) for k in _small_integers() if P.evaluate(drop, k))
    return {"generic_rank": r, "drop_polynomial": P.encode(drop),
            "real_exception_count": len(rational)+len(irrational),
            "rational_exceptions": [{"value": str(x), "rank": rank_at(m, x)} for x in rational],
            "irrational_exceptions": [{"interval": [str(lo), str(hi)]} for lo, hi in irrational],
            "witness": str(witness)}


def _small_integers():
    k = 0
    while True:
        yield k
        if k:
            yield -k
        k += 1


def rank_at(m, t):
    return len(rref(pm_eval(m, t))[1])


def _krylov_stack(a, start, n):
    """[S; S A; ...; S A^(n-1)] for row block S."""
    rows, block = [], start
    for _ in range(n):
        rows.extend(block)
        block = pm_mul(block, a)
    return tuple(rows)


def parametric_analysis(A, B, C):
    """Observable, reachable and minimal (Hankel) orders as functions of t."""
    a, b, c = poly_matrix(A), poly_matrix(B), poly_matrix(C)
    n = len(a)
    if len(a[0]) != n or len(b) != n or len(c[0]) != n:
        raise ValueError("Expected A[n,n], B[n,m], C[p,n]")
    observability = _krylov_stack(a, c, n)
    reachability = pm_transpose(_krylov_stack(pm_transpose(a), pm_transpose(b), n))
    hankel = pm_mul(observability, reachability)
    encode = lambda m: [[P.encode(x) for x in row] for row in m]
    return {"schema": "dynacomp.parametric/1", "arithmetic": "rational polynomial in t",
            "A": encode(a), "B": encode(b), "C": encode(c),
            "observable": rank_profile(observability),
            "reachable": rank_profile(reachability),
            "minimal": rank_profile(hankel),
            "claim": "ranks for real t; exceptional values are exactly the real roots of each drop polynomial"}


def classify(packet, t, quantity="observable"):
    """Is the closure at parameter t structural (generic order) or accidental?"""
    profile = packet[quantity]
    drop = tuple(Fraction(x) for x in profile["drop_polynomial"])
    accidental = P.evaluate(drop, t) == 0
    return {"value": str(Fraction(t)), "quantity": quantity, "accidental": accidental,
            "generic_order": profile["generic_rank"]}


def replay_parametric(packet):
    """Re-derive each profile from the stored model and compare.

    Also evaluates the rank at the generic witness and every rational exception
    directly by rational elimination, independent of the symmetric functions.
    """
    if packet.get("schema") != "dynacomp.parametric/1":
        raise ValueError("Unknown schema")
    fresh = parametric_analysis(packet["A"], packet["B"], packet["C"])
    for key in ("observable", "reachable", "minimal"):
        if fresh[key] != packet[key]:
            raise ValueError(f"{key} profile does not match the model")
    a, b, c = (poly_matrix(packet[k]) for k in ("A", "B", "C"))
    n = len(a)
    stacks = {"observable": _krylov_stack(a, c, n),
              "reachable": pm_transpose(_krylov_stack(pm_transpose(a), pm_transpose(b), n))}
    for key, m in stacks.items():
        profile = packet[key]
        if rank_at(m, Fraction(profile["witness"])) != profile["generic_rank"]:
            raise ValueError(f"{key} witness rank differs from generic rank")
        for x in profile["rational_exceptions"]:
            if rank_at(m, Fraction(x["value"])) >= profile["generic_rank"]:
                raise ValueError(f"{key} exception {x['value']} does not drop rank")
    return True
