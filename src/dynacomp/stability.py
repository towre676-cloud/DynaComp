"""Exact stability decisions and rationally certified quadratic metrics.

Decision: the characteristic polynomial is computed over Q (Faddeev-LeVerrier)
and tested by the Schur-Cohn recursion (discrete time, all eigenvalues in the
open unit disk) or by Schur-Cohn after a Cayley transform (continuous time,
all eigenvalues in the open left half-plane). No eigenvalue is approximated.

Certificate: a rational symmetric P with P > 0 and gamma P - A^T P A > 0 for a
rational gamma < 1 (discrete), or A^T P + P A + 2 sigma P < 0 (continuous).
Floating-point Smith iteration only proposes P; acceptance is an exact
Sylvester (LDL^T pivot) test. A stable matrix always has such a P, but a
proposal can fail at finite precision; failure is reported, never hidden.
"""
from fractions import Fraction
from math import isqrt
from .upoly import primitive
from .rational import matrix, mul, transpose, identity, sub, q, encoded, inverse


def charpoly(A):
    """Coefficients of det(zI - A), lowest degree first, by Faddeev-LeVerrier."""
    a = matrix(A)
    n = len(a)
    if len(a[0]) != n:
        raise ValueError("Expected a square matrix")
    coefficients = [Fraction(1)]
    m = tuple(tuple(Fraction(0) for _ in range(n)) for _ in range(n))
    for k in range(1, n+1):
        shifted = tuple(tuple(m[i][j] + (coefficients[-1] if i == j else 0) for j in range(n)) for i in range(n))
        m = mul(a, shifted)
        coefficients.append(-sum(m[i][i] for i in range(n))/k)
    # coefficients[k] multiplies z^(n-k).
    return tuple(reversed(coefficients))


def schur_cohn(p):
    """Decide whether every root of real p lies in the open unit disk.

    Schur-Cohn: p is Schur stable iff |p_0| < |p_n| and
    (p_n p(z) - p_0 p*(z)) / z is Schur stable, where p* reverses p.
    Returns (stable, steps), with the failing step if unstable.
    """
    p = [q(c) for c in p]
    while p and not p[-1]:
        p.pop()
    if not p:
        raise ValueError("Zero polynomial")
    steps = []
    while len(p) > 1:
        p = list(primitive(p))  # scaling preserves roots and keeps integers small
        a0, an = p[0], p[-1]
        ok = abs(a0) < abs(an)
        steps.append({"degree": len(p)-1, "abs_constant_below_abs_leading": ok})
        if not ok:
            return False, steps
        reduced = [an*x - a0*y for x, y in zip(p, reversed(p))]
        p = reduced[1:]
    return True, steps


def schur_stable(A):
    stable, steps = schur_cohn(charpoly(A))
    return {"stable": stable, "criterion": "Schur-Cohn on the exact characteristic polynomial",
            "charpoly": [str(c) for c in charpoly(A)], "steps": steps}


def cayley(A):
    """(I - A)^(-1) (I + A); None when 1 is an eigenvalue (then A is not Hurwitz)."""
    a = matrix(A)
    n = len(a)
    try:
        return mul(inverse(sub(identity(n), a)), tuple(tuple(x + (i == j) for j, x in enumerate(row)) for i, row in enumerate(a)))
    except ValueError:
        return None


def hurwitz_stable(A):
    """Continuous time: Re(lambda) < 0 iff |(1+lambda)/(1-lambda)| < 1."""
    c = cayley(A)
    if c is None:
        return {"stable": False, "criterion": "eigenvalue 1 makes I - A singular", "steps": []}
    stable, steps = schur_cohn(charpoly(c))
    return {"stable": stable, "criterion": "Schur-Cohn after the Cayley transform", "steps": steps}


def positive_definite(M):
    """Exact Sylvester test by symmetric elimination: every pivot must be positive."""
    a = [list(row) for row in matrix(M)]
    n = len(a)
    if len(a[0]) != n or any(a[i][j] != a[j][i] for i in range(n) for j in range(i)):
        raise ValueError("Expected a symmetric matrix")
    for k in range(n):
        pivot = a[k][k]
        if pivot <= 0:
            return False
        for i in range(k+1, n):
            if a[i][k]:
                f = a[i][k]/pivot
                for j in range(k+1, n):
                    a[i][j] -= f*a[k][j]
    return True


def _float(a):
    return [[float(x) for x in row] for row in a]


def _fmul(a, b):
    bt = list(zip(*b))
    return [[sum(x*y for x, y in zip(row, col)) for col in bt] for row in a]


def spectral_radius_estimate(a, doublings=40):
    """Float estimate from ||A^(2^j)||^(2^-j); a proposal aid only."""
    f = _float(a)
    estimate = None
    for j in range(doublings):
        norm = max(sum(abs(x) for x in row) for row in f)
        if norm == 0:
            return 0.0
        estimate = norm**(1/2**j)
        if norm > 1e150 or norm < 1e-150:
            break
        f = _fmul(f, f)
    return estimate


def _smith(a, continuous=False, doublings=48, scale=1.0):
    """Float proposal for sum_k (A^T/s)^k (A/s)^k; a scale s just above the
    spectral radius yields A^T P A <= s^2 (P - I), so gamma near s^2."""
    n = len(a)
    f = [[x/scale for x in row] for row in _float(a)]
    p = [[float(i == j) for j in range(n)] for i in range(n)]
    power = f
    for _ in range(doublings):
        pt = [list(r) for r in zip(*power)]
        p = [[x+y for x, y in zip(r, s)] for r, s in zip(p, _fmul(_fmul(pt, p), power))]
        power = _fmul(power, power)
        if max(abs(x) for row in power for x in row) < 1e-30:
            break
        if max(abs(x) for row in p for x in row) > 1e200:
            return None
    return p


def _rationalize(p, denominator):
    n = len(p)
    r = [[Fraction((p[i][j]+p[j][i])/2).limit_denominator(denominator) for j in range(n)] for i in range(n)]
    for i in range(n):
        for j in range(i):
            r[i][j] = r[j][i]
    return tuple(tuple(row) for row in r)


def _sym(m):
    return tuple(tuple((m[i][j]+m[j][i])/2 for j in range(len(m))) for i in range(len(m)))


def contraction_holds(A, P, gamma):
    """gamma P - A^T P A > 0 exactly."""
    a, p = matrix(A), matrix(P)
    apa = mul(mul(transpose(a), p), a)
    return positive_definite(tuple(tuple(gamma*x - y for x, y in zip(r, s)) for r, s in zip(p, apa)))


def _best_gamma(a, p, bits=24):
    """Smallest dyadic gamma in (0, 1] found by bisection with an exact acceptance test."""
    if not contraction_holds(a, p, Fraction(1)):
        return None
    lo, hi = Fraction(0), Fraction(1)
    for _ in range(bits):
        mid = (lo+hi)/2
        if contraction_holds(a, p, mid):
            hi = mid
        else:
            lo = mid
    return hi


def quadratic_certificate(A, denominators=(2**12, 2**20, 2**32, 2**48), bits=24):
    """Discover and exactly accept P > 0, A^T P A < gamma P with gamma < 1."""
    a = matrix(A)
    decision = schur_stable(a)
    if not decision["stable"]:
        return {"schema": "dynacomp.quadratic-metric/1", "status": "not_schur_stable", "decision": decision}
    rho = spectral_radius_estimate(a) or 0.0
    best = None
    for scale in sorted({1.0, min(1.0, rho + (1-rho)/2), min(1.0, rho + (1-rho)/8), min(1.0, rho + (1-rho)/32)}):
        if scale <= 0:
            continue
        proposal = _smith(a, scale=scale)
        for d in denominators if proposal is not None else ():
            p = _rationalize(proposal, d)
            if not positive_definite(p):
                continue
            gamma = _best_gamma(a, p, bits)
            if gamma is not None and gamma < 1:
                if best is None or gamma < best[1]:
                    best = (p, gamma)
                break
    if best is not None:
        p, gamma = best
        return {"schema": "dynacomp.quadratic-metric/1", "status": "certified", "A": encoded(a),
                "P": encoded(p), "gamma": str(gamma), "decision": decision,
                "claim": "x^T P x decreases by at least the factor gamma per step"}
    return {"schema": "dynacomp.quadratic-metric/1", "status": "proposal_failed", "decision": decision,
            "claim": "A is Schur stable, but no rational proposal was accepted at the tried precisions"}


def replay_quadratic(packet):
    if packet.get("schema") != "dynacomp.quadratic-metric/1" or packet.get("status") != "certified":
        raise ValueError("Need a certified quadratic metric packet")
    a, p, gamma = matrix(packet["A"]), matrix(packet["P"]), q(packet["gamma"])
    if not 0 <= gamma < 1:
        raise ValueError("Contraction factor must lie in [0, 1)")
    if not positive_definite(p):
        raise ValueError("P is not positive definite")
    if not contraction_holds(a, p, gamma):
        raise ValueError("gamma P - A^T P A is not positive definite")
    return True


def common_quadratic_certificate(modes, denominators=(2**12, 2**20, 2**32), bits=20):
    """One P contracting every mode: stability under arbitrary switching.

    The proposal sums the per-mode Smith matrices; a common P need not exist
    even when every mode is stable, so failure is not proof of instability.
    """
    mats = {str(k): matrix(v) for k, v in modes.items()}
    total = None
    for name, a in mats.items():
        if not schur_stable(a)["stable"]:
            return {"schema": "dynacomp.common-quadratic/1", "status": "mode_unstable", "mode": name}
        s = _smith(a)
        if s is None:
            return {"schema": "dynacomp.common-quadratic/1", "status": "proposal_failed"}
        total = s if total is None else [[x+y for x, y in zip(r, t)] for r, t in zip(total, s)]
    for d in denominators:
        p = _rationalize(total, d)
        if not positive_definite(p):
            continue
        gammas = [_best_gamma(a, p, bits) for a in mats.values()]
        if all(g is not None and g < 1 for g in gammas):
            return {"schema": "dynacomp.common-quadratic/1", "status": "certified",
                    "modes": {k: encoded(a) for k, a in mats.items()}, "P": encoded(p),
                    "gamma": str(max(gammas)),
                    "claim": "x^T P x contracts by gamma under every mode, hence under every switching sequence"}
    return {"schema": "dynacomp.common-quadratic/1", "status": "proposal_failed",
            "claim": "no common quadratic function accepted; switching stability is undecided"}


def replay_common_quadratic(packet):
    if packet.get("schema") != "dynacomp.common-quadratic/1" or packet.get("status") != "certified":
        raise ValueError("Need a certified common quadratic packet")
    p, gamma = matrix(packet["P"]), q(packet["gamma"])
    if not 0 <= gamma < 1 or not positive_definite(p):
        raise ValueError("Invalid metric or factor")
    for name, a in packet["modes"].items():
        if not contraction_holds(matrix(a), p, gamma):
            raise ValueError(f"Mode {name} does not contract")
    return True


def sqrt_up(x, bits=48):
    """Rational upper bound of sqrt(x) for rational x >= 0, within about 2^-bits relative."""
    x = q(x)
    if x < 0:
        raise ValueError("Negative argument")
    if not x:
        return Fraction(0)
    scale = 4**bits
    n = x.numerator*x.denominator*scale
    r = isqrt(n)
    if r*r < n:
        r += 1
    return Fraction(r, x.denominator*2**bits)


def continuous_certificate(A, **options):
    """P > 0 with A^T P + P A < 0, proposed through the Cayley transform.

    If C = (I - A)^(-1) (I + A) satisfies C^T P C < P then A^T P + P A < 0;
    acceptance nevertheless checks the continuous inequality directly.
    """
    a = matrix(A)
    decision = hurwitz_stable(a)
    if not decision["stable"]:
        return {"schema": "dynacomp.continuous-metric/1", "status": "not_hurwitz_stable", "decision": decision}
    discrete = quadratic_certificate(cayley(a), **options)
    if discrete["status"] == "certified":
        p = matrix(discrete["P"])
        if _negative_lyapunov(a, p):
            return {"schema": "dynacomp.continuous-metric/1", "status": "certified", "A": encoded(a),
                    "P": encoded(p), "decision": decision,
                    "claim": "d/dt x^T P x < 0 for every nonzero x, so x' = A x is exponentially stable"}
    return {"schema": "dynacomp.continuous-metric/1", "status": "proposal_failed", "decision": decision}


def _negative_lyapunov(a, p):
    m = mul(transpose(a), p)
    s = tuple(tuple(-(m[i][j] + m[j][i]) for j in range(len(a))) for i in range(len(a)))
    return positive_definite(p) and positive_definite(s)


def replay_continuous(packet):
    if packet.get("schema") != "dynacomp.continuous-metric/1" or packet.get("status") != "certified":
        raise ValueError("Need a certified continuous metric packet")
    if not _negative_lyapunov(matrix(packet["A"]), matrix(packet["P"])):
        raise ValueError("P > 0 and A^T P + P A < 0 do not both hold")
    return True
