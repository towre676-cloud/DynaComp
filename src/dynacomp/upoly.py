"""Univariate polynomials over Q: coefficient tuples, lowest degree first.

Includes Sturm sequences for exact real-root counting and isolation, and a
rational-root finder that needs no integer factorization.
"""
from fractions import Fraction
from math import lcm, gcd
from .rational import q


def trim(p):
    p = [q(x) for x in p]
    while p and not p[-1]:
        p.pop()
    return tuple(p)


def const(c):
    return trim([c])


def add(p, r):
    n = max(len(p), len(r))
    return trim([(p[i] if i < len(p) else 0) + (r[i] if i < len(r) else 0) for i in range(n)])


def neg(p):
    return tuple(-x for x in p)


def sub(p, r):
    return add(p, neg(r))


def scale(p, c):
    return trim([x*c for x in p])


def mul(p, r):
    if not p or not r:
        return ()
    out = [Fraction(0)]*(len(p)+len(r)-1)
    for i, a in enumerate(p):
        if a:
            for j, b in enumerate(r):
                out[i+j] += a*b
    return trim(out)


def degree(p):
    return len(p)-1


def divmod_(p, d):
    if not d:
        raise ZeroDivisionError("Polynomial division by zero")
    p, quotient = list(p), [Fraction(0)]*max(len(p)-len(d)+1, 1)
    while len(p) >= len(d) and p:
        c = p[-1]/d[-1]
        shift = len(p)-len(d)
        quotient[shift] = c
        for i, x in enumerate(d):
            p[shift+i] -= c*x
        p = list(trim(p))
    return trim(quotient), trim(p)


def monic(p):
    return scale(p, 1/p[-1]) if p else ()


def pgcd(p, r):
    while r:
        p, r = r, divmod_(p, r)[1]
    return monic(p)


def derivative(p):
    return trim([i*x for i, x in enumerate(p)][1:])


def evaluate(p, x):
    x, total = q(x), Fraction(0)
    for c in reversed(p):
        total = total*x + c
    return total


def primitive(p):
    """Integer coefficients with gcd 1 and positive leading coefficient."""
    if not p:
        return ()
    m = lcm(*(c.denominator for c in p))
    ints = [int(c*m) for c in p]
    g = gcd(*ints)
    sign = 1 if ints[-1] > 0 else -1
    return tuple(Fraction(sign*c//g) for c in ints)


def squarefree(p):
    if degree(p) < 1:
        return primitive(p)
    return primitive(divmod_(p, pgcd(p, derivative(p)))[0])


def sturm(p):
    seq = [p, derivative(p)]
    while seq[-1]:
        r = neg(divmod_(seq[-2], seq[-1])[1])
        if not r:
            break
        seq.append(r)
    return seq


def _variations(seq, x):
    signs = [s for s in (evaluate(p, x) for p in seq) if s]
    return sum((a > 0) != (b > 0) for a, b in zip(signs, signs[1:]))


def root_bound(p):
    """Cauchy bound: every real root lies strictly inside (-b, b)."""
    return 1 + max(abs(c/p[-1]) for c in p[:-1]) if degree(p) > 0 else Fraction(1)


def count_roots(p, lo, hi, seq=None):
    """Distinct real roots in (lo, hi] of a squarefree p."""
    seq = seq or sturm(p)
    return _variations(seq, lo) - _variations(seq, hi)


def isolate(p, width=None):
    """Disjoint (lo, hi] intervals each holding one real root of squarefree p."""
    if degree(p) < 1:
        return []
    seq, b = sturm(p), root_bound(p)
    stack, out = [(-b, b)], []
    while stack:
        lo, hi = stack.pop()
        k = count_roots(p, lo, hi, seq)
        if k == 0:
            continue
        if k == 1 and (width is None or hi-lo <= width):
            out.append((lo, hi))
            continue
        mid = (lo+hi)/2
        stack.extend([(mid, hi), (lo, mid)])
    return sorted(out)


def real_roots(p):
    """Exact rational roots and isolating intervals for irrational real roots.

    A rational root N/L of an integer polynomial with leading coefficient L
    makes L*x an integer, so refining each interval below width 1/L leaves at
    most one integer candidate. No factorization of coefficients is needed.
    """
    p = squarefree(p)
    if degree(p) < 1:
        return [], []
    lead = abs(int(p[-1]))
    rational, irrational = [], []
    for lo, hi in isolate(p, width=Fraction(1, 2*lead)):
        candidates = [Fraction(N, lead) for N in range(int((lo*lead)//1), int((hi*lead)//1)+2)]
        hit = [x for x in candidates if lo < x <= hi and not evaluate(p, x)]
        if hit:
            rational.append(hit[0])
        else:
            irrational.append((lo, hi))
    return rational, irrational


def encode(p):
    return [str(c) for c in p]


def parse(entry):
    """A rational scalar, or a list of coefficients from lowest degree."""
    if isinstance(entry, (list, tuple)):
        return trim(entry)
    return const(entry)
