"""Two-cusp real evaluation of R(exp(-t)) with outward rational endpoints.

Uses the classical product and reciprocity cited in spectral.PAPER. This is
an evaluation algorithm, not a proof of irrationality of odd zeta values.
"""
from fractions import Fraction
from .rational import q
from .enclosures import exp_negative, pi_interval, sqrt_integer, precision, outward
from .spectral import character5, PAPER


def _product(tlo, thi, bits):
    qlo = exp_negative(thi, bits)[0]
    qhi = exp_negative(tlo, bits)[1]
    if qhi >= 1:
        raise ValueError("Cannot separate q from one at this precision")
    count = 1
    while True:
        tail = qhi**(count+1)/((1-qhi)*(1-qhi**(count+1)))
        if tail <= Fraction(1, 1 << bits):
            break
        count *= 2
        if count > 4096:
            raise ValueError("Product budget exhausted")
    # Retain a small dyadic upper bound, rather than serializing a huge power.
    tail = outward(tail, tail, bits)[1]
    lo, hi = exp_negative(thi/5, bits)[0], exp_negative(tlo/5, bits)[1]
    for j in range(1, count+1):
        char = character5(j)
        if char == 1:
            lo *= 1-qhi**j
            hi *= 1-qlo**j
        elif char == -1:
            lo /= 1-qlo**j
            hi /= 1-qhi**j
        lo, hi = outward(lo, hi, bits)
    factor = exp_negative(tail, bits)[0]
    lo, hi = outward(lo*factor, hi/factor, bits)
    return lo, hi, count, tail


def ramanujan_value(t, bits=80):
    """Evaluate a positive rational t using the faster cusp expansion.

    t<6 uses reciprocity; t>=6 uses the direct product. The threshold only
    selects a representation, never extrapolates an unbounded remainder.
    """
    t = q(t)
    precision(bits)
    if t <= 0:
        raise ValueError("Positive t required")
    if t >= 6:
        lo, hi, count, tail = _product(t, t, bits)
        branch = "direct"
    else:
        pilo, pihi = pi_interval(bits)
        rlo, rhi, count, tail = _product(4*pilo*pilo/t, 4*pihi*pihi/t, bits)
        root = sqrt_integer(5, bits)
        philo, phihi = (1+root[0])/2, (1+root[1])/2
        # Reciprocity rearranged to avoid subtracting two nearby constants.
        lo = max(Fraction(0), (1-phihi*rhi)/(phihi+rhi))
        hi = (1-philo*rlo)/(philo+rlo)
        branch = "reciprocal"
    lo, hi = outward(lo, hi, bits)
    return {"schema": "dynacomp.ramanujan-value/1", "t": str(t), "bits": bits,
            "interval": [str(lo), str(hi)], "branch": branch, "product_terms": count,
            "log_product_tail": str(tail), "source_reference": PAPER,
            "claim": "R(exp(-t)) in rational interval, conditional on classical product and reciprocal identities; no floating evaluation error included."}


def replay_ramanujan_value(packet):
    if packet.get("schema") != "dynacomp.ramanujan-value/1":
        raise ValueError("Unknown Ramanujan value schema")
    if ramanujan_value(packet["t"], packet["bits"]) != packet:
        raise ValueError("Two-cusp evaluation or tail differs")
    return True
