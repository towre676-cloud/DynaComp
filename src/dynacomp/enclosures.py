"""Small outward rational enclosures; no floating special-function oracle.

These routines bound real constants, not runtime floating-point evaluation.
Taylor remainders and integer square roots supply the endpoint inequalities.
"""
from fractions import Fraction
from math import isqrt
from .rational import q


def precision(bits):
    if type(bits) is not int or not 16 <= bits <= 512:
        raise ValueError("Precision must be an integer in [16,512]")
    return 1 << bits


def outward(lo, hi, bits):
    scale = precision(bits)
    lo, hi = q(lo), q(hi)
    if lo > hi:
        raise ValueError("Reversed interval")
    return (Fraction((lo*scale).__floor__(), scale),
            Fraction((hi*scale).__ceil__(), scale))


def sqrt_integer(n, bits=80):
    scale = precision(bits)
    if type(n) is not int or n < 0:
        raise ValueError("Nonnegative integer radicand required")
    root = isqrt(n*scale*scale)
    return Fraction(root, scale), Fraction(root+int(root*root != n*scale*scale), scale)


def log_rational(x, bits=80):
    """Positive rational logarithm, reduced to atanh arguments <= 1/3."""
    precision(bits)
    x = q(x)
    if x <= 0:
        raise ValueError("Positive logarithm argument required")
    if x < 1:
        lo, hi = log_rational(1/x, bits)
        return -hi, -lo
    exponent = x.numerator.bit_length()-x.denominator.bit_length()
    if x < 2**exponent:
        exponent -= 1
    reduced = x / (2**exponent)

    def series(value):
        z = (value-1)/(value+1)
        power = z
        total = Fraction(0)
        for j in range(bits+8):
            total += 2*power/(2*j+1)
            power *= z*z
            remainder = 2*power/((2*j+3)*(1-z*z))
            if remainder <= Fraction(1, 1 << (bits+8)):
                break
        return total, total+remainder

    a, b = series(Fraction(2))
    c, d = series(reduced)
    return outward(exponent*a+c, exponent*b+d, bits)


def pi_interval(bits=80):
    """Machin's identity and the alternating arctangent remainder."""
    precision(bits)

    def atan_inverse(n):
        x = Fraction(1, n)
        total = Fraction(0)
        power = x
        for j in range(bits+8):
            total += (-1)**j*power/(2*j+1)
            power *= x*x
            next_term = (-1)**(j+1)*power/(2*j+3)
            if abs(next_term) <= Fraction(1, 1 << (bits+8)):
                break
        return min(total, total+next_term), max(total, total+next_term)

    a, b = atan_inverse(5)
    c, d = atan_inverse(239)
    return outward(16*a-4*d, 16*b-4*c, bits)


def exp_negative(x, bits=80):
    """Enclose exp(-x), x >= 0, by alternating Taylor and squaring."""
    precision(bits)
    x = q(x)
    if x < 0:
        raise ValueError("Nonnegative exponential argument required")
    if not x:
        return Fraction(1), Fraction(1)
    squarings = 0
    while x > Fraction(1, 2):
        x /= 2
        squarings += 1
    # Extra guard bits control outward errors through repeated squaring.
    work = bits+squarings.bit_length()+8
    # Public precision limit is independent of this internal guard budget.
    scale = 1 << work
    total = Fraction(1)
    term = Fraction(1)
    for j in range(1, bits+16):
        term *= -x/j
        total += term
        next_term = term*(-x)/(j+1)
        if abs(next_term) <= Fraction(1, scale):
            break
    lo, hi = min(total, total+next_term), max(total, total+next_term)
    for _ in range(squarings):
        lo = Fraction((lo*lo*scale).__floor__(), scale)
        hi = Fraction((hi*hi*scale).__ceil__(), scale)
    return outward(max(Fraction(0), lo), min(Fraction(1), hi), bits)
