"""Portable binary64 bits, exact integer accumulation and round-to-nearest-even.

No arithmetic operation here depends on the host floating point rounding mode.
Finite inputs only for arithmetic; NaN payloads survive bits/from_bits transport.
Exact zero sums round to +0. Negative nonzero underflow rounds to -0.
"""
import struct
from fractions import Fraction

SIGN = 1 << 63
FRAC = (1 << 52) - 1
INF = 0x7ff0000000000000
MAX_FINITE = INF - 1


def _word(word):
    if type(word) is not int or not 0 <= word < 1 << 64:
        raise ValueError('Expected an unsigned 64-bit integer')
    return word


def bits(value):
    if type(value) is not float:
        raise TypeError('Expected a Python binary64 float')
    return int.from_bytes(struct.pack('>d', value), 'big')


def from_bits(word):
    return struct.unpack('>d', _word(word).to_bytes(8, 'big'))[0]


def classify(word):
    word = _word(word)
    e, f = (word >> 52) & 2047, word & FRAC
    return ('nan' if f else 'infinity') if e == 2047 else (
        ('subnormal' if f else 'zero') if e == 0 else 'normal')


def dyadic(word):
    """Signed integer m and exponent e with exact stored value m*2**e."""
    word = _word(word)
    exponent, fraction = (word >> 52) & 2047, word & FRAC
    if exponent == 2047:
        raise ValueError('Nonfinite operand')
    m, e = ((1 << 52) | fraction, exponent - 1075) if exponent else (fraction, -1074)
    if word & SIGN:
        m = -m
    if m:
        shift = (abs(m) & -abs(m)).bit_length() - 1
        m >>= shift
        e += shift
    return m, e


def exact(word):
    m, e = dyadic(word)
    return Fraction(m << e) if e >= 0 else Fraction(m, 1 << -e)


def _even_ratio(n, d):
    a, r = divmod(n, d)
    return a + int(2*r > d or (2*r == d and a & 1))


def round_bits(value):
    """Correctly round an exact rational, including subnormals and overflow."""
    if isinstance(value, float):
        raise TypeError('Use exact rationals, or bits() for existing floats')
    x = Fraction(value)
    sign = SIGN if x < 0 else 0
    n, d = abs(x.numerator), x.denominator
    if not n:
        return 0
    e = n.bit_length() - d.bit_length()
    if (n < d << e) if e >= 0 else (n << -e < d):
        e -= 1
    if e > 1023:
        return sign | INF
    unit = max(-1074, e - 52)
    m = _even_ratio(n, d << unit) if unit >= 0 else _even_ratio(n << -unit, d)
    if m == 1 << 53:
        m >>= 1
        unit += 1
    if unit > 971:
        return sign | INF
    if m < 1 << 52:
        return sign | m
    return sign | ((unit + 1075) << 52) | (m - (1 << 52))


def round_dyadic(m, e):
    if type(m) is not int or type(e) is not int:
        raise TypeError('Integer significand and exponent required')
    return round_bits(Fraction(m << e) if e >= 0 else Fraction(m, 1 << -e))


def exact_dot(left, right):
    """Integer superaccumulator: exact sum of products of stored operands."""
    left, right = tuple(left), tuple(right)
    if len(left) != len(right):
        raise ValueError('Dot dimensions differ')
    terms = []
    for a, b in zip(left, right):
        m, e = dyadic(a)
        n, f = dyadic(b)
        if m and n:
            terms.append((m*n, e+f))
    if not terms:
        return Fraction(0)
    exponent = min(e for _, e in terms)
    integer = sum(m << (e-exponent) for m, e in terms)
    return Fraction(integer << exponent) if exponent >= 0 else Fraction(integer, 1 << -exponent)


def dot_bits(left, right):
    """One final rounding, not a BLAS, sequential-dot, or FMA contract."""
    return round_bits(exact_dot(left, right))


def order_key(word):
    word = _word(word)
    if classify(word) == 'nan':
        raise ValueError('NaNs have no numeric order')
    return ((~word) & ((1 << 64)-1)) if word & SIGN else word | SIGN


def next_up(word):
    word = _word(word)
    if classify(word) == 'nan':
        raise ValueError('NaN has no next number')
    if word == INF:
        return word
    if word == SIGN:
        return 1
    return word - 1 if word & SIGN else word + 1


def next_down(word):
    return next_up(_word(word) ^ SIGN) ^ SIGN


def _root(n, degree):
    low, high = 0, 1 << ((n.bit_length()+degree-1)//degree)
    while low < high:
        mid = (low+high+1)//2
        if mid**degree <= n:
            low = mid
        else:
            high = mid-1
    return low


def power_certificate(word, degree):
    """Exact rational d-th power of the stored value; no decimal-intent claim."""
    word = _word(word)
    if type(degree) is not int or not 2 <= degree <= 4096:
        raise ValueError('Power degree must be between 2 and 4096')
    m, e = dyadic(word)
    root = _root(abs(m), degree)
    yes = not m or ((m > 0 or degree % 2) and e % degree == 0 and root**degree == abs(m))
    result = None
    if yes:
        root = -root if m < 0 else root
        k = e//degree if m else 0
        result = Fraction(root << k) if k >= 0 else Fraction(root, 1 << -k)
    return {'schema': 'exact.binary64-power/1', 'word': word, 'hex': f'{word:016x}',
            'degree': degree, 'classification': classify(word), 'odd_significand': m,
            'exponent': e, 'is_power': bool(yes), 'root': None if result is None else str(result),
            'scope': 'exact rational value stored in the finite binary64 word'}


def replay_power(packet):
    if packet != power_certificate(packet['word'], packet['degree']):
        raise ValueError('Binary64 power certificate does not reconstruct')
    return True
