"""Exact rational OGFs and a supported algebraic-series remainder calculus.

Coefficients are low degree first. A realization means a_n=C A^n x0,
including n=0; it is not an input/output transfer-function indexing convention.
The same stdlib-only file is distributed in PerfectPower.
"""
from dataclasses import dataclass
from fractions import Fraction as Q


def _q(x):
    if isinstance(x, float):
        raise TypeError('Exact rational input required')
    value=Q(x)
    if max(abs(value.numerator).bit_length(),value.denominator.bit_length()) > 16384:
        raise ValueError('Exact generating-function arithmetic exceeds 16384 bits')
    return value


def _poly(xs):
    xs = list(map(_q, xs))
    while xs and not xs[-1]:
        xs.pop()
    return tuple(xs) or (Q(0),)


def _add(a, b):
    return _poly([(a[i] if i < len(a) else 0)+(b[i] if i < len(b) else 0)
                  for i in range(max(len(a), len(b)))])


def _mul(a, b):
    c = [Q(0)]*(len(a)+len(b)-1)
    for i, x in enumerate(a):
        for j, y in enumerate(b):
            c[i+j] += x*y
    return _poly(c)


def _div(a, b):
    a, b = list(_poly(a)), _poly(b)
    if b == (0,):
        raise ValueError('Zero denominator')
    c = [Q(0)]*max(1, len(a)-len(b)+1)
    while a and len(a) >= len(b):
        k, v = len(a)-len(b), a[-1]/b[-1]
        c[k] = v
        for j, x in enumerate(b):
            a[j+k] -= v*x
        while a and not a[-1]:
            a.pop()
    return _poly(c), _poly(a)


def _gcd(a, b):
    while b != (0,):
        a, b = b, _div(a, b)[1]
    return tuple(x/a[-1] for x in a)


def _matrix(rows):
    rows = tuple(tuple(_q(x) for x in row) for row in rows)
    if not rows or not rows[0] or any(len(row) != len(rows[0]) for row in rows):
        raise ValueError('Nonempty rectangular matrix required')
    return rows


def _mm(a, b):
    return tuple(tuple(_q(sum(x*y for x, y in zip(row, col))) for col in zip(*b)) for row in a)


def _mv(a, x):
    return tuple(_q(sum(v*w for v, w in zip(row, x))) for row in a)


def _rpow(x, n):
    result=Q(1)
    while n:
        if n & 1:
            result=_q(result*x)
        n//=2
        if n:
            x=_q(x*x)
    return result


def _identity(n):
    return tuple(tuple(Q(i == j) for j in range(n)) for i in range(n))


def _power(a, n):
    result = _identity(len(a))
    while n:
        if n & 1:
            result = _mm(result, a)
        a, n = _mm(a, a), n//2
    return result


def _det(a):
    a = [list(row) for row in a]
    result = Q(1)
    for j in range(len(a)):
        pivot = next((k for k in range(j, len(a)) if a[k][j]), None)
        if pivot is None:
            return Q(0)
        if pivot != j:
            a[j], a[pivot] = a[pivot], a[j]
            result = -result
        v = a[j][j]
        result *= v
        for k in range(j+1, len(a)):
            t = a[k][j]/v
            for i in range(j+1, len(a)):
                a[k][i] -= t*a[j][i]
    return result


def _count(n, limit=10000):
    if type(n) is not int or not 0 <= n <= limit:
        raise ValueError(f'Integer count in 0..{limit} required')
    return n


@dataclass(frozen=True)
class RationalGF:
    numerator: tuple
    denominator: tuple = (Q(1),)

    def __post_init__(self):
        p, q = _poly(self.numerator), _poly(self.denominator)
        if not q[0]:
            raise ValueError('An ordinary generating function needs Q(0) != 0')
        g = _gcd(p, q)
        p, q = _div(p, g)[0], _div(q, g)[0]
        factor = q[0]
        object.__setattr__(self, 'numerator', tuple(x/factor for x in p))
        object.__setattr__(self, 'denominator', tuple(x/factor for x in q))

    @property
    def dimension(self):
        return 0 if self.numerator == (0,) else max(len(self.denominator)-1, len(self.numerator))

    def coefficients(self, count):
        out = []
        for n in range(_count(count)):
            v = self.numerator[n] if n < len(self.numerator) else Q(0)
            out.append(_q(v-sum(self.denominator[j]*out[n-j]
                             for j in range(1, min(n+1, len(self.denominator))))))
        return tuple(out)

    def __add__(self, other):
        return RationalGF(_add(_mul(self.numerator, other.denominator),
                               _mul(other.numerator, self.denominator)),
                          _mul(self.denominator, other.denominator))

    def __mul__(self, other):
        return RationalGF(_mul(self.numerator, other.numerator),
                          _mul(self.denominator, other.denominator))

    def weighted_indices(self):
        """OGF of n*a_n: z*d/dz, retaining index zero."""
        p, q = self.numerator, self.denominator
        dp, dq = _poly([i*p[i] for i in range(1, len(p))]), _poly([i*q[i] for i in range(1, len(q))])
        return RationalGF((0,)+_add(_mul(dp, q), tuple(-x for x in _mul(p, dq))), _mul(q, q))

    def realize(self):
        r = self.dimension
        if r > 128:
            raise ValueError('Realization dimension exceeds 128')
        if not r:
            return (), (), ()
        a = [[Q(0)]*r for _ in range(r)]
        for i in range(r-1):
            a[i][i+1] = Q(1)
        a[-1] = [-self.denominator[r-i] if r-i < len(self.denominator) else Q(0) for i in range(r)]
        return tuple(map(tuple, a)), self.coefficients(r), (Q(1),)+(Q(0),)*(r-1)

    def nth(self, n):
        if type(n) is not int or n < 0:
            raise ValueError('Nonnegative integer index required')
        a, x, c = self.realize()
        return sum(v*w for v, w in zip(c, _mv(_power(a, n), x))) if a else Q(0)

    def subsequence(self, stride, offset=0):
        if type(stride) is not int or stride < 1 or type(offset) is not int or offset < 0:
            raise ValueError('Positive integer stride and nonnegative offset required')
        a, x, c = self.realize()
        if not a:
            return self
        return from_state(_power(a, stride), _mv(_power(a, offset), x), c)

    def packet(self):
        a, x, c = self.realize()
        r = len(a)
        values = self.coefficients(2*r-1) if r else ()
        hankel = [[values[i+j] for j in range(r)] for i in range(r)]
        determinant = _det(hankel)
        if r and not determinant:
            raise ArithmeticError('Minimal realization Hankel witness is singular')
        return {'schema': 'exact.rational-gf/1', 'numerator': list(map(str, self.numerator)),
                'denominator': list(map(str, self.denominator)), 'dimension': r,
                'A': [list(map(str, row)) for row in a], 'initial': list(map(str, x)),
                'C': list(map(str, c)), 'hankel_determinant': str(determinant),
                'indexing': 'a_n = C A^n initial, n >= 0; no input feedthrough'}

    def tail(self, argument, radius, count):
        """Cauchy majorant on |z|=radius; certified absolute series tail."""
        z, radius, count = abs(_q(argument)), _q(radius), _count(count)
        if not 0 <= z < radius:
            raise ValueError('Need |argument| < radius')
        gap = 1-sum(abs(x)*radius**i for i, x in enumerate(self.denominator) if i)
        if gap <= 0:
            raise ValueError('Triangle bound cannot certify a pole-free disk at this radius')
        majorant = sum(abs(x)*radius**i for i, x in enumerate(self.numerator))/gap
        bound = _q(majorant*_rpow(z/radius,count)/(1-z/radius))
        return {'schema': 'exact.rational-gf-tail/1', 'source': self.packet(),
                'argument': str(_q(argument)), 'radius': str(radius), 'count': count,
                'majorant': str(majorant), 'absolute_tail_upper': str(bound)}


def from_recurrence(coefficients, initial):
    c, initial = tuple(map(_q, coefficients)), tuple(map(_q, initial))
    if not c or len(c) != len(initial):
        raise ValueError('Matching positive order and initial prefix required')
    q = (Q(1),)+tuple(-x for x in reversed(c))
    p = [sum(q[j]*initial[i-j] for j in range(i+1)) for i in range(len(c))]
    return RationalGF(tuple(p), q)


def from_state(A, initial, C):
    """Faddeev-LeVerrier exact denominator, followed by full cancellation."""
    a, x, c = _matrix(A), tuple(map(_q, initial)), tuple(map(_q, C))
    n = len(a)
    if n > 128 or len(a[0]) != n or len(x) != n or len(c) != n:
        raise ValueError('Square state matrix, matching vectors, dimension <= 128 required')
    b, denominator = _identity(n), [Q(1)]
    for k in range(1, n+1):
        b = _mm(a, b)
        v = -sum(b[i][i] for i in range(n))/k
        denominator.append(v)
        b = tuple(tuple(b[i][j]+(v if i == j else 0) for j in range(n)) for i in range(n))
    values = []
    for _ in range(n):
        values.append(sum(v*w for v, w in zip(c, x)))
        x = _mv(a, x)
    numerator = [sum(denominator[j]*values[i-j] for j in range(i+1)) for i in range(n)]
    return RationalGF(tuple(numerator), tuple(denominator))


def periodic_gf(prefix, cycle):
    prefix, cycle = tuple(map(_q, prefix)), tuple(map(_q, cycle))
    if not cycle:
        raise ValueError('Nonempty periodic cycle required')
    return RationalGF(prefix or (0,)) + RationalGF((Q(0),)*len(prefix)+cycle,
                                                  (Q(1),)+(Q(0),)*(len(cycle)-1)+(Q(-1),))


def binomial_tail(alpha, argument, count):
    """Enclose (1-z)^(-alpha) for rational alpha>0 and 0<=z<1.

    Bound each future coefficient ratio, not just its asymptotic growth.
    Count is the number of retained coefficients, including the constant.
    """
    a, z, count = _q(alpha), _q(argument), _count(count)
    if a <= 0 or not 0 <= z < 1:
        raise ValueError('Positive alpha and argument in [0,1) required')
    coefficient, partial, power = Q(1), Q(0), Q(1)
    for n in range(count):
        partial = _q(partial+coefficient*power)
        power = _q(power*z)
        coefficient = _q(coefficient*(a+n)/(n+1))
    ratio = z*max(Q(1), (count+a)/(count+1))
    if ratio >= 1:
        raise ValueError('Retain more coefficients to make the geometric majorant contract')
    tail = _q(coefficient*power/(1-ratio))
    return {'schema': 'exact.binomial-tail/1', 'alpha': str(a), 'argument': str(z),
            'count': count, 'lower': str(partial), 'upper': str(partial+tail),
            'tail_upper': str(tail), 'future_ratio_upper': str(ratio),
            'scope': 'positive real binomial branch (1-z)^(-alpha)'}


def replay_generating(packet):
    schema = packet.get('schema')
    if schema == 'exact.rational-gf/1':
        expected = RationalGF(packet['numerator'], packet['denominator']).packet()
    elif schema == 'exact.rational-gf-tail/1':
        source = packet['source']
        replay_generating(source)
        expected = RationalGF(source['numerator'], source['denominator']).tail(
            packet['argument'], packet['radius'], packet['count'])
    elif schema == 'exact.binomial-tail/1':
        expected = binomial_tail(packet['alpha'], packet['argument'], packet['count'])
    else:
        raise ValueError('Unknown generating function schema')
    if packet != expected:
        raise ValueError('Generating function certificate does not reconstruct')
    return True
