"""Small exact linear algebra over Q. No numerical rank decisions."""
from fractions import Fraction


def q(value):
    if isinstance(value, float):
        raise TypeError("Use rational strings; binary floats are not exact inputs")
    return Fraction(value)


def matrix(rows):
    result = tuple(tuple(q(x) for x in row) for row in rows)
    if not result or not result[0] or any(len(r) != len(result[0]) for r in result):
        raise ValueError("Expected a nonempty rectangular matrix")
    return result


def identity(n):
    return matrix([[int(i == j) for j in range(n)] for i in range(n)])


def zeros(m, n):
    return matrix([[0] * n for _ in range(m)])


def transpose(a):
    return tuple(zip(*a))


def mul(a, b):
    if len(a[0]) != len(b):
        raise ValueError("Matrix dimensions do not match")
    return tuple(tuple(sum(x*y for x, y in zip(row, col)) for col in transpose(b)) for row in a)


def sub(a, b):
    if len(a) != len(b) or len(a[0]) != len(b[0]):
        raise ValueError("Matrix dimensions do not match")
    return tuple(tuple(x-y for x, y in zip(r, s)) for r, s in zip(a, b))


def mv(a, x):
    if len(a[0]) != len(x):
        raise ValueError("Vector dimension does not match")
    return tuple(sum(v*q(w) for v, w in zip(row, x)) for row in a)


def addv(x, y):
    if len(x) != len(y):
        raise ValueError("Vector dimensions do not match")
    return tuple(a+b for a, b in zip(x, y))


def rref(rows):
    a = [list(row) for row in matrix(rows)]
    m, n = len(a), len(a[0])
    pivots = []
    for j in range(n):
        k = next((i for i in range(len(pivots), m) if a[i][j]), None)
        if k is None:
            continue
        i = len(pivots)
        a[i], a[k] = a[k], a[i]
        divisor = a[i][j]
        a[i] = [v/divisor for v in a[i]]
        for k in range(m):
            if k != i:
                coefficient = a[k][j]
                a[k] = [v-coefficient*w for v, w in zip(a[k], a[i])]
        pivots.append(j)
        if len(pivots) == m:
            break
    return tuple(tuple(row) for row in a), tuple(pivots)


def row_basis(rows):
    reduced, pivots = rref(rows)
    return reduced[:len(pivots)]


def nullspace(a):
    reduced, pivots = rref(a)
    vectors = []
    for j in range(len(a[0])):
        if j not in pivots:
            v = [q(0)] * len(a[0])
            v[j] = q(1)
            for i, k in enumerate(pivots):
                v[k] = -reduced[i][j]
            vectors.append(tuple(v))
    return tuple(vectors)


def coordinates(row, basis):
    """Express a row in a full-row-rank basis, or return None."""
    reduced, pivots = rref(basis)
    if len(pivots) != len(basis):
        raise ValueError("Basis must have independent rows")
    # Solve B^T c = row^T using an augmented system.
    augmented = tuple(tuple(col) + (q(target),) for col, target in zip(transpose(basis), row))
    result, _ = rref(augmented)
    d = len(basis)
    if any(not any(r[:d]) and r[d] for r in result):
        return None
    answer = [q(0)] * d
    for r in result:
        pivot = next((j for j in range(d) if r[j]), None)
        if pivot is not None:
            answer[pivot] = r[d]
    return tuple(answer)


def norm_inf(a):
    return max(sum(abs(v) for v in row) for row in a)


def vector_norm(x):
    return max((abs(q(v)) for v in x), default=q(0))


def encoded(a):
    return [[str(v) for v in row] for row in a]


def inverse(a):
    n = len(a)
    if len(a[0]) != n:
        raise ValueError("Expected square matrix")
    augmented = tuple(tuple(row)+tuple(int(i == j) for j in range(n)) for i, row in enumerate(a))
    reduced, pivots = rref(augmented)
    if pivots[:n] != tuple(range(n)):
        raise ValueError("Singular matrix")
    return tuple(row[n:] for row in reduced)


def pivot_coordinates(row, basis, pivots):
    """Coordinates in a reduced row echelon basis are its pivot entries; None if outside the span."""
    c = tuple(q(row[p]) for p in pivots)
    rebuilt = [sum(ci*b[j] for ci, b in zip(c, basis) if ci) for j in range(len(row))]
    return c if all(q(x) == y for x, y in zip(row, rebuilt)) else None
