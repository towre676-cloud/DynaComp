"""Labeled breadth-first Krylov spans over Q.

Every retained vector remembers the word that produced it, so a basis is also
an experiment design: row c_i A_w is "observe output i after word w".
"""
from collections import deque
from .rational import q


class Echelon:
    """Incremental exact independence test; rows kept in insertion order."""

    def __init__(self):
        self.rows = []

    def reduce(self, vector):
        v = [q(x) for x in vector]
        for pivot, row in self.rows:
            c = v[pivot]
            if c:
                v = [a-c*b for a, b in zip(v, row)]
        return v

    def add(self, vector):
        v = self.reduce(vector)
        pivot = next((j for j, x in enumerate(v) if x), None)
        if pivot is None:
            return False
        c = v[pivot]
        self.rows.append((pivot, [x/c for x in v]))
        return True

    def __len__(self):
        return len(self.rows)


def vecmat(v, a):
    """Row vector times matrix."""
    if len(v) != len(a):
        raise ValueError("Vector dimension does not match")
    return tuple(sum(v[i]*a[i][j] for i in range(len(v)) if v[i]) for j in range(len(a[0])))


def matvec(a, v):
    if len(a[0]) != len(v):
        raise ValueError("Vector dimension does not match")
    return tuple(sum(x*y for x, y in zip(row, v) if y) for row in a)


def krylov(starts, maps, limit=None, visit=None):
    """Span of {s . w} over start labels s and words w, in shortlex order.

    starts: [(label, vector)], maps: {letter: vector -> vector} (dict order is
    the alphabet order). Returns (vectors, labels) with labels (start, word).
    ``visit(label, vector)`` is called on every enqueued word before the
    independence test and may return a non-None value to stop early.
    """
    echelon, vectors, labels = Echelon(), [], []
    queue = deque((label, (), tuple(q(x) for x in v)) for label, v in starts)
    while queue:
        start, word, v = queue.popleft()
        if visit is not None:
            stop = visit((start, word), v)
            if stop is not None:
                return stop
        if not echelon.add(v):
            continue
        vectors.append(v)
        labels.append((start, word))
        if limit is not None and len(vectors) > limit:
            raise OverflowError("Krylov dimension exceeded limit")
        for letter, apply in maps.items():
            queue.append((start, word+(letter,), apply(v)))
    return tuple(vectors), tuple(labels)
