"""Observable compilation without materializing a dense transition matrix."""
from dataclasses import dataclass
from .rational import q, matrix, row_basis, rref, pivot_coordinates, mul, encoded, mv, addv
from .krylov import krylov


@dataclass(frozen=True)
class SparseOperator:
    n: int
    entries: tuple

    @classmethod
    def from_entries(cls, n, entries):
        if type(n) is not int or n < 1:
            raise ValueError("Positive operator dimension required")
        values = {}
        for i, j, value in entries:
            if type(i) is not int or type(j) is not int or not 0 <= i < n or not 0 <= j < n:
                raise ValueError("Invalid sparse index")
            values[i, j] = values.get((i, j), q(0)) + q(value)
        return cls(n, tuple((i, j, value) for (i, j), value in sorted(values.items()) if value))

    def pullback(self, rows):
        if any(len(row) != self.n for row in rows):
            raise ValueError("Wrong row width")
        out = []
        for row in rows:
            result = [q(0)]*self.n
            for i, j, value in self.entries:
                result[j] += row[i]*value
            out.append(tuple(result))
        return tuple(out)

    def step(self, x):
        if len(x) != self.n:
            raise ValueError("Wrong state size")
        out = [q(0)]*self.n
        for i, j, value in self.entries:
            out[i] += value*q(x[j])
        return tuple(out)


def compile_sparse(operator, B, C, max_rank=128):
    """Observable closure by labeled Krylov search, then one final reduction.

    The incremental echelon tests each pullback once instead of re-reducing
    the whole span every round; coordinates in the final reduced basis are
    read off its pivot columns.
    """
    b, c = matrix(B), matrix(C)
    if len(b) != operator.n or len(c[0]) != operator.n or max_rank < 1:
        raise ValueError("Invalid input/target dimensions or rank limit")
    try:
        rows, labels = krylov([(i, row) for i, row in enumerate(c)],
                              {"A": lambda v: operator.pullback((v,))[0]}, limit=max_rank)
    except OverflowError:
        raise OverflowError("Observable rank exceeded sparse compilation limit") from None
    if not rows:
        raise ValueError("Need a nonzero target")
    e, pivots = rref(rows)
    ea = operator.pullback(e)
    ranks = [len(row_basis(c))]
    for depth in range(1, max(len(w) for _, w in labels)+1):
        ranks.append(sum(1 for _, w in labels if len(w) <= depth))
    return {"schema": "dynacomp.sparse/1", "n": operator.n,
            "entries": [[i, j, str(v)] for i, j, v in operator.entries],
            "B": encoded(b), "C": encoded(c), "E": encoded(e),
            "G": encoded([pivot_coordinates(row, e, pivots) for row in ea]),
            "H": encoded(mul(e, b)), "D": encoded([pivot_coordinates(row, e, pivots) for row in c]),
            "ranks": ranks}


def replay_sparse(packet):
    if packet.get("schema") != "dynacomp.sparse/1":
        raise ValueError("Unknown schema")
    op = SparseOperator.from_entries(packet["n"], packet["entries"])
    b, c, e, g, h, d = (matrix(packet[k]) for k in ("B", "C", "E", "G", "H", "D"))
    if len(b) != op.n or len(e[0]) != op.n or len(row_basis(e)) != len(e):
        raise ValueError("Invalid dimensions or dependent encoder")
    if op.pullback(e) != mul(g, e) or mul(e, b) != h or mul(d, e) != c:
        raise ValueError("Sparse reduction identities fail")
    return True
