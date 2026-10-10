"""Observable compilation without materializing a dense transition matrix."""
from dataclasses import dataclass
from .rational import q, matrix, row_basis, coordinates, mul, encoded, mv, addv


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
    b, c = matrix(B), matrix(C)
    if len(b) != operator.n or len(c[0]) != operator.n or max_rank < 1:
        raise ValueError("Invalid input/target dimensions or rank limit")
    e = row_basis(c)
    if not e:
        raise ValueError("Need a nonzero target")
    ranks = [len(e)]
    while len(e) <= max_rank:
        ea = operator.pullback(e)
        expanded = row_basis(e+ea)
        if len(expanded) == len(e):
            return {"schema": "dynacomp.sparse/1", "n": operator.n,
                    "entries": [[i, j, str(v)] for i, j, v in operator.entries],
                    "B": encoded(b), "C": encoded(c), "E": encoded(e),
                    "G": encoded([coordinates(row, e) for row in ea]),
                    "H": encoded(mul(e, b)), "D": encoded([coordinates(row, e) for row in c]),
                    "ranks": ranks}
        e = expanded
        ranks.append(len(e))
    raise OverflowError("Observable rank exceeded sparse compilation limit")


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
