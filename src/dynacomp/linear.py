"""Compile minimal target-preserving rational linear dynamics."""
from dataclasses import dataclass
from .rational import matrix, row_basis, coordinates, mul, mv, nullspace, encoded, q, addv


@dataclass(frozen=True)
class LinearModel:
    A: tuple
    B: tuple
    C: tuple

    @classmethod
    def from_rows(cls, A, B, C):
        a, b, c = matrix(A), matrix(B), matrix(C)
        n = len(a)
        if len(a[0]) != n or len(b) != n or len(c[0]) != n:
            raise ValueError("Expected A[n,n], B[n,m], C[p,n]")
        return cls(a, b, c)


@dataclass(frozen=True)
class LinearReduction:
    model: LinearModel
    E: tuple
    G: tuple
    H: tuple
    D: tuple
    rounds: tuple

    def packet(self):
        return {"schema": "dynacomp.linear/1", "arithmetic": "rational",
                **{k: encoded(getattr(self.model, k)) for k in ("A", "B", "C")},
                **{k: encoded(getattr(self, k)) for k in ("E", "G", "H", "D")},
                "rounds": list(self.rounds)}

    def simulate(self, x0, inputs):
        z = mv(self.E, tuple(q(v) for v in x0))
        outputs = [mv(self.D, z)]
        for u in inputs:
            z = addv(mv(self.G, z), mv(self.H, u))
            outputs.append(mv(self.D, z))
        return tuple(outputs)


def compile_linear(model):
    """Smallest A-invariant row space containing C, over all initial states."""
    basis = row_basis(model.C)
    if not basis:
        raise ValueError("Zero-only targets need no state; request a nonzero target")
    rounds = []
    while True:
        next_rows = mul(basis, model.A)
        enlarged = row_basis(basis + next_rows)
        if len(enlarged) == len(basis):
            break
        witness = closure_witness(model.A, basis)
        rounds.append({"before": len(basis), "after": len(enlarged), "witness": witness})
        basis = enlarged
    g = tuple(coordinates(row, basis) for row in mul(basis, model.A))
    d = tuple(coordinates(row, basis) for row in model.C)
    return LinearReduction(model, basis, g, mul(basis, model.B), d, tuple(rounds))


def closure_witness(A, E):
    """Same encoded state, different next encoded state: 0 and v."""
    a, e = matrix(A), matrix(E)
    for v in nullspace(e):
        difference = mv(mul(e, a), v)
        if any(difference):
            return {"v": [str(x) for x in v], "Ev": [str(x) for x in mv(e, v)],
                    "EAv": [str(x) for x in difference]}
    return None


def distinguishing_horizon(model, v):
    """Earliest output difference; None proves all-future equality by Cayley–Hamilton."""
    x = tuple(q(w) for w in v)
    for t in range(len(model.A)):
        y = mv(model.C, x)
        if any(y):
            return {"horizon": t, "difference": [str(w) for w in y]}
        x = mv(model.A, x)
    return None


def replay_linear(packet):
    """Independent acceptance by three exact identities, no compiler call."""
    if packet.get("schema") != "dynacomp.linear/1":
        raise ValueError("Unknown certificate schema")
    model = LinearModel.from_rows(packet["A"], packet["B"], packet["C"])
    e, g, h, d = (matrix(packet[k]) for k in ("E", "G", "H", "D"))
    if len(e[0]) != len(model.A) or len(row_basis(e)) != len(e):
        raise ValueError("E must have full row rank and correct width")
    if mul(e, model.A) != mul(g, e):
        raise ValueError("Transition identity EA=GE fails")
    if mul(e, model.B) != h:
        raise ValueError("Input identity EB=H fails")
    if mul(d, e) != model.C:
        raise ValueError("Target identity DE=C fails")
    return True
