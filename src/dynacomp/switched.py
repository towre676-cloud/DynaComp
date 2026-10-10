"""Multi-mode exact closure: switched linear systems and weighted automata over Q.

A switched system applies one of several updates A_s at each step. A reduction
must survive every mode sequence, so its row space must be invariant under all
A_s simultaneously. A weighted automaton (alpha, {A_s}, eta) assigns the word
w = s1...sk the weight alpha A_s1 ... A_sk eta; probabilistic automata, hidden
Markov models and rational series are special cases. The same span closure
yields Schutzenberger minimization and Tzeng-style equivalence checking with
a shortest counterexample word.
"""
from .krylov import krylov, vecmat, matvec
from .rational import matrix, mul, coordinates, encoded, rref, q


def _modes(modes, n=None):
    if not modes:
        raise ValueError("Need at least one mode")
    out = {}
    for name, a in modes.items():
        a = matrix(a)
        if len(a) != len(a[0]) or (n is not None and len(a) != n):
            raise ValueError("Mode matrices must be square and share a dimension")
        n = len(a)
        out[str(name)] = a
    return out, n


def compile_switched(modes, C, inputs=None):
    """Smallest row space containing C and invariant under every mode.

    Each encoder row is c_i A_w for a recorded word w: the experiment it
    represents. Identities: E A_s = G_s E, E B_s = H_s, D E = C.
    """
    a, n = _modes(modes)
    c = matrix(C)
    if len(c[0]) != n:
        raise ValueError("Target width must match state dimension")
    b = {}
    if inputs is not None:
        for name in a:
            b[name] = matrix(inputs[name])
            if len(b[name]) != n:
                raise ValueError("Input matrix height must match state dimension")
    rows, labels = krylov([(i, row) for i, row in enumerate(c)],
                          {s: (lambda v, m=m: vecmat(v, m)) for s, m in a.items()})
    if not rows:
        raise ValueError("Zero-only targets need no state; request a nonzero target")
    G = {s: [list(coordinates(vecmat(r, m), rows)) for r in rows] for s, m in a.items()}
    packet = {"schema": "dynacomp.switched/1", "arithmetic": "rational",
              "modes": {s: encoded(m) for s, m in a.items()}, "C": encoded(c),
              "E": encoded(rows), "G": {s: encoded(g) for s, g in G.items()},
              "D": encoded([coordinates(r, rows) for r in c]),
              "experiments": [{"output": i, "word": list(w)} for i, w in labels]}
    if b:
        packet["inputs"] = {s: encoded(m) for s, m in b.items()}
        packet["H"] = {s: encoded(mul(rows, m)) for s, m in b.items()}
    return packet


def replay_switched(packet):
    if packet.get("schema") != "dynacomp.switched/1":
        raise ValueError("Unknown schema")
    a, n = _modes(packet["modes"])
    c, e, d = matrix(packet["C"]), matrix(packet["E"]), matrix(packet["D"])
    if len(e[0]) != n or len(rref(e)[1]) != len(e):
        raise ValueError("E must have full row rank and correct width")
    if set(packet["G"]) != set(a):
        raise ValueError("Need one reduced update per mode")
    for s, m in a.items():
        if mul(e, m) != mul(matrix(packet["G"][s]), e):
            raise ValueError(f"Transition identity fails for mode {s}")
    if "inputs" in packet:
        if set(packet["inputs"]) != set(a) or set(packet["H"]) != set(a):
            raise ValueError("Need one input matrix per mode")
        for s in a:
            if mul(e, matrix(packet["inputs"][s])) != matrix(packet["H"][s]):
                raise ValueError(f"Input identity fails for mode {s}")
    if mul(d, e) != c:
        raise ValueError("Target identity DE=C fails")
    return True


class WeightedAutomaton:
    def __init__(self, alpha, modes, eta):
        self.modes, n = _modes(modes)
        self.alpha = tuple(q(x) for x in alpha)
        self.eta = tuple(q(x) for x in eta)
        if len(self.alpha) != n or len(self.eta) != n:
            raise ValueError("Initial and final vectors must match dimension")
        self.n = n

    @classmethod
    def from_spec(cls, spec):
        return cls(spec["alpha"], spec["modes"], spec["eta"])

    def weight(self, word):
        v = self.alpha
        for s in word:
            v = vecmat(v, self.modes[s])
        return sum(x*y for x, y in zip(v, self.eta))

    def spec(self):
        return {"alpha": [str(x) for x in self.alpha], "eta": [str(x) for x in self.eta],
                "modes": {s: encoded(m) for s, m in self.modes.items()}}


def minimize_automaton(aut):
    """Forward (reachable) then backward (observable) quotient.

    Result order equals the Hankel rank of the series, so it is minimal among
    all weighted automata over Q with the same word weights. The packet carries
    prefix and suffix words whose Hankel block is nonsingular.
    """
    fwd, _ = krylov([(0, aut.alpha)], {s: (lambda v, m=m: vecmat(v, m)) for s, m in aut.modes.items()})
    if not fwd:
        return {"schema": "dynacomp.automaton/1", "status": "zero_series", "automaton": aut.spec()}
    M = {s: tuple(coordinates(vecmat(v, m), fwd) for v in fwd) for s, m in aut.modes.items()}
    a1 = coordinates(aut.alpha, fwd)
    eta1 = matvec(fwd, aut.eta)
    bwd, _ = krylov([(0, eta1)], {s: (lambda v, m=m: matvec(m, v)) for s, m in M.items()})
    if not bwd:
        return {"schema": "dynacomp.automaton/1", "status": "zero_series", "automaton": aut.spec()}
    # Columns P = bwd^T: M_s P = P N_s, eta1 = P e.
    N = {s: tuple(zip(*(coordinates(matvec(m, v), bwd) for v in bwd))) for s, m in M.items()}
    e = coordinates(eta1, bwd)
    alpha_min = tuple(sum(a1[i]*bwd[k][i] for i in range(len(a1))) for k in range(len(bwd)))
    minimal = WeightedAutomaton(alpha_min, N, e)
    prefixes = krylov([(0, minimal.alpha)], {s: (lambda v, m=m: vecmat(v, m)) for s, m in N.items()})[1]
    suffixes = krylov([(0, minimal.eta)], {s: (lambda v, m=m: matvec(m, v)) for s, m in N.items()})[1]
    return {"schema": "dynacomp.automaton/1", "status": "minimized", "arithmetic": "rational",
            "automaton": aut.spec(), "forward_basis": encoded(fwd),
            "forward_modes": {s: encoded(m) for s, m in M.items()},
            "forward_alpha": [str(x) for x in a1], "forward_eta": [str(x) for x in eta1],
            "backward_basis": encoded(bwd), "minimal": minimal.spec(),
            "order": len(bwd), "original_order": aut.n,
            "hankel_prefixes": [list(w) for _, w in prefixes],
            "hankel_suffixes": [list(reversed(w)) for _, w in suffixes]}


def replay_automaton(packet):
    """Check both quotient identities and the Hankel lower bound on the original."""
    if packet.get("schema") != "dynacomp.automaton/1" or packet.get("status") != "minimized":
        raise ValueError("Need a minimized automaton packet")
    aut = WeightedAutomaton.from_spec(packet["automaton"])
    small = WeightedAutomaton.from_spec(packet["minimal"])
    fwd, bwd = matrix(packet["forward_basis"]), matrix(packet["backward_basis"])
    M, _ = _modes(packet["forward_modes"])
    a1 = tuple(q(x) for x in packet["forward_alpha"])
    eta1 = tuple(q(x) for x in packet["forward_eta"])
    r, s = len(fwd), len(bwd)
    if set(M) != set(aut.modes) or set(small.modes) != set(aut.modes):
        raise ValueError("Alphabets differ")
    if len(fwd[0]) != aut.n or len(bwd[0]) != r or len(a1) != r or len(eta1) != r or small.n != s:
        raise ValueError("Certificate dimensions do not match")
    if vecmat(a1, fwd) != aut.alpha or matvec(fwd, aut.eta) != eta1:
        raise ValueError("Forward boundary identities fail")
    for x, m in aut.modes.items():
        if mul(fwd, m) != mul(M[x], fwd):
            raise ValueError("Forward identity fails")
    P = tuple(zip(*bwd))
    for x, m in M.items():
        if mul(m, P) != mul(P, small.modes[x]):
            raise ValueError("Backward identity fails")
    if matvec(P, small.eta) != eta1 or vecmat(a1, P) != small.alpha:
        raise ValueError("Backward boundary identities fail")
    pre = [tuple(w) for w in packet["hankel_prefixes"]]
    suf = [tuple(w) for w in packet["hankel_suffixes"]]
    if len(pre) != s or len(suf) != s or any(x not in aut.modes for w in pre+suf for x in w):
        raise ValueError("Hankel words must be square of minimal order over the alphabet")
    hankel = tuple(tuple(aut.weight(u+v) for v in suf) for u in pre)
    if len(rref(hankel)[1]) != s:
        raise ValueError("Hankel block is singular: minimality not established")
    return True


def check_equivalence(first, second):
    """Decide whether two weighted automata assign every word the same weight.

    Explores the joint row span of (alpha1 A1_w, alpha2 A2_w) in shortlex order,
    testing every visited word. The first difference found is a shortest
    distinguishing word. Otherwise the span itself is an equivalence certificate.
    """
    if set(first.modes) != set(second.modes):
        raise ValueError("Automata must share an alphabet")
    letters = list(first.modes)
    n1 = first.n
    final = first.eta + tuple(-x for x in second.eta)

    def step(v, s):
        return vecmat(v[:n1], first.modes[s]) + vecmat(v[n1:], second.modes[s])

    def visit(label, v):
        if sum(x*y for x, y in zip(v, final)):
            return label
        return None

    result = krylov([(0, first.alpha + second.alpha)],
                    {s: (lambda v, s=s: step(v, s)) for s in letters}, visit=visit)
    pair = {"first": first.spec(), "second": second.spec()}
    if isinstance(result[0], int):
        word = list(result[1])
        return {"schema": "dynacomp.equivalence/1", "equivalent": False, **pair, "word": word,
                "weights": [str(first.weight(word)), str(second.weight(word))]}
    basis = result[0]
    if not basis:
        return {"schema": "dynacomp.equivalence/1", "equivalent": True, **pair,
                "basis": [], "G": {s: [] for s in letters}, "start": []}
    G = {s: encoded([coordinates(step(v, s), basis) for v in basis]) for s in letters}
    return {"schema": "dynacomp.equivalence/1", "equivalent": True, **pair,
            "basis": encoded(basis), "G": G,
            "start": [str(x) for x in coordinates(first.alpha + second.alpha, basis)]}


def replay_equivalence(packet):
    if packet.get("schema") != "dynacomp.equivalence/1":
        raise ValueError("Unknown schema")
    first = WeightedAutomaton.from_spec(packet["first"])
    second = WeightedAutomaton.from_spec(packet["second"])
    if set(first.modes) != set(second.modes):
        raise ValueError("Automata must share an alphabet")
    if not packet["equivalent"]:
        word = packet["word"]
        if any(s not in first.modes for s in word) or first.weight(word) == second.weight(word):
            raise ValueError("Counterexample word does not distinguish")
        return True
    if not packet["basis"]:
        if any(first.alpha) or any(second.alpha):
            raise ValueError("Empty span certifies only two zero initial vectors")
        return True
    basis = matrix(packet["basis"])
    n1, n = first.n, first.n + second.n
    if len(basis[0]) != n:
        raise ValueError("Basis width must be the joint dimension")
    start = tuple(q(x) for x in packet["start"])
    if vecmat(start, basis) != first.alpha + second.alpha:
        raise ValueError("Start vector is not in the certified span")
    final = first.eta + tuple(-x for x in second.eta)
    if any(matvec(basis, final)):
        raise ValueError("Certified span does not annihilate the weight difference")
    for s in first.modes:
        g = matrix(packet["G"][s])
        for row, coefficients in zip(basis, g):
            image = vecmat(row[:n1], first.modes[s]) + vecmat(row[n1:], second.modes[s])
            if vecmat(coefficients, basis) != image:
                raise ValueError(f"Span is not invariant under letter {s}")
    return True
