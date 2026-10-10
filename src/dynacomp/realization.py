"""Minimal input-output realization with an independent Hankel-rank certificate.

The observable quotient keeps every distinction any initial state could make.
When the system starts at rest (x0 = 0) only the reachable part matters, and
the order can drop further: to the rank of the block Hankel matrix of Markov
parameters C A^k B, which no realization of the same input-output map can beat.
"""
from .krylov import krylov, vecmat, matvec
from .linear import LinearModel
from .rational import matrix, mul, coordinates, encoded, rref, transpose, q, mv, addv


def _block_rank(rows):
    return len(rref(rows)[1]) if rows else 0


def minimal_realization(model):
    """Kalman-minimal (G, H, D) with identical Markov parameters.

    Stage 1 keeps the reachable column space P: A P = P F, B = P K.
    Stage 2 takes the observable row space of (F, K, CP): E F = G E, E K = H,
    D E = C P. Then C A^k B = D G^k H for every k.
    """
    a, b, c = model.A, model.B, model.C
    m, p = len(b[0]), len(c)
    columns = [(j, tuple(row[j] for row in b)) for j in range(m)]
    reach, reach_labels = krylov(columns, {"A": lambda v: matvec(a, v)})
    if not reach:
        raise ValueError("No input reaches the state; the input-output map is zero")
    P = transpose(reach)
    F = transpose(tuple(coordinates(matvec(a, v), reach) for v in reach))
    K = transpose(tuple(coordinates(col, reach) for _, col in columns))
    cp = mul(c, P)
    rows = [(i, cp[i]) for i in range(p)]
    obs, obs_labels = krylov(rows, {"A": lambda v: vecmat(v, F)})
    if not obs:
        raise ValueError("No output observes a reachable state; the input-output map is zero")
    E = obs
    G = tuple(coordinates(vecmat(row, F), E) for row in E)
    H = mul(E, K)
    D = tuple(coordinates(row, E) for row in cp)
    s = len(E)
    # Hankel certificate: experiments from the minimal system's own Krylov spans.
    hrow_vectors, hrow_labels = krylov([(i, D[i]) for i in range(p)], {"A": lambda v: vecmat(v, G)})
    hcol_vectors, hcol_labels = krylov([(j, tuple(row[j] for row in H)) for j in range(m)],
                                       {"A": lambda v: matvec(G, v)})
    assert len(hrow_vectors) == len(hcol_vectors) == s
    hankel_rows = [[i, len(w)] for i, w in hrow_labels]
    hankel_cols = [[j, len(w)] for j, w in hcol_labels]
    return {"schema": "dynacomp.realization/1", "arithmetic": "rational",
            "A": encoded(a), "B": encoded(b), "C": encoded(c),
            "P": encoded(P), "F": encoded(F), "K": encoded(K),
            "E": encoded(E), "G": encoded(G), "H": encoded(H), "D": encoded(D),
            "reachable_dimension": len(reach), "minimal_order": s,
            "hankel_rows": hankel_rows, "hankel_columns": hankel_cols,
            "claim": "zero-initial-state input-output map preserved; order equals Hankel rank"}


def markov_parameter(a, b, c, k):
    """C A^k B, computed by k vector-matrix products per output row."""
    rows = c
    for _ in range(k):
        rows = tuple(vecmat(r, a) for r in rows)
    return mul(rows, b)


def hankel_submatrix(a, b, c, row_labels, col_labels):
    """Entries C_i A^(k+l) B_j, building C A^t once for t up to the largest k+l."""
    top = max(k for _, k in row_labels) + max(l for _, l in col_labels)
    rows, markov = c, []
    for _ in range(top+1):
        markov.append(mul(rows, b))
        rows = tuple(vecmat(r, a) for r in rows)
    return tuple(tuple(markov[k+l][i][j] for j, l in col_labels) for i, k in row_labels)


def replay_realization(packet):
    """Accept identities, then prove minimality from the ORIGINAL Markov parameters."""
    if packet.get("schema") != "dynacomp.realization/1":
        raise ValueError("Unknown schema")
    model = LinearModel.from_rows(packet["A"], packet["B"], packet["C"])
    P, F, K = (matrix(packet[k]) for k in ("P", "F", "K"))
    E, G, H, D = (matrix(packet[k]) for k in ("E", "G", "H", "D"))
    n, m, p = len(model.A), len(model.B[0]), len(model.C)
    r, s = len(F), len(G)
    if len(P) != n or len(P[0]) != r or len(F[0]) != r or len(K) != r or len(K[0]) != m:
        raise ValueError("Reachability block dimensions do not match")
    if len(E) != s or len(E[0]) != r or len(G[0]) != s or len(H) != s or len(D) != p or len(D[0]) != s:
        raise ValueError("Observability block dimensions do not match")
    if mul(model.A, P) != mul(P, F) or mul(P, K) != model.B:
        raise ValueError("Reachability identities AP=PF, B=PK fail")
    if mul(E, F) != mul(G, E) or mul(E, K) != H or mul(D, E) != mul(model.C, P):
        raise ValueError("Observability identities EF=GE, EK=H, DE=CP fail")
    rows = [tuple(x) for x in packet["hankel_rows"]]
    cols = [tuple(x) for x in packet["hankel_columns"]]
    if len(rows) != s or len(cols) != s:
        raise ValueError("Hankel certificate must be square of minimal order")
    for i, k in rows:
        if type(i) is not int or type(k) is not int or not 0 <= i < p or k < 0:
            raise ValueError("Invalid Hankel row label")
    for j, l in cols:
        if type(j) is not int or type(l) is not int or not 0 <= j < m or l < 0:
            raise ValueError("Invalid Hankel column label")
    if _block_rank(hankel_submatrix(model.A, model.B, model.C, rows, cols)) != s:
        raise ValueError("Hankel submatrix is singular: minimality not established")
    return True


def simulate_zero_state(G, H, D, inputs):
    z = tuple(q(0) for _ in G)
    outputs = [mv(D, z)]
    for u in inputs:
        z = addv(mv(G, z), mv(H, u))
        outputs.append(mv(D, z))
    return tuple(outputs)
