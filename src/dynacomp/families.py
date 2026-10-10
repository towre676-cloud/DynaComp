"""Common exact quotients for coefficient families and switching systems."""
from .linear import LinearModel, compile_linear
from .rational import matrix, row_basis, mul, coordinates, encoded, q, transpose, inverse, sub


def compile_family(models):
    models = tuple(models)
    if not models:
        raise ValueError("Nonempty coefficient family required")
    first = models[0]
    n, m = len(first.A), len(first.B[0])
    if any(len(model.A) != n or len(model.B[0]) != m or model.C != first.C for model in models):
        raise ValueError("Models must share dimensions and target")
    e = row_basis(first.C)
    if not e:
        raise ValueError("Need nonzero target")
    ranks = [len(e)]
    while True:
        rows = e + tuple(row for model in models for row in mul(e, model.A))
        enlarged = row_basis(rows)
        if len(enlarged) == len(e):
            break
        e = enlarged
        ranks.append(len(e))
    return {"schema": "dynacomp.family/1", "E": encoded(e),
            "D": encoded([coordinates(row, e) for row in first.C]), "C": encoded(first.C),
            "A": [encoded(model.A) for model in models], "B": [encoded(model.B) for model in models],
            "G": [encoded([coordinates(row, e) for row in mul(e, model.A)]) for model in models],
            "H": [encoded(mul(e, model.B)) for model in models], "ranks": ranks,
            "coverage": "arbitrary linear coefficient combinations, including time-varying coefficients"}


def replay_family(packet):
    if packet.get("schema") != "dynacomp.family/1":
        raise ValueError("Unknown schema")
    from .linear import replay_linear
    count = len(packet["A"])
    if not count or any(len(packet[k]) != count for k in ("B", "G", "H")):
        raise ValueError("Family counts disagree")
    for i in range(count):
        replay_linear({"schema": "dynacomp.linear/1", "A": packet["A"][i], "B": packet["B"][i],
                       "C": packet["C"], "E": packet["E"], "D": packet["D"],
                       "G": packet["G"][i], "H": packet["H"][i]})
    return True


def closure_conditions(A_coefficients, E):
    """Exact affine equations for where a fixed row encoding closes.

    A(theta)=A0+theta1*A1+... . Each returned equation is a coefficient
    list [constant,theta1,...] required to vanish. No sampling is involved.
    """
    coeffs = [matrix(a) for a in A_coefficients]
    e = matrix(E)
    if not coeffs or any(len(a) != len(e[0]) or len(a[0]) != len(e[0]) for a in coeffs):
        raise ValueError("Square coefficient operators with matching encoder required")
    et = transpose(e)
    v = mul(et, inverse(mul(e, et)))
    defects = []
    for a in coeffs:
        ea = mul(e, a)
        defects.append(sub(ea, mul(mul(ea, v), e)))
    equations = []
    for i in range(len(e)):
        for j in range(len(e[0])):
            equation = tuple(d[i][j] for d in defects)
            if any(equation) and equation not in equations:
                equations.append(equation)
    return {"equations": [[str(v) for v in row] for row in equations],
            "parameters": len(coeffs)-1, "meaning": "every affine equation must vanish for closure"}


def specialize(coefficients, parameters):
    matrices = [matrix(a) for a in coefficients]
    values = [q(1)]+[q(v) for v in parameters]
    if not matrices or len(values) != len(matrices):
        raise ValueError("Coefficient count must be one plus parameter count")
    shape = len(matrices[0]), len(matrices[0][0])
    if any((len(a), len(a[0])) != shape for a in matrices):
        raise ValueError("Coefficient shapes differ")
    return tuple(tuple(sum(value*a[i][j] for value, a in zip(values, matrices))
                       for j in range(shape[1])) for i in range(shape[0]))
