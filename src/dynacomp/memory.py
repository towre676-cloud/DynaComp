"""Derive exact scalar output recurrences with finite input memory."""
from .rational import matrix, mul, coordinates, mv, q, encoded, row_basis


def compile_output_memory(model):
    """Minimal scalar output recurrence for all initial states and input words."""
    if len(model.C) != 1 or not any(model.C[0]):
        raise ValueError("Need a nonzero scalar target")
    rows = [model.C[0]]
    while True:
        next_row = mul((rows[-1],), model.A)[0]
        coefficients = coordinates(next_row, matrix(rows))
        if coefficients is not None:
            break
        rows.append(next_row)
    order = len(rows)
    markov = [mul((row,), model.B)[0] for row in rows]
    kernels = []
    for j in range(order):
        kernel = list(markov[order-1-j])
        for i in range(j+1, order):
            for action in range(len(kernel)):
                kernel[action] -= coefficients[i]*markov[i-1-j][action]
        kernels.append(tuple(kernel))
    return {"schema": "dynacomp.memory/1", "order": order, "A": encoded(model.A),
            "B": encoded(model.B), "C": encoded(model.C),
            "output_coefficients": [str(v) for v in coefficients],
            "input_kernels": encoded(kernels),
            "meaning": "y[t+q]=sum_i a[i]y[t+i]+sum_j K[j]u[t+j]"}


def replay_memory(packet):
    """Check output and forcing polynomial identities, not training traces."""
    from .linear import LinearModel
    if packet.get("schema") != "dynacomp.memory/1":
        raise ValueError("Unknown schema")
    model = LinearModel.from_rows(packet["A"], packet["B"], packet["C"])
    order = packet["order"]
    if type(order) is not int or not 1 <= order <= len(model.A) or len(model.C) != 1:
        raise ValueError("Invalid scalar recurrence order")
    a = tuple(q(v) for v in packet["output_coefficients"])
    kernels = matrix(packet["input_kernels"])
    if len(a) != order or len(kernels) != order or len(kernels[0]) != len(model.B[0]):
        raise ValueError("Invalid recurrence dimensions")
    rows = [model.C[0]]
    for _ in range(order):
        rows.append(mul((rows[-1],), model.A)[0])
    if len(row_basis(rows[:-1])) != order:
        raise ValueError("Recurrence is not of minimal observable order")
    expected = tuple(sum(a[i]*rows[i][j] for i in range(order)) for j in range(len(model.A)))
    if expected != rows[-1]:
        raise ValueError("Output recurrence identity fails")
    markov = [mul((row,), model.B)[0] for row in rows[:-1]]
    for j in range(order):
        expected = tuple(markov[order-1-j][k]-sum(a[i]*markov[i-1-j][k] for i in range(j+1, order))
                         for k in range(len(model.B[0])))
        if expected != kernels[j]:
            raise ValueError("Input-memory identity fails")
    return True


def simulate_memory(packet, initial_outputs, inputs):
    """Initial q outputs must be consistent with the chosen initial input prefix."""
    replay_memory(packet)
    order = packet["order"]
    ys = [q(v) for v in initial_outputs]
    us = [tuple(q(v) for v in u) for u in inputs]
    a = tuple(q(v) for v in packet["output_coefficients"])
    kernels = matrix(packet["input_kernels"])
    if len(ys) != order or len(us) < order-1 or any(len(u) != len(kernels[0]) for u in us):
        raise ValueError("Need q initial outputs and at least q-1 matching inputs")
    for t in range(len(us)-order+1):
        output = sum(a[i]*ys[t+i] for i in range(order))
        output += sum(sum(kernels[j][k]*us[t+j][k] for k in range(len(kernels[0]))) for j in range(order))
        ys.append(output)
    return tuple(ys)
