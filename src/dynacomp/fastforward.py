"""Exact logarithmic composition of reduced affine update blocks."""
from .rational import matrix, identity, mul, mv, q, addv


def compose(after, before):
    """(A,b) after (C,d) is (AC,Ad+b)."""
    a, b = after
    c, d = before
    return mul(a, c), addv(mv(a, d), b)


def repeat_affine(A, b, repeats):
    a = matrix(A)
    if type(repeats) is not int or repeats < 0 or len(a) != len(a[0]) or len(b) != len(a):
        raise ValueError("Need square A, matching b and a nonnegative integer repeat count")
    result = identity(len(a)), tuple(q(0) for _ in a)
    power = a, tuple(q(v) for v in b)
    while repeats:
        if repeats & 1:
            result = compose(power, result)
        repeats >>= 1
        if repeats:
            power = compose(power, power)
    return result


def repeat_input_block(reduction, x0, inputs, repeats):
    """Repeat a fixed control word without replaying every individual step."""
    r = len(reduction.G)
    block = identity(r), tuple(q(0) for _ in range(r))
    for u in inputs:
        step = reduction.G, mv(reduction.H, u)
        block = compose(step, block)
    power, offset = repeat_affine(*block, repeats)
    final = addv(mv(power, mv(reduction.E, x0)), offset)
    return mv(reduction.D, final)
