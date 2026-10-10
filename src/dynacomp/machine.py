"""Uniform error bounds for a specified integer-accumulated binary64 runtime.

Every row sums exact products of stored operands, then rounds once to nearest
even. Encoder, transition+forcing, and output each follow that contract. This
does not certify NumPy/BLAS or arbitrary sequential floating point execution.
"""
from fractions import Fraction as Q
from .rational import matrix, q, encoded
from .binary64 import bits, from_bits, exact, round_bits, dot_bits, classify, MAX_FINITE
from .validated import replay_validated

U = Q(1, 1 << 53)
ETA = Q(1, 1 << 1075)
LIMIT = exact(MAX_FINITE)


def _up(x):
    scale = 1 << 96
    return Q((x*scale).__ceil__(), scale)


def _vec(values, n):
    result = tuple(map(q, values))
    if len(result) != n:
        raise ValueError('Domain dimensions differ')
    return result


def _abs_mv(a, x):
    return tuple(sum(abs(v)*w for v, w in zip(row, x)) for row in a)


def _round_budget(magnitude):
    if magnitude > LIMIT:
        raise ValueError('Cannot exclude binary64 overflow on this domain')
    return _up(U*magnitude+ETA) if magnitude else Q(0)


def _certify_execution(E, G, H, D, initial_center, initial_radius, input_radii,
                       tolerance, input_centers, budget, contract, schema):
    e, g, h, d = map(matrix, (E, G, H, D))
    r, n, m = len(g), len(e[0]), len(h[0])
    if len(e) != r or len(g[0]) != r or len(h) != r or len(d[0]) != r:
        raise ValueError('Expected E[r,n], G[r,r], H[r,m], D[p,r]')
    center, radius = _vec(initial_center, n), _vec(initial_radius, n)
    radii = tuple(_vec(row, m) for row in input_radii)
    means = tuple(_vec(row, m) for row in input_centers) if input_centers is not None else ((Q(0),)*m,)*len(radii)
    if len(means) != len(radii) or len(radii) > 4096:
        raise ValueError('Matching input domain and horizon <= 4096 required')
    if any(x < 0 for x in radius+sum(radii, ())):
        raise ValueError('Nonnegative radii required')
    tol = q(tolerance)
    if tol < 0:
        raise ValueError('Nonnegative tolerance required')
    words, loaded, differences = {}, {}, {}
    for name, a in zip(('E', 'G', 'H', 'D'), (e, g, h, d)):
        w = tuple(tuple(round_bits(v) for v in row) for row in a)
        if any(classify(v) == 'infinity' for row in w for v in row):
            raise ValueError('Matrix loading overflows binary64')
        loaded[name] = tuple(tuple(exact(v) for v in row) for row in w)
        differences[name] = tuple(tuple(abs(x-y) for x, y in zip(row, other))
                                  for row, other in zip(a, loaded[name]))
        words[name] = [list(row) for row in w]
    eh, gh, hh, dh = (loaded[k] for k in ('E', 'G', 'H', 'D'))
    de, dg, dh_input, dd = (differences[k] for k in ('E', 'G', 'H', 'D'))
    xbound = tuple(abs(c)+b for c, b in zip(center, radius))
    if any(v > LIMIT for v in xbound):
        raise ValueError('Initial domain exceeds finite binary64 range')
    z = tuple(_up(v) for v in _abs_mv(e, xbound))
    error = tuple(_up(a+budget(b, n)) for a, b in zip(_abs_mv(de, xbound), _abs_mv(eh, xbound)))
    outputs, states = [], []
    for t in range(len(radii)+1):
        actual_bound = tuple(a+b for a, b in zip(z, error))
        output_error = tuple(_up(a+b+budget(c, r)) for a, b, c in zip(
            _abs_mv(dh, error), _abs_mv(dd, z), _abs_mv(dh, actual_bound)))
        states.append(list(map(str, error)))
        outputs.append(list(map(str, output_error)))
        if t == len(radii):
            break
        ubound = tuple(abs(c)+b for c, b in zip(means[t], radii[t]))
        if any(v > LIMIT for v in ubound):
            raise ValueError('Input domain exceeds finite binary64 range')
        next_error = tuple(_up(a+b+c+budget(v+w, r+m)) for a, b, c, v, w in zip(
            _abs_mv(gh, error), _abs_mv(dg, z), _abs_mv(dh_input, ubound),
            _abs_mv(gh, actual_bound), _abs_mv(hh, ubound)))
        z = tuple(_up(v+w) for v, w in zip(_abs_mv(g, z), _abs_mv(h, ubound)))
        error = next_error
    worst = max(q(v) for row in outputs for v in row)
    return {'schema': schema,
            **{k: encoded(a) for k, a in zip(('E', 'G', 'H', 'D'), (e, g, h, d))},
            'initial_center': list(map(str, center)), 'initial_radius': list(map(str, radius)),
            'input_centers': [list(map(str, row)) for row in means],
            'input_radii': [list(map(str, row)) for row in radii], 'loaded_words': words,
            'state_error_upper': states, 'output_error_upper': outputs,
            'worst_execution_upper': str(worst), 'tolerance': str(tol),
            'status': 'tolerance_met' if worst <= tol else 'not_certified',
            'contract': contract,
            'scope': 'uniform deviation from the supplied rational reduced dynamics; finite horizon; includes coefficient loading'}


def certify_execution(E, G, H, D, initial_center, initial_radius, input_radii,
                      tolerance, input_centers=None):
    return _certify_execution(E, G, H, D, initial_center, initial_radius,
                             input_radii, tolerance, input_centers,
                             lambda magnitude, count: _round_budget(magnitude),
                             'finite stored binary64 inputs; integer exact dot; one RN-even rounding per row; +0 for exact zero',
                             'dynacomp.binary64-execution/1')


def replay_execution(packet):
    if packet.get('schema') != 'dynacomp.binary64-execution/1':
        raise ValueError('Unknown execution schema')
    fresh = certify_execution(*(packet[k] for k in ('E', 'G', 'H', 'D')),
                              packet['initial_center'], packet['initial_radius'],
                              packet['input_radii'], packet['tolerance'], packet['input_centers'])
    if packet != fresh:
        raise ValueError('Execution certificate does not reconstruct')
    return True


def compose_execution(source, tolerance):
    """Combine replayed finite model/source guarantees with machine execution."""
    if source.get('schema') == 'dynacomp.validated-box/1':
        replay_validated(source)
        finite, source_error = source, Q(0)
    elif source.get('schema') == 'dynacomp.spectral-composition/1':
        from .spectral import replay_spectral_composition
        replay_spectral_composition(source)
        finite = source['finite']
        source_error = q(source['combined_upper'])-q(finite['certificate']['worst_upper'])
    else:
        raise ValueError('Need a validated-box or spectral-composition source')
    tol = q(tolerance)
    runtime = certify_execution(*(finite[k] for k in ('E', 'G', 'H', 'D')),
                                finite['initial_center'], finite['initial_radius'],
                                finite['input_radii'], tol, finite['input_centers'])
    combined = [[str(source_error+q(a)+q(b)) for a, b in zip(model, execution)]
                for model, execution in zip(finite['certificate']['output_upper'], runtime['output_error_upper'])]
    worst = max(q(x) for row in combined for x in row)
    stable = finite['reduced_contraction_certified'] is not False
    return {'schema': 'dynacomp.machine-composition/1', 'source': source, 'execution': runtime,
            'tolerance': str(tol), 'combined_output_upper': combined, 'combined_upper': str(worst),
            'status': 'tolerance_met' if worst <= tol and stable else 'not_certified',
            'scope': 'source/model error plus actual specified binary64 execution; finite horizon'}


def replay_machine_composition(packet):
    if packet.get('schema') != 'dynacomp.machine-composition/1' or packet != compose_execution(packet['source'], packet['tolerance']):
        raise ValueError('Machine composition does not reconstruct')
    return True


class Binary64Runtime:
    """Execute precisely the machine protocol certified above; no BLAS calls."""
    def __init__(self, packet, require_accepted=True):
        import copy
        self._packet = copy.deepcopy(packet)
        if packet.get('schema') == 'dynacomp.machine-composition/1':
            replay_machine_composition(self._packet)
            self._execution = self._packet['execution']
        else:
            replay_execution(self._packet)
            self._execution = self._packet
        if require_accepted and self._packet['status'] != 'tolerance_met':
            raise ValueError('Tolerance not certified')
        self._words = {name: tuple(tuple(row) for row in rows) for name, rows in self._execution['loaded_words'].items()}

    @property
    def packet(self):
        import copy
        return copy.deepcopy(self._packet)

    def _domain(self, words, center, radius):
        if len(words) != len(center):
            raise ValueError('Input dimensions differ')
        for word, c, r in zip(words, center, radius):
            if abs(exact(word)-q(c)) > q(r):
                raise ValueError('Stored value is outside the certified domain')

    def simulate_bits(self, initial, inputs):
        p = self._execution
        initial, inputs = tuple(initial), tuple(tuple(row) for row in inputs)
        if len(inputs) > len(p['input_radii']):
            raise ValueError('Execution exceeds the certified horizon')
        self._domain(initial, p['initial_center'], p['initial_radius'])
        for t, row in enumerate(inputs):
            self._domain(row, p['input_centers'][t], p['input_radii'][t])
        z = tuple(dot_bits(row, initial) for row in self._words['E'])
        output = []
        for t in range(len(inputs)+1):
            output.append(tuple(dot_bits(row, z) for row in self._words['D']))
            if t < len(inputs):
                z = tuple(dot_bits(tuple(g)+tuple(h), z+inputs[t])
                          for g, h in zip(self._words['G'], self._words['H']))
        return tuple(output)

    def simulate(self, initial, inputs):
        result = self.simulate_bits(tuple(map(bits, initial)), [tuple(map(bits, row)) for row in inputs])
        return tuple(tuple(from_bits(word) for word in row) for row in result)
