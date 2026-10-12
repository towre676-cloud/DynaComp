"""Plain-language scenario adapter for the web/mobile explorer.

Real DynaComp rational compilation and replay drive the comparison. Illustrative
room/supply models are explicitly labelled; they are not fitted customer models.
"""
from fractions import Fraction
from .linear import LinearModel, compile_linear, replay_linear
from .adaptive import adaptive_coordinates, replay_adaptive
from .rational import q, mv, addv, matrix, encoded


def _number(value, lo, hi):
    if isinstance(value, bool):
        raise ValueError('A number is required')
    value = Fraction(str(value))
    if not lo <= value <= hi:
        raise ValueError('The setting is outside its supported range')
    return value


def scenario(request):
    kind = request.get('example', 'warmth')
    if kind not in ('warmth', 'supplies', 'custom'):
        raise ValueError('Choose warmth, supplies, or your own model')
    step_value = _number(request.get('steps', 16), 1, 40)
    if step_value.denominator != 1:
        raise ValueError('Choose a whole number of steps')
    steps = int(step_value)
    start = _number(request.get('start', 20), -100, 1000)
    strength = _number(request.get('strength', 2), -20, 20)
    tolerance = _number(request.get('tolerance', 0), 0, 100)
    output = request.get('output', 'whole')
    if output not in ('whole', 'one'):
        raise ValueError('Choose the whole system or one place')
    feedthrough = None
    descriptor = None
    if kind == 'custom':
        workload = request['model']
        supplied = workload.get('model', workload)
        if any(key in supplied for key in ('time', 'sample_period')):
            raise ValueError('Open a prepared sampled model here. Continuous-time preparation uses the full downloadable DynaComp package.')
        n = len(supplied['A'])
        if not 1 <= n <= 12:
            raise ValueError('The interactive workspace supports 1–12 model values')
        a, b = supplied['A'], supplied['B']
        if 'E' in supplied:
            from .engineering import descriptor_transform
            a, b, descriptor = descriptor_transform(matrix(a), matrix(b), matrix(supplied['E']))
            if descriptor['status'] != 'exact':
                raise ValueError('This interactive workspace requires exact descriptor elimination')
        model = LinearModel.from_rows(a, b, supplied['C'])
        f = workload.get('feedthrough', supplied.get('D'))
        if f is not None:
            feedthrough = matrix(f)
            if len(feedthrough) != len(model.C) or len(feedthrough[0]) != len(model.B[0]):
                raise ValueError('The direct-output inputs do not match the model')
        if len(model.B[0]) > 4 or len(model.C) > 4:
            raise ValueError('Use at most four inputs and outputs in this workspace')
        initial = tuple(q(x) for x in supplied.get('initial', workload.get('initial_center', [str(start)] * n)))
        inputs = [tuple(q(x) for x in u) for u in supplied.get('inputs', [[str(strength)] * len(model.B[0])] * steps)]
        if any(len(u) != len(model.B[0]) for u in inputs):
            raise ValueError('The input values do not match the model')
        if len(inputs) > 40:
            raise ValueError('Use at most forty updates in this workspace')
        title, unit = 'Your model', supplied.get('unit', 'output units')
    else:
        n = 6
        a = [[Fraction(0) for _ in range(n)] for _ in range(n)]
        for i in range(n):
            a[i][i] = Fraction(3, 4) if kind == 'warmth' else Fraction(7, 8)
            a[i][(i - 1) % n] = Fraction(1, 16)
            a[i][(i + 1) % n] = Fraction(1, 16)
        c = [[Fraction(1, n)] * n] if output == 'whole' else [[int(i == 0) for i in range(n)]]
        model = LinearModel.from_rows(a, [[Fraction(1, 4)] for _ in range(n)], c)
        initial = tuple(start + Fraction(i - 2, 4) for i in range(n))
        inputs = [(strength if t < steps // 2 else strength / 2,) for t in range(steps)]
        title, unit = ('Warmth across six connected rooms', 'illustrative temperature units') if kind == 'warmth' else ('Supplies across six connected places', 'illustrative stock units')
    if len(initial) != n:
        raise ValueError('The starting values do not match the model')
    exact = compile_linear(model)
    packet = exact.packet()
    replay_linear(packet)
    reduced = exact.simulate(initial, inputs)
    bound = [Fraction(0)] * (len(inputs) + 1)
    rank, method, accepted = len(exact.E), 'exact output-preserving reduction', True
    if tolerance > 0:
        # Bound a box around the supplied initial values and all inputs up to
        # each displayed magnitude. The nominal run is one member of that box.
        adaptive = adaptive_coordinates(model, initial, [0] * n,
                    [[str(abs(v)) for v in u] for u in inputs], str(tolerance),
                    initial_coordinates=[max(range(n), key=lambda i: sum(abs(row[i]) for row in model.C))])
        replay_adaptive(adaptive)
        if adaptive['status'] == 'tolerance_met' and len(adaptive['candidate']['selected']) < rank:
            packet = adaptive
            candidate = adaptive['candidate']
            z = mv(matrix(candidate['E']), initial)
            reduced = [mv(matrix(candidate['D']), z)]
            for u in inputs:
                z = addv(mv(matrix(candidate['G']), z), mv(matrix(candidate['H']), u))
                reduced.append(mv(matrix(candidate['D']), z))
            rank = len(candidate['selected'])
            bound = [q(x) for x in adaptive['box_bound']['bound']]
            method = 'bounded approximate reduction'
    x = initial
    full = [mv(model.C, x)]
    for u in inputs:
        x = addv(mv(model.A, x), mv(model.B, u))
        full.append(mv(model.C, x))
    if feedthrough is not None:
        zero = (q(0),) * len(model.B[0])
        # Same right-continuous terminal-hold convention as the engineering
        # importer: y_t=C x_t + F u_t, holding the final input at the last sample.
        direct = [mv(feedthrough, inputs[t] if t < len(inputs) else inputs[-1] if inputs else zero) for t in range(len(full))]
        full = [addv(y, d) for y, d in zip(full, direct)]
        reduced = [addv(y, d) for y, d in zip(reduced, direct)]
    differences = [max(abs(a - b) for a, b in zip(f, r)) for f, r in zip(full, reduced)]
    if any(d > b for d, b in zip(differences, bound)):
        raise AssertionError('The actual trajectory escaped its error bound')
    return {'schema': 'dynacomp.explorer-run/1', 'title': title, 'unit': unit,
            'original_values': n, 'retained_values': rank, 'method': method,
            'tolerance': str(tolerance), 'max_error_bound': str(max(bound)),
            'max_actual_difference': str(max(differences)), 'accepted': accepted,
            'scope': 'Exact rational computation for the supplied sampled model and displayed finite horizon. Illustrative scenarios are not customer measurements.',
            'request': request, 'full': [[str(v) for v in row] for row in full],
            'reduced': [[str(v) for v in row] for row in reduced],
            'error_bound': [str(v) for v in bound], 'certificate': packet,
            'descriptor': descriptor, 'feedthrough': encoded(feedthrough) if feedthrough is not None else None,
            'output_input_semantics': 'right-continuous-terminal-hold',
            'measurement_and_continuous_sampling_errors': 'not covered by the displayed model-reduction bound'}
