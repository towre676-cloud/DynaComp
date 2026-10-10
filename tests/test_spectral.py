"""Independent formal-series, Decimal, tail, domain and composition checks."""
import copy
import json
import random
import unittest
from decimal import Decimal, localcontext
from fractions import Fraction as F

from dynacomp.enclosures import exp_negative, log_rational, pi_interval, sqrt_integer
from dynacomp.spectral import (compile_spectral, replay_spectral, search_spectral,
    ramanujan_coefficients, spectral_error_at, compose_spectral,
    replay_spectral_composition, compress_spectral_modes, replay_spectral_subset)
from dynacomp.spectral import simulate_spectral
from dynacomp.ramanujan import ramanujan_value, replay_ramanujan_value
from dynacomp.rational import q

PI = Decimal('3.14159265358979323846264338327950288419716939937510582097494459230781640628620899862803482534211706798214808651328230664709384460955058223172535940812848111745028410270193852110555964462294895493038196')


def dec(x):
    x = q(x)
    return Decimal(x.numerator)/Decimal(x.denominator)


def add(a, b):
    return a[0]+b[0], a[1]+b[1]


def product(a, b):
    return a[0]*b[0]+5*a[1]*b[1], a[0]*b[1]+a[1]*b[0]


def cf_value(t, depth=1000):
    """Direct finite continued fraction, independent of product/divisor formulas."""
    with localcontext() as ctx:
        ctx.prec = 100
        t = dec(t)
        z = (-t).exp()
        denominator = Decimal(1)
        for k in range(depth, 0, -1):
            denominator = 1+z**k/denominator
        return (-t/5).exp()/denominator


class EnclosureTests(unittest.TestCase):
    def test_sqrt_endpoints_are_exact_integer_inequalities(self):
        for bits in [16, 64, 128]:
            for n in [0, 1, 2, 5, 121, 1000003]:
                lo, hi = sqrt_integer(n, bits)
                self.assertLessEqual(lo*lo, n)
                self.assertGreaterEqual(hi*hi, n)
                self.assertLessEqual(hi-lo, F(1, 1 << bits))

    def test_log_against_high_precision_decimal(self):
        with localcontext() as ctx:
            ctx.prec = 180
            for x in [F(1, 100), F(1), F(3, 2), F(2), F(4096), F(1000001, 7)]:
                for bits in [16, 80, 128]:
                    lo, hi = log_rational(x, bits)
                    value = dec(x).ln()
                    self.assertLessEqual(dec(lo), value)
                    self.assertGreaterEqual(dec(hi), value)
                    self.assertLessEqual(hi-lo, F(2, 1 << bits))

    def test_machin_pi_contains_independent_digits(self):
        with localcontext() as ctx:
            ctx.prec = 180
            for bits in [16, 80, 128]:
                lo, hi = pi_interval(bits)
                self.assertLess(dec(lo), PI)
                self.assertGreater(dec(hi), PI)

    def test_exp_with_scaling_and_extreme_decay(self):
        with localcontext() as ctx:
            ctx.prec = 160
            for x in [0, F(1, 1000000), F(1, 2), F(19, 7), 100, 10000]:
                for bits in [16, 80, 128]:
                    lo, hi = exp_negative(x, bits)
                    exact = (-dec(x)).exp()
                    self.assertLessEqual(dec(lo), exact)
                    self.assertGreaterEqual(dec(hi), exact)
                    self.assertLess(hi-lo, F(10, 1 << bits))

    def test_invalid_constant_domains(self):
        for func, arg in [(log_rational, 0), (log_rational, -1), (exp_negative, -1), (sqrt_integer, -1)]:
            with self.assertRaises(ValueError):
                func(arg)
        for bits in [True, 8, 513]:
            with self.assertRaises(ValueError):
                pi_interval(bits)


class SpectralTests(unittest.TestCase):
    def test_divisor_coefficients_equal_independent_formal_product_log(self):
        # F(q) = product (1-q^j)^chi(j), expanded by polynomial multiplication.
        count = 50
        degree = count//5
        f = [1]+[0]*degree
        char = {0: 0, 1: 1, 2: -1, 3: -1, 4: 1}
        for j in range(1, degree+1):
            sign = char[j % 5]
            if sign == 1:
                f = [f[n]-(f[n-j] if n >= j else 0) for n in range(degree+1)]
            elif sign == -1:
                old = f[:]
                f = [sum(old[n-a*j] for a in range(n//j+1)) for n in range(degree+1)]
        w = [(F(0), F(0))]*(count+1)
        w[0] = (F(1), F(0))
        for m, value in enumerate(f):
            if 5*m+1 <= count:
                w[5*m+1] = (-F(value, 2), F(value, 2))
        coeff = [(F(0), F(0))]*(count+1)
        for n in range(1, count+1):
            v = (n*w[n][0], n*w[n][1])
            for j in range(1, n):
                term = product((j*coeff[j][0], j*coeff[j][1]), w[n-j])
                v = (v[0]-term[0], v[1]-term[1])
            coeff[n] = v[0]/n, v[1]/n
        self.assertEqual(tuple(coeff[1:]), ramanujan_coefficients(count))

    def test_positive_quadratic_part_and_log_growth(self):
        with localcontext() as ctx:
            ctx.prec = 100
            for n, (a, b) in enumerate(ramanujan_coefficients(200), 1):
                self.assertGreater(b, 0)
                val = dec(a)+dec(b)*Decimal(5).sqrt()
                bound = (Decimal(5).sqrt()+3)/2*(1+Decimal(n).ln())
                self.assertLessEqual(abs(val), bound)

    def test_power_family_tail_against_large_independent_sum(self):
        packet = compile_spectral({'kind': 'power_moments', 'power': 2, 'amplitude': '-3/2'}, 8, '1/100')
        tail = sum(F(3, 2*k**4) for k in range(9, 1000))
        self.assertLess(tail, q(packet['all_time_upper']))
        self.assertEqual(packet['status'], 'tolerance_met')
        self.assertTrue(replay_spectral(json.loads(json.dumps(packet))))

    def test_ramanujan_tail_against_extra_exact_coefficients(self):
        packet = compile_spectral({'kind': 'ramanujan_log', 'start': 2}, 16, '1/100')
        with localcontext() as ctx:
            ctx.prec = 100
            root = Decimal(5).sqrt()
            pairs = ramanujan_coefficients(300)
            for step in [0, 1, 4]:
                remainder = sum(abs(dec(a)+dec(b)*root)/Decimal(k)**(4+2*step)
                                for k, (a, b) in enumerate(pairs, 1) if k > 16)
                self.assertLess(remainder, dec(spectral_error_at(packet, step)))

    def test_search_exhaustion_is_never_acceptance(self):
        packet, history = search_spectral({'kind': 'ramanujan_log'}, '1/1000000', max_modes=7)
        self.assertEqual([row['modes'] for row in history], [1, 2, 4, 7])
        self.assertEqual(packet['status'], 'not_certified')
        good, _ = search_spectral({'kind': 'ramanujan_log', 'start': 2}, '1/100', max_modes=32)
        self.assertEqual(good['status'], 'tolerance_met')

    def test_heat_initial_field_and_steady_state_against_analytic_values(self):
        with localcontext() as ctx:
            ctx.prec = 100
            for modes in [1, 4, 16]:
                packet = compile_spectral({'kind': 'heat_rod'}, modes, '1/100')
                initial = sum(dec(v) for v in packet['model']['C'][0])
                steady = sum(dec(v)/Decimal((2*i+1)**2) for i, v in enumerate(packet['model']['C'][0]))
                self.assertLess(abs(1-initial), dec(packet['all_time_upper']))
                self.assertLess(abs(PI**2/12-steady), dec(packet['all_time_upper']))

    def test_heat_all_time_budget_with_independent_decimal_zoh_simulation(self):
        packet = compile_spectral({'kind': 'heat_rod', 'dt': '1/10'}, 8, '1/10', bits=64)
        rng = random.Random(712)
        with localcontext() as ctx:
            ctx.prec = 80
            modes = 300
            poles = [(-Decimal((2*i+1)**2)/10).exp() for i in range(modes)]
            weights = [8/(PI**2*Decimal((2*i+1)**2)) for i in range(modes)]
            exact = [Decimal('-0.73')]*modes
            retained = [Decimal('-0.73')]*8
            a = [dec(row[2]) for row in packet['model']['entries']]
            b = [dec(row[0]) for row in packet['model']['B']]
            c = [dec(v) for v in packet['model']['C'][0]]
            for t in range(41):
                y = sum(w*x for w, x in zip(weights, exact))
                yr = sum(w*x for w, x in zip(c, retained))
                # Charge the independent 300-mode reference's own omitted tail.
                reference_tail = Decimal(8)/(PI**2)*(Decimal(1)/(4*modes)+Decimal(1)/(48*modes**3))
                self.assertLess(abs(y-yr)+reference_tail, dec(packet['all_time_upper']))
                u = Decimal(rng.randrange(-100, 101))/100
                exact = [p*x+(1-p)*u/Decimal((2*i+1)**2) for i, (p, x) in enumerate(zip(poles, exact))]
                retained = [p*x+gain*u for p, x, gain in zip(a, retained, b)]

    def test_heat_precision_error_is_accounted_and_shrinks(self):
        coarse = compile_spectral({'kind': 'heat_rod', 'dt': '1/100'}, 8, 1, bits=16)
        fine = compile_spectral({'kind': 'heat_rod', 'dt': '1/100'}, 8, 1, bits=80)
        key = 'retained_coefficient_and_sampling'
        self.assertGreater(q(coarse['error_budget'][key]), q(fine['error_budget'][key]))
        self.assertGreater(q(coarse['all_time_upper']), q(fine['all_time_upper']))
        with self.assertRaises(ValueError):
            compile_spectral({'kind': 'heat_rod', 'dt': '1/1000000000'}, 1, 1, bits=16)

    def test_packet_tampering_rejected(self):
        p = compile_spectral({'kind': 'ramanujan_log'}, 4, 1)
        for change in [lambda x: x.update(all_time_upper='0'),
                       lambda x: x['model']['C'][0].__setitem__(0, '0'),
                       lambda x: x['domain']['initial_radius'].__setitem__(0, '1'),
                       lambda x: x['modal_enclosures'][0]['weight_interval'].__setitem__(0, '0'),
                       lambda x: x.update(status='refuted')]:
            forged = copy.deepcopy(p)
            change(forged)
            with self.assertRaises(ValueError):
                replay_spectral(forged)

    def test_invalid_source_parameters(self):
        for family in [{'kind': 'unknown'}, {'kind': 'heat_rod', 'dt': 0},
                       {'kind': 'heat_rod', 'input_bound': -1}, {'kind': 'ramanujan_log', 'start': 0},
                       {'kind': 'power_moments', 'power': -1}, {'kind': 'heat_rod', 'mass_matrix': 1}]:
            with self.assertRaises(ValueError):
                compile_spectral(family, 4, 1)
        for count in [0, True, -1]:
            with self.assertRaises(ValueError):
                compile_spectral({'kind': 'power_moments'}, count, 1)

    def test_zero_power_amplitude_exact(self):
        p = compile_spectral({'kind': 'power_moments', 'amplitude': 0}, 2, 0)
        self.assertEqual(p['all_time_upper'], '0')
        self.assertEqual(p['status'], 'tolerance_met')

    def test_subset_all_time_and_arbitrary_order(self):
        p = compile_spectral({'kind': 'heat_rod', 'initial_amplitude': 0}, 32, '1/50')
        out = compress_spectral_modes(p, [0], '1/50')
        self.assertEqual(out['status'], 'tolerance_met')
        self.assertLess(float(q(out['combined_all_time_upper'])), 0.012)
        self.assertTrue(replay_spectral_subset(out))
        r = compress_spectral_modes(p, [2, 0, 1], '1/50')
        self.assertEqual(r['G'][0][0], p['model']['entries'][2][2])
        for keep in [[], [0, 0], [32], [True]]:
            with self.assertRaises(ValueError):
                compress_spectral_modes(p, keep, 1)

    def test_subset_ramanujan_exact_evolution_covered_forever(self):
        p = compile_spectral({'kind': 'ramanujan_log', 'start': 2}, 32, '1/50')
        out = compress_spectral_modes(p, [0], '1/50')
        self.assertEqual(out['status'], 'tolerance_met')
        for t in [0, 1, 10, 100]:
            dropped = sum(abs(q(p['model']['C'][0][i]))*q(p['domain']['initial_center'][i])*
                          q(p['model']['entries'][i][2])**t for i in range(1, 32))
            self.assertLessEqual(dropped, q(out['removed_modes_upper']))
        forged = copy.deepcopy(out)
        forged['G'][0][0] = '2'
        with self.assertRaises(ValueError):
            replay_spectral_subset(forged)

    def test_arbitrary_reduction_composes_and_refuses_exhausted_budget(self):
        p = compile_spectral({'kind': 'power_moments', 'power': 2}, 4, '1/100')
        args = ([[1, 0, 0, 0]], [[1]], [[0]], [[1]])
        out = compose_spectral(p, *args, 4, '1/10')
        self.assertEqual(out['status'], 'tolerance_met')
        self.assertTrue(replay_spectral_composition(json.loads(json.dumps(out))))
        fail = compose_spectral(p, *args, 4, 0)
        self.assertEqual(fail['status'], 'not_certified')
        self.assertTrue(replay_spectral_composition(fail))

    def test_composition_cross_source_domain_tampering_rejected(self):
        p = compile_spectral({'kind': 'power_moments', 'power': 2}, 4, 1)
        args = ([[1, 0, 0, 0]], [[1]], [[0]], [[1]])
        out = compose_spectral(p, *args, 4, 1)
        for change in [lambda x: x.update(combined_upper='0'), lambda x: x.update(horizon=5),
                       lambda x: x['finite']['initial_center'].__setitem__(0, '0')]:
            forged = copy.deepcopy(out)
            change(forged)
            with self.assertRaises(ValueError):
                replay_spectral_composition(forged)

    def test_exact_execution_and_source_input_guards(self):
        p = compile_spectral({'kind': 'power_moments', 'power': 2}, 4, 1)
        values = simulate_spectral(p, [0]*4)
        for t, row in enumerate(values):
            self.assertEqual(row[0], sum(F(1, k**(4+2*t)) for k in range(1, 5)))
        for options in [{'inputs': [1]}, {'inputs': [0], 'initial_amplitude': 1}]:
            with self.assertRaises(ValueError):
                simulate_spectral(p, **options)

    def test_heat_execution_common_field_and_outside_domain(self):
        p = compile_spectral({'kind': 'heat_rod'}, 16, 1)
        outputs = simulate_spectral(p, [1, -1], initial_amplitude='1/2')
        self.assertEqual(outputs[0][0], sum(q(v) for v in p['model']['C'][0])/2)
        for inputs, initial in [([2], 0), ([0], 2)]:
            with self.assertRaises(ValueError):
                simulate_spectral(p, inputs, initial)

    def test_execution_refuses_unaccepted_and_over_horizon(self):
        p = compile_spectral({'kind': 'power_moments'}, 2, 0)
        with self.assertRaises(ValueError):
            simulate_spectral(p, [0])
        self.assertEqual(len(simulate_spectral(p, [0], allow_unaccepted=True)), 2)
        packet = compose_spectral(p, [[1,0]], [[1]], [[0]], [[1]], 2, 2)
        self.assertEqual(len(simulate_spectral(packet, [0,0])), 3)
        with self.assertRaises(ValueError):
            simulate_spectral(packet, [0,0,0])

    def test_subset_simulation_matches_direct_selected_modes(self):
        p = compile_spectral({'kind': 'heat_rod'}, 4, 1)
        selected = compress_spectral_modes(p, [2,0], 1)
        values = simulate_spectral(selected, [1, 0], initial_amplitude='-1/2')
        state = [F(-1,2)]*2
        for t, u in enumerate([1,0]):
            self.assertEqual(values[t][0], sum(q(p['model']['C'][0][i])*v for i,v in zip([2,0],state)))
            state = [q(p['model']['entries'][i][2])*v+q(p['model']['B'][i][0])*u for i,v in zip([2,0],state)]

    def test_general_direction_projection_exact_simulation_and_combined_bound(self):
        try:
            import numpy
            import scipy
        except ImportError:
            self.skipTest('Numerical proposal extras not installed')
        from dynacomp import SparseLinearModel, StableProposalContext
        source = compile_spectral({'kind':'heat_rod'}, 16, '1/50')
        proposal = StableProposalContext(SparseLinearModel.from_spec(source['model']), 20).propose(3)
        packet = compose_spectral(source, *(proposal[k] for k in ('E','G','H','D')), 20, '1/50', require_contraction=True)
        self.assertEqual(packet['status'], 'tolerance_met')
        self.assertTrue(packet['finite']['reduced_contraction_certified'])
        inputs = [(-1)**j for j in range(20)]
        full = simulate_spectral(source, inputs, '3/4')
        reduced = simulate_spectral(packet, inputs, '3/4')
        finite_error = q(packet['finite']['certificate']['worst_upper'])
        for a, b in zip(full, reduced):
            self.assertLessEqual(abs(a[0]-b[0]), finite_error)
        self.assertGreaterEqual(q(packet['combined_upper']), finite_error+q(source['all_time_upper']))


class TwoCuspTests(unittest.TestCase):
    def test_direct_and_reciprocal_against_continued_fraction(self):
        with localcontext() as ctx:
            ctx.prec = 100
            for t in ['1/10', '1', '59/10', '6', '20']:
                p = ramanujan_value(t, bits=80)
                value = cf_value(t)
                self.assertLess(dec(p['interval'][0]), value)
                self.assertGreater(dec(p['interval'][1]), value)
                self.assertTrue(replay_ramanujan_value(p))

    def test_extreme_cusp_uses_bounded_work(self):
        p = ramanujan_value('1/1000000', bits=80)
        self.assertEqual(p['branch'], 'reciprocal')
        self.assertEqual(p['product_terms'], 1)
        with localcontext() as ctx:
            ctx.prec = 100
            limit = (Decimal(5).sqrt()-1)/2
            self.assertLess(dec(p['interval'][0]), limit)
            self.assertGreater(dec(p['interval'][1]), limit)

    def test_high_precision_interval_stays_serializable(self):
        p = ramanujan_value(6, bits=512)
        self.assertLess(len(json.dumps(p)), 10000)
        self.assertTrue(replay_ramanujan_value(p))

    def test_tampered_value_and_invalid_t_rejected(self):
        p = ramanujan_value(1)
        p['interval'][0] = '0'
        with self.assertRaises(ValueError):
            replay_ramanujan_value(p)
        for t in [0, -1]:
            with self.assertRaises(ValueError):
                ramanujan_value(t)
