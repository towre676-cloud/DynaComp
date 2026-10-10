import copy
import math
import random
import unittest
from fractions import Fraction as Q
from dynacomp import binary64 as b
from dynacomp.generating import (RationalGF, from_recurrence, from_state,
                                periodic_gf, binomial_tail, replay_generating)


class BinaryTests(unittest.TestCase):
    def test_fields_and_roundtrip_all_boundaries(self):
        words = [0, b.SIGN, 1, b.FRAC, 1 << 52, b.MAX_FINITE,
                 b.INF, b.INF | b.SIGN, b.INF | 12345]
        self.assertEqual([b.classify(x) for x in words],
                         ['zero','zero','subnormal','subnormal','normal','normal','infinity','infinity','nan'])
        for w in words:
            self.assertEqual(b.bits(b.from_bits(w)), w)
        rng = random.Random(61)
        for _ in range(3000):
            word = rng.randrange(1 << 64)
            if b.classify(word) in ('nan','infinity'):
                continue
            value = b.exact(word)
            self.assertEqual(value, Q.from_float(b.from_bits(word)))
            self.assertEqual(b.round_bits(value), word if value else 0)

    def test_ties_subnormals_and_overflow(self):
        self.assertEqual(b.round_bits(Q(1)+Q(1, 1 << 53)), b.bits(1.0))
        self.assertEqual(b.round_bits(Q(1)+Q(3, 1 << 53)), b.bits(1.0)+2)
        self.assertEqual(b.round_bits(Q(1, 1 << 1075)), 0)
        self.assertEqual(b.round_bits(-Q(1, 1 << 1075)), b.SIGN)
        self.assertEqual(b.round_bits(Q(3, 1 << 1075)), 2)
        halfway = b.exact(b.MAX_FINITE)+Q(1 << 970)
        self.assertEqual(b.round_bits(halfway), b.INF)
        self.assertEqual(b.round_bits(halfway-1), b.MAX_FINITE)
        boundary = (b.exact(b.FRAC)+b.exact(1 << 52))/2
        self.assertEqual(b.round_bits(boundary), 1 << 52)

    def test_rational_rounding_against_host_reference(self):
        rng = random.Random(14)
        for _ in range(1500):
            x = Q(rng.randrange(-(1 << 70), 1 << 70), rng.randrange(1, 1 << 60))
            x *= Q(2)**rng.randrange(-1100, 1100)
            try:
                expected = b.bits(float(x))
            except OverflowError:
                expected = b.INF | (b.SIGN if x < 0 else 0)
            self.assertEqual(b.round_bits(x), expected)

    def test_integer_accumulator_cancellation_and_overflowing_products(self):
        xs = list(map(b.bits, [1e16, 1.0, -1e16]))
        ones = [b.bits(1.0)]*3
        self.assertEqual(b.from_bits(b.dot_bits(xs, ones)), 1.0)
        self.assertEqual((1e16+1.0)-1e16, 0.0)
        self.assertEqual(b.dot_bits([b.MAX_FINITE,b.MAX_FINITE|b.SIGN], [b.bits(2.0)]*2), 0)
        rng = random.Random(101)
        for _ in range(200):
            a = [b.bits(rng.uniform(-100,100)) for _ in range(17)]
            c = [b.bits(rng.uniform(-100,100)) for _ in range(17)]
            exact = sum((Q.from_float(b.from_bits(x))*Q.from_float(b.from_bits(y)) for x,y in zip(a,c)), Q(0))
            self.assertEqual(b.exact_dot(a,c), exact)
            self.assertEqual(b.dot_bits(a,c), b.bits(float(exact)))

    def test_neighbors_and_order(self):
        rng = random.Random(99)
        for w in [0,b.SIGN,1,b.FRAC,1 << 52,b.MAX_FINITE,b.INF,b.SIGN|b.INF]+[rng.randrange(b.INF) for _ in range(500)]:
            value = b.from_bits(w)
            self.assertEqual(b.next_up(w), b.bits(math.nextafter(value, math.inf)))
            self.assertEqual(b.next_down(w), b.bits(math.nextafter(value, -math.inf)))
        words = list(map(b.bits, [-math.inf,-100.,-1.,-0.,0.,1.,100.,math.inf]))
        self.assertEqual(sorted(words,key=b.order_key), words)

    def test_exact_power_and_decimal_trap(self):
        cases = [(4.,2,True),(0.25,2,True),(-0.125,3,True),(-4.,2,False),
                 (0.1,2,False),(2.,2,False),(0.,17,True)]
        for value,degree,expected in cases:
            p = b.power_certificate(b.bits(value),degree)
            self.assertEqual(p['is_power'],expected)
            self.assertTrue(b.replay_power(p))
            if expected:
                self.assertEqual(Q(p['root'])**degree,Q.from_float(value))
        self.assertTrue(b.power_certificate(1,2)['is_power'])
        self.assertFalse(b.power_certificate(1,5)['is_power'])
        p=b.power_certificate(b.bits(4.),2);p['root']='3'
        with self.assertRaises(ValueError):b.replay_power(p)

    def test_invalid_words_and_operands(self):
        for word in [-1,1 << 64,True,'1']:
            with self.assertRaises(ValueError):b.exact(word)
        for word in [b.INF,b.INF|1]:
            with self.assertRaises(ValueError):b.exact(word)
        with self.assertRaises(ValueError):b.dot_bits([0],[])
        with self.assertRaises(ValueError):b.order_key(b.INF|1)
        with self.assertRaises(TypeError):b.round_bits(0.1)


class GeneratingTests(unittest.TestCase):
    def test_cancellation_transients_and_minimality(self):
        for p,q,r in [([1,-1],[1,-1],1),([0,1],[1,-1],2),([0,0,1],[1],3),
                      ([0],[1,-1],0),([0,1],[1,-1,-1],2)]:
            gf=RationalGF(p,q);packet=gf.packet()
            self.assertEqual(gf.dimension,r)
            self.assertTrue(replay_generating(packet))
            if r:
                self.assertNotEqual(Q(packet['hankel_determinant']),0)
                self.assertEqual(from_state(*gf.realize()),gf)
        gf=RationalGF([0,0,1],[1]);self.assertEqual(gf.coefficients(6),(0,0,1,0,0,0))

    def test_random_recurrence_state_and_gf_agreement(self):
        rng=random.Random(572)
        for _ in range(80):
            order=rng.randrange(1,6)
            c=[Q(rng.randrange(-3,4),3) for _ in range(order)]
            initial=[Q(rng.randrange(-3,4),2) for _ in range(order)]
            gf=from_recurrence(c,initial);out=initial[:]
            for n in range(20-order):out.append(sum(x*y for x,y in zip(c,out[-order:])))
            self.assertEqual(gf.coefficients(20),tuple(out))
            self.assertTrue(replay_generating(gf.packet()))
            for n in [0,3,19]:self.assertEqual(gf.nth(n),out[n])
            if gf.dimension:self.assertEqual(from_state(*gf.realize()),gf)

    def test_arbitrary_state_extraction_cancels_invisible_modes(self):
        gf=from_state([['1/2',0,0],[0,'1/2',0],[0,0,'1/3']],[1,1,1],[1,1,0])
        self.assertEqual(gf,RationalGF([2],[1,'-1/2']))
        self.assertEqual(gf.dimension,1)
        self.assertTrue(replay_generating(gf.packet()))
        with self.assertRaises(ValueError):from_recurrence([1,1],[0,1]).nth(10**12)

    def test_arithmetic_filters_and_weighted_counts(self):
        gf=from_recurrence([1,1],[0,1])
        out=gf.coefficients(50)
        self.assertEqual(gf.subsequence(3,2).coefficients(10),out[2:32:3])
        self.assertEqual(gf.weighted_indices().coefficients(20),tuple(i*x for i,x in enumerate(out[:20])))
        periodic=periodic_gf([1,0],[1,0,1])
        self.assertEqual(periodic.nth(10**12),(1,0,1)[(10**12-2)%3])
        self.assertEqual((gf*gf).coefficients(20),tuple(sum(out[j]*out[i-j] for j in range(i+1)) for i in range(20)))

    def test_rational_cauchy_bound(self):
        gf=RationalGF([1,2],[1,'-1/3'])
        p=gf.tail('1/4','1/2',20)
        true_value=Q(3,2)/(1-Q(1,12))
        partial=sum(a*Q(1,4)**i for i,a in enumerate(gf.coefficients(20)))
        self.assertLessEqual(abs(true_value-partial),Q(p['absolute_tail_upper']))
        self.assertTrue(replay_generating(p))
        with self.assertRaises(ValueError):gf.tail('1/4',4,20)

    def test_algebraic_binomial_enclosures(self):
        for z in [Q(0),Q(1,10),Q(1,2),Q(9,10)]:
            p=binomial_tail('1/2',z,50)
            self.assertLessEqual(Q(p['lower'])**2,1/(1-z))
            self.assertGreaterEqual(Q(p['upper'])**2,1/(1-z))
            self.assertTrue(replay_generating(p))
        for alpha in [1,2,3]:
            p=binomial_tail(alpha,'2/3',30);truth=Q(3)**alpha
            self.assertLessEqual(Q(p['lower']),truth);self.assertGreaterEqual(Q(p['upper']),truth)
        with self.assertRaises(ValueError):binomial_tail(10,'9/10',0)

    def test_tamper_and_bad_definitions(self):
        p=RationalGF([0,1],[1,-1,-1]).packet()
        for key in ['dimension','initial','hankel_determinant','A','indexing']:
            bad=copy.deepcopy(p);bad[key]=None
            with self.assertRaises(ValueError):replay_generating(bad)
        with self.assertRaises(ValueError):RationalGF([1],[0,1])
        with self.assertRaises(TypeError):RationalGF([0.1],[1])
        with self.assertRaises(ValueError):from_recurrence([1],[])


if __name__ == '__main__':unittest.main()
