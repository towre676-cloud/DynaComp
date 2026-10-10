import copy
import unittest
from fractions import Fraction as Q
from dynacomp.polynomial import compile_polynomial, replay_polynomial, unpack, evaluate, compose
from dynacomp.rational import matrix, mv


class PolynomialTests(unittest.TestCase):
    def test_nonlinear_feature_closes_exactly(self):
        update = [{(1, 0): 1, (0, 2): 1}, {(0, 1): 1}]
        packet = compile_polynomial(update, [{(1, 0): 1}])
        self.assertEqual(packet["status"], "exact")
        self.assertEqual(len(packet["features"]), 2)
        self.assertTrue(replay_polynomial(packet))
        features = [unpack(p, 2) for p in packet["features"]]
        g, d = matrix(packet["G"]), matrix(packet["D"])
        for x0 in ((2, 3), (-1, -2), (Q(1, 7), Q(2, 3))):
            x = x0
            z = tuple(evaluate(p, x) for p in features)
            for _ in range(12):
                self.assertEqual(mv(d, z), (x[0],))
                x = tuple(evaluate(p, x) for p in update)
                z = mv(g, z)

    def test_unbounded_degree_reports_limit_not_impossibility(self):
        p = compile_polynomial([{(2,): 1}], [{(1,): 1}], max_degree=8)
        self.assertEqual(p["status"], "resource_limit")
        self.assertNotIn("G", p)

    def test_constant_and_affine_features(self):
        p = compile_polynomial([{(1,): 1, (0,): 1}], [{(1,): 1}])
        self.assertEqual(len(p["features"]), 2)
        self.assertTrue(replay_polynomial(p))

    def test_mutated_polynomial_packet_rejected(self):
        p = compile_polynomial([{(1,): 1, (0,): 1}], [{(1,): 1}])
        p["G"][0][0] = str(Q(p["G"][0][0])+1)
        with self.assertRaises(ValueError):
            replay_polynomial(p)

    def test_dependent_features_rejected(self):
        p = compile_polynomial([{(1,): 1}], [{(1,): 1}])
        p["features"].append(copy.deepcopy(p["features"][0]))
        with self.assertRaises(ValueError):
            replay_polynomial(p)

    def test_invalid_exponents(self):
        with self.assertRaises(ValueError):
            compile_polynomial([{(-1,): 1}], [{(1,): 1}])

    def test_round_limit(self):
        p = compile_polynomial([{(1, 0): 1, (0, 1): 1}, {(0, 1): 1}], [{(1, 0): 1}], max_rounds=0)
        self.assertEqual(p["status"], "resource_limit")
