import copy
import random
import unittest
from fractions import Fraction as Q
from dynacomp import LinearModel
from dynacomp.enrich import greedy_enrichment, replay_enrichment
from dynacomp.polynomial import (compile_controlled_polynomial, replay_controlled_polynomial, controlled_step,
                                 evaluate, unpack)
from dynacomp.rational import mv, addv, matrix, vector_norm


def chain(n):
    A = [[0]*n for _ in range(n)]
    for i in range(n):
        A[i][i] = "1/2"
        if i+1 < n:
            A[i][i+1], A[i+1][i] = "1/4", "1/8"
    return LinearModel.from_rows(A, [[1]]+[[0]]*(n-1), [[0]*(n-1)+[1]])


class EnrichmentTests(unittest.TestCase):
    def test_certified_and_true_error_below_bound(self):
        m = chain(6)
        scenarios = [([0]*6, [[1]]*8), ([1, 0, -1, 0, 1, 0], [[0]]*8)]
        p = greedy_enrichment(m, scenarios, "1/20")
        self.assertEqual(p["status"], "certified")
        self.assertLess(p["reduced_dimension"], 6)
        self.assertTrue(replay_enrichment(p))
        v, w, g, h = (matrix(p[k]) for k in ("V", "W", "G", "H"))
        for (x0, inputs), bound in zip(scenarios, p["certificates"]):
            x = tuple(Q(a) for a in x0)
            z = mv(w, x)
            for t, b in enumerate(bound):
                self.assertLessEqual(vector_norm(mv(m.C, x)) - vector_norm(mv(m.C, mv(v, z))), Q(b))
                self.assertLessEqual(vector_norm(tuple(a-c for a, c in zip(mv(m.C, x), mv(m.C, mv(v, z))))), Q(b))
                if t < len(inputs):
                    x = addv(mv(m.A, x), mv(m.B, inputs[t]))
                    z = addv(mv(g, z), mv(h, inputs[t]))

    def test_zero_tolerance_reaches_exactness_and_budget_stops(self):
        m = chain(4)
        p = greedy_enrichment(m, [([1, 1, 1, 1], [[1]]*5)], 0)
        self.assertEqual(p["status"], "certified")
        self.assertTrue(replay_enrichment(p))
        q = greedy_enrichment(m, [([1, 1, 1, 1], [[1]]*5)], 0, max_dimension=1)
        self.assertEqual(q["status"], "budget_exhausted")

    def test_tampered_certificate_rejected(self):
        p = greedy_enrichment(chain(4), [([1, 0, 0, 0], [[1]]*4)], "1/10")
        p["certificates"][0][1] = str(Q(p["certificates"][0][1]) + 1)
        with self.assertRaises(ValueError):
            replay_enrichment(p)


class ControlledPolynomialTests(unittest.TestCase):
    def test_bilinear_lift_matches_trajectories(self):
        # x' = x + u y^2, y' = y/2 + v : variables (x, y, u, v)
        update = [{(1, 0, 0, 0): 1, (0, 2, 1, 0): 1}, {(0, 1, 0, 0): "1/2", (0, 0, 0, 1): 1}]
        p = compile_controlled_polynomial(update, [{(1, 0): 1}], 2)
        self.assertEqual(p["status"], "exact")
        self.assertTrue(replay_controlled_polynomial(p))
        features = [unpack(f, 4) for f in p["features"]]
        d = matrix(p["D"])
        rng = random.Random(3)
        for _ in range(10):
            x = (Q(rng.randrange(-3, 4)), Q(rng.randrange(-3, 4)))
            z = tuple(evaluate(f, x+(0, 0)) for f in features)
            for _ in range(5):
                u = (Q(rng.randrange(-2, 3)), Q(rng.randrange(-2, 3)))
                x = (x[0] + u[0]*x[1]**2, x[1]/2 + u[1])
                z = controlled_step(p, z, u)
                self.assertEqual(mv(d, z)[0], x[0])

    def test_resource_limit_and_tamper(self):
        update = [{(1, 1, 0): 1, (0, 0, 1): 1}, {(0, 1, 0): 1}]
        self.assertEqual(compile_controlled_polynomial(update, [{(1, 0): 1}], 1)["status"], "resource_limit")
        update = [{(1, 0, 0): 1, (0, 2, 1): 1}, {(0, 1, 0): "1/2"}]
        p = compile_controlled_polynomial(update, [{(1, 0): 1}], 1)
        bad = copy.deepcopy(p)
        bad["blocks"][1]["G"][1][0] = "2"
        with self.assertRaises(ValueError):
            replay_controlled_polynomial(bad)
